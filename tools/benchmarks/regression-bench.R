# Regression benchmark for the instenginer matrix and tabular APIs.
#
# Run with:
#   Rscript tools/benchmarks/regression-bench.R
#
# The script measures wall-clock time for representative scenarios and
# writes a CSV to tools/benchmarks/results/regression-bench-<timestamp>.csv.
# When tools/benchmarks/baseline.csv exists, the script also prints a
# side-by-side comparison so CI can flag regressions.
#
# Scenarios cover:
#   - n / p / minority_frac sweeps for ADASYN, NearMiss-1 and ADANEAR.
#   - FNN exact (kd_tree + brute) and native BLAS brute paths.
#   - RcppHNSW path for cosine metric and for large high-dim bases.
#   - Tabular API overhead vs matrix API.
#
# Each measurement is taken as the median of three timed runs after a
# single warmup run.

suppressPackageStartupMessages({
  library(instenginer)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_data <- function(n, p, minority_frac, seed = 42L){
  set.seed(seed)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  storage.mode(x) <- "double"
  n_min <- max(2L, as.integer(round(n * minority_frac)))
  y <- factor(
    c(rep("min", n_min), rep("maj", n - n_min)),
    levels = c("min", "maj")
  )
  dat <- as.data.frame(x)
  names(dat) <- paste0("x", seq_len(p))
  dat$TARGET <- y
  list(x = x, y = y, dat = dat, n = n, p = p, n_min = n_min)
}

# Time a thunk. We accept a zero-arg function (not a pre-evaluated
# expression) because R semantics evaluate function arguments eagerly,
# which would make the system.time() call measure essentially nothing.
time_call <- function(thunk, repeats = 3L){
  # One warmup run (not timed).
  invisible(thunk())
  timings <- numeric(repeats)
  for(i in seq_len(repeats)){
    gc(verbose = FALSE)
    t <- system.time(thunk())[["elapsed"]]
    timings[i] <- t
  }
  list(
    median = stats::median(timings),
    min = min(timings),
    max = max(timings)
  )
}

bench_row <- function(scenario, fn, args){
  thunk <- function() do.call(fn, args)
  res <- time_call(thunk, repeats = 3L)
  data.frame(
    scenario = scenario,
    fn = fn,
    median_s = round(res$median, 4L),
    min_s = round(res$min, 4L),
    max_s = round(res$max, 4L),
    stringsAsFactors = FALSE
  )
}

run_scenarios <- function(){
  scenarios <- list(
    list(label = "small_n1k_p10", n = 1000L, p = 10L, frac = 0.05),
    list(label = "med_n5k_p20",  n = 5000L, p = 20L, frac = 0.05),
    list(label = "med_n10k_p50", n = 10000L, p = 50L, frac = 0.02),
    list(label = "highp_n5k_p200", n = 5000L, p = 200L, frac = 0.05)
  )

  rows <- list()
  for(sc in scenarios){
    cat(sprintf("\n=== Scenario: %s (n=%d, p=%d, frac=%.3f) ===\n",
                sc$label, sc$n, sc$p, sc$frac))
    d <- make_data(sc$n, sc$p, sc$frac)

    # ADASYN matrix - FNN kd_tree (small p) or FNN brute (large p)
    rows[[length(rows) + 1L]] <- bench_row(
      sc$label, "sby_adasyn_matrix",
      list(d$x, d$y, sby_over_ratio = 0.5, sby_seed = 1L,
           sby_knn_engine = "FNN")
    )

    # NearMiss matrix - FNN
    rows[[length(rows) + 1L]] <- bench_row(
      sc$label, "sby_nearmiss_matrix",
      list(d$x, d$y, sby_under_ratio = 0.5, sby_seed = 1L,
           sby_knn_engine = "FNN")
    )

    # ADANEAR matrix - FNN
    rows[[length(rows) + 1L]] <- bench_row(
      sc$label, "sby_adanear_matrix",
      list(d$x, d$y, sby_over_ratio = 0.3, sby_under_ratio = 0.5,
           sby_seed = 1L, sby_knn_engine = "FNN")
    )

    # ADASYN tabular wrapper (overhead measurement)
    rows[[length(rows) + 1L]] <- bench_row(
      paste0(sc$label, "_tabular"), "sby_adasyn",
      list(TARGET ~ ., d$dat, sby_over_ratio = 0.5, sby_seed = 1L,
           sby_knn_engine = "FNN")
    )

    # ADASYN with native BLAS brute KNN
    options(instenginer.sby_use_native_brute = TRUE)
    rows[[length(rows) + 1L]] <- bench_row(
      paste0(sc$label, "_native_brute"), "sby_adasyn_matrix",
      list(d$x, d$y, sby_over_ratio = 0.5, sby_seed = 1L,
           sby_knn_engine = "FNN", sby_knn_algorithm = "brute")
    )

    # ADASYN with FNN brute (no native BLAS) for comparison
    options(instenginer.sby_use_native_brute = FALSE)
    rows[[length(rows) + 1L]] <- bench_row(
      paste0(sc$label, "_fnn_brute"), "sby_adasyn_matrix",
      list(d$x, d$y, sby_over_ratio = 0.5, sby_seed = 1L,
           sby_knn_engine = "FNN", sby_knn_algorithm = "brute")
    )
    options(instenginer.sby_use_native_brute = TRUE)
  }

  # HNSW cosine: separate scenario because it needs n large enough to matter.
  cat("\n=== Scenario: hnsw_cosine_n5k_p50 ===\n")
  d <- make_data(5000L, 50L, 0.05)
  rows[[length(rows) + 1L]] <- bench_row(
    "hnsw_cosine_n5k_p50", "sby_adasyn_matrix",
    list(d$x, d$y, sby_over_ratio = 0.5, sby_seed = 1L,
         sby_knn_engine = "RcppHNSW", sby_knn_distance_metric = "cosine")
  )

  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Run and report
# ---------------------------------------------------------------------------

results <- run_scenarios()

# Print formatted table.
cat("\n\n=== Results ===\n")
print(results, row.names = FALSE)

# Write CSV with timestamp.
out_dir <- file.path("tools", "benchmarks", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_file <- file.path(out_dir, sprintf("regression-bench-%s.csv", ts))
write.csv(results, out_file, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_file))

# Compare against baseline if it exists.
baseline_file <- file.path("tools", "benchmarks", "baseline.csv")
if(file.exists(baseline_file)){
  baseline <- read.csv(baseline_file, stringsAsFactors = FALSE)
  cat("\n=== Comparison against baseline ===\n")
  cmp <- merge(
    results[, c("scenario", "fn", "median_s")],
    baseline[, c("scenario", "fn", "median_s")],
    by = c("scenario", "fn"),
    suffixes = c("_new", "_baseline"),
    all = FALSE
  )
  cmp$ratio <- round(cmp$median_s_new / cmp$median_s_baseline, 3L)
  cmp$delta_pct <- round(100 * (cmp$median_s_new - cmp$median_s_baseline) /
                         cmp$median_s_baseline, 1L)
  print(cmp, row.names = FALSE)
  regressions <- cmp[!is.na(cmp$ratio) & cmp$ratio > 1.10, ]
  if(nrow(regressions) > 0L){
    cat("\nREGRESSIONS DETECTED (>10% slower):\n")
    print(regressions, row.names = FALSE)
  }else{
    cat("\nNo regressions vs baseline (threshold 10%).\n")
  }
}else{
  cat(sprintf("\nNo baseline at %s; consider saving this run as baseline.csv.\n",
              baseline_file))
}
