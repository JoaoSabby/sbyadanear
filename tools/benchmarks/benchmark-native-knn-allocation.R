#!/usr/bin/env Rscript
# Benchmark focused on native KNN allocation pressure in brute paths.

suppressPackageStartupMessages({
  library(sbyadanear)
})

make_case <- function(n, p, minority_ratio, seed){
  set.seed(seed)
  minority_n <- max(2L, floor(n * minority_ratio))
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  y <- factor(c(rep("minor", minority_n), rep("major", n - minority_n)))
  list(x = x, y = y)
}

run_one <- function(name, n, p, minority_ratio, workers){
  dat <- make_case(n, p, minority_ratio, seed = 100 + n + p)
  gc()
  start <- proc.time()["elapsed"]
  result <- sby_nearmiss_matrix(
    sby_x_matrix = dat$x,
    sby_y_vector = dat$y,
    sby_nearmiss_ratio = 1,
    sby_knn_under_k = 5L,
    sby_seed = 1L,
    sby_knn_algorithm = "brute",
    sby_knn_engine = "FNN",
    sby_knn_workers = workers,
    sby_knn_query_chunk_size = 500L
  )
  elapsed <- proc.time()["elapsed"] - start
  data.frame(
    case = name,
    n = n,
    p = p,
    minority_ratio = minority_ratio,
    workers = workers,
    elapsed_sec = as.numeric(elapsed),
    output_rows = nrow(result$sby_x_matrix),
    minor_out = as.integer(result$sby_output_class_distribution["minor"]),
    major_out = as.integer(result$sby_output_class_distribution["major"])
  )
}

cases <- list(
  small_low_dim = list(n = 1000L, p = 10L, minority_ratio = 0.10),
  medium_low_dim = list(n = 100000L, p = 10L, minority_ratio = 0.05),
  small_high_dim = list(n = 1000L, p = 200L, minority_ratio = 0.10)
)

workers <- unique(c(1L, max(1L, min(2L, parallel::detectCores(logical = FALSE)))))
out <- do.call(rbind, lapply(names(cases), function(case_name){
  spec <- cases[[case_name]]
  do.call(rbind, lapply(workers, function(worker_count){
    run_one(case_name, spec$n, spec$p, spec$minority_ratio, worker_count)
  }))
}))

print(out)
