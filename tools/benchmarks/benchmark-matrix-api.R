# Benchmarks manuais da API matricial industrial do sbyadanear.
# Execute com: Rscript tools/benchmarks/benchmark-matrix-api.R

suppressPackageStartupMessages({
  library(sbyadanear)
})

make_data <- function(n = 2000L, p = 20L, minority = 0.15){
  set.seed(42)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  storage.mode(x) <- "double"
  y <- factor(ifelse(seq_len(n) <= ceiling(n * minority), "minor", "major"), levels = c("minor", "major"))
  dat <- as.data.frame(x)
  names(dat) <- paste0("x", seq_len(p))
  dat$TARGET <- y
  list(x = x, y = y, dat = dat)
}

time_call <- function(label, expr){
  gc()
  t <- system.time(out <- force(expr))
  cat(label, "elapsed=", t[["elapsed"]], "rows=", nrow(out$sby_x_matrix %||% out), "\n")
  invisible(out)
}

`%||%` <- function(x, y) if(is.null(x)) y else x

d <- make_data()

time_call("adasyn_matrix", sby_adasyn_matrix(d$x, d$y, sby_seed = 1, sby_knn_engine = "FNN"))
time_call("nearmiss_matrix", sby_nearmiss_matrix(d$x, d$y, sby_seed = 1, sby_knn_engine = "FNN"))
time_call("adanear_matrix", sby_adanear_matrix(d$x, d$y, sby_seed = 1, sby_knn_engine = "FNN"))

gc()
print(system.time(sbyadanear:::sby_get_knnx(d$x, d$x[1:500, , drop = FALSE], 5L, "kd_tree", "FNN", "euclidean", 1L, 16L, 200L, sby_knn_return = "both")))
print(system.time(sbyadanear:::sby_get_knnx(d$x, d$x[1:500, , drop = FALSE], 5L, "kd_tree", "FNN", "euclidean", 1L, 16L, 200L, sby_knn_return = "index")))
print(system.time(sbyadanear:::sby_get_knnx(d$x, d$x[1:500, , drop = FALSE], 5L, "kd_tree", "FNN", "euclidean", 1L, 16L, 200L, sby_knn_return = "dist")))
