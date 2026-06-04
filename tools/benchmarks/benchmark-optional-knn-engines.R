#!/usr/bin/env Rscript
# Benchmark comparativo dos engines KNN opcionais na API interna comum.

suppressPackageStartupMessages({
  library(sbyadanear)
})

make_case <- function(n, p, query_n, seed){
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  query <- x[seq_len(query_n), , drop = FALSE]
  list(x = x, query = query)
}

run_engine <- function(engine, x, query, k, workers){
  gc()
  start <- proc.time()["elapsed"]
  result <- sbyadanear:::sby_get_knnx(
    x,
    query,
    k,
    "auto",
    engine,
    "euclidean",
    workers,
    16L,
    200L,
    sby_knn_query_chunk_size = 500L,
    sby_knn_return = "both"
  )
  elapsed <- proc.time()["elapsed"] - start
  data.frame(
    engine = engine,
    elapsed_sec = as.numeric(elapsed),
    query_rows = nrow(query),
    reference_rows = nrow(x),
    columns = ncol(x),
    k = k,
    workers = workers,
    index_rows = nrow(result$nn.index),
    index_cols = ncol(result$nn.index)
  )
}

case <- make_case(n = 5000L, p = 50L, query_n = 1000L, seed = 20260604L)
engines <- c("FNN", "KernelKnn", "bigKNN")
available <- vapply(engines, function(engine){
  if(identical(engine, "FNN")){
    return(requireNamespace("FNN", quietly = TRUE))
  }
  if(identical(engine, "KernelKnn")){
    return(requireNamespace("KernelKnn", quietly = TRUE))
  }
  requireNamespace("bigKNN", quietly = TRUE) && requireNamespace("bigmemory", quietly = TRUE)
}, logical(1L))
engines <- engines[available]

if(!length(engines)){
  stop("Nenhum engine opcional disponivel para benchmark")
}

out <- do.call(rbind, lapply(engines, run_engine, x = case$x, query = case$query, k = 5L, workers = 1L))
print(out)
