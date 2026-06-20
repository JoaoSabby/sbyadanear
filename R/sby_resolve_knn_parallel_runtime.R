#' Resolver runtime efetivo do backend de paralelismo KNN
#'
#' @param sby_knn_parallel_backend Backend de paralelismo KNN validado
#'
#' @return String descritiva do runtime efetivo
#'
#' @noRd
sby_resolve_knn_parallel_runtime <- function(sby_knn_parallel_backend){
  sby_knn_parallel_backend <- sby_validate_knn_parallel_backend(
    sby_knn_parallel_backend = sby_knn_parallel_backend
  )

  if(identical(sby_knn_parallel_backend, "parallel")){
    return("parallel")
  }

  sby_uses_tbb <- isTRUE(sby_call_native("rcpp_parallel_uses_tbb_c"))
  if(isTRUE(sby_uses_tbb)){
    return("RcppParallel::TBB")
  }

  return("RcppParallel::TinyThread")
}
