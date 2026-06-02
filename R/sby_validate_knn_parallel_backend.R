#' Validar backend de paralelismo KNN
#'
#' @param sby_knn_parallel_backend Backend solicitado para paralelismo KNN
#'
#' @return String validada com o backend de paralelismo
#' @noRd
sby_validate_knn_parallel_backend <- function(sby_knn_parallel_backend){
  sby_knn_parallel_backend <- match.arg(
    arg = sby_knn_parallel_backend,
    choices = c("parallel", "RcppParallel")
  )

  if(identical(sby_knn_parallel_backend, "RcppParallel") &&
     !requireNamespace(package = "RcppParallel", quietly = TRUE)){
    sby_adanear_abort(
      sby_message = "'sby_knn_parallel_backend = RcppParallel' requer o pacote RcppParallel"
    )
  }

  return(sby_knn_parallel_backend)
}
