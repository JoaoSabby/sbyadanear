#' Resolve Intel oneAPI MKL configuration
#'
#' @details
#' Uses automatic defaults to minimize user-facing tuning.
#'
#' @return List with enable flag and suggested thread count.
#' @noRd
sby_resolve_oneapi_mkl <- function(){
  # Enable automatically when MKL runtime appears available.
  sby_has_mkl_runtime <- nzchar(Sys.getenv("MKLROOT", unset = "")) ||
    nzchar(Sys.getenv("ONEAPI_ROOT", unset = ""))

  # Keep only one optional override for advanced users.
  sby_mode <- getOption("sbyadanear.perf_mode", "auto")

  sby_enabled <- if(identical(sby_mode, "manual")){
    FALSE
  }else{
    sby_has_mkl_runtime
  }

  list(enabled = sby_enabled, threads = NA_integer_)
}

#' Apply BLAS/MKL thread hints
#'
#' @param sby_workers Number of KNN workers in R.
#' @return Previous MKL_NUM_THREADS and OMP_NUM_THREADS values.
#' @noRd
sby_configure_blas_threads <- function(sby_workers){
  sby_cfg <- sby_resolve_oneapi_mkl()
  sby_previous <- sby_hpc_capture_env()

  if(!isTRUE(sby_cfg$enabled)){
    return(sby_previous)
  }

  sby_target_threads <- if(isTRUE(sby_workers > 1L)) 1L else
    max(1L, as.integer(parallel::detectCores(logical = FALSE) %||% 1L))

  Sys.setenv(
    MKL_NUM_THREADS = as.character(sby_target_threads),
    OMP_NUM_THREADS = as.character(sby_target_threads)
  )

  return(sby_previous)
}

#' @title Restaurar threads BLAS e MKL
#' @usage sby_restore_blas_threads(sby_previous)
#' @description
#' Restaura o ambiente de paralelismo capturado antes de uma chamada nativa que
#' ajustou limites de BLAS, oneMKL ou OpenMP.
#' @details
#' Esta função interna delega o restauro para `sby_hpc_restore_env`, mantendo uma
#' API curta para pares `configure` e `restore`.
#' @param sby_previous Estado anterior produzido por `sby_configure_blas_threads`.
#' @return Resultado invisível de `sby_hpc_restore_env`.
#' @seealso sby_configure_blas_threads, sby_hpc_restore_env
#' @examples
#' \dontrun{
#' env <- sbyadanear:::sby_configure_blas_threads(1L)
#' sbyadanear:::sby_restore_blas_threads(env)
#' }
#' @keywords internal
sby_restore_blas_threads <- function(sby_previous){
  sby_hpc_restore_env(sby_previous)
}

`%||%` <- function(x, y) if(is.null(x) || is.na(x)) y else x
