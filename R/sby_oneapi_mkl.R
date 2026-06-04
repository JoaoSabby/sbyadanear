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
#' @return Invisible TRUE when completed successfully.
#' @noRd
sby_configure_blas_threads <- function(sby_workers){
  sby_cfg <- sby_resolve_oneapi_mkl()

  if(!isTRUE(sby_cfg$enabled)){
    return(invisible(TRUE))
  }

  sby_target_threads <- if(isTRUE(sby_workers > 1L)) 1L else
    max(1L, as.integer(parallel::detectCores(logical = FALSE) %||% 1L))

  Sys.setenv(
    MKL_NUM_THREADS = as.character(sby_target_threads),
    OMP_NUM_THREADS = as.character(sby_target_threads),
    MKL_DYNAMIC = "FALSE"
  )

  invisible(TRUE)
}

`%||%` <- function(x, y) if(is.null(x) || is.na(x)) y else x
