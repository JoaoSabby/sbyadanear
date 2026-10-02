#' Intel oneAPI runtime diagnostics
#' @description Report resolved libraries, oneMKL version, CPU affinity,
#'   physical/logical cores and container quota. No AVX-512 requirement.
#' @return A list containing runtime and the automatic thread plan. Configured
#'   MKL threads do not measure a BLAS team; use MKL_VERBOSE=1 to observe NThr.
#' @export
sby_hpc_cpu_report <- function() {
  runtime <- sby_native("sby_runtime_cpp")
  c(list(runtime=runtime), sby_thread_plan(-1L,runtime))
}
