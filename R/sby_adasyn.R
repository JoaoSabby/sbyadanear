#' ADASYN resampling
#'
#' @description
#' Apply ADASYN through the common exact Intel oneAPI engine.
#'
#' @param sby_formula Formula selecting the outcome and existing numeric predictor columns.
#'
#' @param sby_data Data frame containing the outcome and plain numeric predictors.
#'
#' @param sby_adasyn_ratio Nonnegative expansion relative to the original minority: G = floor(n_min * ratio). Zero disables generation. This is a reparameterization, not the beta used in the original paper. Inactive when d_th prevents ADASYN or beta is supplied instead.
#'
#' @param sby_knn_over_k Positive integer ADASYN neighbor count; equivalent to sby_adasyn_k in HPC interfaces.
#'
#' @param sby_seed Integer seed from 0 to .Machine$integer.max; default sample.int(10e7, 1). ADASYN uses a scoped Mersenne-Twister/Inversion/Rejection RNG and restores RNGkind and .Random.seed, including after errors. Identical input, seed, parameters and numerical environment reproduce the result. Evaluating the default sample.int consumes the caller RNG; supply a seed explicitly to avoid this. NearMiss itself is deterministic and does not draw random numbers.
#'
#' @param sby_audit FALSE retains the always-present sbyaudit and sby attributes. TRUE also attaches the detailed audit. HPC interfaces always return a tibble; classic tabular interfaces return a list with sby_balanced_data when audited; matrix interfaces return lists. For recipes, the step bake method returns data with attributes; final recipes::bake(recipe) may drop them. The last audit remains in `prepared_recipe$steps[[i]]$audit_log$last`.
#'
#' @param sby_return_scaled Include an additional standardized matrix. Tabular primary output remains in the original scale; matrix results use sby_x_scaled.
#'
#' @param sby_restore_types Must be TRUE. Binary synthetics use threshold >= 0.5; integer domains use round() with ties to even; every synthetic predictor is clamped to its original minimum and maximum. Original rows are copied without rounding.
#'
#' @param sby_knn_algorithm Compatibility selector: only auto or brute. Every route uses the common exact neighbor engine; alternative trees are not implemented.
#'
#' @param sby_knn_engine Compatibility selector: only auto or native, both using Intel oneAPI. Approximate or external engines are rejected.
#'
#' @param sby_knn_distance_metric Only euclidean is supported by this scientific contract. Other metrics are rejected.
#'
#' @param sby_knn_workers Classic-interface thread ceiling when sby_config_max_threads = -1L; default 1L. Supply -1L or a positive integer. Validated even when an explicit sby_config_max_threads takes precedence. HPC uses its own sby_config_max_threads.
#'
#' @param sby_knn_parallel_backend Validated legacy selector: parallel or RcppParallel. Both map to Intel OpenMP/oneMKL, without fork or TBB.
#'
#' @param sby_knn_hnsw_m Legacy compatibility parameter. Only its default 16L is accepted; HNSW is not executed.
#'
#' @param sby_knn_hnsw_ef Legacy compatibility parameter. Only its default 200L is accepted; HNSW is not executed.
#'
#' @param sby_knn_query_chunk_size Positive integer requested query tile size, capped at 128 to bound the distance buffer. Does not change neighbor geometry.
#'
#' @param sby_config_max_threads Positive integer per-call thread ceiling, or -1L for detection. Capped by physical cores, CPU affinity, container quota and the hard OpenMP thread limit. Both Intel OpenMP and oneMKL receive this resolved ceiling. Local controls are restored on success or error. BLAS runs outside OpenMP regions; small BLAS calls may use fewer threads. HPC interfaces use this parameter directly. No AVX-512 requirement.
#'
#' @param sby_adasyn_beta Optional number from 0 to 1: G = floor((n_maj - n_min) * beta), using the original paper parameterization. Cannot be supplied together with an explicitly supplied sby_adasyn_ratio. NULL uses ratio.
#'
#' @param sby_adasyn_d_th Threshold from 0 to 1. ADASYN executes only when n_min / n_maj < d_th. Default 1; zero disables generation. Class roles are determined on the original data and remain fixed.
#'
#' @param sby_adasyn_zero_difficulty Policy when all difficulty values are zero: "error" (default) stops because the paper normalization is undefined; "uniform" explicitly requests the documented uniform-quota extension. Consult the audit for the resolved policy and whether fallback was used.
#'
#' @inherit sby_adanear_hpc details references
#'
#' @return A tibble, or a list with sby_balanced_data and sby_scaling_info when sby_audit = TRUE or sby_return_scaled = TRUE. sbyaudit and sby are always attached; audit follows sby_audit. See Details.
#' @export
sby_adasyn <- function(
  sby_formula,
  sby_data,
  sby_adasyn_ratio = 0.2,
  sby_knn_over_k = 5L,
  sby_seed = sample.int(10e7, 1),
  sby_audit = FALSE,
  sby_return_scaled = FALSE,
  sby_restore_types = TRUE,
  sby_knn_algorithm = "auto",
  sby_knn_engine = "auto",
  sby_knn_distance_metric = "euclidean",
  sby_knn_workers = 1L,
  sby_knn_parallel_backend = "parallel",
  sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L,
  sby_knn_query_chunk_size = 1000L,
  sby_config_max_threads = -1L,
  sby_adasyn_beta = NULL,
  sby_adasyn_d_th = 1,
  sby_adasyn_zero_difficulty = c("error", "uniform")
) {
  parameters <- mget(names(formals(sys.function())), envir=environment())
  sby_dispatch("adasyn", parameters, "sby_adasyn", !missing(sby_adasyn_ratio))
}
