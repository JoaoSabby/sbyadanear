#' ADASYN and NearMiss resampling
#'
#' @description
#' Apply ADASYN and NearMiss through the common exact Intel oneAPI engine.
#'
#' @param sby_formula Formula selecting the outcome and existing numeric predictor columns.
#'
#' @param sby_data Data frame containing the outcome and plain numeric predictors.
#'
#' @param sby_adasyn_ratio Nonnegative expansion relative to the original minority: G = floor(n_min * ratio). Zero disables generation. This is a reparameterization, not the beta used in the original paper. Inactive when d_th prevents ADASYN or beta is supplied instead.
#'
#' @param sby_nearmiss_ratio Nonnegative majority retention relative to the expanded minority: min(n_maj, floor((n_min + G) * ratio)). Zero disables NearMiss. A positive ratio rounding to zero produces an informative error. If the target retains the entire majority, neighbor scoring is skipped.
#'
#' @param sby_knn_over_k Positive integer ADASYN neighbor count; equivalent to sby_adasyn_k in HPC interfaces.
#'
#' @param sby_knn_under_k Positive integer NearMiss score neighbor count; equivalent to sby_nearmiss_k in HPC interfaces.
#'
#' @param sby_seed Integer seed from 0 to .Machine$integer.max; default sample.int(10e7, 1). ADASYN uses a scoped Mersenne-Twister/Inversion/Rejection RNG and restores RNGkind and .Random.seed, including after errors. Identical input, seed, parameters and numerical environment reproduce the result. Evaluating the default sample.int consumes the caller RNG; supply a seed explicitly to avoid this. NearMiss itself is deterministic and does not draw random numbers.
#'
#' @param sby_audit FALSE retains the always-present sbyaudit and sby attributes. TRUE also attaches the detailed audit. HPC interfaces always return a tibble; classic tabular interfaces return a list with sby_balanced_data when audited; matrix interfaces return lists. For recipes, the step bake method returns data with attributes; final recipes::bake(recipe) may drop them. The last audit remains in `prepared_recipe$steps[[i]]$audit_log$last`.
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
#' @param nearmiss_model Integer 1L, 2L or 3L; default 3L. NearMiss-1 retains majority rows with the SMALLEST mean distance to their K NEAREST rare rows. NearMiss-2 retains rows with the SMALLEST mean distance to their K FARTHEST rare rows. NearMiss-3 first takes the union of the M nearest majority rows of EACH rare row, then retains candidates with the LARGEST mean distance to their K nearest rare rows. M = sby_nearmiss_m and K = sby_nearmiss_k (or sby_knn_under_k) are independent. Continuous synthetic rare rows participate in the combined pipeline before domain restoration. Candidate shortage retains all candidates and warns; it does not fill from non-candidates. Ties use increasing original row index. Other values produce an informative error. The model is validated even when NearMiss is inactive.
#'
#' @param sby_nearmiss_m Positive integer preselection count M for NearMiss-3; default 3L, capped at majority size. Does not alter K. Validated for every model, but used in neighbor selection only for model 3 when NearMiss executes.
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
sby_adanear <- function(
  sby_formula,
  sby_data,
  sby_adasyn_ratio = 0.2,
  sby_nearmiss_ratio = 1,
  sby_knn_over_k = 5L,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10e7, 1),
  sby_audit = FALSE,
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
  nearmiss_model = 3L,
  sby_nearmiss_m = 3L,
  sby_adasyn_beta = NULL,
  sby_adasyn_d_th = 1,
  sby_adasyn_zero_difficulty = c("error", "uniform")
) {
  parameters <- mget(names(formals(sys.function())), envir=environment())
  sby_dispatch("adanear", parameters, "sby_adanear", !missing(sby_adasyn_ratio))
}
