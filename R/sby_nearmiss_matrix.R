#' NearMiss resampling
#'
#' @description
#' Apply NearMiss through the common exact Intel oneAPI engine.
#'
#' @param sby_x_matrix Dense numeric matrix. Double columns containing only whole numbers are recognized as integer domains.
#'
#' @param sby_y_vector Outcome with exactly two observed classes and no missing values. Factor levels and ordering are preserved.
#'
#' @param sby_nearmiss_ratio Nonnegative majority retention relative to the expanded minority: min(n_maj, floor((n_min + G) * ratio)). Zero disables NearMiss. A positive ratio rounding to zero produces an informative error. If the target retains the entire majority, neighbor scoring is skipped.
#'
#' @param sby_knn_under_k Positive integer NearMiss score neighbor count; equivalent to sby_nearmiss_k in HPC interfaces.
#'
#' @param sby_seed Integer seed from 0 to .Machine$integer.max; default sample.int(10e7, 1). ADASYN uses a scoped Mersenne-Twister/Inversion/Rejection RNG and restores RNGkind and .Random.seed, including after errors. Identical input, seed, parameters and numerical environment reproduce the result. Evaluating the default sample.int consumes the caller RNG; supply a seed explicitly to avoid this. NearMiss itself is deterministic and does not draw random numbers.
#'
#' @param sby_audit FALSE retains the always-present sbyaudit and sby attributes. TRUE also attaches the detailed audit. HPC interfaces always return a tibble; classic tabular interfaces return a list with sby_balanced_data when audited; matrix interfaces return lists. For recipes, the step bake method returns data with attributes; final recipes::bake(recipe) may drop them. The last audit remains in `prepared_recipe$steps[[i]]$audit_log$last`.
#'
#' @param sby_audit_level Compatibility selector: none, light or full. light and full enable the same detailed audit in this version.
#'
#' @param sby_return_index Include sby_retained_index for retained original input rows; synthetic rows do not have an original row index.
#'
#' @param sby_return_scaled Include an additional standardized matrix. Tabular primary output remains in the original scale; matrix results use sby_x_scaled.
#'
#' @param sby_return_original_scale TRUE returns sby_x_matrix in the original scale; FALSE returns it standardized. Integer-domain restrictions apply before this optional presentation transformation.
#'
#' @param sby_scaling_info Optional list with finite numeric sby_center and positive sby_scale vectors, one value per predictor in input-column order. Aliases centers/scales and means/sds are accepted. Overrides automatic scaling.
#'
#' @param sby_input_already_scaled TRUE reconstructs original-scale values using supplied scaling information before processing. Such reconstructed originals cannot be guaranteed bit-identical to originals not supplied to the function.
#'
#' @param sby_fixed_minority_label Optional observed label defining the preserved class in classic NearMiss interfaces.
#'
#' @param sby_fixed_majority_label Optional observed label identifying the other class. Must agree with the preserved-class choice.
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
#' @param sby_memory_guard Whether to enforce the conservative dense-buffer estimate. This estimate is not an exact total-process memory ceiling.
#'
#' @param sby_config_max_threads Positive integer per-call thread ceiling, or -1L for detection. Capped by physical cores, CPU affinity, container quota and the hard OpenMP thread limit. Both Intel OpenMP and oneMKL receive this resolved ceiling. Local controls are restored on success or error. BLAS runs outside OpenMP regions; small BLAS calls may use fewer threads. HPC interfaces use this parameter directly. No AVX-512 requirement.
#'
#' @param nearmiss_model Integer 1L, 2L or 3L; default 3L. NearMiss-1 retains majority rows with the SMALLEST mean distance to their K NEAREST rare rows. NearMiss-2 retains rows with the SMALLEST mean distance to their K FARTHEST rare rows. NearMiss-3 first takes the union of the M nearest majority rows of EACH rare row, then retains candidates with the LARGEST mean distance to their K nearest rare rows. M = sby_nearmiss_m and K = sby_nearmiss_k (or sby_knn_under_k) are independent. Continuous synthetic rare rows participate in the combined pipeline before domain restoration. Candidate shortage retains all candidates and warns; it does not fill from non-candidates. Ties use increasing original row index. Other values produce an informative error. The model is validated even when NearMiss is inactive.
#'
#' @param sby_nearmiss_m Positive integer preselection count M for NearMiss-3; default 3L, capped at majority size. Does not alter K. Validated for every model, but used in neighbor selection only for model 3 when NearMiss executes.
#'
#' @inherit sby_adanear_hpc details references
#'
#' @return A list with sby_x_matrix, sby_y_vector, sby_scaling_info, class distributions/ratios and diagnostics. Additional indices/standardized output follows its flags. Attributes follow Details.
#' @export
sby_nearmiss_matrix <- function(
  sby_x_matrix,
  sby_y_vector,
  sby_nearmiss_ratio = 1,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10e7, 1),
  sby_audit = FALSE,
  sby_audit_level = c("none", "light", "full"),
  sby_return_index = TRUE,
  sby_return_scaled = FALSE,
  sby_return_original_scale = TRUE,
  sby_scaling_info = NULL,
  sby_input_already_scaled = FALSE,
  sby_fixed_minority_label = NULL,
  sby_fixed_majority_label = NULL,
  sby_knn_algorithm = "auto",
  sby_knn_engine = "auto",
  sby_knn_distance_metric = "euclidean",
  sby_knn_workers = 1L,
  sby_knn_parallel_backend = "parallel",
  sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L,
  sby_knn_query_chunk_size = 1000L,
  sby_memory_guard = TRUE,
  sby_config_max_threads = -1L,
  nearmiss_model = 3L,
  sby_nearmiss_m = 3L
) {
  parameters <- mget(names(formals(sys.function())), envir=environment())
  sby_dispatch("nearmiss", parameters, "sby_nearmiss_matrix", FALSE)
}
