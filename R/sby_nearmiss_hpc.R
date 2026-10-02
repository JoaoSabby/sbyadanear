#' NearMiss resampling
#'
#' @description
#' Apply NearMiss through the common exact Intel oneAPI engine.
#'
#' @param .data Data frame or tibble containing the outcome and plain numeric predictors.
#'
#' @param formula Formula outcome ~ predictors. Select existing columns only; transformations and interactions are rejected.
#'
#' @param sby_nearmiss_k Positive integer number K of Euclidean distances d averaged for NearMiss, capped at the expanded rare-class size. Used only when majority retention actually reduces the data.
#'
#' @param sby_nearmiss_ratio Nonnegative majority retention relative to the expanded minority: min(n_maj, floor((n_min + G) * ratio)). Zero disables NearMiss. A positive ratio rounding to zero produces an informative error. If the target retains the entire majority, neighbor scoring is skipped.
#'
#' @param sby_config_max_threads Positive integer per-call thread ceiling, or -1L for detection. Capped by physical cores, CPU affinity, container quota and the hard OpenMP thread limit. Both Intel OpenMP and oneMKL receive this resolved ceiling. Local controls are restored on success or error. BLAS runs outside OpenMP regions; small BLAS calls may use fewer threads. HPC interfaces use this parameter directly. No AVX-512 requirement.
#'
#' @param sby_seed Integer seed from 0 to .Machine$integer.max; default sample.int(10e7, 1). ADASYN uses a scoped Mersenne-Twister/Inversion/Rejection RNG and restores RNGkind and .Random.seed, including after errors. Identical input, seed, parameters and numerical environment reproduce the result. Evaluating the default sample.int consumes the caller RNG; supply a seed explicitly to avoid this. NearMiss itself is deterministic and does not draw random numbers.
#'
#' @param sby_audit FALSE retains the always-present sbyaudit and sby attributes. TRUE also attaches the detailed audit. HPC interfaces always return a tibble; classic tabular interfaces return a list with sby_balanced_data when audited; matrix interfaces return lists. For recipes, the step bake method returns data with attributes; final recipes::bake(recipe) may drop them. The last audit remains in `prepared_recipe$steps[[i]]$audit_log$last`.
#'
#' @param nearmiss_model Integer 1L, 2L or 3L; default 3L. NearMiss-1 retains majority rows with the SMALLEST mean distance to their K NEAREST rare rows. NearMiss-2 retains rows with the SMALLEST mean distance to their K FARTHEST rare rows. NearMiss-3 first takes the union of the M nearest majority rows of EACH rare row, then retains candidates with the LARGEST mean distance to their K nearest rare rows. M = sby_nearmiss_m and K = sby_nearmiss_k (or sby_knn_under_k) are independent. Continuous synthetic rare rows participate in the combined pipeline before domain restoration. Candidate shortage retains all candidates and warns; it does not fill from non-candidates. Ties use increasing original row index. Other values produce an informative error. The model is validated even when NearMiss is inactive.
#'
#' @param sby_nearmiss_m Positive integer preselection count M for NearMiss-3; default 3L, capped at majority size. Does not alter K. Validated for every model, but used in neighbor selection only for model 3 when NearMiss executes.
#'
#' @inherit sby_adanear_hpc details references
#'
#' @return A balanced tibble with sbyaudit and sby attributes; sby_audit = TRUE additionally attaches audit. See Details for indices and telemetry.
#' @export
sby_nearmiss_hpc <- function(
  .data,
  formula,
  sby_nearmiss_k = 7,
  sby_nearmiss_ratio = 1,
  sby_config_max_threads = -1,
  sby_seed = sample.int(10e7, 1),
  sby_audit = FALSE,
  nearmiss_model = 3L,
  sby_nearmiss_m = 3L
) {
  parameters <- mget(names(formals(sys.function())), envir=environment())
  sby_dispatch("nearmiss", parameters, "sby_nearmiss_hpc", FALSE)
}
