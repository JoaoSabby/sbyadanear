#' ADASYN resampling
#'
#' @description
#' Apply ADASYN through the common exact Intel oneAPI engine.
#'
#' @param .data Data frame or tibble containing the outcome and plain numeric predictors.
#'
#' @param formula Formula outcome ~ predictors. Select existing columns only; transformations and interactions are rejected.
#'
#' @param sby_adasyn_k Positive integer number of ADASYN neighbors. Difficulty uses min(k, n - 1); interpolation uses min(k, n_min - 1). These searches run only when G > 0.
#'
#' @param sby_adasyn_ratio Nonnegative expansion relative to the original minority: G = floor(n_min * ratio). Zero disables generation. This is a reparameterization, not the beta used in the original paper. Inactive when d_th prevents ADASYN or beta is supplied instead.
#'
#' @param sby_config_max_threads Positive integer per-call thread ceiling, or -1L for detection. Capped by physical cores, CPU affinity, container quota and the hard OpenMP thread limit. Both Intel OpenMP and oneMKL receive this resolved ceiling. Local controls are restored on success or error. BLAS runs outside OpenMP regions; small BLAS calls may use fewer threads. HPC interfaces use this parameter directly. No AVX-512 requirement.
#'
#' @param sby_seed Integer seed from 0 to .Machine$integer.max; default sample.int(10e7, 1). ADASYN uses a scoped Mersenne-Twister/Inversion/Rejection RNG and restores RNGkind and .Random.seed, including after errors. Identical input, seed, parameters and numerical environment reproduce the result. Evaluating the default sample.int consumes the caller RNG; supply a seed explicitly to avoid this. NearMiss itself is deterministic and does not draw random numbers.
#'
#' @param sby_audit FALSE retains the always-present sbyaudit and sby attributes. TRUE also attaches the detailed audit. HPC interfaces always return a tibble; classic tabular interfaces return a list with sby_balanced_data when audited; matrix interfaces return lists. For recipes, the step bake method returns data with attributes; final recipes::bake(recipe) may drop them. The last audit remains in `prepared_recipe$steps[[i]]$audit_log$last`.
#'
#' @param sby_adasyn_beta Optional number from 0 to 1: G = floor((n_maj - n_min) * beta), using the original paper parameterization. Cannot be supplied together with an explicitly supplied sby_adasyn_ratio. NULL uses ratio.
#'
#' @param sby_adasyn_d_th Threshold from 0 to 1. ADASYN executes only when n_min / n_maj < d_th. Default 1; zero disables generation. Class roles are determined on the original data and remain fixed.
#'
#' @param sby_adasyn_zero_difficulty Policy when all difficulty values are zero: "error" (default) stops because the paper normalization is undefined; "uniform" explicitly requests the documented uniform-quota extension. Consult the audit for the resolved policy and whether fallback was used.
#'
#' @inherit sby_adanear_hpc details references
#'
#' @return A balanced tibble with sbyaudit and sby attributes; sby_audit = TRUE additionally attaches audit. See Details for indices and telemetry.
#' @export
sby_adasyn_hpc <- function(
  .data,
  formula,
  sby_adasyn_k = 3,
  sby_adasyn_ratio = 0.2,
  sby_config_max_threads = -1,
  sby_seed = sample.int(10e7, 1),
  sby_audit = FALSE,
  sby_adasyn_beta = NULL,
  sby_adasyn_d_th = 1,
  sby_adasyn_zero_difficulty = c("error", "uniform")
) {
  parameters <- mget(names(formals(sys.function())), envir=environment())
  sby_dispatch("adasyn", parameters, "sby_adasyn_hpc", !missing(sby_adasyn_ratio))
}
