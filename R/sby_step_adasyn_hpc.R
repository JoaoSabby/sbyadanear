#' ADASYN recipe step
#'
#' @description
#' Apply ADASYN through the common exact Intel oneAPI engine.
#'
#' @param recipe recipes::recipe object to which the step is added.
#'
#' @param ... Recipes selectors resolving to exactly one outcome.
#'
#' @param role Recipes role assigned to the step.
#'
#' @param trained Internal training indicator. Use prep() to train the step.
#'
#' @param columns Outcome column resolved during prep().
#'
#' @param sby_adasyn_ratio Nonnegative expansion relative to the original minority: G = floor(n_min * ratio). Zero disables generation. This is a reparameterization, not the beta used in the original paper. Inactive when d_th prevents ADASYN or beta is supplied instead.
#'
#' @param sby_adasyn_k Positive integer number of ADASYN neighbors. Difficulty uses min(k, n - 1); interpolation uses min(k, n_min - 1). These searches run only when G > 0.
#'
#' @param sby_config_max_threads Positive integer per-call thread ceiling, or -1L for detection. Capped by physical cores, CPU affinity, container quota and the hard OpenMP thread limit. Both Intel OpenMP and oneMKL receive this resolved ceiling. Local controls are restored on success or error. BLAS runs outside OpenMP regions; small BLAS calls may use fewer threads. HPC interfaces use this parameter directly. No AVX-512 requirement.
#'
#' @param sby_seed Integer seed from 0 to .Machine$integer.max; default sample.int(10e7, 1). ADASYN uses a scoped Mersenne-Twister/Inversion/Rejection RNG and restores RNGkind and .Random.seed, including after errors. Identical input, seed, parameters and numerical environment reproduce the result. Evaluating the default sample.int consumes the caller RNG; supply a seed explicitly to avoid this. NearMiss itself is deterministic and does not draw random numbers.
#'
#' @param sby_audit FALSE retains the always-present sbyaudit and sby attributes. TRUE also attaches the detailed audit. HPC interfaces always return a tibble; classic tabular interfaces return a list with sby_balanced_data when audited; matrix interfaces return lists. For recipes, the step bake method returns data with attributes; final recipes::bake(recipe) may drop them. The last audit remains in `prepared_recipe$steps[[i]]$audit_log$last`.
#'
#' @param sby_restore_types Must be TRUE. Binary synthetics use threshold >= 0.5; integer domains use round() with ties to even; every synthetic predictor is clamped to its original minimum and maximum. Original rows are copied without rounding.
#'
#' @param skip TRUE skips resampling new data in bake(); prep() still applies the step to training data.
#'
#' @param id Recipes step identifier.
#'
#' @param sby_adasyn_beta Optional number from 0 to 1: G = floor((n_maj - n_min) * beta), using the original paper parameterization. Cannot be supplied together with an explicitly supplied sby_adasyn_ratio. NULL uses ratio.
#'
#' @param sby_adasyn_d_th Threshold from 0 to 1. ADASYN executes only when n_min / n_maj < d_th. Default 1; zero disables generation. Class roles are determined on the original data and remain fixed.
#'
#' @param sby_adasyn_zero_difficulty Policy when all difficulty values are zero: "error" (default) stops because the paper normalization is undefined; "uniform" explicitly requests the documented uniform-quota extension. Consult the audit for the resolved policy and whether fallback was used.
#'
#' @inherit sby_adanear_hpc details references
#'
#' @return A recipe with the step added. The step bake method returns data. See Details for auditing after recipes::bake(recipe).
#' @export
sby_step_adasyn_hpc <- function(
  recipe,
  ...,
  role = NA,
  trained = FALSE,
  columns = NULL,
  sby_adasyn_ratio = 0.2,
  sby_adasyn_k = 3L,
  sby_config_max_threads = -1L,
  sby_seed = sample.int(10e7, 1),
  sby_audit = FALSE,
  sby_restore_types = TRUE,
  skip = TRUE,
  id = recipes::rand_id("adasyn_hpc"),
  sby_adasyn_beta = NULL,
  sby_adasyn_d_th = 1,
  sby_adasyn_zero_difficulty = c("error", "uniform")
) {
  terms <- rlang::enquos(...)
  parameters <- mget(setdiff(names(formals(sys.function())), c("recipe", "...", "role", "trained", "columns", "skip", "id")), envir=environment())
  sby_step_add(recipe, terms, parameters, "sby_step_adasyn_hpc", "adasyn", role, trained, columns, skip, id, !missing(sby_adasyn_ratio))
}
