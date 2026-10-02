#' ADASYN and NearMiss resampling
#'
#' @description
#' Apply ADASYN and NearMiss through the common exact Intel oneAPI engine.
#'
#' @param .data Data frame or tibble containing the outcome and plain numeric predictors.
#'
#' @param formula Formula outcome ~ predictors. Select existing columns only; transformations and interactions are rejected.
#'
#' @param sby_adasyn_k Positive integer number of ADASYN neighbors. Difficulty uses min(k, n - 1); interpolation uses min(k, n_min - 1). These searches run only when G > 0.
#'
#' @param sby_nearmiss_k Positive integer number K of Euclidean distances d averaged for NearMiss, capped at the expanded rare-class size. Used only when majority retention actually reduces the data.
#'
#' @param sby_config_max_threads Positive integer per-call thread ceiling, or -1L for detection. Capped by physical cores, CPU affinity, container quota and the hard OpenMP thread limit. Both Intel OpenMP and oneMKL receive this resolved ceiling. Local controls are restored on success or error. BLAS runs outside OpenMP regions; small BLAS calls may use fewer threads. HPC interfaces use this parameter directly. No AVX-512 requirement.
#'
#' @param sby_seed Integer seed from 0 to .Machine$integer.max; default sample.int(10e7, 1). ADASYN uses a scoped Mersenne-Twister/Inversion/Rejection RNG and restores RNGkind and .Random.seed, including after errors. Identical input, seed, parameters and numerical environment reproduce the result. Evaluating the default sample.int consumes the caller RNG; supply a seed explicitly to avoid this. NearMiss itself is deterministic and does not draw random numbers.
#'
#' @param sby_adasyn_ratio Nonnegative expansion relative to the original minority: G = floor(n_min * ratio). Zero disables generation. This is a reparameterization, not the beta used in the original paper. Inactive when d_th prevents ADASYN or beta is supplied instead.
#'
#' @param sby_nearmiss_ratio Nonnegative majority retention relative to the expanded minority: min(n_maj, floor((n_min + G) * ratio)). Zero disables NearMiss. A positive ratio rounding to zero produces an informative error. If the target retains the entire majority, neighbor scoring is skipped.
#'
#' @param sby_audit FALSE retains the always-present sbyaudit and sby attributes. TRUE also attaches the detailed audit. HPC interfaces always return a tibble; classic tabular interfaces return a list with sby_balanced_data when audited; matrix interfaces return lists. For recipes, the step bake method returns data with attributes; final recipes::bake(recipe) may drop them. The last audit remains in `prepared_recipe$steps[[i]]$audit_log$last`.
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
#' @details
#' The combined pipeline computes one initial population z-scale on all
#' original rows; constant columns have scale 1. ADASYN difficulty is the
#' proportion of majority neighbors among the K mixed original neighbors,
#' excluding only the exact self index. It normalizes these difficulties and
#' allocates G synthetics using largest remainders: floor(G * w_i), then one
#' unit per largest fractional remainder, breaking ties by original row index.
#' This integer policy explicitly resolves a choice left open by the paper.
#' The quotas sum to exactly G. Zero difficulty requires the explicit uniform
#' extension or produces an error. Each synthetic chooses one ORIGINAL rare
#' neighbor and one uniform lambda in [0,1), shared by all columns.
#' Synthetic rows are never reused as parents.
#'
#' NearMiss averages Euclidean d, not squared d, in the same initial z-space.
#' In the combined pipeline, it uses rare originals plus continuous synthetics
#' BEFORE domain restoration. NearMiss-only routes do not execute ADASYN;
#' ADASYN-only routes do not execute NearMiss. NearMiss-3 uses the two-stage
#' operational definition documented by imbalanced-learn: the Mani and Zhang
#' prefilter followed by ranking candidates by their largest mean nearest-rare
#' distance. This latter selection must be distinguished from the prefilter.
#'
#' The affine inverse is evaluated stably using original parents (algebraically
#' equivalent to z * sigma + mu), avoiding cancellation when the global center
#' is large relative to rare values. Only synthetics undergo domain restoration:
#' observed binary 0/1 columns use threshold >= 0.5; observed integer columns,
#' including whole-valued doubles, use round() with ties to even; all columns
#' are clamped to their original minima and maxima. Arbitrary categorical codes
#' and decimal grids are not inferred. Domain restoration is an additional
#' postprocessing policy, not the continuous ADASYN rule in the paper.
#'
#' Retained originals are copied in input order, followed by restored synthetics.
#' Tabular routes return only formula-selected predictors and the outcome,
#' in their input-column order. Factor levels and ordering are preserved.
#' On equal class counts, the first factor level or first observed label is
#' the preserved class. If G = 0, difficulty is not calculated (NA in audit).
#' If retention includes the entire majority, scoring is skipped (NA).
#' audit$adasyn$executed and audit$nearmiss$executed record these conditions.
#'
#' The sbyaudit attribute is always present. adasyn and nearmiss contain
#' 1-based OUTPUT positions of synthetics and retained majority rows.
#' original_minority_indices and original_majority_indices refer to INPUT rows.
#' The element named after the public function contains evaluated parameters,
#' input dimensions and effective controls. Input records are not copied into
#' this attribute. sby$synthetic_rows is retained for compatibility (0L when
#' there are no synthetics). With sby_audit = TRUE, audit includes initial_scale,
#' domains, original_indices, parameters, counts, nearmiss, adasyn, synthetics
#' (parent, neighbor, lambda), synthetic_continuous, synthetic_standardized,
#' threads and a per-stage telemetry tibble.
#'
#' Telemetry records elapsed time, current RSS, cumulative PROCESS peak RSS,
#' physical/logical cores, quota/affinity, configured ceilings and observed
#' OpenMP team sizes. Unavailable metrics are NA. mkl_configured_threads is
#' configuration, whereas mkl_observed_threads is NA: use MKL_VERBOSE=1 to
#' observe DGEMM NThr. Buffer estimates and output-object sizes are not RSS.
#' Detailed telemetry and extra lineage matrices are collected only when
#' sby_audit = TRUE. Both audited and unaudited calls produce identical values.
#'
#' Scaling, exact-neighbor searches, interpolation and NearMiss score averaging
#' use Intel OpenMP. DGEMM uses threaded oneMKL outside OpenMP teams, avoiding
#' nested teams. RNG and R assembly remain serial. Query tiles are bounded;
#' row norms are cached per call, distance buffers reused and each worker
#' retains only K neighbors in an exact heap. NearMiss sorts only the retained
#' prefix. Double precision and cancellation checks are preserved; extreme
#' squared distances are square-rooted in extended precision when required.
#' Reproducibility requires the same numerical environment; cross-version or
#' cross-platform floating-point identity is not promised.
#'
#' The recipes step bake method returns data with these attributes, but the
#' final column selection in recipes::bake(recipe) may drop them. Consult
#' `prepared_recipe$steps[[i]]$audit_log$last` for the last sbyaudit and audit.
#' Each call replaces this entry; it is not a complete execution history.
#' Installation requires Linux x86_64, icpx, oneMKL LP64 and libiomp5. There is
#' no alternative-BLAS or serial native backend, and no AVX-512 requirement.
#'
#' @references
#' He, H., Bai, Y., Garcia, E. A., and Li, S. (2008). ADASYN: Adaptive
#' synthetic sampling approach for imbalanced learning. IJCNN, 1322-1328.
#' doi:10.1109/IJCNN.2008.4633969.
#'
#' Mani, I., and Zhang, I. (2003). kNN approach to unbalanced data distributions:
#' a case study involving information extraction. ICML Workshop on Learning
#' from Imbalanced Data Sets.
#'
#' Two-stage NearMiss operational definition:
#' https://imbalanced-learn.org/stable/under_sampling.html#near-miss
#'
#' @return A balanced tibble with sbyaudit and sby attributes; sby_audit = TRUE additionally attaches audit. See Details for indices and telemetry.
#' @export
sby_adanear_hpc <- function(
  .data,
  formula,
  sby_adasyn_k = 3,
  sby_nearmiss_k = 7,
  sby_config_max_threads = -1,
  sby_seed = sample.int(10e7, 1),
  sby_adasyn_ratio = 0.2,
  sby_nearmiss_ratio = 1,
  sby_audit = FALSE,
  nearmiss_model = 3L,
  sby_nearmiss_m = 3L,
  sby_adasyn_beta = NULL,
  sby_adasyn_d_th = 1,
  sby_adasyn_zero_difficulty = c("error", "uniform")
) {
  parameters <- mget(names(formals(sys.function())), envir=environment())
  sby_dispatch("adanear", parameters, "sby_adanear_hpc", !missing(sby_adasyn_ratio))
}
