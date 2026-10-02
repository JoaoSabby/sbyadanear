# sbyadanear

ADASYN and NearMiss-1, NearMiss-2 and NearMiss-3 for binary outcomes with
numeric predictors. Version 0.4.0 supports Linux x86_64 with Intel oneAPI:
one exact double-precision engine, Intel OpenMP and threaded oneMKL BLAS.
Installation fails without that toolchain. AVX-512 is not required.

## Installation

```sh
source /opt/intel/oneapi/setvars.sh
Rscript -e 'install.packages(c("Rcpp", "tibble", "recipes", "rlang", "generics", "testthat", "roxygen2"))'
R CMD INSTALL .
```

`MKLROOT` must identify the oneMKL installation and `icpx` must be on PATH.
Source setvars.sh before each R session unless Intel library paths are
registered with the system dynamic loader. Linkage uses LP64, `mkl_intel_thread`, `mkl_core` and `libiomp5`. The runtime
diagnostic checks the OpenMP library actually resolved. The canonical
Dockerfile in `docker/oraclelinux97-r453/` builds Oracle Linux 9 with R and
oneAPI. Its name is retained for existing paths; the image tracks Oracle
Linux 9, and R_VERSION controls the R source release (default 4.5.3).

```sh
docker build -f docker/oraclelinux97-r453/Dockerfile -t sbyadanear-oneapi .
docker run --rm --cpus=2 -e MKL_VERBOSE=1 -e MKL_DYNAMIC=FALSE \
  -e SBY_REQUIRE_INTEL=true -v "$PWD:/workspace/sbyadanear" \
  sbyadanear-oneapi bash tools/ci-check.sh
```

The workflow checks installation, linkage, scientific tests, synchronized
help, observed OpenMP teams and DGEMM `NThr` from MKL_VERBOSE. Passing local
algorithm checks does not establish successful Linux/Intel compilation or
observed oneMKL parallelism: these require the target workflow/container.

## Usage

```r
library(sbyadanear)
result <- sby_adanear_hpc(
  .data = data, formula = outcome ~ .,
  sby_adasyn_ratio = 0.5, sby_nearmiss_ratio = 1,
  sby_adasyn_k = 5L, sby_nearmiss_k = 3L,
  nearmiss_model = 3L, sby_nearmiss_m = 3L,
  sby_config_max_threads = 4L,
  sby_seed = 123L, sby_audit = TRUE
)
attr(result, "sbyaudit")
attr(result, "audit")$telemetry
attr(result, "audit")$synthetics
sby_hpc_cpu_report()
```

Use an explicit seed for repeatability. The default `sample.int(10e7, 1)`
consumes one caller RNG draw when evaluated. The sampling call fixes a
scoped RNG and restores the previous state even after errors. Repeatability
requires identical input, seed, parameters and numerical environment;
cross-platform or library-version floating-point identity is not promised.
NearMiss itself makes no random draws.

## Scientific contract and pipeline

1. Determine class roles from the originals. On equal counts, preserve the
   first factor level or first observed non-factor label. Numeric labels
   are matched numerically, without lossy conversion to text.
2. Keep originals and compute one population z-scale over all original
   rows: mean and standard deviation with denominator n. Constant columns
   use scale 1.
3. ADASYN estimates difficulty from mixed original neighbors, excluding
   only the exact self index. It normalizes the majority-neighbor fractions,
   chooses another ORIGINAL rare neighbor, and draws one lambda in [0,1)
   shared by every predictor. Synthetics are never reused as parents.
4. NearMiss selects original majority rows using original rare rows plus
   continuous synthetic rare rows in the SAME initial z-space. Scores
   average Euclidean **d**, not squared distance.
5. Evaluate the affine inverse stably from the original parents, avoiding
   cancellation. Restore only synthetic domains: observed binary 0/1
   columns use threshold `>= 0.5`; integer-valued columns, including doubles,
   use `round()` with ties to even. Clamp EVERY synthetic predictor to its
   original minimum and maximum. Arbitrary categorical codes and decimal
   grids are not inferred.
6. Return retained originals in input order, followed by restored
   synthetics. Original rare and retained majority rows are copied without
   inverse scaling or rounding. Formula-selected predictors and the outcome
   retain their input-column order; factor levels/ordering are preserved.

ADASYN uses `G = floor(n_min * sby_adasyn_ratio)`. Alternatively,
`sby_adasyn_beta` uses `G = floor((n_maj - n_min) * beta)`, the paper's
parameterization. Do not explicitly supply both. Ratios beyond balance imply beta > 1 and
extend the paper parameter range; the beta parameter itself remains [0,1].
ADASYN runs only if
`n_min / n_maj < sby_adasyn_d_th` (default 1); zero disables it. Neighbor
counts are capped at available rows. G = 0 skips difficulty calculations,
recorded as NA with `audit$adasyn$executed = FALSE`.

Largest-remainder integer quotas sum to exactly G, with original-index
breaks for ties. This resolves a choice not fully specified by the paper.
Zero difficulty stops by default because normalization is undefined;
`sby_adasyn_zero_difficulty = "uniform"` explicitly requests uniform quotas.
Domain restoration is also additional postprocessing of continuous ADASYN.
Different quotas, scales, ties or fallback policies can change results
relative to other implementations.

NearMiss retention is `min(n_maj, floor((n_min + G) * ratio))`. Zero disables
it. The cap precedes integer conversion to prevent overflow. A positive
ratio rounding to zero produces an informative error. A target covering the
entire majority skips scoring, recorded as NA with `executed = FALSE`.

## NearMiss variants

| Model | Score neighbors | Retention |
|---|---|---|
| 1L | K nearest rare rows | Smallest mean distance |
| 2L | K farthest rare rows | Smallest mean distance |
| 3L (default) | K nearest rare rows after preselection | Largest candidate mean distance |

NearMiss-3 preselects the union of the M nearest majority rows of EACH rare
row. `sby_nearmiss_m` defines M independently of K, and only affects model 3.
If there are too few candidates, retain all candidates, warn and record
shortage; never fill from non-candidates. All ties use increasing original
row index. Models outside 1, 2 and 3 produce an informative error. The
two-stage variant follows the imbalanced-learn operational definition:
the Mani/Zhang prefilter and subsequent candidate ranking are distinguished.

## Audit

`sbyaudit` is always present. `adasyn` and `nearmiss` are 1-based OUTPUT row
positions of synthetics and retained majority rows. The original indices
refer to INPUT rows. Package, function and evaluated/effective parameters
are recorded; the full input dataset is not duplicated in this attribute.
`sby$synthetic_rows` remains available (0L without synthetics).

With `sby_audit = TRUE`, `audit` also provides the initial scale, domains,
indices, counts, NearMiss candidates/scores, ADASYN difficulty/quotas,
parents/neighbors/lambda, continuous and standardized synthetic matrices,
thread diagnostics and per-stage `telemetry` as a tibble. Apply documented
domain restoration to the parents' interpolation to reconstruct final
synthetics. Complete audit requires additional memory.

Telemetry distinguishes elapsed time, current RSS, cumulative PROCESS peak
RSS, buffer estimates and output-object size. Peak RSS is not an exclusive
stage peak. It includes cores, quota/affinity, resolved ceilings and observed
OpenMP team sizes. Configured MKL threads do not prove actual DGEMM teams:
`mkl_observed_threads` is NA; use `MKL_VERBOSE=1` to observe `NThr`.
Unavailable metrics are NA. Detailed RSS telemetry is skipped when audit is
disabled; sampling values are unchanged.

## Parallelism and bounded memory

The HPC ceiling comes from **sby_config_max_threads**, capped by physical
cores, CPU affinity, container quota and the hard OpenMP limit. Scaling,
exact-neighbor searches, interpolation and score averaging use Intel
OpenMP. DGEMM uses threaded oneMKL outside those teams, preventing nested
teams. RNG and R assembly remain serial. The call restores local MKL and
OpenMP controls on success or error without changing environment variables.
Small BLAS calls may use fewer threads than configured.

Per-call row norms are cached; query distance tiles are reused. Each worker
keeps only K neighbors in an exact heap, reducing neighbor scratch from
O(threads * reference_rows) to O(threads * K). Tie ordering and distances are
preserved. NearMiss sorts only the retained prefix. Cancellation checks
remain; extended-precision square roots prevent loss of finite distances
when their squared values overflow or underflow double storage. No full
all-pairs distance matrix is constructed, although exact search still has
quadratic time complexity in the compared row counts.

## Compatibility in 0.4.0

All tabular, matrix, index-selector and recipes interfaces delegate to the
same engine. HPC always returns a tibble; audited classic tabular calls or
calls requesting additional z output return a list with sby_balanced_data.
Matrix/index interfaces return lists. Legacy engine/algorithm selectors
accept auto/native and auto/brute; other engines and non-Euclidean metrics
are rejected. Legacy parallel/RcppParallel selectors map to Intel OpenMP,
without fork/TBB. HNSW is inactive and accepts only its defaults.
`sby_restore_types` must be TRUE; external sby_type_info is rejected.

Classic sby_knn_workers is used when sby_config_max_threads = -1L; an
explicit ceiling takes precedence. Deprecated selectors are validated
compatibility controls, not independent computation backends. The default
NearMiss is now 3L; request 1L to keep the former variant. Uniform ADASYN
fallback now requires explicit selection.

Classic interfaces accept external scaling and already-standardized input.
Reconstructed originals cannot be guaranteed bit-identical to originals not
provided. Standardized matrix output is a presentation option applied AFTER
domain restoration. HPC always computes its scale from original input.

Recipes step bake returns data with audit attributes. Final
recipes::bake(recipe) may drop them; consult
`prepared_recipe$steps[[i]]$audit_log$last` (sbyaudit and audit). Each call
replaces this entry; it is not a full history. Default skip = TRUE prevents
resampling new data while prep still resamples training data.

## References

- He et al. (2008), [ADASYN DOI](https://doi.org/10.1109/IJCNN.2008.4633969).
- Mani and Zhang (2003), *kNN approach to unbalanced data distributions*,
  ICML Workshop. Original PDFs are retained in references/.
- [NearMiss operational definitions](https://imbalanced-learn.org/stable/under_sampling.html#near-miss).
- [Intel OpenMP libraries](https://www.intel.com/content/www/us/en/docs/dpcpp-cpp-compiler/developer-guide-reference/2025-2/use-the-openmp-libraries.html).
- [Intel oneMKL linkage](https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2024-0/selecting-libraries-to-link-with.html).
