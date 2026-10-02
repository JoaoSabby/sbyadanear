# Technical review: sbyadanear 0.4.0

Review date: 2026-10-02. Scope: complete sby_adanear_hpc execution, all public
sampling/recipes interfaces, scientific definitions, original-row provenance,
precision, domains, parallel controls, build configuration and documentation.

## Outcome and validation boundary

The reviewed implementation uses one shared exact double-precision engine.
All 14 sby_adanear_hpc parameters have an identified execution or validation
path. Parameters controlling a disabled stage remain validated; their
conditional use is explicitly documented. The returned data copies all
original rare rows and retained original majority rows, followed by restored
rare synthetics. Class roles stay fixed from the input.

Local validation passed 23 scientific/regression cases with 231 expectations,
zero failures and zero warnings. It also passed after installing the R
namespace and exercising the exported interfaces and recipes S3 methods.
Roxygen/NAMESPACE/Rd synchronization and static R usage checks passed.
A before/after comparison preserved all numerical sampling outputs, lineage,
quotas, candidates and scores in 24 scenarios across models 1/2/3, K=1/7,
thread ceilings 1/2 and query sets crossing the 128-row tile boundary.
The only normalized comparison metadata was the intentionally corrected
resolved zero-difficulty policy, formerly stored as its default choice vector.

These local native tests used an excluded GCC/OpenMP/CBLAS diagnostic build
on R 4.6.0. It explicitly reports that Intel was NOT tested. It is neither a
supported Windows implementation nor a production fallback. The production
package requires Linux x86_64, icpx, oneMKL and libiomp5. Full production
compilation/linkage, observed Intel OpenMP teams and public-pipeline oneMKL
NThr remain UNVERIFIED until the Oracle Linux workflow/container is run.
No local diagnostic result establishes Intel runtime success.

## Parameter trace for sby_adanear_hpc

| Parameter | Execution path and condition |
|---|---|
| .data | Formula extraction; validate numeric predictors/outcome; preserve original rows and predictor storage. |
| formula | Resolve existing outcome/predictor columns; reject transformations/interactions; preserve selected input-column order. |
| sby_adasyn_k | Mixed-neighbor difficulty and rare-original interpolation K; capped separately; searches run only when G > 0. |
| sby_nearmiss_k | K nearest/farthest rare distances used by the selected NearMiss score; includes continuous synthetics; scoring only when majority rows are removed. |
| sby_config_max_threads | Integer ceiling or -1; resolve physical cores, affinity, current/ancestor CPU quota and hard OpenMP limit; pass effective value to native ThreadScope for both OpenMP and oneMKL; restore prior controls. |
| sby_seed | Validate integer; scoped Mersenne-Twister/Inversion/Rejection; R-thread neighbor choices and lambda; restore caller RNG even on error. Default sampling consumes one RNG draw before this scope. |
| sby_adasyn_ratio | G = floor(original_minority * ratio); inactive with beta or a preventing d_th. Zero disables synthetics. |
| sby_nearmiss_ratio | min(original_majority, floor(expanded_minority * ratio)); cap before integer cast; zero disables selection; a positive target rounding to zero errors. |
| sby_audit | Always keep sbyaudit/sby; TRUE adds scale, indices, lineage, candidates, scores, thread diagnostics and stage telemetry. Detailed RSS/lineage transfer is skipped otherwise. |
| nearmiss_model | Validate 1/2/3, default 3L; choose nearest-smallest, farthest-smallest or two-stage nearest-largest scoring. |
| sby_nearmiss_m | Independently cap and use M in NearMiss-3 candidate preselection; inactive for models 1/2. |
| sby_adasyn_beta | Optional original-paper G = floor((majority-minority)*beta), beta in [0,1]; reject simultaneously explicit ratio. |
| sby_adasyn_d_th | Generate only if original_minority/original_majority < d_th; zero disables ADASYN; equal class counts are preserved without synthesis at default threshold. |
| sby_adasyn_zero_difficulty | Resolve error/uniform; error for undefined normalization by default; uniform is an explicit extension, recorded only as used when all difficulties are zero. |

R/sby_adanear_hpc.R delegates through sby_dispatch and sby_run in R/sby_core.R
into the registered three-argument sby_pipeline_cpp routine. The effective
options list contains the native controls; results are assembled and audited
in R. No separate HPC geometry or hidden secondary thread budget remains.

## Scientific and provenance findings

- ADASYN difficulty excludes the exact self index, not every duplicate row.
  Normalize majority-neighbor counts (algebraically equivalent to normalizing
  r_i = hits_i/K because K is common). Neighbor parents are original rare
  rows. One lambda applies to every predictor of each synthetic.
- Largest-remainder quotas sum to G with deterministic original-index ties.
  This integer policy, explicit uniform fallback, initial z scaling and
  domain restoration are documented choices/extensions. Ratios exceeding
  balance imply beta > 1; they extend the paper's beta range, while the beta
  argument itself is limited to [0,1]. Literal equality to unspecified paper
  choices or other libraries is not claimed.
- NearMiss-1 and NearMiss-2 average Euclidean d rather than d squared.
  NearMiss-3 preselects per rare row and ranks candidates by the largest
  mean nearest-rare distance. Its second stage is identified as the
  imbalanced-learn operational definition, distinct from the original
  Mani/Zhang prefilter. Candidate shortage warns and never fills arbitrarily.
- Initial population mean/SD are computed once on all originals. The same
  z-space is used for difficulty, interpolation and NearMiss. Continuous
  synthetic z values participate before binary/integer restoration.
- The affine inverse uses original parents to avoid cancellation. Binary
  0/1 synthetics use >=0.5; whole-number domains use ties-to-even rounding;
  every synthetic predictor is bounded by original minima/maxima. Original
  rows are copied and never rounded or inverse-transformed.
- Audit positions distinguish output rows from original input indices.
  Optional lineage records parents, neighbor and lambda, plus continuous/z
  matrices before domain restoration. Caller data are not duplicated as
  a full dataset in the parameter summary.

## Bugs and inconsistencies corrected

| Finding | Correction/evidence |
|---|---|
| Squared-distance NearMiss scores could select different rows. | Shared mean of Euclidean d; dedicated d-versus-d-squared oracle fixture. |
| Retention could overflow before capping. | Double retention cap before integer conversion; tests include extreme finite ratios. |
| R integer count multiplication could overflow memory estimates. | Counts promoted to double before products. |
| Distinct close numeric outcomes could collapse through text conversion. | Direct match of outcome values; regression with 1 and 1+2*epsilon. |
| Recipes formulas failed for outcome names containing spaces. | Build formula from the outcome symbol, without parsing an unquoted name. |
| Squared distances could overflow/underflow despite a finite distance. | Extended-precision square root in those cases; reject truly unrepresentable distances on the R calling thread. |
| An unlimited mount-root quota missed restricted current/ancestor cgroups. | Inspect process membership and ancestors for cgroup v2 and conventional v1 CPU mounts; fixtures for fractional and inherited quotas. |
| Invalid external scaling/budget types produced obscure internal errors. | Explicit type/shape validation before numerical use. |
| Legacy workers could escape validation when another ceiling prevailed. | Validate them even when an explicit ceiling takes precedence. |
| Default zero-difficulty choice vector was stored instead of resolved policy. | Record the resolved error/uniform value. |
| Recipes audit attributes could be dropped by final generic processing. | Preserve last-call sbyaudit/audit in the trained step's audit_log; document replacement semantics. |
| Library detection could choose an empty legacy directory or miss symlinked oneMKL roots. | Locate the actual libmkl_core.so file; follow symlinks during container library discovery. |
| Aggregate MKL evidence could miss over-budget or serial public calls. | Label each probe/public call; require NThr evidence for budgets 1/2 and reject any over-budget DGEMM. Parser positive/negative fixtures passed. |

## Performance changes without changing scientific selections

- Cache each immutable matrix's row norms within the call and reuse caches
  for subsets. Preserve the original extended-precision accumulation order.
- Allocate one query tile buffer and reuse it for successive DGEMM calls.
- Maintain an exact K-neighbor heap per worker, retaining Euclidean distances
  and original-index ties. Scratch falls from O(threads*reference_rows) to
  O(threads*K); no approximation or squared-distance ranking is introduced.
- Sort only the retained NearMiss prefix; average candidate scores in
  independent OpenMP iterations.
- Reuse the R input matrix for validation/native transfer; avoid a second
  conversion. Skip domain inference without synthetics unless needed for
  audit, preserving empty-column storage in no-op returns.
- Collect detailed RSS telemetry only when requested. RNG/R assembly remain
  serial because using R APIs from OpenMP workers would be unsafe.
- Allocate worker scratch on the calling thread. Allocation errors and
  invalid-distance errors unwind through Rcpp and restore ThreadScope.
  User interrupts are checked between neighbor tiles and parent batches.

A three-run diagnostic timing comparison on a fixed 3500-row/12-predictor
workload had median 0.17 s before and 0.14 s after. This is a small local
GCC/CBLAS measurement, not a production Intel benchmark or a guaranteed
speedup. Reduced allocations and scratch bounds are the established changes.
Exact search still has quadratic time in the compared row counts. Full
lineage audit increases memory. Hardware capacity and user ceilings may
legitimately resolve to one thread; this is not an alternative serial backend.

## Package-wide cleanup and documentation

Obsolete Fortran/C++ engines, float/SGEMM buffers, optional MKL switching,
fork/TBB/approximate routes, uncalled helpers and their obsolete wrappers,
exports/tests/documentation were replaced by the shared implementation.
All 18 public exports are preserved. Native registration now declares the
three routines actually used; dynamic symbol lookup is disabled. DESCRIPTION
imports and namespace/S3 declarations match current code. Rcpp namespace
loading is retained intentionally for the native bridge.

Legacy engine/backend/HNSW arguments remain only as documented and validated
compatibility selectors; they are not hidden alternate algorithms. Removing
the public arguments would create another unnecessary interface break.
The unrelated session-wide future/database/locale configuration was removed
from the optional development profile. Generated local artifacts and prior
private audit files remain outside the package/Git publication scope.

Roxygen, generated help, README, NEWS, headers/comments and target-build
instructions now describe the same behavior in technical English. Shared
scientific details are inherited from the HPC topic to avoid conflicting
copies. The help documents no-op stages, index conventions, observed versus
configured threads, cumulative memory peaks, RNG default consumption,
external scaling limits and recipes audit behavior.

The target workflow uses Oracle Linux 9 inside a container with Intel compiler,
OpenMP and oneMKL, and R 4.5.3 by default. It verifies documentation, installed
scientific tests, actual linkage, OpenMP teams, per-call DGEMM evidence and
R CMD check. Action references follow the current official [checkout](https://github.com/actions/checkout) and [upload-artifact](https://github.com/actions/upload-artifact) usage, with checkout credential
persistence disabled. Build/test evidence is retained even when CI fails.
The workflow accompanies this change. Production verification remains pending until its Oracle Linux execution succeeds; workflow logs provide the authoritative target-runtime evidence.

## Remaining target verification

Run the documented container command or the repository workflow on Linux.
Success must include icpx compilation, Intel libraries actually resolved,
observed OpenMP teams and DGEMM NThr=2 for the public pipeline's two-thread
case. The local machine has no configured Linux/Intel execution environment.
Custom cgroup mount locations, very large CPU masks beyond CPU_SETSIZE and
cross-library numerical reproducibility remain portability limits rather
than locally verified configurations. No production performance claim is
made before this evidence exists.
