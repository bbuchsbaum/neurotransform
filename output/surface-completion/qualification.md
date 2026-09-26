# Surface engine qualification

Implementation revision: `9ab40ffd9b3b9fd46ac2e805826cbc8fde60b429`.
[Source binding](source-binding.json) verifies every installed R/C++ source,
DESCRIPTION and NAMESPACE against that revision and records installed R
serialization and DLL hashes. This is a native geometry implementation with
known strict Workbench parity failures, not a fully admitted production route.

| Phase | Result | Evidence |
|---|---|---|
| S1 shared geometry | Pass | `s1-before.log`, `s1-final.log`, compiled entry-point and sampler regressions |
| S2 spherical admission | Pass for specified admitted domain | `s2-accepted.log`, degree-two/fold/topology/underflow rejection fixtures |
| S3 indexing and reuse | Pass | `s3-tests.log`, `s3-subnormal-tests.log`, `s8-thread-parity.log`, [engine benchmark](s8-benchmark/receipt.json) |
| S4 ordinary full-template parity | Fail: three of four directions | [Final ordinary comparison](s8-ordinary/receipt.json), [independent failure diagnosis](s4-full/diagnostic-receipt.json) |
| S5 coverage, missingness, labels | Pass scoped policy fixtures | `s5-s6-oracle-tests.log`, independent ROI/label fixtures and deterministic tie contracts |
| S6 adjoint/reverse | Pass scoped algebra and geometry fixtures | `s5-s6-testthat.Rout`, rectangular masked dense-matrix and rebuilt-reverse tests |
| S7 adaptive area | Implemented, experimental; full upsampling parity fails | [Mathematical contract](adaptive-contract.md), `s7-area-measure-tests.log`, [final adaptive comparison](s8-adaptive/receipt.json) |
| S8 package/consumer | Local package and analytic consumer probe pass; production pin held | `s8-00check.log`, `s8-testthat.Rout`, [consumer probe](s8-consumer-probe.json); published-revision CI is recorded by GitHub Actions |

## Full-target comparison

All runs use pinned fsaverage 164k and registered fsLR 32k inputs, both
hemispheres. Tolerances were fixed before evaluation: `1e-12` for controlled
analytic double calculations, `2e-6` for small Workbench basis weights, and
`5e-5` maximum absolute error for the bounded template fields.

| Route | Ordinary maximum error | Adaptive, all source vertices | Adaptive, medial-wall source ROI |
|---|---:|---:|---:|
| L 164k to 32k | 4.40729e-5 | 2.44982e-5 | 2.44982e-5 |
| L 32k to 164k | 1.64334e-4 | 1.64512e-4 | 8.05538e-5 |
| R 164k to 32k | 9.74984e-5 | 2.58036e-5 | 2.58036e-5 |
| R 32k to 164k | 8.09991e-5 | 8.25291e-5 | 8.25291e-5 |

Ordinary comparisons use eleven bounded fields, including ramps, oscillations,
seeded values and localized impulses. Adaptive comparisons add a constant field
and both categorical policies, with pinned anatomical area metrics and source
ROIs. Complete basis responses are checked on synthetic meshes.

All template rows have valid weights and unit mass before ROI application.
Reversing source face order leaves outputs unchanged. Adaptive source-coverage
masks agree with Workbench, and both categorical policies have zero mismatches
on supported outputs in all eight cases. Unsupported label encodings are
reported separately: this harness explicitly requests native key 0, while
Workbench emits key 1. They are not interchangeable with available measurements.
The library preserves supplied label tables and specifies its own deterministic
tie policy; bitwise categorical parity for every possible tie is not claimed.

The ordinary failures involve 17 targets. An independent least-squares interior
projection plus edge/vertex search agrees with native weights within 1.58e-14
at every failing target and the worst target on the passing route. Every face
whose AABB can improve a nearest-vertex upper bound is considered. Impulse
responses show that Workbench selects an edge at these queries and removes a
small positive third weight. The independent face remains closer even after
float32 radius normalization. This diagnoses a comparator disagreement; it
does not convert the original strict parity failure into a pass. No small
weights were pruned to fit the observed failures.

## Performance and package checks

The isolated ordinary-engine benchmark uses all target vertices and eleven
columns at four OpenMP threads. Construction takes 0.424–1.39 seconds;
additional peak resident memory is 334–544 MiB. Both meet the declared limits
of 60 seconds and 1 GiB additional memory. Cached sampling takes
0.023–0.224 seconds after index creation, and agrees with plan application to
2.23e-16. R policy-aware plan application takes 0.488–4.88 seconds in these runs.
These are measured local workloads, not hardware-independent guarantees.
Worst-case indexed search remains linear in triangle count.

The full comparator process additionally retains parsed input/expected data,
multiple plans and categorical results; its peak memory sometimes exceeds
1 GiB. Those receipts remain available. The engine benchmark explicitly records
baseline and additional memory and excludes CSV parsing from the timed engine
work. The first sandboxed timing attempt failed to read macOS counters; the
subsequent measured attempts resolved that instrumentation error.

The version 0.2.0 package passes R CMD check with no errors, 1689 assertions and
69 optional-fixture skips. One pre-existing compiler warning comes from R's
`Boolean.h` requesting a warning group unknown to local Clang. Examples and
vignettes pass. The new vignette's rendered numerical output was inspected.
One-thread and four-thread indexed/exhaustive weights and sampler results are
bitwise identical in the retained concurrency probe.

Cross-platform checks run on macOS, Linux and Windows in the repository's
[GitHub Actions workflow](https://github.com/bbuchsbaum/neurotransform/actions/workflows/R-CMD-check.yaml).
Use the exact published commit's run for its result; local package checks do not
stand in for those jobs.

## Remaining admission boundary

The original neuroatlas analytic engine probe now returns +1/3 for both face
orders. Its old engine pin has not been advanced because the required full
method-specific parity gates remain failed. No production route was activated,
and no unrelated neuroatlas working-tree changes were modified.

The next scientific decision is the intended near-edge agreement contract:
retain the more accurate native closest-point solution with a prospectively
validated comparator-error policy, or develop an explicitly named Workbench
compatibility mode and qualify it separately. The present work does neither
silently. The adaptive method continues to require `experimental=TRUE`.

Spherical admission conservatively rejects uncertain predicates; it is not an
exact-arithmetic acceptance guarantee for every ill-conditioned mesh. Tiny
origin translations within the radius-spread tolerance are not distinguished.
Area inputs must be strictly positive and geometry-bound; zero-area vertices
are rejected. Registration accuracy, new surface families and every possible
label tie remain outside the demonstrated qualification.
