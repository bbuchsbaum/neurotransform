# Barycentric face-selection repair: validation

2026-09-26. Uncommitted changes over
`9e7550d45d6d6b06155b27717057f4137aaf666e`.
`source.json` binds the implementation, regression tests and oracle generators
to SHA-256 hashes. Changes are confined to neurotransform; downstream route
admission and engine pinning have not been changed.

## Behavior

The compiled weight kernel now selects the closest point on a nondegenerate
triangle by squared Euclidean distance, including edge/vertex projections.
Canonical vertex order and an exact-distance tie-break remove face-row and
winding dependence. Degeneracy is checked relative to triangle scale; the old
fixed denominator perturbation and weight-sum score are gone. Positive small
weights are retained instead of being dropped below 1e-10.

Spherical plans use this closest-point rule at a common radius and error if no
valid triangle provides support. They do not substitute nearest vertices for
missing support. Generic morphisms explicitly request interior-only orthogonal
projections, preserving their existing uncovered-query `NA` behavior.

## Evidence

- `probe-before.json`: the original installed engine fails the retained
  neuroatlas compiled admission probe (`+1/3`, `-1/3`; exit 1). Its source
  revision is not asserted; its DLL hash is recorded.
- `probe-after.json`: the rebuilt package passes the same unmodified probe
  (`+1/3`, `+1/3`; exit 0), with its own DLL hash and session information.
- `focused-tests.log`: 92 expectations pass across the geometry, Workbench,
  morphism and surface-plan test files. Geometry tests cover all octants,
  boundary/near-boundary queries, triangle winding, face and vertex permutations,
  scale changes, distance ranking, ties, degeneracy and missing support. An
  independent L1-ball projection oracle checks full weights for 212 queries.
- `inst/extdata/barycentric_oracle/oracle.json`: independent Workbench 2.2.1
  BARYCENTRIC impulse-response fixtures, with three synthetic cases and 66
  targets each. Maximum absolute weight error is `2.081e-7`, below the preset
  `2e-6` tolerance. Both original and reversed face rows pass.
- `template-comparison.json`: 64 seeded target vertices in each direction and
  hemisphere, six nonconstant metrics, using fsaverage 164k and fsLR 32k spheres
  registered in fsaverage space. All four comparisons pass the preset `5e-5`
  absolute tolerance. The receipt records source/target geometry hashes, query
  indices, Workbench commands, versions, compiled DLL hash and R sessions.

| Hemisphere / direction | Maximum absolute error |
| --- | ---: |
| L, 164k to 32k | 8.77549199596039e-6 |
| L, 32k to 164k | 3.07433680163394e-6 |
| R, 164k to 32k | 7.34706291827258e-6 |
| R, 32k to 164k | 1.91484686606902e-6 |

`R CMD build` and `R CMD check --no-manual` complete. The final package check
has zero errors and one compiler warning from R's `R_ext/Boolean.h`: the local
Homebrew clang does not recognize `-Wfixed-enum-extension`. Tests report
`FAIL 0 | WARN 0 | SKIP 69 | PASS 1309`; optional fixture/CRAN skips are listed
in `testthat.Rout`. Full logs are retained as `build.log`, `00check.log`,
`00install.out` and `testthat.Rout`. The source-tree suite also passed before
the final API-default adjustment; the final package check and focused tests
cover that adjustment.

## Limits and reproduction

This validates ordinary closest-point barycentric interpolation, not radial
intersection or ADAP_BARY_AREA. Template comparisons are sampled, not
full-density qualification. Masks, labels, missing metric values and area
correction remain outside this evidence. The kernel still scans all source
faces per query; no acceleration claim is made. The separate surface-sampler
implementation in `src/rcpp_sample.cpp` is unchanged by this repair.

The independent oracle generator and template-comparison script are in `tools/`.
Set `R_LIBS` to a freshly built package library for the template comparison.
Raw GIFTIs, template CSVs and attempt logs are retained under
`/tmp/neurotransform-bary-qa/`; the official Workbench disk image was unmounted
after the comparisons. The original neuroatlas findings and its failing receipt
were left intact.
