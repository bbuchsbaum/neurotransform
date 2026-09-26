# Surface completion execution

S1: the centroid-ranking regression failed twice against the original installed
engine (`s1-before.log`) and passes after sharing projection code. The final
S1 checks (`s1-final.log`) include shape/index regressions and Workbench fixtures.
An intermediate build (`s1-after.log`) failed because two local build commands
overlapped and object files were removed during linking. This is retained as
a build-orchestration failure, not a numerical pass. Subsequent builds ran
sequentially.

S2: `s2-accepted.log` passes all focused tests after fixing an independently
identified subnormal-determinant false-acceptance case and extreme-radius
quality checks. The deterministic generic-ray count establishes radial degree;
area is only a sanity diagnostic. Radius spread is origin-based and tolerant,
so sufficiently small translations are not claimed to be distinguishable.
Uncertain signs or failure to find a generic ray are rejected.

S3: shared indexed/exhaustive search passes `s3-tests.log`, including large-
coordinate cancellation, deterministic ties and serialized sampler rebuilding.
Distance is reconstructed from original convex corner coordinates and clamped
to their AABB so floating-point pruning bounds remain conservative. Full-size
performance and S4-S8 qualification remain open.

S3 performance: the optimized build constructed full-template plans in under
1.4 seconds. The repeat with authorized macOS counters measured 493–621 MiB
peak process RSS, below the 1 GiB budget (`s4-measured/receipt.json`).

S4: all targets in both directions and hemispheres were compared on eleven
bounded fields. Row sums and source-face permutation invariance pass. Three
routes fail the prospectively fixed 5e-5 maximum-error gate (max 1.64334e-4).
The first attempt also failed memory instrumentation because sandboxed
`/usr/bin/time -l` could not read kern.clockrate; the measured repeat resolves
that instrumentation failure, not the numerical failures. Both receipts remain.

Independent local least-squares projection at all 17 threshold-failing targets
(and the worst target on the passing route) agrees with native weights within
1.58e-14. All faces with bounding-box distance below a nearest-vertex upper
bound were considered. Workbench impulse outputs use only an edge at all 18
queries, omitting a small positive third native weight. Repeating the independent
face/edge comparison after float32 radius scaling still finds the interior face
closer. These diagnosed differences do not establish Workbench parity: the
original gate remains failed. No small weights are silently removed.

S5/S6: explicit ROI, missing-value and label policies; frozen Euclidean adjoint;
geometry-verified reverse construction; and deprecated legacy inverse are
implemented. Analytic tests and ordinary Workbench mask/non-tied-label fixtures
pass. `s5-s6-testthat.Rout` records 1521 passes, zero failures and 69 skips for
unavailable optional fixtures. Full R CMD check has zero errors and one existing
R-header/Clang unknown-warning-option warning. Vignettes and examples pass.
The first documentation attempt using roxygen's source loader failed because
its temporary S4 package namespace was unavailable; generating with pkgload and
compile=FALSE succeeds, with both logs retained.

S7: the separately named adaptive method is implemented with explicit
`experimental=TRUE`, geometry-bound positive areas, units and provenance.
Analytic support/ROI/area identities and twelve small independent Workbench
basis/label/ROI cases pass. An independent review found that effective area
needed to follow per-column missing-data omission; the corrected diagnostics
and reproduced regressions pass. Exact source-vertex queries now retain exact
structural zeros, preventing tiny spurious adaptive/categorical contributors.

The first S7 fixture-generation attempt reused one ROI filename; the retained
attempt uses immutable per-case ROI files and hashes. An intermediate focused
test failed in its JSON-to-matrix fixture reader, which was corrected. Those
attempt logs remain and are not counted as numerical passes.

S7/S8 full templates: both adaptive downsampling hemispheres pass the original
5e-5 gate with and without published medial-wall masks. Upsampling still fails
(maximum 1.64512e-4). Supported label outputs match for both policies in all eight
cases; coverage masks agree. The first full label report counted a missing-key
encoding mismatch (native explicit 0 versus Workbench 1) as prediction error.
The final report separates that mapping from supported categorical predictions;
no numerical tolerance was relaxed. The effective-area identity has scaled
error at most 2.28e-14 across these workloads. Ordinary full-template maxima
remain unchanged after the exact-vertex fix, with three of four routes failing.

S8 local package gate: version 0.2.0, 1689 passing assertions, zero failures,
69 optional-fixture skips, zero errors and one existing R-header/Clang warning.
The new vignette renders and its results were inspected: ordinary +1/3 and
edge zero, retained weight 2/3 under omission, matching 2/3 adjoint products,
and matching effective-area integrals 6248.2. One/four-thread results are bitwise
identical. The original neuroatlas compiled admission probe now passes both
face orders at +1/3. Source/artifact hashes are bound to the exact implementation
revision in `source-binding.json`; package checks do not admit production routes.

S8 performance: the full comparison harness includes CSV parsing, two plans,
wide expected/result matrices and label applications. Its final absolute peak
RSS sometimes exceeds 1 GiB. A separate engine benchmark removes CSV parsing
from measurement and records baseline RSS, construction, application, reusable
sampler creation and cached sampling separately. All four ordinary routes take
0.424–1.39 seconds to construct and add 334–544 MiB peak RSS, passing the declared
60-second/1-GiB additional-memory budget. Plan and cached sampler values agree
within 2.23e-16. Full comparison-harness memory receipts are retained as well.

S8 downstream admission is held: neuroatlas's existing engine pin is unchanged,
and its unrelated working-tree changes were preserved. Full method-specific
parity gates remain failed; the repaired analytic probe alone is insufficient
to activate routes. The Workbench image mounted for this task was unmounted.
