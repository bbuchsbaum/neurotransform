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
