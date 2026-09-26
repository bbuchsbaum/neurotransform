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
