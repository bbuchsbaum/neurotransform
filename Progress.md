Progress toward first release
=============================

Status snapshot (tests passing as of last run):
- Core morphisms, affine/warp sampling (linear/cubic), Jacobians: covered by tests.
- Resampling plans: C-backed mapping/interpolation, caching keyed by morphism hash, grids, interpolation, and warp file mtimes. Tests cover plan parity, 4D volumes, and flattened paths.
- Fixtures: AFNI/FSL/ANTs affines and warps validated in tests.
- AFNI, FSL ingestion helpers: covered.
- Surface/barycentric sampling: unskipped and tested.
- Default warp I/O now routes through neuroim2 (RNifti kept only as optional fallback); neuroim2 affine extraction and tests added.

Pending/nice-to-have items:
- Additional cache invalidation cues (e.g., hash underlying files) if workflows involve frequent file changes.
- Benchmark harness exists (`tools/bench_resample.R`); could integrate into CI or docs.
- Documentation/vignette for plan usage and performance characteristics.

Next candidate tasks:
- Optional CI microbenchmarks to track regressions.
- Double-check plan caching when morphism files are regenerated in place.
