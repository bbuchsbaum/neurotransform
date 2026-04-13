# neurotransform news

## Unreleased

- Fixed FSL handedness handling for right-handed affines. `fsl_vox_to_fsl()`
  and the higher-level FSL affine conversion stack now require image
  dimensions when a handedness swap is needed instead of silently computing the
  wrong matrix.
- Added `source_dim` / `target_dim` (and `ref_dim` where applicable) plumbing
  through the exported FSL affine IO helpers:
  `read_linear_transform()`, `write_linear_transform()`,
  `read_linear_transform_array()`, and `write_linear_transform_array()`.
- Corrected `fsl_flirt_to_internal_affine()` for handedness-aware conversion.
  For right-handed identity affines plus a FLIRT x-translation, the resulting
  pullback translation now matches the handedness-correct math rather than the
  previous silently wrong sign.
- `cpp_triplets_to_dgC()` now errors on out-of-bounds indices and aggregates
  duplicate triplets deterministically.
- `SurfToSurfMorphism` barycentric queries outside all faces now return `NA`
  rows instead of `(0,0,0)`.
- `build_affine_matrix(anchor = "centre"/"center")` now errors until an
  explicit anchor point is supplied.
- `detect_fnirt_def_type()` and `warp_from_field()` no longer mutate the
  caller RNG state.
