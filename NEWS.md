# neurotransform 0.1.0

- Fixed NIfTI warp export to write the ANTs 5D vector layout and intent code.
- Reject unsupported H5 components and malformed affine/displacement data
  instead of silently applying an incomplete transform.
- Fixed cubic warp interpolation mixing adjacent vector components.
- Corrected warp Jacobians to use the full grid orientation when converting
  voxel derivatives to physical RAS derivatives.
- Added small, independently generated SimpleITK fixtures for both H5 orders,
  centered affines, oblique vector fields, and image resampling.

- Corrected ANTs/ITK affine direction, ANTs NIfTI/H5 LPS vector conversion,
  H5 displacement parameter layout, and composite transform ordering.
- Corrected AFNI affine axis conversion and the default direction of
  `3dAllineate -1Dmatrix_save` matrices.

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

- Explicit `ants_h5` loading now routes affine-only H5 files through the existing
  ITK affine reader. Unsupported/ambiguous components still fail explicitly.
