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

- Dense FSL relative and absolute fields now decode scaled-voxel coordinates
  using source and reference image geometry, including handedness flips and
  offsets. Evaluating these fields requires `source_affine` and `source_dim`;
  reference geometry defaults to the warp lattice or can be supplied explicitly.
  The normalized RAS displacement is shared by point transforms, resampling
  plans, Jacobians, and exports. Inversion swaps the stored geometries.
- Explicit `ants_h5` loading now routes affine-only H5 files through the existing
  ITK affine reader. Unsupported/ambiguous components still fail explicitly.
- Warp instances now own independent caches, preventing one instance from
  reusing another instance's decoded field after a file is rewritten.
- Added native FSL 5.0.9 fixtures for all four source/reference handedness pairs,
  with differing oblique grids, relative/absolute fields, image and coordinate
  comparisons, flattened resampling, exports, and Jacobian regression checks.
- Dense FSL warps now fail when constructed (`Warp3DMorphism()`,
  `read_transform()`) if `source_affine`/`source_dim` are missing or malformed,
  instead of at first use.
- `detect_transform_type()` classifies NIfTI vector fields from the header
  before the filename: FSL intents 2006 (dense) and 2007--2009 (coefficients),
  and the ITK 5D vector layout (intent 1007), are definitive, and a 4D
  `(X, Y, Z, 3)` field is FSL's layout. Real FSL fields were previously
  auto-detected as ANTs and silently misread. An inferred FSL field without
  source geometry now errors with a message naming the fix.
- `detect_fnirt_def_type()` no longer uses a displacement-magnitude threshold,
  which classified real FNIRT relative fields containing a FLIRT affine as
  absolute. Given `source_affine`/`source_dim` it picks the reading whose
  implied source coordinates fall inside the source image, then falls back to
  a volume-preserving Jacobian test; when neither is decisive it stops instead
  of guessing. `read_transform()` passes the source geometry and no longer
  defaults ambiguous fields to relative. `threshold_mm` is deprecated and
  ignored.
- `invert()` on a dense FSL warp reads only the forward field header to obtain
  the default reference geometry.
- FNIRT spline-coefficient files (`fnirt --cout`, `warp_type = "fsl_coef"`) are
  now decoded with FSL's conventions, established against native FSL 5.0.9:
  knots over FSL's x-flipped reference index, unnormalized cubic (or quadratic)
  B-spline weights, the `--aff` matrix stored in the sform, and
  `source_FSL = inv(A) ref_FSL + d`. The coefficients are decoded to the dense
  field `fnirtfileutils --withaff` would write and share the dense FSL path, so
  transforms, flattened resampling plans, Jacobians, and exports all apply. The
  previous evaluator treated the header as a world grid, normalized the
  weights, and ignored FSL coordinates; it has been removed together with
  `cpp_apply_bspline_coeff_field()`. Coefficient warps require source and
  reference geometry (the reference orientation and origin are not stored in
  the file), and files without FSL spline intents or with a reflecting `--aff`
  matrix are rejected.
- Added native FSL 5.0.9 coefficient fixtures (`inst/extdata/fsl_coef_oracle`)
  for every handedness pair with and without `--aff`, and a script to generate
  real FLIRT/FNIRT validation data locally (`inst/extdata/fsl/register_to_mni.sh`).
- `invert()` on FSL warps takes the inverse's format from the inverse file's
  header, so coefficient and dense inverses round-trip. Coefficient warps
  reject `def_type = "absolute"`.
- The dense FSL reader trusts a header intent of 2006 over a filename that
  mentions coefficients, and checks the header only on a cache miss.
- Single-slice FNIRT coefficient files no longer fail to decode.
- Custom loaders registered for `"fsl_coef"` must now return the raw
  coefficient list of `load_warp_fsl_coef()` (coefficients, knot spacing,
  spline order, reference dimensions and voxel size, and the FLIRT matrix).
