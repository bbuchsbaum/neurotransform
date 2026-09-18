# FSL Test Data

This directory contains test fixtures for validating FSL FLIRT and FNIRT transform handling.

## Current Contents

### Placeholder Data (small fixtures for basic tests)
- `S01_lin_6dof.mat` - Small 6-DOF FLIRT matrix (rotation + translation)
- `S01_lin_12dof.mat` - Small 12-DOF FLIRT matrix (full affine)
- `S01_warp.nii.gz` - Tiny 3×3×3×3 FNIRT-compatible displacement field
- `S01_coef.nii.gz` - Tiny coefficient field

These files are **mock/placeholder data** that exercise code paths without adding significant size to the package. They are NOT validated against real FSL output.

## Generating Real Test Data

Real FLIRT/FNIRT fixtures are generated locally, not committed (they are
large; `.gitignore` and `.Rbuildignore` exclude them). With FSL installed:

```bash
cd inst/extdata/fsl
./register_to_mni.sh
```

Without FSL, the script header gives a Docker command using the same pinned
FSL 5.0.9 image as the other native fixtures (about 7 minutes under
emulation). The source image is `../afni/ss_sub-1001_T1w.nii.gz`; the
reference is FSL's `MNI152_T1_2mm_brain`. FSL's outputs are deterministic:
repeated runs produce byte-identical files.

**Generated files:** the FLIRT matrix and its FLIRT-only resampling (with
`-noresampblur`, since FLIRT otherwise blurs when downsampling), the FNIRT
relative field, spline coefficients, affine-free field, absolute-coordinate
field, FSL Jacobian determinant, FNIRT and `applywarp` outputs (float), and the
`invwarp` inverse field. The script header lists each file.

## Alternative: Download FSL Course Data

For larger, more diverse test data:

```bash
cd _testdata
curl -L -C - -O https://fsl.fmrib.ox.ac.uk/fslcourse/downloads/registration.tar.gz
tar -xzf registration.tar.gz
```

This provides ~1.3GB of registration examples with various transforms.

## Test Coverage

| Test File | Data Used | Validates Against |
|-----------|-----------|-------------------|
| `test_fsl_ingest.R` | Synthetic | Roundtrip correctness |
| `test_fsl_warp.R` | S01_* placeholders | Basic loading, synthetic transforms |
| `test_fsl_dense_oracle.R` | `../fsl_dense_oracle` (shipped) | Native `convertwarp`/`applywarp`, all handedness pairs |
| `test_fsl_fnirt_resample.R` | Real FNIRT data (local) | `applywarp`, `flirt -applyxfm`, `convertwarp`, `fnirtfileutils`, `invwarp` |

## Notes

- FNIRT warps can be "relative" (displacement) or "absolute" (coordinate) - use `detect_fnirt_def_type()` to auto-detect
- Dense FNIRT fields need the source image geometry (`source_affine`, `source_dim`); `fnirt --fout` fields already include the FLIRT affine
- FSL uses "scaled voxel" coordinates for FLIRT matrices - use `fsl_flirt_to_internal_affine()` to convert
- Large test data files should NOT be committed to git - use `.gitignore` or local cache
