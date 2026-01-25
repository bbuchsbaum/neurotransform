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

For comprehensive validation against FSL reference output, run the registration script:

```bash
cd inst/extdata/fsl
chmod +x register_to_mni.sh
./register_to_mni.sh
```

**Prerequisites:**
- FSL installed with `FSLDIR` set
- AFNI test data generated first (the script uses `../afni/ss_sub-1001_T1w.nii.gz`)

**Generated files:**
- `highres2standard.mat` - Real FLIRT affine matrix (12 DOF)
- `highres2standard_warp.nii.gz` - Real FNIRT warp field (displacement)
- `highres2standard_warp_coef.nii.gz` - FNIRT coefficient field
- `highres_in_mni.nii.gz` - FNIRT warped output (reference)
- `highres_in_mni_applywarp.nii.gz` - applywarp output (reference for validation)
- `standard2highres_warp.nii.gz` - Inverse warp field

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
| `test_fsl_fnirt_resample.R` | Real FNIRT data | FSL applywarp reference output (r > 0.8) |

## Notes

- FNIRT warps can be "relative" (displacement) or "absolute" (coordinate) - use `detect_fnirt_def_type()` to auto-detect
- FSL uses "scaled voxel" coordinates for FLIRT matrices - use `fsl_flirt_to_internal_affine()` to convert
- Large test data files should NOT be committed to git - use `.gitignore` or local cache
