# Test Transform Fixtures

This directory contains test data for validating AFNI, FSL, and ANTs transform handling in neurotransform.

## Test Validation Status

| Tool | Real Data | Reference Output | Correlation Test |
|------|-----------|------------------|------------------|
| **AFNI** | Yes (sub-1001 T1w) | Yes (3dQwarp output) | r > 0.8 |
| **ANTs** | Yes (chris_t1) | Yes (antsApplyTransforms) | r > 0.8 |
| **FSL** | Partial* | Partial* | r > 0.8* |

*FSL requires running `fsl/register_to_mni.sh` to generate real test data.

## AFNI (`afni/`)

**Real neuroimaging data with reference validation:**

- `ss_sub-1001_T1w.nii.gz` — Skull-stripped T1w source image
- `mni.nii` — MNI template (target space)
- `sub-1001_T1w_to_mni.aff12.1D` — Linear affine transform
- `sub-1001_T1w_in_mni_WARP.nii.gz` — Forward nonlinear warp (3dQwarp)
- `sub-1001_T1w_in_mni_WARP_inv.nii.gz` — Inverse nonlinear warp
- `sub-1001_T1w_in_mni.nii.gz` — **Reference output** for validation

**Generation:** Run `align_sub-1001_to_mni.sh` (requires AFNI installed)

## ANTs (`chris/ants/`)

**Real neuroimaging data with reference validation:**

- `chris_t1.nii.gz` — T1w source image
- `chris_to_mni_Composite.h5` — Full composite transform (warp + affine)
- `chris_to_mni_InverseComposite.h5` — Inverse composite
- `reg_0GenericAffine.mat` — ITK affine component
- `reg_1Warp.nii.gz` — Forward warp (NIfTI format)
- `reg_1InverseWarp.nii.gz` — Inverse warp
- `chris_in_mni.nii.gz` — **Reference output** for validation

**Generation:** Run `tools/warp_chris_to_mni.R` (requires ANTs/rpyANTs)

## FSL (`fsl/`)

### Placeholder Data (basic tests)

- `S01_lin_6dof.mat`, `S01_lin_12dof.mat` — Small FLIRT matrices
- `S01_warp.nii.gz`, `S01_coef.nii.gz` — Tiny 3×3×3×3 placeholder warps

These exercise code paths but are NOT validated against FSL output.

### Real Data (requires generation)

Run `register_to_mni.sh` to generate:

- `highres2standard.mat` — Real FLIRT affine (12 DOF)
- `highres2standard_warp.nii.gz` — Real FNIRT warp field
- `highres_in_mni.nii.gz` — Reference output from FNIRT
- `highres_in_mni_applywarp.nii.gz` — **Reference output** for validation
- `standard2highres_warp.nii.gz` — Inverse warp

**Requirements:** FSL installed with `FSLDIR` set, AFNI test data generated first.

## Downloading Additional Test Data

### FSL Course Data (~1.3 GB)

```bash
curl -L -O -C - https://fsl.fmrib.ox.ac.uk/fslcourse/downloads/registration.tar.gz
tar -xzf registration.tar.gz
```

### AFNI Demo Data

```bash
# From AFNI tutorials
curl -L -O https://afni.nimh.nih.gov/pub/dist/edu/data/CD.expanded/AFNI_data6.tar.gz
```

## Notes

- The `afni/`, `chris/`, and `demo1/` folders are developer fixtures and are excluded from `R CMD build` via `.Rbuildignore` to keep the package tarball small.
- Tests use `skip_if_not(file.exists(...))` so CI/CRAN-style checks still run without these large files.
- If you want full end-to-end validation locally, keep these directories present (or re-generate them via the scripts referenced above).
