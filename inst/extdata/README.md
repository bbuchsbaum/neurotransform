Test transform fixtures
======================

Small example files for exercising AFNI, FSL and ANTs transform parsers.
Heavyweight course/demo datasets are **not** vendored here; see the
download recipes below if you need the full originals.

AFNI
----
- `afni/anatQQ.FT.aff12.1D` — real affine from AFNI's `@SSwarper` output (downloaded from `https://afni.nimh.nih.gov/pub/dist/edu/data/CD.expanded/AFNI_data6/FT_analysis/Qwarp/anat_warped/`).

FSL
---
- `fsl/S01_lin_6dof.mat`, `fsl/S01_lin_12dof.mat` — small non-identity 4×4 FLIRT matrices (rotation + translation, slight scaling/shear).
- `fsl/S01_warp.nii.gz`, `fsl/S01_coef.nii.gz` — 3×3×3×3 float displacement fields with non-zero gradients (FNIRT-compatible shape).
- Full course data (large, ~1.3 GB) is available via:
  - `curl -L -O -C - https://fsl.fmrib.ox.ac.uk/fslcourse/downloads/registration.tar.gz`
  - Unpack to get `fsl_course_data/registration/` with the real `S01_*` transforms.

ANTs
----
- `ants/sample_ANTs_0GenericAffine.mat` — non-identity ITK affine (5° Z-rotation, small translation).
- `ants/sample_ANTs_1Warp.nii.gz`, `ants/sample_ANTs_1InverseWarp.nii.gz` — paired 3×3×3×3 displacement fields (forward/inverse, small gradients).

Notes and usage
---------------
- NIfTI warps are tiny but non-zero so resampling/Jacobian code paths are exercised without bloating the repo.
- If you need higher-resolution or subject-realistic data, fetch the upstream datasets above (or AFNI QC demos) locally and point tests/benchmarks at them without committing to the repo.
- For CRAN safety, keep large downloads out of the package tarball; use a local cache or `.gitignore` a download directory.  
