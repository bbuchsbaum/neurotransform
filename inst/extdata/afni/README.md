# AFNI Registration Fixtures

Files produced by running `align_sub-1001_to_mni.sh` (AFNI 3dQwarp with `-allineate`):

- `sub-1001_T1w_in_mni.nii.gz`: T1 warped into MNI space (forward).
- `sub-1001_T1w_in_mni_WARP.nii.gz`: forward nonlinear warp (5D NIfTI: x,y,z,t=1,vec=3).
- `sub-1001_T1w_in_mni_WARP_inv.nii.gz`: inverse nonlinear warp (5D NIfTI).
- `sub-1001_T1w_to_mni.aff12.1D`: 3dAllineate affine (DICOM-to-DICOM).
- `sub-1001_T1w_to_mni.inv.aff12.1D`: inverse affine.
- `sub-1001_T1w_to_mni.params.1D`: shift/angle parameters.
- `sub-1001_T1w_in_mni_Allin.aff12.1D`: affine from 3dAllineate stage.
- `anatQQ.FT.aff12.1D`: additional affine fixture.
- Inputs: `sub-1001_T1w.nii.gz`, `ss_sub-1001_T1w.nii.gz`, `mni.nii`.
- Script: `align_sub-1001_to_mni.sh` (cleans old outputs, runs 3dQwarp, inverts warp).

Notes:
- Warps are stored with a singleton time dim; last dim holds 3 displacement components.
- Affines are row-wise DICOM-to-DICOM matrices in .1D format.
