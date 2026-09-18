#!/usr/bin/env bash
# FSL FNIRT test data generation script
# Creates real FNIRT warps and reference outputs for validation testing
#
# Prerequisites:
#   - FSL installed and FSLDIR set
#   - MNI152 template available (comes with FSL)
#   - Source T1w image (we use the same subject as AFNI tests for consistency)
#
# Usage:
#   cd inst/extdata/fsl
#   ./register_to_mni.sh
#
# Or, without a local FSL install, with the pinned FSL 5.0.9 image used for the
# other native fixtures (run from the package root; about 7 minutes emulated):
#   docker run --rm --network none --platform linux/amd64 \
#     -e FSLDIR=/usr/share/fsl/5.0 -e FSLOUTPUTTYPE=NIFTI_GZ \
#     -v "$PWD/inst/extdata:/work" --entrypoint /bin/bash \
#     brainlife/fsl@sha256:fbd262c385e9de22aa58bf7b6311cbd5cd96c7b4eaff151f879191e869bf224e \
#     -c 'export PATH=$FSLDIR/bin:$PATH; cd /work/fsl && bash register_to_mni.sh'
#
# Outputs (large; gitignored and excluded from the package build). Resampled
# references are written as float so tests can compare them tightly:
#   - highres2standard.mat                FLIRT affine (12 DOF)
#   - highres_in_mni_flirt.nii.gz         FLIRT-only resampled source (no resampling blur)
#   - highres2standard_warp.nii.gz        FNIRT relative field (includes the affine)
#   - highres2standard_warp_coef.nii.gz   FNIRT spline coefficients
#   - highres2standard_warp_noaff.nii.gz  Nonlinear field without the affine
#   - highres2standard_warp_abs.nii.gz    The FNIRT mapping as absolute coordinates
#   - highres2standard_jac.nii.gz         FSL Jacobian determinant (with affine)
#   - highres_in_mni.nii.gz               FNIRT's own warped output
#   - highres_in_mni_applywarp.nii.gz     applywarp output (float)
#   - standard2highres_warp.nii.gz        invwarp inverse field

set -euo pipefail

# Check FSL is available
if [ -z "${FSLDIR:-}" ]; then
    echo "Error: FSLDIR not set. Please source FSL setup script first."
    exit 1
fi

# Configuration
MNI_TEMPLATE="${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz"

# Use the same source as AFNI tests for cross-tool comparison
SOURCE="../afni/ss_sub-1001_T1w.nii.gz"

# Output names
FLIRT_MAT="highres2standard.mat"
FNIRT_WARP="highres2standard_warp.nii.gz"
FNIRT_COEF="highres2standard_warp_coef.nii.gz"
WARPED_OUT="highres_in_mni.nii.gz"
INV_WARP="standard2highres_warp.nii.gz"
FNIRT_WARP_NOAFF="highres2standard_warp_noaff.nii.gz"
FNIRT_WARP_ABS="highres2standard_warp_abs.nii.gz"
FNIRT_JAC="highres2standard_jac.nii.gz"
FNIRT_LOG="highres2standard_fnirt.log"
FLIRT_OUT="highres_in_mni_flirt.nii.gz"

# Verify inputs exist
if [ ! -f "$MNI_TEMPLATE" ]; then
    echo "Error: MNI template not found at $MNI_TEMPLATE"
    echo "Trying 1mm template..."
    MNI_TEMPLATE="${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz"
    if [ ! -f "$MNI_TEMPLATE" ]; then
        echo "Error: No MNI template found. Check FSL installation."
        exit 1
    fi
fi

if [ ! -f "$SOURCE" ]; then
    echo "Error: Source image not found at $SOURCE"
    echo "Please ensure AFNI test data exists first (run AFNI alignment script)"
    exit 1
fi

echo "=== FSL Registration Pipeline ==="
echo "Source: $SOURCE"
echo "Template: $MNI_TEMPLATE"
echo ""

# Clean previous outputs
rm -f "$FLIRT_MAT" "$FNIRT_WARP" "$FNIRT_COEF" "$WARPED_OUT" "$INV_WARP" "$FNIRT_WARP_NOAFF" \
      "$FNIRT_WARP_ABS" "$FNIRT_JAC" "$FNIRT_LOG" "$FLIRT_OUT" "${WARPED_OUT%.nii.gz}_applywarp.nii.gz"

# Step 1: Linear registration with FLIRT (12 DOF affine)
echo "Step 1: Running FLIRT (linear registration)..."
flirt \
    -in "$SOURCE" \
    -ref "$MNI_TEMPLATE" \
    -out highres_flirt.nii.gz \
    -omat "$FLIRT_MAT" \
    -dof 12 \
    -interp trilinear

echo "  FLIRT complete: $FLIRT_MAT"

# Step 2: Nonlinear registration with FNIRT
echo "Step 2: Running FNIRT (nonlinear registration)..."
fnirt \
    --in="$SOURCE" \
    --ref="$MNI_TEMPLATE" \
    --aff="$FLIRT_MAT" \
    --cout="$FNIRT_COEF" \
    --fout="$FNIRT_WARP" \
    --iout="$WARPED_OUT" \
    --logout="$FNIRT_LOG" \
    --config=T1_2_MNI152_2mm

echo "  FNIRT complete: $FNIRT_WARP"

# Step 3: Generate inverse warp
echo "Step 3: Generating inverse warp..."
invwarp \
    --ref="$SOURCE" \
    --warp="$FNIRT_WARP" \
    --out="$INV_WARP"

echo "  Inverse warp complete: $INV_WARP"

# Step 4: Verify by applying warp to generate reference output
echo "Step 4: Generating reference output with applywarp..."
applywarp \
    --in="$SOURCE" \
    --ref="$MNI_TEMPLATE" \
    --warp="$FNIRT_WARP" \
    --out="${WARPED_OUT%.nii.gz}_applywarp.nii.gz" \
    --interp=trilinear \
    --datatype=float

# Step 5: Export the nonlinear part alone. The --fout field above already
# contains the FLIRT affine; this field lets tests check FLIRT + FNIRT ordering.
echo "Step 5: Exporting the nonlinear field without the affine..."
fnirtfileutils \
    --in="$FNIRT_COEF" \
    --ref="$MNI_TEMPLATE" \
    --out="$FNIRT_WARP_NOAFF"

# Step 6: The same mapping as absolute source coordinates.
echo "Step 6: Converting the field to absolute coordinates..."
convertwarp \
    --ref="$MNI_TEMPLATE" \
    --warp1="$FNIRT_WARP" \
    --rel \
    --absout \
    --out="$FNIRT_WARP_ABS"

# Step 7: FSL's own Jacobian determinant of the full mapping (with the affine).
echo "Step 7: Computing the FNIRT Jacobian determinant..."
fnirtfileutils \
    --in="$FNIRT_COEF" \
    --ref="$MNI_TEMPLATE" \
    --jac="$FNIRT_JAC" \
    --withaff

# Step 8: FLIRT-only reference. FLIRT blurs the input when resampling to a
# coarser grid unless -noresampblur is given, so its registration -out image is
# not a plain trilinear sample of the matrix; resample again without the blur.
echo "Step 8: Resampling with the FLIRT matrix alone..."
flirt \
    -in "$SOURCE" \
    -ref "$MNI_TEMPLATE" \
    -applyxfm -init "$FLIRT_MAT" \
    -interp trilinear \
    -noresampblur \
    -datatype float \
    -out "$FLIRT_OUT"

# Clean up intermediate files
rm -f highres_flirt.nii.gz

echo ""
echo "=== Registration Complete ==="
echo "Outputs:"
echo "  - $FLIRT_MAT (FLIRT affine matrix)"
echo "  - $FNIRT_WARP (FNIRT warp field)"
echo "  - $FNIRT_COEF (FNIRT coefficient field)"
echo "  - $WARPED_OUT (warped output from FNIRT)"
echo "  - ${WARPED_OUT%.nii.gz}_applywarp.nii.gz (reference from applywarp)"
echo "  - $INV_WARP (inverse warp)"
echo "  - $FNIRT_WARP_NOAFF (nonlinear field without the affine)"
echo "  - $FNIRT_WARP_ABS (absolute-coordinate field)"
echo "  - $FNIRT_JAC (FSL Jacobian determinant, including the affine)"
echo "  - $FLIRT_OUT (FLIRT-only resampled reference, float)"
echo ""
echo "These files can be used to test neurotransform FSL warp handling."
