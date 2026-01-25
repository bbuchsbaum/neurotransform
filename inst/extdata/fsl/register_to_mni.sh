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
# Outputs:
#   - highres2standard.mat          FLIRT affine (12 DOF)
#   - highres2standard_warp.nii.gz  FNIRT warp field
#   - highres_in_mni.nii.gz         Reference output (source warped to MNI)
#   - standard2highres_warp.nii.gz  Inverse warp field (optional)

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
rm -f "$FLIRT_MAT" "$FNIRT_WARP" "$FNIRT_COEF" "$WARPED_OUT" "$INV_WARP"

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
    --interp=trilinear

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
echo ""
echo "These files can be used to test neurotransform FSL warp handling."
