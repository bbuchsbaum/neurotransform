#!/usr/bin/env bash
set -euo pipefail

BASE="mni.nii"
SOURCE="ss_sub-1001_T1w.nii.gz"
OUT_ALIGNED="sub-1001_T1w_in_mni.nii.gz"
WARP="sub-1001_T1w_in_mni_WARP.nii.gz"
WARP_INV="sub-1001_T1w_in_mni_WARP_inv.nii.gz"

# Clear any previous outputs so reruns always overwrite old results.
rm -f "$OUT_ALIGNED" "$WARP" "$WARP_INV"

# Nonlinear registration to MNI using 3dQwarp (includes an initial affine).
3dQwarp \
  -base "$BASE" \
  -source "$SOURCE" \
  -prefix "$OUT_ALIGNED" \
  -allineate \
  -blur 0 3 \
  -useweight \
  -overwrite \
  -verb

# Invert the nonlinear warp for downstream use.
3dNwarpCat -prefix "$WARP_INV" -iwarp "$WARP"
