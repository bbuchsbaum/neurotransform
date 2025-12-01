#!/usr/bin/env bash
set -euo pipefail

# Optional helper to pull larger, real-world transform datasets.
# Nothing here is run automatically; invoke manually when you need
# more than the tiny fixtures in inst/extdata.

DEST="${1:-./_testdata}"
mkdir -p "$DEST"

echo "Download destination: $DEST"

fetch_fsl_registration() {
  local url="https://fsl.fmrib.ox.ac.uk/fslcourse/downloads/registration.tar.gz"
  echo "Fetching FSL registration practical (~1.3 GB)..."
  curl -L -C - -o "$DEST/registration.tar.gz" "$url"
  echo "Unpacking..."
  tar -xzf "$DEST/registration.tar.gz" -C "$DEST"
  echo "FSL course data unpacked to $DEST/fsl_course_data/registration"
}

fetch_afni_affines() {
  local base="https://afni.nimh.nih.gov/pub/dist/edu/data/CD.expanded/AFNI_data6/FT_analysis/Qwarp/anat_warped"
  mkdir -p "$DEST/afni"
  echo "Fetching AFNI anat affine..."
  curl -L -o "$DEST/afni/anatQQ.FT.aff12.1D" "$base/anatQQ.FT.aff12.1D"
  echo "AFNI affine saved to $DEST/afni"
}

case "${2:-}" in
  fsl)
    fetch_fsl_registration
    ;;
  afni)
    fetch_afni_affines
    ;;
  all|"")
    fetch_fsl_registration
    fetch_afni_affines
    ;;
esac

echo "Done. Remember to keep large downloads out of git/CRAN."
