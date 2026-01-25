test_that("AFNI nonlinear warp resample overlaps reference (forward)", {
  src_path <- system.file("extdata/afni/ss_sub-1001_T1w.nii.gz", package = "neurotransform")
  tgt_path <- system.file("extdata/afni/mni.nii", package = "neurotransform")
  # Use FORWARD warp for pullback resampling (AFNI semantics: at each MNI voxel,

  # displacement points to corresponding source location)
  warp_path <- system.file("extdata/afni/sub-1001_T1w_in_mni_WARP.nii.gz", package = "neurotransform")
  ref_path <- system.file("extdata/afni/sub-1001_T1w_in_mni.nii.gz", package = "neurotransform")

  skip_if_not(file.exists(src_path))
  skip_if_not(file.exists(tgt_path))
  skip_if_not(file.exists(warp_path))
  skip_if_not(file.exists(ref_path))

  src <- suppressWarnings(read_image(src_path))
  tgt <- suppressWarnings(read_image(tgt_path))
  ref <- suppressWarnings(read_image(ref_path))

  morph <- Warp3DMorphism("src", "tgt", warp_path = warp_path, warp_type = "afni")

  out <- suppressWarnings(resample_to(src, target = tgt, transform = morph, method = "linear"))
  out_arr <- as.array(out); if (length(dim(out_arr)) == 4) out_arr <- out_arr[, , , 1, drop = TRUE]
  ref_arr <- as.array(ref); if (length(dim(ref_arr)) == 4) ref_arr <- ref_arr[, , , 1, drop = TRUE]

  mask <- is.finite(out_arr) & is.finite(ref_arr) & (ref_arr != 0)
  expect_gt(mean(mask), 0.1)

  r <- suppressWarnings(cor(as.vector(out_arr[mask]), as.vector(ref_arr[mask])))
  expect_gt(r, 0.8)
})

test_that("AFNI nonlinear warp transforms coordinates within expected range", {
  # Test that AFNI warps transform coordinates to reasonable locations
  # and that forward/inverse warps exist and can be loaded

  warp_path <- system.file("extdata/afni/sub-1001_T1w_in_mni_WARP.nii.gz", package = "neurotransform")
  inv_warp_path <- system.file("extdata/afni/sub-1001_T1w_in_mni_WARP_inv.nii.gz", package = "neurotransform")

  skip_if_not(file.exists(warp_path))
  skip_if_not(file.exists(inv_warp_path))

  # Create morphisms
  fwd_morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "afni")
  inv_morph <- Warp3DMorphism("mni", "native", warp_path = inv_warp_path, warp_type = "afni")

  # Test coords in MNI space (where both warps are defined)
  test_coords <- matrix(c(
     0,  0,  0,   # MNI origin
    10, 20, 30,   # Arbitrary point
   -20, 40, 50    # Another point
  ), ncol = 3, byrow = TRUE)

  # Both warps should transform coordinates
  warped_fwd <- suppressWarnings(transform(fwd_morph, test_coords))
  warped_inv <- suppressWarnings(transform(inv_morph, test_coords))

  # Warped coordinates should be finite
  expect_true(all(is.finite(warped_fwd)))
  expect_true(all(is.finite(warped_inv)))

  # Warped coordinates should be within reasonable neuroimaging bounds (±200mm)
  expect_true(all(abs(warped_fwd) < 200))
  expect_true(all(abs(warped_inv) < 200))

  # Displacements should be non-trivial (warps are not identity)
  disp_fwd <- warped_fwd - test_coords
  disp_inv <- warped_inv - test_coords
  expect_gt(max(abs(disp_fwd)), 1)  # At least 1mm displacement somewhere
  expect_gt(max(abs(disp_inv)), 1)
})
