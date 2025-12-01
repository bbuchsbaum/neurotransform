test_that("afni_aff12_to_ras flips Z axis correctly", {
  mat_rai <- diag(4)
  mat_rai[3, 4] <- 5
  mat_ras <- afni_aff12_to_ras(mat_rai)
  expect_equal(mat_ras[3, 4], -5)
  expect_equal(mat_ras[1:2, 1:2], diag(2))
})

test_that("afni_load_affine_morphism reads and converts aff12", {
  path <- system.file("extdata/afni/anatQQ.FT.aff12.1D", package = "neurotransform")
  m <- afni_load_affine_morphism("src", "tgt", path, direction = "source_to_target")
  expect_s4_class(m, "Affine3DMorphism")
  expect_equal(source_of(m), "src")
  expect_equal(target_of(m), "tgt")
  expect_true(is.matrix(m@matrix))
})

test_that("afni_warp_transform_coords flips RAI correctly with sample field", {
  # Build a tiny warp in RAI space: +1 in z (inferior->superior flip should invert sign)
  dimf <- c(2L, 2L, 2L)
  disp <- array(0, dim = c(3, dimf))
  disp[3, , , ] <- 1
  # Save as temporary NIfTI using RNifti if available
  skip_if_not_installed("RNifti")
  tmp <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(RNifti::asNifti(aperm(disp, c(2,3,4,1))), tmp)

  m <- Warp3DMorphism("src", "tgt", tmp, warp_type = "afni")
  coords <- matrix(c(0, 0, 0), ncol = 3)
  out <- afni_warp_transform_coords(m, coords)
  # In RAS, expect -1 in z due to RAI flip before/after
  expect_equal(out[1, 3], -1, tolerance = 1e-8)
})
