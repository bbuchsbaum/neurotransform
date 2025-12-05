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

test_that("afni_warp_transform_coords applies LPS to RAS conversion", {
  # Build a tiny warp with displacement in X direction (LPS convention)
  # LPS: +X=Left, so +1 displacement in X should become -1 in RAS
  dimf <- c(2L, 2L, 2L)
  disp <- array(0, dim = c(3, dimf))
  disp[1, , , ] <- 1  # +1 displacement in X (LPS Left direction)

  # Write using neuroim2 - need 4D array (X, Y, Z, 3)
  tmp <- tempfile(fileext = ".nii.gz")
  arr4d <- aperm(disp, c(2, 3, 4, 1))  # Convert to (X, Y, Z, 3) for NIfTI
  space <- neuroim2::NeuroSpace(dim = dim(arr4d))
  vol <- neuroim2::NeuroVec(arr4d, space)
  neuroim2::write_vec(vol, tmp)

  m <- Warp3DMorphism("src", "tgt", tmp, warp_type = "afni")
  coords <- matrix(c(0, 0, 0), ncol = 3)
  out <- afni_warp_transform_coords(m, coords)
  # LPS->RAS: X is negated, so +1 becomes -1
  expect_equal(out[1, 1], -1, tolerance = 1e-8)
  # Z should be unchanged (same in LPS and RAS)
  expect_equal(out[1, 3], 0, tolerance = 1e-8)
})
