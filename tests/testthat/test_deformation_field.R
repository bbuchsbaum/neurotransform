test_that("deformation_field builds array and jacobian map", {
  warp_path <- system.file("extdata/ants/sample_ANTs_1Warp.nii.gz", package = "neurotransform")
  inv_path <- system.file("extdata/ants/sample_ANTs_1InverseWarp.nii.gz", package = "neurotransform")
  morph <- Warp3DMorphism("src", "tgt", warp_path, warp_type = "ants", inverse_path = inv_path)

  grid <- grid_spec(dims = c(3L, 3L, 3L), affine = diag(4))
  field <- deformation_field(morph, grid, with_jacobian = TRUE)

  expect_s3_class(field, "DeformationField")
  expect_equal(dim(field), c(3L, 3L, 3L, 3L))
  expect_equal(attr(field, "grid")@dims, grid@dims)
  expect_true(!is.null(attr(field, "jacobian_det")))

  det_vol <- jacobian_det_field(morph, grid)
  expect_equal(dim(det_vol), grid@dims)
  expect_equal(det_vol, attr(field, "jacobian_det"), tolerance = 1e-6)
})
