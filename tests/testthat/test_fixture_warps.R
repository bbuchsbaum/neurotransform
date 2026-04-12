test_that("fixture warps load and carry expected gradients", {
  path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  morph <- Warp3DMorphism("a", "b", path, warp_type = "fsl")
  warp <- load_warp_array(morph)

  expect_equal(warp$dim, c(3L, 3L, 3L))
  expect_true(!is.null(warp$vox_to_world))

  arr <- array(warp$array, dim = c(3L, warp$dim))
  center <- arr[, 2, 2, 2]
  expect_equal(center, c(0.5, -0.25, 0.1), tolerance = 1e-6)
})

test_that("warp transform plugs into resample pipeline", {
  path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  morph <- Warp3DMorphism("src", "tgt", path, warp_type = "fsl")

  vol <- array(0, dim = c(3, 3, 3))
  for (x in 0:2) for (y in 0:2) for (z in 0:2) {
    vol[x + 1, y + 1, z + 1] <- x + 10 * y + 100 * z
  }

  sampler <- volume_sampler(vol, affine = diag(4), method = "linear",
                            outside = -99)
  coords <- rbind(c(1, 1, 1), c(0.5, 0.5, 0.5))

  src_coords <- transform(morph, coords)
  expected <- cpp_sample_volume(vol, src_coords, diag(4), "linear", outside = -99)

  res <- resample(sampler, morphism = morph, coords = coords, modulate = "none")
  expect_equal(res, expected, tolerance = 1e-6)
})

test_that("warp jacobian uses fixture gradients", {
  path <- system.file("extdata/ants/sample_ANTs_1Warp.nii.gz", package = "neurotransform")
  morph <- Warp3DMorphism("src", "tgt", path, warp_type = "ants",
                          inverse_path = system.file("extdata/ants/sample_ANTs_1InverseWarp.nii.gz",
                                                     package = "neurotransform"))

  coords <- matrix(c(1, 1, 1), ncol = 3)
  J <- jacobian(morph, coords)
  Jmat <- J@values[1, , ]

  expect_equal(diag(Jmat), c(1.5, 0.75, 1.1), tolerance = 1e-2)
  dets <- jacobian_det(morph, coords)
  expect_equal(as.numeric(dets), 1.2375, tolerance = 1e-2)
})

test_that("warp_method selects cubic sampler and jacobian path", {
  warp_path <- system.file("extdata/ants/sample_ANTs_1Warp.nii.gz", package = "neurotransform")
  inv_path <- system.file("extdata/ants/sample_ANTs_1InverseWarp.nii.gz", package = "neurotransform")

  morph_lin <- Warp3DMorphism("src", "tgt", warp_path, warp_type = "ants", inverse_path = inv_path,
                              warp_method = "linear")
  morph_cub <- Warp3DMorphism("src", "tgt", warp_path, warp_type = "ants", inverse_path = inv_path,
                              warp_method = "cubic")

  coords <- matrix(c(1.2, 1.1, 0.7), ncol = 3)
  warp <- load_warp_array(morph_lin)

  lin_out <- transform(morph_lin, coords)
  cub_out <- transform(morph_cub, coords)

  expect_equal(lin_out,
               cpp_apply_warp_field(coords, warp$array, warp$dim, warp$world_to_vox),
               tolerance = 1e-8)
  expect_equal(cub_out,
               cpp_apply_warp_field_cubic(coords, warp$array, warp$dim, warp$world_to_vox),
               tolerance = 1e-8)

  lin_Jdet <- as.numeric(jacobian_det(morph_lin, coords))
  cub_Jdet <- as.numeric(jacobian_det(morph_cub, coords))

  expect_equal(lin_Jdet,
               cpp_warp_jacobian_det(coords, warp$array, warp$dim, warp$world_to_vox, warp$vox_to_world),
               tolerance = 1e-8)
  expect_equal(cub_Jdet,
               cpp_warp_jacobian_det_cubic(coords, warp$array, warp$dim, warp$world_to_vox, warp$vox_to_world),
               tolerance = 1e-8)
})
