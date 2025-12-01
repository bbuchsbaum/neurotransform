test_that("resampling plan matches resample_volume output", {
  warp_path <- system.file("extdata/ants/sample_ANTs_1Warp.nii.gz", package = "neurotransform")
  morph <- Warp3DMorphism("src", "tgt", warp_path, warp_type = "ants")

  src_grid <- grid_spec(dims = c(3L, 3L, 3L), affine = diag(4), domain = "src")
  tgt_grid <- grid_spec(dims = c(3L, 3L, 3L), affine = diag(4), domain = "tgt")

  vol <- array(runif(27), dim = c(3, 3, 3))

  plan <- make_resampling_plan(morph, source_grid = src_grid, target_grid = tgt_grid,
                               interpolation = "linear", reuse_count = 2L)

  via_plan <- apply_resampling_plan(plan, vol, outside = NA_real_)
  direct <- resample_volume(vol, morphism = morph, target = tgt_grid,
                            method = "linear", modulate = "none")

  expect_equal(via_plan, direct, tolerance = 1e-8)
})

test_that("resampling plan caches and handles 4D volumes", {
  warp_path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  morph <- Warp3DMorphism("src", "tgt", warp_path, warp_type = "fsl")

  src_grid <- grid_spec(dims = c(3L, 3L, 3L), affine = diag(4), domain = "src")
  tgt_grid <- grid_spec(dims = c(3L, 3L, 3L), affine = diag(4), domain = "tgt")

  vol4d <- array(runif(3*3*3*2), dim = c(3, 3, 3, 2))

  plan1 <- make_resampling_plan(morph, src_grid, tgt_grid, interpolation = "linear", reuse_count = 3L)
  plan2 <- make_resampling_plan(morph, src_grid, tgt_grid, interpolation = "linear", reuse_count = 3L)

  expect_identical(plan1$key, plan2$key)
  expect_true(!is.null(plan1$handle))

  out <- apply_resampling_plan(plan1, vol4d, outside = 0)
  expect_equal(dim(out), c(tgt_grid@dims, 2))
})
