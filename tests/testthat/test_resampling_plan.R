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

test_that("flattened resampling plan matches direct resample for affine-warp-affine paths", {
  warp_path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")

  aff1 <- diag(4)
  aff1[1, 4] <- -0.1
  aff2 <- diag(4)
  aff2[2, 4] <- 0.2

  path <- compose(
    compose(
      Affine3DMorphism("src", "mid1", aff1),
      Warp3DMorphism("mid1", "mid2", warp_path, warp_type = "fsl")
    ),
    Affine3DMorphism("mid2", "tgt", aff2)
  )

  src_grid <- grid_spec(dims = c(4L, 4L, 4L), affine = diag(4), domain = "src")
  tgt_affine <- diag(4)
  tgt_affine[1:3, 4] <- c(0.5, 0.5, 0.5)
  tgt_grid <- grid_spec(dims = c(2L, 2L, 2L), affine = tgt_affine, domain = "tgt")

  vol <- array(seq_len(prod(src_grid@dims)), dim = src_grid@dims)

  plan <- make_resampling_plan(path, src_grid, tgt_grid,
                               interpolation = "linear", reuse_count = 3L)
  via_plan <- apply_resampling_plan(plan, vol, outside = -99)

  sampler <- volume_sampler(vol, affine = src_grid@affine, method = "linear",
                            outside = -99, domain = "src")
  direct <- array(resample(sampler, path, tgt_grid, modulate = "none"),
                  dim = tgt_grid@dims)

  expect_true(!is.null(plan$handle))
  expect_equal(via_plan, direct, tolerance = 1e-8)
})

test_that("reuse_count leaves resampling output invariant", {
  warp_path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  morph <- Warp3DMorphism("src", "tgt", warp_path, warp_type = "fsl")

  src_grid <- grid_spec(dims = c(4L, 4L, 4L), affine = diag(4), domain = "src")
  tgt_affine <- diag(4)
  tgt_affine[1:3, 4] <- c(0.5, 0.5, 0.5)
  tgt_grid <- grid_spec(dims = c(2L, 2L, 2L), affine = tgt_affine, domain = "tgt")

  vol <- array(seq_len(prod(src_grid@dims)), dim = src_grid@dims)

  plan_direct <- make_resampling_plan(morph, src_grid, tgt_grid,
                                      interpolation = "linear", reuse_count = 1L)
  plan_cached <- make_resampling_plan(morph, src_grid, tgt_grid,
                                      interpolation = "linear", reuse_count = 3L)

  expect_null(plan_direct$handle)
  expect_true(!is.null(plan_cached$handle))

  out_direct <- apply_resampling_plan(plan_direct, vol, outside = -99)
  out_cached <- apply_resampling_plan(plan_cached, vol, outside = -99)

  expect_equal(out_cached, out_direct, tolerance = 1e-8)
})

test_that("resampling plan modulation fallback matches direct resample for 4D data", {
  scale <- diag(4)
  scale[1:3, 1:3] <- diag(c(1.5, 0.5, 1.2))
  morph <- Affine3DMorphism("src", "tgt", scale)

  src_grid <- grid_spec(dims = c(6L, 6L, 6L), affine = diag(4), domain = "src")
  tgt_grid <- grid_spec(dims = c(2L, 2L, 2L), affine = diag(4), domain = "tgt")

  vol4d <- array(seq_len(prod(c(src_grid@dims, 2L))), dim = c(src_grid@dims, 2L))

  plan <- make_resampling_plan(morph, src_grid, tgt_grid,
                               interpolation = "linear", reuse_count = 3L)
  via_plan <- apply_resampling_plan(plan, vol4d, outside = -99,
                                    modulate = "sqrt_jacobian")

  sampler <- volume_sampler(vol4d, affine = src_grid@affine, method = "linear",
                            outside = -99, domain = "src")
  direct <- array(resample(sampler, morph, tgt_grid, modulate = "sqrt_jacobian"),
                  dim = c(tgt_grid@dims, 2L))

  expect_equal(via_plan, direct, tolerance = 1e-8)
})

test_that("resampling plan cache key changes when warp file mtime changes", {
  src_warp <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  tmp_warp <- tempfile(fileext = ".nii.gz")
  expect_true(file.copy(src_warp, tmp_warp, overwrite = TRUE))
  on.exit(unlink(tmp_warp), add = TRUE)

  morph <- Warp3DMorphism("src", "tgt", tmp_warp, warp_type = "fsl")
  src_grid <- grid_spec(dims = c(3L, 3L, 3L), affine = diag(4), domain = "src")
  tgt_grid <- grid_spec(dims = c(3L, 3L, 3L), affine = diag(4), domain = "tgt")

  plan1 <- make_resampling_plan(morph, src_grid, tgt_grid,
                                interpolation = "linear", reuse_count = 3L)
  old_mtime <- unclass(file.info(tmp_warp)$mtime)

  Sys.setFileTime(tmp_warp, as.POSIXct(old_mtime + 10, origin = "1970-01-01", tz = "UTC"))
  new_mtime <- unclass(file.info(tmp_warp)$mtime)

  plan2 <- make_resampling_plan(morph, src_grid, tgt_grid,
                                interpolation = "linear", reuse_count = 3L)

  expect_gt(new_mtime, old_mtime)
  expect_false(identical(plan1$key, plan2$key))
})
