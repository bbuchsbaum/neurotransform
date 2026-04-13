test_that("fsl spacing and scaling matrices round-trip", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 3; aff[3, 3] <- 4
  sp <- fsl_spacing_from_affine(aff)
  expect_equal(sp, c(2, 3, 4))
  dims <- c(5L, 6L, 7L)
  vox2fsl <- fsl_vox_to_fsl(aff, dim = dims)
  fsl2vox <- fsl_fsl_to_vox(aff, dim = dims)
  expect_equal(vox2fsl %*% fsl2vox, diag(4))
})

test_that("fsl world<->fsl transforms invert", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 3; aff[3, 3] <- 4
  dims <- c(5L, 6L, 7L)
  w2f <- fsl_world_to_fsl(aff, dim = dims)
  f2w <- fsl_fsl_to_world(aff, dim = dims)
  expect_equal(w2f %*% f2w, diag(4), tolerance = 1e-8)
})

test_that("fsl_flirt_to_internal_affine matches known conversion", {
  # Simple FLIRT matrix: translate +1 in x in FSL vox
  flirt <- diag(4); flirt[1, 4] <- 1
  src_aff <- diag(4)
  ref_aff <- diag(4)
  dims <- c(5L, 5L, 5L)
  internal <- fsl_flirt_to_internal_affine(
    flirt, src_aff, ref_aff,
    source_dim = dims, ref_dim = dims
  )
  # With handedness-aware FSL scaling, this pullback lands at +1 mm in x.
  expect_equal(internal[1, 4], 1, tolerance = 1e-8)
})

test_that("fsl_vox_to_fsl errors for right-handed affine without dims", {
  expect_error(
    fsl_vox_to_fsl(diag(4)),
    "dim must be supplied"
  )
})

test_that("detect_fnirt_def_type preserves RNG state", {
  skip_if_not_installed("neuroim2")
  dimf <- c(2, 2, 2, 3)
  arr <- array(0, dim = dimf)
  tmp <- tempfile(fileext = ".nii.gz")
  space <- neuroim2::NeuroSpace(dimf, trans = diag(4))
  neuroim2::write_vec(neuroim2::DenseNeuroVec(arr, space), tmp, format = "nifti")

  set.seed(42)
  before <- .Random.seed
  detect_fnirt_def_type(tmp, sample_n = 10, threshold_mm = 1)
  after <- .Random.seed

  expect_identical(after, before)
})

test_that("detect_fnirt_def_type handles tiny synthetic warp", {
  skip_if_not_installed("neuroim2")
  dimf <- c(2, 2, 2, 3)
  arr <- array(0, dim = dimf)
  tmp <- tempfile(fileext = ".nii.gz")
  space <- neuroim2::NeuroSpace(dimf, trans = diag(4))
  neuroim2::write_vec(neuroim2::DenseNeuroVec(arr, space), tmp, format = "nifti")
  def_type <- detect_fnirt_def_type(tmp, sample_n = 10, threshold_mm = 1)
  expect_equal(def_type, "relative")
})
