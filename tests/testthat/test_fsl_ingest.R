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

test_that("FLIRT conversion matches the scaled-voxel point equation", {
  source_affine <- diag(4)
  source_affine[1:3, 1:3] <- diag(c(2, 3, 4))
  source_affine[1:3, 4] <- c(10, -20, 5)
  target_affine <- diag(4)
  target_affine[1:3, 1:3] <- diag(c(-1.5, 2.5, 3.5))
  target_affine[1:3, 4] <- c(30, 40, -10)
  source_dim <- c(7L, 8L, 9L)
  target_dim <- c(10L, 11L, 12L)
  flirt <- diag(4)
  flirt[1:3, 1:3] <- matrix(c(
    1.0, 0.1, 0.0,
    0.0, 1.0, 0.2,
    0.0, 0.0, 0.9
  ), 3, byrow = TRUE)
  flirt[1:3, 4] <- c(5, -3, 2)

  # Independent statement of FSL's convention: scaled voxel coordinates, with
  # an x swap only for a right-handed voxel-to-world matrix.
  scaled_vox <- function(affine, dims) {
    spacing <- sqrt(colSums(affine[1:3, 1:3, drop = FALSE]^2))
    out <- diag(c(spacing, 1))
    if (det(affine[1:3, 1:3]) > 0) {
      swap <- diag(4)
      swap[1, 1] <- -1
      swap[1, 4] <- (dims[1] - 1) * spacing[1]
      out <- swap %*% out
    }
    out
  }
  src_vox_to_fsl <- scaled_vox(source_affine, source_dim)
  tgt_vox_to_fsl <- scaled_vox(target_affine, target_dim)
  expected <- source_affine %*% solve(src_vox_to_fsl) %*%
    solve(flirt) %*% tgt_vox_to_fsl %*% solve(target_affine)

  actual <- fsl_flirt_to_internal_affine(
    flirt, source_affine, target_affine,
    source_dim = source_dim, ref_dim = target_dim
  )
  target_points <- rbind(c(30, 40, -10), c(15, 52, 4))
  apply_affine <- function(mat, points) {
    (cbind(points, 1) %*% t(mat))[, 1:3, drop = FALSE]
  }

  expect_equal(actual, expected, tolerance = 1e-10)
  expect_equal(
    apply_affine(actual, target_points),
    apply_affine(expected, target_points),
    tolerance = 1e-10
  )
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
