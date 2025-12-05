test_that("identity affine resampling preserves volume", {
  src <- array(seq_len(27), dim = c(3, 3, 3))
  grid <- grid_spec(c(3L, 3L, 3L), diag(4))
  morph <- Affine3DMorphism("src", "tgt", diag(4))

  out <- resample_volume(src, morphism = morph, target = grid,
                         method = "nearest", modulate = "none")
  expect_equal(out, src)
})

test_that("zero displacement warp preserves volume", {
  skip_if_not_installed("neuroim2")

  src <- array(runif(27), dim = c(3, 3, 3))
  grid <- grid_spec(c(3L, 3L, 3L), diag(4))

  # zero displacement field (components in 4th dimension)
  disp <- array(0, dim = c(3, 3, 3, 3))
  warp_path <- tempfile(fileext = ".nii.gz")
  space <- neuroim2::NeuroSpace(dim(disp), trans = diag(4))
  neuroim2::write_vec(neuroim2::DenseNeuroVec(disp, space), warp_path, format = "nifti")

  morph <- Warp3DMorphism("src", "tgt", warp_path = warp_path, warp_type = "ants")

  out <- resample_volume(src, morphism = morph, target = grid,
                         method = "nearest", modulate = "none")
  expect_equal(out, src)
})

test_that("pure translation affine shifts volume correctly", {
  # Source grid: 4x4x4, 1 mm spacing, origin (0,0,0)
  src <- array(0, dim = c(5, 5, 5))
  src[3, 3, 3] <- 1  # impulse at world (2,2,2)
  src_aff <- diag(4)
  src_grid <- grid_spec(dim(src), src_aff)

  # target grid: same size/affine
  tgt_grid <- src_grid

  # Morphism: translate -1 mm in each axis (target -> source) so target coords map one voxel back.
  T <- diag(4)
  T[1:3, 4] <- c(-1, -1, -1)
  morph <- Affine3DMorphism("src", "tgt", matrix = T)

  out <- resample_volume(src, morphism = morph, target = tgt_grid,
                         method = "nearest", modulate = "none")

  # Target voxel at (3,3,3) (1-based index 4) should sample source at (3,3,3)
  expect_equal(out[4, 4, 4], 1)
})

test_that("translation + scale affine works", {
  src <- array(0, dim = c(4, 4, 4))
  src[2, 2, 2] <- 1

  # Source affine: 1 mm voxels
  src_aff <- diag(4)
  src_grid <- grid_spec(dim(src), src_aff)

  # Target affine: 2 mm voxels, shifted origin
  tgt_aff <- diag(4)
  tgt_aff[1:3, 1:3] <- diag(c(2, 2, 2))
  tgt_aff[1:3, 4] <- c(1, 1, 1)
  tgt_grid <- grid_spec(dim(src), tgt_aff)

  # Morphism: map target->source; here, since target voxels are coarser and shifted,
  # use identity to let world mapping handle it.
  morph <- Affine3DMorphism("src", "tgt", matrix = diag(4))

  out <- resample_volume(src, morphism = morph, target = tgt_grid,
                         method = "nearest", modulate = "none")

  # Compute where the impulse lands: source impulse at world (1,1,1)
  # Target grid voxel (1,1,1) has world (1 + 2*(i-1), etc.). Expect impulse near target index (1,1,1).
  expect_equal(out[1, 1, 1], 1)
})
