# FSL Warp Transform Tests
#
# FSL FNIRT produces warp fields that can be either:
# - "relative" (displacement vectors): coord + disp = warped_coord
# - "absolute" (coordinate fields): the warp directly stores target coordinates
#
# FSL uses a special "FSL space" coordinate system that differs from world coords.

# The test data in extdata/fsl is placeholder/mock data, not real warps.

test_that("FSL warp morphism can be created", {
  warp_path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  # Test that the morphism can be created
  morph <- Warp3DMorphism("native", "standard", warp_path = warp_path, warp_type = "fsl")
  expect_s4_class(morph, "Warp3DMorphism")
  expect_equal(morph@warp_type, "fsl")
})

test_that("detect_fnirt_def_type works on synthetic relative warp",
{
  skip_if_not_installed("neuroim2")

  # Create a synthetic relative displacement field (small displacements)
  dimf <- c(10, 10, 10, 3)
  arr <- array(rnorm(prod(dimf), mean = 0, sd = 2), dim = dimf)  # Small displacements
  tmp <- tempfile(fileext = ".nii.gz")
  space <- neuroim2::NeuroSpace(dimf, trans = diag(4))
  neuroim2::write_vec(neuroim2::DenseNeuroVec(arr, space), tmp, format = "nifti")

  def_type <- detect_fnirt_def_type(tmp, sample_n = 100, threshold_mm = 50)
  expect_equal(def_type, "relative")
  unlink(tmp)
})

test_that("detect_fnirt_def_type works on synthetic absolute warp", {
  skip_if_not_installed("neuroim2")

  # Create a synthetic absolute coordinate field (large values = world coords)
  dimf <- c(10, 10, 10, 3)
  # Create coordinate values that look like world coordinates
  arr <- array(0, dim = dimf)
  for (i in 1:dimf[1]) {
    for (j in 1:dimf[2]) {
      for (k in 1:dimf[3]) {
        arr[i, j, k, 1] <- i * 2  # X coord scaled by voxel size
        arr[i, j, k, 2] <- j * 2  # Y coord
        arr[i, j, k, 3] <- k * 2  # Z coord
      }
    }
  }

  tmp <- tempfile(fileext = ".nii.gz")
  space <- neuroim2::NeuroSpace(dimf, trans = diag(c(2, 2, 2, 1)))  # 2mm voxels
  neuroim2::write_vec(neuroim2::DenseNeuroVec(arr, space), tmp, format = "nifti")

  def_type <- detect_fnirt_def_type(tmp, sample_n = 100, threshold_mm = 50)
  expect_equal(def_type, "absolute")
  unlink(tmp)
})

test_that("FSL synthetic warp transform produces finite results", {
  skip_if_not_installed("RNifti")

  # Create a simple identity-like warp (small displacements)
  # With identity transform, warp coords = voxel coords
  # So we need to test at integer voxel coordinates
  dimf <- c(20, 20, 20, 3)
  arr <- array(0, dim = dimf)  # Zero displacement = identity

  tmp <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(arr, tmp)

  morph <- Warp3DMorphism("src", "tgt", warp_path = tmp, warp_type = "fsl",
                          def_type = "relative")

  # With identity transform, world coords = voxel coords
  # Valid voxel range is 0 to 19 (20 voxels)
  test_coords <- matrix(c(
    5, 5, 5,     # Well inside
    10, 10, 10,  # Center
    15, 15, 15   # Still inside
  ), ncol = 3, byrow = TRUE)

  warped <- transform(morph, test_coords)

  # Should be finite
  expect_true(all(is.finite(warped)))

  # For zero-displacement warp, output should equal input
  expect_equal(warped, test_coords, tolerance = 0.1)

  unlink(tmp)
})

test_that("FSL warp with small displacement transforms correctly", {
  skip_if_not_installed("RNifti")

  # Create warp with known small displacement (+1mm in X direction)
  dimf <- c(20, 20, 20, 3)
  arr <- array(0, dim = dimf)
  arr[, , , 1] <- 1  # +1mm displacement in X

  tmp <- tempfile(fileext = ".nii.gz")

  nii <- RNifti::asNifti(arr)
  xform_mat <- diag(c(2, 2, 2, 1))
  xform_mat[1:3, 4] <- c(0, 0, 0)
  RNifti::sform(nii) <- xform_mat
  RNifti::writeNifti(nii, tmp)

  morph <- Warp3DMorphism("src", "tgt", warp_path = tmp, warp_type = "fsl",
                          def_type = "relative")

  # Test at a point within the warp volume (10mm from origin)
  test_coords <- matrix(c(10, 10, 10), ncol = 3)
  warped <- transform(morph, test_coords)

  # Should have +1mm in X (10 + 1 = 11)
  expect_equal(warped[1, 1], 11, tolerance = 0.1)
  expect_equal(warped[1, 2], 10, tolerance = 0.1)
  expect_equal(warped[1, 3], 10, tolerance = 0.1)

  unlink(tmp)
})

test_that("fsl_flirt_to_internal_affine handles non-identity transforms", {
  # Test with a rotation + translation
  flirt <- matrix(c(
    0.9962, -0.0872, 0, 5,
    0.0872,  0.9962, 0, -3,
    0,       0,      1, 2,
    0,       0,      0, 1
  ), nrow = 4, byrow = TRUE)

  src_aff <- diag(c(2, 2, 2, 1))  # 2mm voxels
  ref_aff <- diag(c(1, 1, 1, 1))  # 1mm voxels

  internal <- fsl_flirt_to_internal_affine(
    flirt, src_aff, ref_aff,
    source_dim = c(5L, 5L, 5L), ref_dim = c(5L, 5L, 5L)
  )

  # Result should be a valid 4x4 affine
  expect_equal(dim(internal), c(4, 4))
  expect_equal(internal[4, ], c(0, 0, 0, 1))

  # Should be invertible
  expect_true(abs(det(internal)) > 0.01)
})
