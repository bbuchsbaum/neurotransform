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
  morph <- Warp3DMorphism("native", "standard", warp_path = warp_path, warp_type = "fsl",
                          source_affine = diag(4), source_dim = c(3L, 3L, 3L))
  expect_s4_class(morph, "Warp3DMorphism")
  expect_equal(morph@warp_type, "fsl")
})

test_that("FSL warp construction fails without usable image geometry", {
  warp_path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  expect_error(Warp3DMorphism("native", "standard", warp_path, warp_type = "fsl"),
               "require source_affine and source_dim")
  expect_error(Warp3DMorphism("native", "standard", warp_path, warp_type = "fsl",
                              source_affine = diag(4)),
               "Supply both source_affine and source_dim")
  expect_error(Warp3DMorphism("native", "standard", warp_path, warp_type = "fsl",
                              source_affine = diag(4), source_dim = c(3, 3, 3, 1)),
               "dim\\(image\\)\\[1:3\\]")
  expect_error(Warp3DMorphism("native", "standard", warp_path, warp_type = "fsl",
                              source_affine = diag(4), source_dim = c(3, 3, 3),
                              target_dim = c(3, 3, 3)),
               "Supply both target_affine and target_dim")
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

# Synthetic dense FSL fields on a left-handed 2 mm reference. The mapping to
# source FSL coordinates is L %*% ref_fsl + t; the source image is a
# left-handed 1 mm grid (FSL coordinates = voxel indices) that contains it.
synthetic_fsl_field <- function(L, absolute) {
  dims <- c(20L, 24L, 20L)
  ref_aff <- diag(c(-2, 2, 2, 1))
  ref_aff[1:3, 4] <- c(38, -46, -38)
  vox <- as.matrix(expand.grid(lapply(dims, function(n) 0:(n - 1))))
  ref_fsl <- (cbind(vox, 1) %*% t(fsl_vox_to_fsl(ref_aff, dims)))[, 1:3]
  lin <- ref_fsl %*% t(L)
  src <- sweep(lin, 2, 5 - apply(lin, 2, min), "+")
  img <- RNifti::asNifti(array(if (absolute) src else src - ref_fsl, c(dims, 3L)))
  RNifti::sform(img) <- structure(ref_aff, code = 1L)
  RNifti::qform(img) <- structure(ref_aff, code = 1L)
  path <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(img, path, datatype = "float")
  list(path = path, source_affine = diag(c(-1, 1, 1, 1)),
       source_dim = as.integer(ceiling(apply(src, 2, max) + 5)) + 1L)
}

test_that("detect_fnirt_def_type identifies FSL absolute and relative fields", {
  skip_if_not_installed("RNifti")
  cyc <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), 3, byrow = TRUE)
  rz90 <- matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 1), 3, byrow = TRUE)
  cases <- list(
    list(diag(3), TRUE), list(diag(3), FALSE),
    # Scaled or axis-permuted mappings (small or sagittally stored sources)
    # leave the Jacobian test ambiguous; the field-of-view test decides.
    list(0.85 * cyc, TRUE), list(0.7 * rz90, TRUE), list(0.6 * diag(3), TRUE),
    list(0.75 * cyc, FALSE), list(0.7 * rz90, FALSE), list(1.3 * diag(3), FALSE)
  )
  for (case in cases) {
    f <- synthetic_fsl_field(case[[1]], case[[2]])
    expected <- if (case[[2]]) "absolute" else "relative"
    expect_equal(detect_fnirt_def_type(f$path, source_affine = f$source_affine,
                                       source_dim = f$source_dim), expected)
    m <- read_transform(f$path, source_affine = f$source_affine, source_dim = f$source_dim)
    expect_equal(m@params$def_type, expected)
  }
})

test_that("detect_fnirt_def_type refuses to guess ambiguous fields", {
  skip_if_not_installed("RNifti")
  cyc <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), 3, byrow = TRUE)
  f <- synthetic_fsl_field(0.85 * cyc, TRUE)
  expect_error(detect_fnirt_def_type(f$path), "pass def_type explicitly")
  expect_equal(detect_fnirt_def_type(synthetic_fsl_field(diag(3), TRUE)$path), "absolute")
})

test_that("FSL synthetic warp transform produces finite results", {
  skip_if_not_installed("RNifti")

  # Create a simple identity-like warp (small displacements)
  # With identity transform, warp coords = voxel coords
  # So we need to test at integer voxel coordinates
  dimf <- c(20, 20, 20, 3)
  arr <- array(0, dim = dimf)  # Zero displacement = identity

  tmp <- tempfile(fileext = ".nii.gz")
  nii <- RNifti::asNifti(arr)
  RNifti::qform(nii) <- structure(diag(4), code = 1L)
  RNifti::sform(nii) <- structure(diag(4), code = 1L)
  RNifti::writeNifti(nii, tmp)

  morph <- Warp3DMorphism("src", "tgt", warp_path = tmp, warp_type = "fsl",
                          def_type = "relative",
                          source_affine = diag(4), source_dim = dimf[1:3])

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
  RNifti::pixdim(nii) <- c(2, 2, 2, 1)
  RNifti::qform(nii) <- structure(xform_mat, code = 1L)
  RNifti::sform(nii) <- structure(xform_mat, code = 1L)
  RNifti::writeNifti(nii, tmp)

  morph <- Warp3DMorphism("src", "tgt", warp_path = tmp, warp_type = "fsl",
                          def_type = "relative",
                          source_affine = xform_mat, source_dim = dimf[1:3])

  # Test at a point within the warp volume (10mm from origin)
  test_coords <- matrix(c(10, 10, 10), ncol = 3)
  warped <- transform(morph, test_coords)

  # FSL positive X is negative RAS X on this right-handed source grid.
  expect_equal(warped[1, 1], 9, tolerance = 0.1)
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
