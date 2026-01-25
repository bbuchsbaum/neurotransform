# Coverage Tests for sampler.R
#
# Target: Improve coverage from 76.22% to >80%

# ==============================================================================
# SAMPLER CLASS TESTS
# ==============================================================================

test_that("Sampler show method works", {
  vol <- array(1:27, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4))

  expect_output(show(sampler), "Sampler")
  expect_output(show(sampler), "dims=")
  expect_output(show(sampler), "method=")
})

test_that("Sampler prototype has default values", {
  s <- new("Sampler")

  expect_equal(s@domain, "")
  expect_equal(s@method, "linear")
  expect_equal(length(s@dims), 0)
})

# ==============================================================================
# GRID CLASS TESTS
# ==============================================================================

test_that("Grid show method works", {
  grid <- grid_spec(c(10, 20, 30), diag(c(2, 2, 2, 1)))

  expect_output(show(grid), "Grid")
  expect_output(show(grid), "dims=")
  expect_output(show(grid), "spacing=")
})

test_that("grid_spec creates valid Grid", {
  grid <- grid_spec(c(5, 6, 7), diag(4))

  expect_s4_class(grid, "Grid")
  expect_equal(grid@dims, c(5L, 6L, 7L))
  expect_equal(grid@affine, diag(4))
  expect_true(nzchar(grid@domain))
})

test_that("grid_spec accepts custom domain", {
  grid <- grid_spec(c(5, 5, 5), diag(4), domain = "custom_domain")
  expect_equal(grid@domain, "custom_domain")
})

test_that("grid_coords returns correct world coordinates", {
  grid <- grid_spec(c(2, 2, 2), diag(4))
  coords <- grid_coords(grid)

  expect_equal(nrow(coords), 8)  # 2^3
  expect_equal(ncol(coords), 3)
})

test_that("grid_coords respects affine", {
  # 2mm spacing with offset
  aff <- diag(c(2, 2, 2, 1))
  aff[1:3, 4] <- c(10, 20, 30)
  grid <- grid_spec(c(2, 2, 2), aff)

  coords <- grid_coords(grid)

  # First voxel (0,0,0) should map to (10,20,30)
  first_coord <- coords[1, ]
  expect_equal(first_coord, c(10, 20, 30))
})

test_that("grid_from_data creates grid from array", {
  arr <- array(1, dim = c(5, 6, 7))
  grid <- grid_from_data(arr)

  expect_s4_class(grid, "Grid")
  expect_equal(grid@dims, c(5L, 6L, 7L))
})

# ==============================================================================
# EXTRACT AFFINE TESTS
# ==============================================================================

test_that("extract_affine returns identity for array", {
  arr <- array(1, dim = c(3, 3, 3))
  aff <- extract_affine(arr)

  expect_equal(aff, diag(4))
})

test_that("extract_affine returns matrix as-is if 4x4", {
  mat <- diag(c(2, 2, 2, 1))
  result <- extract_affine(mat)

  expect_equal(result, mat)
})

test_that("extract_affine returns identity for non-4x4 matrix", {
  mat <- matrix(1:9, 3, 3)
  result <- extract_affine(mat)

  expect_equal(result, diag(4))
})

test_that("extract_affine uses RNifti for niftiImage", {
  skip_if_not_installed("RNifti")

  arr <- array(1, dim = c(3, 3, 3))
  nii <- RNifti::asNifti(arr)

  # extract_affine should return a 4x4 matrix
  result <- extract_affine(nii)
  expect_equal(dim(result), c(4, 4))
  expect_true(is.matrix(result))
})

test_that("extract_affine uses neuroim2 for DenseNeuroVol", {
  skip_if_not_installed("neuroim2")

  arr <- array(1, dim = c(3, 3, 3))
  aff <- diag(c(3, 3, 3, 1))
  space <- neuroim2::NeuroSpace(c(3, 3, 3), trans = aff)
  vol <- neuroim2::DenseNeuroVol(arr, space)

  result <- extract_affine(vol)
  expect_equal(result, aff)
})

# ==============================================================================
# VOLUME SAMPLER TESTS
# ==============================================================================

test_that("volume_sampler creates valid sampler", {
  vol <- array(1:27, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), method = "linear")

  expect_s4_class(sampler, "Sampler")
  expect_equal(sampler@dims, c(3L, 3L, 3L))
  expect_equal(sampler@vdim, 1L)
  expect_equal(sampler@method, "linear")
})

test_that("volume_sampler extracts affine if not provided", {
  skip_if_not_installed("RNifti")

  arr <- array(1, dim = c(3, 3, 3))
  nii <- RNifti::asNifti(arr)
  RNifti::sform(nii) <- diag(c(2, 2, 2, 1))

  sampler <- volume_sampler(nii)
  expect_s4_class(sampler, "Sampler")
})

test_that("volume_sampler supports nearest interpolation", {
  vol <- array(1:27, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), method = "nearest")

  coords <- matrix(c(0.5, 0.5, 0.5), ncol = 3)
  result <- sampler@evaluate(coords)

  expect_true(is.numeric(result))
})

test_that("volume_sampler supports cubic interpolation", {
  vol <- array(rnorm(64), dim = c(4, 4, 4))
  sampler <- volume_sampler(vol, affine = diag(4), method = "cubic")

  coords <- matrix(c(1.5, 1.5, 1.5), ncol = 3)
  result <- sampler@evaluate(coords)

  expect_true(is.numeric(result))
})

test_that("volume_sampler handles 4D data", {
  vol <- array(rnorm(54), dim = c(3, 3, 3, 2))
  sampler <- volume_sampler(vol, affine = diag(4))

  expect_equal(sampler@vdim, 2L)
})

test_that("volume_sampler respects outside value", {
  vol <- array(1:27, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), outside = -999)

  # Query outside bounds
  coords <- matrix(c(100, 100, 100), ncol = 3)
  result <- sampler@evaluate(coords)

  # Should get outside value
  expect_equal(as.numeric(result), -999)
})

test_that("volume_sampler computes correct bounds", {
  vol <- array(1, dim = c(10, 10, 10))
  aff <- diag(c(2, 2, 2, 1))
  sampler <- volume_sampler(vol, affine = aff)

  # Bounds should span 0 to 18mm (0-9 voxels * 2mm)
  expect_equal(sampler@bounds$min, c(0, 0, 0))
  expect_equal(sampler@bounds$max, c(18, 18, 18))
})

# ==============================================================================
# SURFACE SAMPLER TESTS
# ==============================================================================

test_that("surface_sampler creates valid sampler with nearest", {
  vertices <- matrix(c(0, 0, 0,
                       1, 0, 0,
                       0, 1, 0,
                       0, 0, 1), ncol = 3, byrow = TRUE)
  data <- c(1, 2, 3, 4)

  sampler <- surface_sampler(vertices, data, method = "nearest")

  expect_s4_class(sampler, "Sampler")
  expect_equal(sampler@method, "nearest")
})

test_that("surface_sampler nearest finds closest vertex", {
  vertices <- matrix(c(0, 0, 0,
                       10, 0, 0,
                       0, 10, 0), ncol = 3, byrow = TRUE)
  data <- c(100, 200, 300)

  sampler <- surface_sampler(vertices, data, method = "nearest")

  # Query near first vertex
  coords <- matrix(c(0.1, 0.1, 0.1), ncol = 3)
  result <- sampler@evaluate(coords)

  expect_equal(result, 100)
})

test_that("surface_sampler errors on barycentric without faces", {
  vertices <- matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0), ncol = 3, byrow = TRUE)
  data <- c(1, 2, 3)

  expect_error(surface_sampler(vertices, data, method = "barycentric"),
               "faces")
})

test_that("surface_sampler handles matrix data", {
  vertices <- matrix(c(0, 0, 0,
                       1, 0, 0,
                       0, 1, 0), ncol = 3, byrow = TRUE)
  data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)

  sampler <- surface_sampler(vertices, data, method = "nearest")

  expect_equal(sampler@vdim, 2L)
})

# ==============================================================================
# RESAMPLE FUNCTION TESTS
# ==============================================================================

test_that("resample works with identity morphism", {
  vol <- array(1:27, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), domain = "test")
  morph <- IdentityMorphism("test")

  coords <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE)
  result <- resample(sampler, morph, coords)

  expect_equal(length(result), 2)
})

test_that("resample works with Grid input", {
  vol <- array(1:27, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), domain = "test")
  morph <- IdentityMorphism("test")
  grid <- grid_spec(c(2, 2, 2), diag(4), domain = "test")

  result <- resample(sampler, morph, grid)

  expect_equal(length(result), 8)  # 2^3
})

test_that("resample with jacobian modulation", {
  vol <- array(rep(1, 27), dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), domain = "a")

  # Scaling morphism: det = 8
  mat <- diag(c(2, 2, 2, 1))
  morph <- Affine3DMorphism("a", "b", mat)

  coords <- matrix(c(1, 1, 1), ncol = 3)
  result_none <- resample(sampler, morph, coords, modulate = "none")
  result_jac <- resample(sampler, morph, coords, modulate = "jacobian")

  # Jacobian modulated should be scaled by det
  expect_equal(result_jac, result_none * 8)
})

test_that("resample with sqrt_jacobian modulation", {
  vol <- array(rep(1, 27), dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), domain = "a")

  mat <- diag(c(2, 2, 2, 1))
  morph <- Affine3DMorphism("a", "b", mat)

  coords <- matrix(c(1, 1, 1), ncol = 3)
  result_sqrt <- resample(sampler, morph, coords, modulate = "sqrt_jacobian")

  # sqrt(8) = 2.83...
  expect_equal(result_sqrt, sqrt(8), tolerance = 1e-6)
})

test_that("resample warns on domain mismatch", {
  vol <- array(1:27, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), domain = "domain_a")
  morph <- IdentityMorphism("domain_b")

  coords <- matrix(c(1, 1, 1), ncol = 3)

  expect_warning(resample(sampler, morph, coords), "domain")
})

# ==============================================================================
# RESAMPLE VOLUME TESTS
# ==============================================================================

test_that("resample_volume works with identity", {
  vol <- array(1:27, dim = c(3, 3, 3))
  morph <- IdentityMorphism("test")
  target <- grid_spec(c(2, 2, 2), diag(4))

  result <- resample_volume(vol, morph, target)

  expect_equal(dim(result), c(2, 2, 2))
})

test_that("resample_volume accepts volume as target", {
  skip_if_not_installed("neuroim2")

  vol <- array(1:27, dim = c(3, 3, 3))
  morph <- IdentityMorphism("test")

  # Target as array
  target <- array(0, dim = c(4, 4, 4))
  result <- resample_volume(vol, morph, target)

  expect_equal(dim(result), c(4, 4, 4))
})

# ==============================================================================
# SAMPLE VOLUME ON SURFACE TESTS
# ==============================================================================

test_that("sample_volume_on_surface errors on non-VolToSurf morphism", {
  vol <- array(1:27, dim = c(3, 3, 3))
  morph <- IdentityMorphism("test")

  expect_error(sample_volume_on_surface(vol, morph), "VolToSurfMorphism")
})

test_that("sample_volume_on_surface errors without coords", {
  vol <- array(1:27, dim = c(3, 3, 3))
  morph <- VolToSurfMorphism("vol", "surf", method = "trilinear")

  expect_error(sample_volume_on_surface(vol, morph), "coordinates")
})

test_that("sample_volume_on_surface works with provided coords", {
  vol <- array(1:125, dim = c(5, 5, 5))
  morph <- VolToSurfMorphism("vol", "surf", method = "trilinear")

  # Surface coords (in world space, identity affine)
  surface_coords <- matrix(c(2, 2, 2,
                             3, 3, 3), ncol = 3, byrow = TRUE)

  result <- sample_volume_on_surface(vol, morph, surface_coords)

  expect_equal(length(result), 2)
})
