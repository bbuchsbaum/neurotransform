# Coverage Tests for io_helpers.R
#
# Target: Improve coverage from 51.19% to >80%

# ==============================================================================
# READ_IMAGE TESTS
# ==============================================================================

test_that("read_image reads existing volume", {
  path <- system.file("extdata/chris/chris_t1.nii.gz", package = "neurotransform")
  skip_if_not(file.exists(path))

  vol <- read_image(path)
  expect_true(inherits(vol, "DenseNeuroVol") || inherits(vol, "DenseNeuroVec"))
})

test_that("read_image errors on missing file", {
  expect_error(read_image("/nonexistent/path.nii.gz"), "not found")
})

# ==============================================================================
# WRITE_IMAGE TESTS
# ==============================================================================

test_that("write_image writes DenseNeuroVol", {
  skip_if_not_installed("neuroim2")

  # Create simple volume
  arr <- array(rnorm(1000), dim = c(10, 10, 10))
  space <- neuroim2::NeuroSpace(c(10, 10, 10), trans = diag(4))
  vol <- neuroim2::DenseNeuroVol(arr, space)

  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_image(vol, tmp)
  expect_true(file.exists(tmp))
})

test_that("write_image writes array with default affine", {
  skip_if_not_installed("neuroim2")

  arr <- array(rnorm(1000), dim = c(10, 10, 10))

  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_image(arr, tmp)
  expect_true(file.exists(tmp))
})

test_that("write_image writes 4D array", {
  skip_if_not_installed("neuroim2")

  arr <- array(rnorm(2000), dim = c(10, 10, 10, 2))

  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_image(arr, tmp)
  expect_true(file.exists(tmp))
})

# ==============================================================================
# GRID_OF TESTS
# ==============================================================================

test_that("grid_of returns Grid unchanged", {
  grid <- grid_spec(c(10, 10, 10), diag(4))
  result <- grid_of(grid)
  expect_identical(result, grid)
})

test_that("grid_of extracts grid from DenseNeuroVol", {
  skip_if_not_installed("neuroim2")

  arr <- array(1, dim = c(10, 10, 10))
  aff <- diag(c(2, 2, 2, 1))
  space <- neuroim2::NeuroSpace(c(10, 10, 10), trans = aff)
  vol <- neuroim2::DenseNeuroVol(arr, space)

  grid <- grid_of(vol)
  expect_s4_class(grid, "Grid")
  expect_equal(grid@dims, c(10, 10, 10))
})

test_that("grid_of extracts grid from plain array", {
  arr <- array(1, dim = c(5, 6, 7))
  grid <- grid_of(arr)

  expect_s4_class(grid, "Grid")
  expect_equal(grid@dims, c(5, 6, 7))
  expect_equal(grid@affine, diag(4))
})

test_that("grid_of errors on unsupported type", {
  expect_error(grid_of("string"), "Unsupported type")
  expect_error(grid_of(42), "Unsupported type")
})

# ==============================================================================
# READ_TRANSFORM TESTS
# ==============================================================================

test_that("read_transform passes through Morphism", {
  aff <- Affine3DMorphism("a", "b", diag(4))
  result <- read_transform(aff)
  expect_identical(result, aff)
})

test_that("read_transform errors on invalid path type", {
  expect_error(read_transform(c("a", "b")), "must be a file path")
  expect_error(read_transform(42), "must be a file path")
})

test_that("read_transform errors on missing file", {
  expect_error(read_transform("/nonexistent/file.nii.gz"), "not found")
})

test_that("read_transform auto-detects H5 type", {
  h5_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                         package = "neurotransform")
  skip_if_not(file.exists(h5_path))
  skip_if_not_installed("hdf5r")

  result <- read_transform(h5_path)
  expect_true(is(result, "Morphism") || is(result, "MorphismPath"))
})

test_that("read_transform auto-detects NIfTI as ANTs type", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  result <- read_transform(warp_path)
  expect_s4_class(result, "Warp3DMorphism")
})

test_that("read_transform respects explicit type", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  result <- read_transform(warp_path, type = "fsl")
  expect_s4_class(result, "Warp3DMorphism")
  expect_equal(result@warp_type, "fsl")
})

test_that("read_transform errors on unsupported type", {
  tmp <- tempfile(fileext = ".xyz")
  file.create(tmp)
  on.exit(unlink(tmp))

  expect_error(read_transform(tmp, type = "unsupported"), "Unsupported")
})

# ==============================================================================
# AS_MORPHISM_PATH TESTS
# ==============================================================================

test_that("as_morphism_path creates path from morphisms", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))

  path <- as_morphism_path(f, g)

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)
  expect_equal(source_of(path), "a")
  expect_equal(target_of(path), "c")
})

test_that("as_morphism_path handles list input", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))

  path <- as_morphism_path(list(f, g))

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)
})

test_that("as_morphism_path loads from file paths", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  path <- as_morphism_path(warp_path)
  expect_s4_class(path, "MorphismPath")
})

# ==============================================================================
# RESAMPLE_TO TESTS
# ==============================================================================

test_that("resample_to works with file path transform", {
  src_path <- system.file("extdata/chris/chris_t1.nii.gz", package = "neurotransform")
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")

  skip_if_not(file.exists(src_path))
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("neuroim2")

  src <- neuroim2::read_vol(src_path)

  # Create target grid
  target_grid <- grid_spec(c(20, 20, 20), diag(c(5, 5, 5, 1)))

  # Should work with path string
  result <- resample_to(src, target_grid, warp_path)

  expect_true(inherits(result, "DenseNeuroVol"))
  expect_equal(dim(result)[1:3], c(20, 20, 20))
})

test_that("resample_to returns 3D for single-volume 4D output", {
  skip_if_not_installed("neuroim2")

  # Create simple 3D volume
  arr <- array(rnorm(1000), dim = c(10, 10, 10))
  space <- neuroim2::NeuroSpace(c(10, 10, 10), trans = diag(4))
  src <- neuroim2::DenseNeuroVol(arr, space)

  # Identity transform
  morph <- IdentityMorphism("test")
  target_grid <- grid_spec(c(5, 5, 5), diag(c(2, 2, 2, 1)))

  result <- resample_to(src, target_grid, morph)

  expect_true(inherits(result, "DenseNeuroVol"))
  expect_equal(length(dim(result)), 3)
})

# ==============================================================================
# ANTS_H5_MORPHISM TESTS
# ==============================================================================

test_that("ants_h5_morphism with apply_affine=FALSE returns single warp", {
  h5_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                         package = "neurotransform")
  skip_if_not(file.exists(h5_path))
  skip_if_not_installed("hdf5r")

  result <- ants_h5_morphism(h5_path, apply_affine = FALSE)

  expect_s4_class(result, "Warp3DMorphism")
})

test_that("ants_h5_morphism uses custom source/target", {
  h5_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                         package = "neurotransform")
  skip_if_not(file.exists(h5_path))
  skip_if_not_installed("hdf5r")

  result <- ants_h5_morphism(h5_path, source = "native_space", target = "mni_space",
                             apply_affine = FALSE)

  expect_equal(source_of(result), "native_space")
  expect_equal(target_of(result), "mni_space")
})
