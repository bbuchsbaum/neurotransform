# Coverage Tests for warp_transform.R
#
# Target: Improve coverage from 51.35% to >80%

# ==============================================================================
# WARP TRANSFORM COORDS TESTS
# ==============================================================================

test_that("warp_transform_coords handles non-warp morphism", {
  # Should fall back to transform() for non-warp morphisms
  aff <- Affine3DMorphism("a", "b", diag(4))
  coords <- matrix(c(1, 2, 3), ncol = 3)

  # This should work via fallback
  result <- neurotransform:::warp_transform_coords(aff, coords)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
})

test_that("warp_transform_coords errors on missing warp_path", {
  # Create a warp morphism with empty path (edge case)
  morph <- methods::new("Warp3DMorphism",
                        source = "a",
                        target = "b",
                        warp_path = "",
                        warp_type = "ants",
                        params = list(),
                        inverse_type = "none",
                        cost = 1,
                        method_tag = "test",
                        cache = new.env(),
                        hash = "test")

  coords <- matrix(c(1, 2, 3), ncol = 3)

  expect_error(neurotransform:::warp_transform_coords(morph, coords), "warp_path")
})

test_that("warp_transform_coords handles ANTs warp", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  # Get valid coordinates
  warp <- neurotransform:::load_warp_array(morph)
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]

  coords <- matrix(center_world, nrow = 1)
  result <- neurotransform:::warp_transform_coords(morph, coords)

  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("warp_transform_coords handles cubic interpolation", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants",
                          warp_method = "cubic")

  warp <- neurotransform:::load_warp_array(morph)
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]

  coords <- matrix(center_world, nrow = 1)
  result <- neurotransform:::warp_transform_coords(morph, coords)

  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("warp_transform_coords handles AFNI warp via dispatch", {
  warp_path <- system.file("extdata/afni/sub-1001_T1w_in_mni_WARP.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "afni")

  # Test coords in MNI space
  coords <- matrix(c(0, 0, 0), ncol = 3)
  result <- neurotransform:::warp_transform_coords(morph, coords)

  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# COMPOSE WARPS TESTS
# ==============================================================================

test_that("compose_warps combines two warp fields", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  warpA <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")
  warpB <- Warp3DMorphism("b", "c", warp_path = warp_path, warp_type = "ants")

  result <- neurotransform:::compose_warps(warpB, warpA)

  expect_true(is.list(result))
  expect_true("array" %in% names(result))
  expect_true("dim" %in% names(result))
  expect_true("world_to_vox" %in% names(result))
})

# ==============================================================================
# ABSOLUTE VS RELATIVE DEFORMATION TESTS
# ==============================================================================

test_that("warp_transform_coords handles absolute deformation type", {
  skip_if_not_installed("RNifti")

  # Create synthetic absolute coordinate field
  dimf <- c(10, 10, 10, 3)
  arr <- array(0, dim = dimf)

  # Fill with absolute coordinates (voxel = world for identity affine)
  for (i in 1:dimf[1]) {
    for (j in 1:dimf[2]) {
      for (k in 1:dimf[3]) {
        arr[i, j, k, 1] <- i - 1  # X coordinate
        arr[i, j, k, 2] <- j - 1  # Y coordinate
        arr[i, j, k, 3] <- k - 1  # Z coordinate
      }
    }
  }

  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  nii <- RNifti::asNifti(arr)
  RNifti::sform(nii) <- diag(4)
  RNifti::writeNifti(nii, tmp)

  morph <- Warp3DMorphism("a", "b", warp_path = tmp, warp_type = "fsl",
                          def_type = "absolute")

  # Test at center
  coords <- matrix(c(5, 5, 5), ncol = 3)
  result <- transform(morph, coords)

  expect_true(all(is.finite(result)))
})

test_that("warp_transform_coords caches converted absolute displacement", {
  skip_if_not_installed("RNifti")

  # Create synthetic absolute field
  dimf <- c(5, 5, 5, 3)
  arr <- array(0, dim = dimf)
  for (i in 1:dimf[1]) {
    for (j in 1:dimf[2]) {
      for (k in 1:dimf[3]) {
        arr[i, j, k, 1] <- i - 1
        arr[i, j, k, 2] <- j - 1
        arr[i, j, k, 3] <- k - 1
      }
    }
  }

  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  nii <- RNifti::asNifti(arr)
  RNifti::sform(nii) <- diag(4)
  RNifti::writeNifti(nii, tmp)

  morph <- Warp3DMorphism("a", "b", warp_path = tmp, warp_type = "fsl",
                          def_type = "absolute")

  coords <- matrix(c(2, 2, 2), ncol = 3)

  # First call
  result1 <- transform(morph, coords)

  # Second call should use cache
  result2 <- transform(morph, coords)

  expect_equal(result1, result2)
})

test_that("warp_transform_coords applies absolute deformations exactly on non-identity grids", {
  aff <- diag(4)
  aff[1:3, 1:3] <- diag(c(2, 3, 4))
  aff[1:3, 4] <- c(10, -6, 4)
  grid <- grid_spec(dims = c(6L, 6L, 6L), affine = aff)

  delta <- c(1.25, -0.75, 0.5)
  world <- grid_coords(grid)
  field <- array(0, dim = c(grid@dims, 3L))
  for (k in 1:3) {
    field[, , , k] <- array(world[, k] + delta[k], dim = grid@dims)
  }

  morph <- warp_from_field("src", "tgt", field, grid = grid,
                           representation = "deformations")

  coords <- matrix(c(
    14, 0, 8,
    18, 6, 16
  ), ncol = 3, byrow = TRUE)

  result <- transform(morph, coords)

  expect_equal(result,
               sweep(coords, 2, delta, "+"),
               tolerance = 1e-8)

  cache_key <- paste0(morph@warp_path, "::relative")
  expect_true(exists(cache_key, envir = morph@cache, inherits = FALSE))
})

# ==============================================================================
# EMBEDDED AFFINE HANDLING
# ==============================================================================

test_that("warp_transform_coords does not apply embedded affine", {
  h5_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                         package = "neurotransform")
  skip_if_not(file.exists(h5_path))
  skip_if_not_installed("hdf5r")

  # Direct warp without affine application
  morph <- Warp3DMorphism("a", "b", warp_path = h5_path, warp_type = "ants_h5")

  warp <- neurotransform:::load_warp_array(morph)
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]

  coords <- matrix(center_world, nrow = 1)

  # This should work without applying embedded affine
  result <- transform(morph, coords)

  expect_true(all(is.finite(result)))
})

# ==============================================================================
# CACHING TESTS
# ==============================================================================

test_that("warp_transform_coords uses morphism cache", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  warp <- neurotransform:::load_warp_array(morph)
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]

  coords <- matrix(center_world, nrow = 1)

  # First transform loads and caches
  result1 <- transform(morph, coords)

  # Second transform should use cache
  result2 <- transform(morph, coords)

  expect_equal(result1, result2)
})

# ==============================================================================
# EDGE CASES
# ==============================================================================

test_that("warp_transform handles multiple coordinates", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  warp <- neurotransform:::load_warp_array(morph)
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]

  # Multiple coordinates
  offsets <- matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, byrow = TRUE)
  coords <- sweep(offsets, 2, center_world, "+")

  result <- transform(morph, coords)

  expect_equal(nrow(result), 3)
  expect_true(all(is.finite(result)))
})

test_that("warp_transform handles out-of-bounds coordinates gracefully", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  # Way outside the warp field
  coords <- matrix(c(10000, 10000, 10000), ncol = 3)

  # Should not error, may return NA or extrapolated values
  result <- transform(morph, coords)
  expect_equal(nrow(result), 1)
})
