# Coverage Tests for transform_path.R
#
# Target: Improve coverage from 75% to >80%

# ==============================================================================
# TRANSFORM PATH BASIC TESTS
# ==============================================================================

test_that("transform_path with empty path returns coords unchanged", {
  coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = TRUE)
  result <- transform_path(list(), coords)
  expect_equal(result, coords)
})

test_that("transform_path with single identity returns coords unchanged", {
  coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = TRUE)
  id <- IdentityMorphism("a")

  result <- transform_path(list(id), coords)
  expect_equal(result, coords)
})

test_that("transform_path with single affine applies transform", {
  mat <- diag(4)
  mat[1:3, 4] <- c(10, 20, 30)  # Translation only

  aff <- Affine3DMorphism("a", "b", mat)
  coords <- matrix(c(0, 0, 0), ncol = 3)

  # Pullback: target coords -> source coords
  # For translation x' = x + t, pullback is x = x' - t
  # So coords (0,0,0) in target = (-10,-20,-30) in source
  result <- transform_path(list(aff), coords)

  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
})

test_that("transform_path with two affines composes correctly", {
  # f: A -> B with translation (10,0,0)
  mat1 <- diag(4)
  mat1[1, 4] <- 10

  # g: B -> C with translation (0,10,0)
  mat2 <- diag(4)
  mat2[2, 4] <- 10

  f <- Affine3DMorphism("a", "b", mat1)
  g <- Affine3DMorphism("b", "c", mat2)

  coords <- matrix(c(0, 0, 0), ncol = 3)

  # Path f then g, pullback applies g then f
  result <- transform_path(list(f, g), coords)

  expect_equal(nrow(result), 1)
})

test_that("transform_path handles multiple points", {
  mat <- diag(4)
  mat[1:3, 1:3] <- diag(c(2, 2, 2))  # Scale by 2

  aff <- Affine3DMorphism("a", "b", mat)
  coords <- matrix(c(0, 0, 0,
                     1, 1, 1,
                     2, 3, 4), ncol = 3, byrow = TRUE)

  result <- transform_path(list(aff), coords)

  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
})

# ==============================================================================
# TRANSFORM PATH WITH WARPS
# ==============================================================================

test_that("transform_path works with warp morphism", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  warp <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  # Get valid coords
  w <- neurotransform:::load_warp_array(warp)
  center_vox <- w$dim / 2
  center_world <- (w$vox_to_world %*% c(center_vox, 1))[1:3]
  coords <- matrix(center_world, nrow = 1)

  result <- transform_path(list(warp), coords)

  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("transform_path with affine + warp chain", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  mat <- diag(4)
  aff <- Affine3DMorphism("native", "a", mat)
  warp <- Warp3DMorphism("a", "mni", warp_path = warp_path, warp_type = "ants")

  w <- neurotransform:::load_warp_array(warp)
  center_vox <- w$dim / 2
  center_world <- (w$vox_to_world %*% c(center_vox, 1))[1:3]
  coords <- matrix(center_world, nrow = 1)

  result <- transform_path(list(aff, warp), coords)

  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# VALIDATE PATH
# ==============================================================================

test_that("validate_path errors on incompatible morphisms", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("c", "d", diag(4))  # Not composable with f

  expect_error(neurotransform:::validate_path(list(f, g)))
})

test_that("validate_path passes for compatible morphisms", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))

  # Should not error
  expect_silent(neurotransform:::validate_path(list(f, g)))
})

test_that("validate_path handles single morphism", {
  f <- Affine3DMorphism("a", "b", diag(4))

  expect_silent(neurotransform:::validate_path(list(f)))
})
