# Coverage Tests for jacobian.R
#
# Target: Improve coverage from 34.21% to >80%
# Focus: JacobianField class methods, MorphismPath Jacobian, edge cases

# ==============================================================================
# JACOBIAN FIELD CLASS TESTS
# ==============================================================================

test_that("JacobianField show method works", {
  # Create a simple JacobianField
  values <- array(0, dim = c(2, 3, 3))
  values[1, , ] <- diag(3)
  values[2, , ] <- diag(3) * 2

  coords <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE)

  jf <- new("JacobianField",
            values = values,
            coords = coords,
            mode = "pullback",
            morphism_hash = "test")

  # Show should not error
  expect_output(show(jf), "JacobianField")
  expect_output(show(jf), "n=2")
  expect_output(show(jf), "pullback")
})

test_that("JacobianField show method works with empty field", {
  jf <- new("JacobianField",
            values = array(NA_real_, dim = c(0, 3, 3)),
            coords = matrix(NA_real_, nrow = 0, ncol = 3),
            mode = "pullback",
            morphism_hash = "")

  # Should not error on empty field
  expect_output(show(jf), "n=0")
})

test_that("JacobianField subsetting works", {
  values <- array(0, dim = c(3, 3, 3))
  values[1, , ] <- diag(3)
  values[2, , ] <- diag(3) * 2
  values[3, , ] <- diag(3) * 3

  coords <- matrix(1:9, ncol = 3, byrow = TRUE)

  jf <- new("JacobianField",
            values = values,
            coords = coords,
            mode = "pullback",
            morphism_hash = "test")

  # Single element extraction (drop = TRUE)
  single <- jf[1]
  expect_true(is.matrix(single))
  expect_equal(dim(single), c(3, 3))
  expect_equal(single, diag(3))

  # Multiple element extraction
  subset <- jf[1:2]
  expect_s4_class(subset, "JacobianField")
  expect_equal(length(subset), 2)

  # Subsetting with j parameter
  elem <- jf[1, 1, 1]
  expect_equal(elem, 1)
})

test_that("JacobianField length method works", {
  values <- array(0, dim = c(5, 3, 3))
  coords <- matrix(1:15, ncol = 3)

  jf <- new("JacobianField",
            values = values,
            coords = coords,
            mode = "pullback",
            morphism_hash = "test")

  expect_equal(length(jf), 5)
})

test_that("JacobianField det method works", {
  values <- array(0, dim = c(2, 3, 3))
  values[1, , ] <- diag(3)       # det = 1
  values[2, , ] <- diag(3) * 2   # det = 8

  coords <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE)

  jf <- new("JacobianField",
            values = values,
            coords = coords,
            mode = "pullback",
            morphism_hash = "test")

  dets <- det(jf)
  expect_equal(dets, c(1, 8))
})

test_that("JacobianField solve method inverts and switches mode", {
  values <- array(0, dim = c(2, 3, 3))
  values[1, , ] <- diag(3)
  values[2, , ] <- diag(c(2, 3, 4))

  coords <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE)

  jf <- new("JacobianField",
            values = values,
            coords = coords,
            mode = "pullback",
            morphism_hash = "test")

  inv_jf <- solve(jf)

  # Mode should switch
  expect_equal(inv_jf@mode, "pushforward")

  # Values should be inverted
  expect_equal(inv_jf[1], diag(3))
  expect_equal(inv_jf[2], diag(c(0.5, 1/3, 0.25)), tolerance = 1e-10)

  # Double inversion should recover original mode
  double_inv <- solve(inv_jf)
  expect_equal(double_inv@mode, "pullback")
})

# ==============================================================================
# IDENTITY MORPHISM JACOBIAN TESTS
# ==============================================================================

test_that("IdentityMorphism jacobian returns identity matrices", {
  id <- IdentityMorphism("domain")
  coords <- matrix(c(0, 0, 0, 10, 20, 30), ncol = 3, byrow = TRUE)

  jac <- jacobian(id, coords)

  expect_s4_class(jac, "JacobianField")
  expect_equal(length(jac), 2)
  expect_equal(jac[1], diag(3))
  expect_equal(jac[2], diag(3))
})

test_that("IdentityMorphism jacobian_det returns 1 or log(1)=0", {
  id <- IdentityMorphism("domain")
  coords <- matrix(c(0, 0, 0, 10, 20, 30), ncol = 3, byrow = TRUE)

  dets <- jacobian_det(id, coords)
  expect_equal(dets, c(1, 1))

  log_dets <- jacobian_det(id, coords, log = TRUE)
  expect_equal(log_dets, c(0, 0))
})

test_that("IdentityMorphism jacobian works in pushforward mode", {
  id <- IdentityMorphism("domain")
  coords <- matrix(c(0, 0, 0), ncol = 3)

  jac <- jacobian(id, coords, mode = "pushforward")
  expect_equal(jac@mode, "pushforward")
  expect_equal(jac[1], diag(3))
})

# ==============================================================================
# AFFINE MORPHISM JACOBIAN TESTS
# ==============================================================================

test_that("Affine3DMorphism jacobian returns constant linear part", {
  mat <- diag(4)
  mat[1:3, 1:3] <- matrix(c(2, 0, 0, 0, 3, 0, 0, 0, 4), nrow = 3)
  mat[1:3, 4] <- c(10, 20, 30)  # translation doesn't affect Jacobian

  aff <- Affine3DMorphism("a", "b", mat)
  coords <- matrix(c(0, 0, 0, 100, 100, 100), ncol = 3, byrow = TRUE)

  jac <- jacobian(aff, coords)

  # Jacobian should be the linear part
  expect_equal(jac[1], mat[1:3, 1:3])
  expect_equal(jac[2], mat[1:3, 1:3])  # Same at all points
})

test_that("Affine3DMorphism jacobian_det is determinant of linear part", {
  mat <- diag(4)
  mat[1:3, 1:3] <- diag(c(2, 3, 4))

  aff <- Affine3DMorphism("a", "b", mat)
  coords <- matrix(c(0, 0, 0), ncol = 3)

  det_val <- jacobian_det(aff, coords)
  expect_equal(det_val, 24)  # 2 * 3 * 4

  log_det <- jacobian_det(aff, coords, log = TRUE)
  expect_equal(log_det, log(24))
})

test_that("Affine3DMorphism jacobian pushforward mode inverts", {
  mat <- diag(4)
  mat[1:3, 1:3] <- diag(c(2, 3, 4))

  aff <- Affine3DMorphism("a", "b", mat)
  coords <- matrix(c(0, 0, 0), ncol = 3)

  jac_pull <- jacobian(aff, coords, mode = "pullback")
  jac_push <- jacobian(aff, coords, mode = "pushforward")

  # Pushforward should be inverse of pullback
  expect_equal(jac_push[1], solve(jac_pull[1]), tolerance = 1e-10)
})

test_that("Affine3DMorphism jacobian_det pushforward mode inverts", {
  mat <- diag(4)
  mat[1:3, 1:3] <- diag(c(2, 3, 4))

  aff <- Affine3DMorphism("a", "b", mat)
  coords <- matrix(c(0, 0, 0), ncol = 3)

  det_pull <- jacobian_det(aff, coords, mode = "pullback")
  det_push <- jacobian_det(aff, coords, mode = "pushforward")

  expect_equal(det_push, 1 / det_pull)
})

# ==============================================================================
# MORPHISM PATH JACOBIAN (CHAIN RULE)
# ==============================================================================

test_that("MorphismPath jacobian applies chain rule correctly", {
  # Two affine transforms: f then g
  mat1 <- diag(4)
  mat1[1:3, 1:3] <- diag(c(2, 1, 1))  # Scale x by 2

  mat2 <- diag(4)
  mat2[1:3, 1:3] <- diag(c(1, 3, 1))  # Scale y by 3

  f <- Affine3DMorphism("a", "b", mat1)
  g <- Affine3DMorphism("b", "c", mat2)

  path <- compose(f, g)
  coords <- matrix(c(0, 0, 0), ncol = 3)

  jac <- jacobian(path, coords)

  # Chain rule: J(g . f) = J(f) %*% J(g) evaluated at appropriate coords
  # For affines, this should be the product of the linear parts
  expected <- mat1[1:3, 1:3] %*% mat2[1:3, 1:3]
  expect_equal(jac[1], expected)
})

test_that("MorphismPath jacobian_det matches det of jacobian", {
  mat1 <- diag(4)
  mat1[1:3, 1:3] <- diag(c(2, 2, 2))

  mat2 <- diag(4)
  mat2[1:3, 1:3] <- diag(c(3, 3, 3))

  f <- Affine3DMorphism("a", "b", mat1)
  g <- Affine3DMorphism("b", "c", mat2)

  path <- compose(f, g)
  coords <- matrix(c(0, 0, 0), ncol = 3)

  jac <- jacobian(path, coords)
  det_from_jac <- det(jac)

  det_direct <- jacobian_det(path, coords)

  expect_equal(det_from_jac, det_direct)
  expect_equal(det_direct, 8 * 27)  # (2^3) * (3^3)
})

test_that("MorphismPath jacobian pushforward mode works", {
  mat <- diag(4)
  mat[1:3, 1:3] <- diag(c(2, 3, 4))

  f <- Affine3DMorphism("a", "b", mat)
  g <- Affine3DMorphism("b", "c", mat)

  path <- compose(f, g)
  coords <- matrix(c(0, 0, 0), ncol = 3)

  jac_push <- jacobian(path, coords, mode = "pushforward")

  expect_equal(jac_push@mode, "pushforward")
})

test_that("Empty MorphismPath jacobian returns identity", {
  # Create empty path (though this is an edge case)
  empty_path <- methods::new("MorphismPath",
                             morphisms = list(),
                             source = "a",
                             target = "a")

  coords <- matrix(c(0, 0, 0), ncol = 3)
  jac <- jacobian(empty_path, coords)

  expect_equal(jac[1], diag(3))
})

# ==============================================================================
# SURFACE MORPHISM JACOBIAN ERRORS
# ==============================================================================

test_that("VolToSurfMorphism jacobian throws appropriate error", {
  morph <- VolToSurfMorphism("vol", "surf", method = "trilinear")
  coords <- matrix(c(0, 0, 0), ncol = 3)

  expect_error(jacobian(morph, coords), "not defined")
})

test_that("SurfToSurfMorphism jacobian throws appropriate error", {
  morph <- SurfToSurfMorphism("surf1", "surf2", method = "sphere")
  coords <- matrix(c(0, 0, 0), ncol = 3)

  expect_error(jacobian(morph, coords), "not implemented")
})

# ==============================================================================
# WARP MORPHISM JACOBIAN TESTS
# ==============================================================================

test_that("Warp3DMorphism jacobian_det log mode works", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  warp <- suppressWarnings(neurotransform:::load_warp_array(morph))
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]

  coords <- matrix(center_world, nrow = 1)

  det_normal <- jacobian_det(morph, coords, log = FALSE)
  det_log <- jacobian_det(morph, coords, log = TRUE)

  expect_equal(det_log, log(abs(det_normal)), tolerance = 1e-6)
})

test_that("Warp3DMorphism jacobian pushforward requires invertible morphism",
{
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  # Create non-invertible warp (no inverse_path)
  morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  warp <- suppressWarnings(neurotransform:::load_warp_array(morph))
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]

  coords <- matrix(center_world, nrow = 1)

  # Pushforward should error without inverse
  expect_error(jacobian(morph, coords, mode = "pushforward"), "invertible")
})
