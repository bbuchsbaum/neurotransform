# Test coordinate conversion functions

# =============================================================================
# LPS <-> RAS CONVERSION TESTS
# =============================================================================

test_that("lps_to_ras flips X and Y for vector", {
  lps <- c(-10, -20, 30)
  ras <- lps_to_ras(lps)

  expect_equal(ras[1], 10)   # X flipped

  expect_equal(ras[2], 20)   # Y flipped
  expect_equal(ras[3], 30)   # Z unchanged
})

test_that("lps_to_ras flips X and Y for matrix", {
  lps <- matrix(c(
    -10, -20, 30,
    5, -15, 25,
    0, 0, 0
  ), ncol = 3, byrow = TRUE)

  ras <- lps_to_ras(lps)

  expect_equal(ras[1, 1], 10)
  expect_equal(ras[1, 2], 20)
  expect_equal(ras[1, 3], 30)
  expect_equal(ras[2, 1], -5)
  expect_equal(ras[2, 2], 15)
})

test_that("ras_to_lps flips X and Y for vector", {
  ras <- c(10, 20, 30)
  lps <- ras_to_lps(ras)

  expect_equal(lps[1], -10)  # X flipped
  expect_equal(lps[2], -20)  # Y flipped
  expect_equal(lps[3], 30)   # Z unchanged
})

test_that("ras_to_lps flips X and Y for matrix", {
  ras <- matrix(c(
    10, 20, 30,
    -5, 15, 25
  ), ncol = 3, byrow = TRUE)

  lps <- ras_to_lps(ras)

  expect_equal(lps[1, 1], -10)
  expect_equal(lps[1, 2], -20)
  expect_equal(lps[2, 1], 5)
  expect_equal(lps[2, 2], -15)
})

test_that("lps_to_ras round-trip preserves coordinates", {
  original <- c(10, 20, 30)

  # RAS -> LPS -> RAS
  lps <- ras_to_lps(original)
  back <- lps_to_ras(lps)

  expect_equal(back, original)
})

test_that("lps_to_ras round-trip preserves matrix coordinates", {
  original <- matrix(c(
    10, 20, 30,
    -5, -15, 25,
    0, 0, 0
  ), ncol = 3, byrow = TRUE)

  back <- lps_to_ras(ras_to_lps(original))

  expect_equal(back, original)
})

test_that("lps_to_ras errors on wrong vector length", {
  expect_error(lps_to_ras(c(1, 2)), "length 3")
  expect_error(lps_to_ras(c(1, 2, 3, 4)), "length 3")
})

test_that("lps_to_ras errors on wrong matrix columns", {
  expect_error(lps_to_ras(matrix(1:6, ncol = 2)), "3 columns")
})

test_that("lps_to_ras errors on non-vector/matrix input", {
  # list(1,2,3) passes is.vector but fails on numeric operations
  # Data frames should fail
  expect_error(lps_to_ras(data.frame(x = 1, y = 2, z = 3)), "must be a vector or matrix")
})

# =============================================================================
# TKRAS <-> RAS CONVERSION TESTS
# =============================================================================

test_that("tkras_to_ras adds c_ras offset for vector", {
  tkras <- c(0, 0, 0)
  c_ras <- c(0.5, -2.3, 1.1)

  ras <- tkras_to_ras(tkras, c_ras)

  expect_equal(ras, c(0.5, -2.3, 1.1))
})

test_that("tkras_to_ras adds c_ras offset for matrix", {
  tkras <- matrix(c(
    0, 0, 0,
    10, 20, 30
  ), ncol = 3, byrow = TRUE)
  c_ras <- c(1, 2, 3)

  ras <- tkras_to_ras(tkras, c_ras)

  expect_equal(ras[1, ], c(1, 2, 3))
  expect_equal(ras[2, ], c(11, 22, 33))
})

test_that("ras_to_tkras subtracts c_ras offset for vector", {
  ras <- c(0.5, -2.3, 1.1)
  c_ras <- c(0.5, -2.3, 1.1)

  tkras <- ras_to_tkras(ras, c_ras)

  expect_equal(tkras, c(0, 0, 0))
})

test_that("ras_to_tkras subtracts c_ras offset for matrix", {
  ras <- matrix(c(
    1, 2, 3,
    11, 22, 33
  ), ncol = 3, byrow = TRUE)
  c_ras <- c(1, 2, 3)

  tkras <- ras_to_tkras(ras, c_ras)

  expect_equal(tkras[1, ], c(0, 0, 0))
  expect_equal(tkras[2, ], c(10, 20, 30))
})

test_that("tkras_to_ras round-trip preserves coordinates", {
  original <- c(10, 20, 30)
  c_ras <- c(0.5, -2.3, 1.1)

  # RAS -> tkRAS -> RAS
  tkras <- ras_to_tkras(original, c_ras)
  back <- tkras_to_ras(tkras, c_ras)

  expect_equal(back, original)
})

test_that("tkras_to_ras errors on wrong c_ras length", {
  expect_error(tkras_to_ras(c(0, 0, 0), c(1, 2)), "length 3")
})

test_that("ras_to_tkras errors on wrong c_ras length", {
  expect_error(ras_to_tkras(c(0, 0, 0), c(1, 2)), "length 3")
})

test_that("tkras_to_ras errors on non-vector/matrix input", {
  # The function validates tkras first, so expect the validation error
  expect_error(tkras_to_ras(list(1, 2, 3), c(0, 0, 0)))
})

# =============================================================================
# APPLY_AFFINE TESTS
# =============================================================================

test_that("apply_affine with identity preserves coordinates", {
  coords <- c(10, 20, 30)
  result <- apply_affine(coords, diag(4))

  expect_equal(result, coords)
})

test_that("apply_affine applies translation", {
  coords <- c(0, 0, 0)
  aff <- diag(4)
  aff[1:3, 4] <- c(5, 10, 15)

  result <- apply_affine(coords, aff)

  expect_equal(result, c(5, 10, 15))
})

test_that("apply_affine applies scaling", {
  coords <- c(1, 2, 3)
  aff <- diag(c(2, 3, 4, 1))

  result <- apply_affine(coords, aff)

  expect_equal(result, c(2, 6, 12))
})

test_that("apply_affine applies rotation", {
  # 90 degree rotation around Z axis
  # R fills matrices by column, so we need to be careful about order
  theta <- pi / 2
  aff <- diag(4)
  # Standard rotation matrix: [[cos, -sin], [sin, cos]]
  # In R's column-major order: c(cos, sin, -sin, cos)
  aff[1:2, 1:2] <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)

  coords <- c(1, 0, 0)
  result <- apply_affine(coords, aff)

  # (1, 0, 0) rotated 90 degrees CCW around Z -> (0, 1, 0)
  expect_equal(result[1], 0, tolerance = 1e-10)
  expect_equal(result[2], 1, tolerance = 1e-10)
  expect_equal(result[3], 0, tolerance = 1e-10)
})

test_that("apply_affine works with matrix input", {
  coords <- matrix(c(
    0, 0, 0,
    1, 1, 1,
    2, 2, 2
  ), ncol = 3, byrow = TRUE)

  aff <- diag(4)
  aff[1:3, 4] <- c(10, 20, 30)

  result <- apply_affine(coords, aff)

  expect_equal(result[1, ], c(10, 20, 30))
  expect_equal(result[2, ], c(11, 21, 31))
  expect_equal(result[3, ], c(12, 22, 32))
})

test_that("apply_affine errors on non-4x4 affine", {
  expect_error(apply_affine(c(1, 2, 3), diag(3)), "4x4")
})

test_that("apply_affine errors on wrong vector length", {
  expect_error(apply_affine(c(1, 2), diag(4)), "length 3")
})

test_that("apply_affine errors on wrong matrix columns", {
  expect_error(apply_affine(matrix(1:6, ncol = 2), diag(4)), "3 columns")
})

# =============================================================================
# COMPOSE_AFFINES TESTS
# =============================================================================

test_that("compose_affines computes B * A", {
  A <- diag(4)
  A[1:3, 4] <- c(1, 0, 0)  # Translate by (1, 0, 0)

  B <- diag(4)
  B[1:3, 4] <- c(0, 2, 0)  # Translate by (0, 2, 0)

  composed <- compose_affines(A, B)

  # B(A(x)) = x + (1,0,0) + (0,2,0) = x + (1,2,0)
  expect_equal(composed[1:3, 4], c(1, 2, 0))
})

test_that("compose_affines with identity yields other matrix", {
  A <- diag(4)
  A[1:3, 4] <- c(5, 10, 15)

  # A * I = A
  expect_equal(compose_affines(diag(4), A), A)

  # I * A = A
  expect_equal(compose_affines(A, diag(4)), A)
})

test_that("compose_affines scale then translate", {
  scale <- diag(c(2, 2, 2, 1))
  translate <- diag(4)
  translate[1:3, 4] <- c(10, 10, 10)

  # Scale first, then translate
  composed <- compose_affines(scale, translate)

  # Point (1,1,1): scale -> (2,2,2), translate -> (12,12,12)
  result <- apply_affine(c(1, 1, 1), composed)
  expect_equal(result, c(12, 12, 12))
})

test_that("compose_affines validates inputs", {
  expect_error(compose_affines(diag(3), diag(4)), "4x4")
  expect_error(compose_affines(diag(4), diag(3)), "4x4")
})

# =============================================================================
# INVERT_AFFINE TESTS
# =============================================================================

test_that("invert_affine of identity is identity", {
  inv <- invert_affine(diag(4))
  expect_equal(inv, diag(4))
})

test_that("invert_affine of translation negates translation", {
  aff <- diag(4)
  aff[1:3, 4] <- c(5, 10, 15)

  inv <- invert_affine(aff)

  expect_equal(inv[1:3, 4], c(-5, -10, -15))
})

test_that("invert_affine of scale inverts scale factors", {
  aff <- diag(c(2, 4, 8, 1))

  inv <- invert_affine(aff)

  expect_equal(diag(inv), c(0.5, 0.25, 0.125, 1))
})

test_that("invert_affine round-trip yields identity transformation", {
  # Complex affine: rotation + scale + translation
  theta <- pi / 4
  aff <- diag(c(2, 2, 2, 1))  # Scale
  aff[1:2, 1:2] <- 2 * matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  aff[1:3, 4] <- c(10, 20, 30)  # Translation

  inv <- invert_affine(aff)

  # A * A^-1 should be identity
  product <- aff %*% inv
  expect_equal(product, diag(4), tolerance = 1e-10)
})

test_that("invert_affine validates input", {
  expect_error(invert_affine(diag(3)), "4x4")
})

# =============================================================================
# VOXELS_TO_WORLD / WORLD_TO_VOXELS TESTS
# =============================================================================

test_that("voxels_to_world with identity returns voxel coordinates", {
  voxels <- matrix(c(0, 0, 0, 1, 2, 3), ncol = 3, byrow = TRUE)
  world <- voxels_to_world(voxels, diag(4))

  expect_equal(world, voxels)
})

test_that("voxels_to_world applies affine correctly", {
  voxels <- matrix(c(0, 0, 0), ncol = 3)

  # 2mm voxels with (10, 20, 30) origin
  aff <- diag(c(2, 2, 2, 1))
  aff[1:3, 4] <- c(10, 20, 30)

  world <- voxels_to_world(voxels, aff)

  expect_equal(world[1, ], c(10, 20, 30))
})

test_that("world_to_voxels inverts voxels_to_world", {
  original_voxels <- matrix(c(5, 10, 15, 0, 0, 0), ncol = 3, byrow = TRUE)

  aff <- diag(c(2, 2, 2, 1))
  aff[1:3, 4] <- c(-10, -20, -30)

  world <- voxels_to_world(original_voxels, aff)
  back <- world_to_voxels(world, aff)

  expect_equal(back, original_voxels, tolerance = 1e-10)
})

# =============================================================================
# COORDINATE CONVERSION INTEGRATION TESTS
# =============================================================================

test_that("LPS coordinates can be converted and used with affine", {
  # Create a point in LPS space
  lps_point <- c(-45, -30, 20)

  # Convert to RAS
  ras_point <- lps_to_ras(lps_point)

  # Apply an affine in RAS space
  aff <- diag(4)
  aff[1:3, 4] <- c(10, 10, 10)
  transformed <- apply_affine(ras_point, aff)

  # Convert back to LPS
  lps_result <- ras_to_lps(transformed)

  # Should be original + offset in LPS terms
  expect_equal(lps_result, c(-55, -40, 30))
})

test_that("tkRAS conversion preserves relative positions", {
  # Create two points in tkRAS
  p1_tkras <- c(0, 0, 0)
  p2_tkras <- c(10, 10, 10)
  c_ras <- c(5, -3, 2)

  # Convert to scanner RAS
  p1_ras <- tkras_to_ras(p1_tkras, c_ras)
  p2_ras <- tkras_to_ras(p2_tkras, c_ras)

  # Distance should be preserved
  dist_tkras <- sqrt(sum((p2_tkras - p1_tkras)^2))
  dist_ras <- sqrt(sum((p2_ras - p1_ras)^2))

  expect_equal(dist_tkras, dist_ras)
})

test_that("coordinate systems handle multiple points correctly", {
  # Matrix of LPS points
  lps_points <- matrix(c(
    -10, -20, 30,
    0, 0, 0,
    50, -50, 100
  ), ncol = 3, byrow = TRUE)

  # Convert to RAS and back
  ras_points <- lps_to_ras(lps_points)
  back_to_lps <- ras_to_lps(ras_points)

  expect_equal(back_to_lps, lps_points)
})

# =============================================================================
# RAI <-> RAS CONVERSION TESTS (AFNI)
# =============================================================================
# RAI (AFNI): +X=Right, +Y=Anterior, +Z=Inferior
# RAS (NIfTI): +X=Right, +Y=Anterior, +Z=Superior
# X and Y are identical; only Z differs.

test_that("RAI to RAS conversion only flips Z", {
  # In RAI, +Z points Inferior. In RAS, +Z points Superior.
  # So RAI Z = -RAS Z
  rai <- c(10, 20, 30)  # Point in RAI

  # Convert to RAS: only Z should be negated
  ras <- rai
  ras[3] <- -rai[3]

  expect_equal(ras, c(10, 20, -30))
})

test_that("RAS to RAI conversion only flips Z", {
  ras <- c(10, 20, 30)  # Point in RAS

  # Convert to RAI: only Z should be negated
  rai <- ras
  rai[3] <- -ras[3]

  expect_equal(rai, c(10, 20, -30))
})

test_that("RAI round-trip preserves coordinates", {
  original <- c(10, 20, 30)

  # RAS -> RAI -> RAS (flip Z twice)
  rai <- original
  rai[3] <- -original[3]
  back <- rai
  back[3] <- -rai[3]

  expect_equal(back, original)
})

test_that("AFNI coordinate system differs from LPS/DICOM", {
  # This test documents that RAI != LPS (DICOM)
  # RAI: +X=Right, +Y=Anterior, +Z=Inferior
  # LPS: +X=Left, +Y=Posterior, +Z=Superior

  ras <- c(10, 20, 30)

  # RAI conversion: negate Z only
  rai <- c(ras[1], ras[2], -ras[3])

 # LPS conversion: negate X and Y
  lps <- c(-ras[1], -ras[2], ras[3])

  # They should NOT be equal
  expect_false(all(rai == lps))

  # Verify the expected values
  expect_equal(rai, c(10, 20, -30))
  expect_equal(lps, c(-10, -20, 30))
})
