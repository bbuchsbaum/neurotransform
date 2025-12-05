# Test morphism invert() methods comprehensively

test_that("invert(IdentityMorphism) returns self", {
  id <- IdentityMorphism("test")
  inv <- invert(id)

  expect_identical(inv, id)
})

test_that("invert(IdentityMorphism) maintains domain", {
  id <- IdentityMorphism("my_domain")
  inv <- invert(id)

  expect_equal(source_of(inv), "my_domain")
  expect_equal(target_of(inv), "my_domain")
})

test_that("invert(Affine3DMorphism) swaps source and target", {
  aff <- Affine3DMorphism("domain_a", "domain_b", diag(4))
  inv <- invert(aff)

  expect_equal(source_of(inv), "domain_b")
  expect_equal(target_of(inv), "domain_a")
})

test_that("invert(Affine3DMorphism) computes correct inverse matrix", {
  # Translation matrix
  A <- diag(4)
  A[1:3, 4] <- c(10, 20, 30)

  aff <- Affine3DMorphism("a", "b", A)
  inv <- invert(aff)

  # Inverse of translation (10,20,30) is translation (-10,-20,-30)
  expect_equal(inv@matrix[1, 4], -10, tolerance = 1e-10)
  expect_equal(inv@matrix[2, 4], -20, tolerance = 1e-10)
  expect_equal(inv@matrix[3, 4], -30, tolerance = 1e-10)
})

test_that("invert(Affine3DMorphism) preserves cost", {
  aff <- Affine3DMorphism("a", "b", diag(4), cost = 2.5)
  inv <- invert(aff)

  expect_equal(inv@cost, 2.5)
})

test_that("invert(Affine3DMorphism) preserves method_tag", {
  aff <- Affine3DMorphism("a", "b", diag(4), method_tag = "functional")
  inv <- invert(aff)

  expect_equal(inv@method_tag, "functional")
})

test_that("invert(Affine3DMorphism) double inversion returns equivalent matrix", {
  # Scale by 2
  A <- diag(c(2, 2, 2, 1))

  aff <- Affine3DMorphism("a", "b", A)
  inv_inv <- invert(invert(aff))

  # Matrix should be equivalent to original
  expect_equal(inv_inv@matrix, aff@matrix, tolerance = 1e-10)
})

test_that("invert(Affine3DMorphism) compose with inverse yields identity transform", {
  # Complex affine: scale + translation
  A <- diag(4)
  A[1, 1] <- 2
  A[2, 2] <- 0.5
  A[1:3, 4] <- c(5, -3, 7)

  aff <- Affine3DMorphism("a", "b", A)
  inv <- invert(aff)

  coords <- matrix(c(
    0, 0, 0,
    10, 20, 30,
    -5, -5, -5,
    1.5, 2.5, 3.5
  ), ncol = 3, byrow = TRUE)

  # Apply aff then inv should return original coords
  result <- transform(inv, transform(aff, coords))
  expect_equal(result, coords, tolerance = 1e-10)
})

test_that("invert(Affine3DMorphism) handles rotation matrix", {
  # 90 degree rotation around Z axis
  theta <- pi / 2
  A <- diag(4)
  A[1:2, 1:2] <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)

  aff <- Affine3DMorphism("a", "b", A)
  inv <- invert(aff)

  # Inverse of 90 degree rotation is -90 degree rotation
  coords <- matrix(c(1, 0, 0), ncol = 3)
  rotated <- transform(aff, coords)
  back <- transform(inv, rotated)

  expect_equal(back, coords, tolerance = 1e-10)
})

test_that("invert(Warp3DMorphism) errors without inverse_path", {
  w <- Warp3DMorphism("a", "b", "fwd.nii", warp_type = "ants")

  expect_error(invert(w), "inverse")
})

test_that("invert(Warp3DMorphism) swaps paths when inverse available", {
  w <- Warp3DMorphism("a", "b", "fwd.nii", warp_type = "ants", inverse_path = "inv.nii")
  inv <- invert(w)

  # Source/target swapped
  expect_equal(source_of(inv), "b")
  expect_equal(target_of(inv), "a")

  # Paths swapped
  expect_equal(inv@warp_path, "inv.nii")
  expect_equal(inv@inverse_path, "fwd.nii")
})

test_that("invert(Warp3DMorphism) preserves cost", {
  w <- Warp3DMorphism("a", "b", "fwd.nii", warp_type = "ants",
                      inverse_path = "inv.nii", cost = 3.0)
  inv <- invert(w)

  expect_equal(inv@cost, 3.0)
})

test_that("invert(Warp3DMorphism) preserves method_tag", {
  w <- Warp3DMorphism("a", "b", "fwd.nii", warp_type = "ants",
                      inverse_path = "inv.nii", method_tag = "nonlinear_reg")
  inv <- invert(w)

  expect_equal(inv@method_tag, "nonlinear_reg")
})

test_that("invert(Warp3DMorphism) preserves warp_type", {
  w <- Warp3DMorphism("a", "b", "fwd.nii", warp_type = "fsl",
                      inverse_path = "inv.nii")
  inv <- invert(w)

  expect_equal(inv@warp_type, "fsl")
})

test_that("invert(Warp3DMorphism) double inversion restores paths", {
  w <- Warp3DMorphism("a", "b", "fwd.nii", warp_type = "ants",
                      inverse_path = "inv.nii")
  inv_inv <- invert(invert(w))

  expect_equal(inv_inv@warp_path, "fwd.nii")
  expect_equal(inv_inv@inverse_path, "inv.nii")
  expect_equal(source_of(inv_inv), "a")
  expect_equal(target_of(inv_inv), "b")
})

test_that("invert(VolToSurfMorphism) errors (adjoint only)", {
  v2s <- VolToSurfMorphism("vol", "surf", method = "trilinear")

  expect_error(invert(v2s), "inverse")
})

test_that("invert base Morphism errors for non-invertible types", {
  # The base class invert method should error
  v2s <- VolToSurfMorphism("a", "b", method = "trilinear")

  expect_error(invert(v2s), "inverse")
})

test_that("invert preserves hash uniqueness", {
  aff <- Affine3DMorphism("a", "b", diag(4))
  inv <- invert(aff)

  # Different morphisms should have different hashes
  expect_false(morphism_hash(aff) == morphism_hash(inv))
})

test_that("invert creates new object (not in-place modification)", {
  A <- diag(4)
  A[1, 4] <- 10

  aff <- Affine3DMorphism("a", "b", A)
  original_matrix <- aff@matrix

  inv <- invert(aff)

  # Original should be unchanged
  expect_equal(aff@matrix, original_matrix)
  expect_equal(aff@source, "a")
  expect_equal(aff@target, "b")
})

test_that("invert handles scaling affine correctly", {
  # Non-uniform scale
  A <- diag(c(2, 3, 4, 1))

  aff <- Affine3DMorphism("a", "b", A)
  inv <- invert(aff)

  coords <- matrix(c(10, 10, 10), ncol = 3)
  scaled <- transform(aff, coords)
  unscaled <- transform(inv, scaled)

  expect_equal(unscaled, coords, tolerance = 1e-10)
})

test_that("invert handles shear affine correctly", {
  # Shear matrix
  A <- diag(4)
  A[1, 2] <- 0.5  # X sheared by Y

  aff <- Affine3DMorphism("a", "b", A)
  inv <- invert(aff)

  coords <- matrix(c(5, 10, 15), ncol = 3)
  sheared <- transform(aff, coords)
  unsheared <- transform(inv, sheared)

  expect_equal(unsheared, coords, tolerance = 1e-10)
})
