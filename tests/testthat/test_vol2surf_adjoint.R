test_that("backproject_surface_to_volume matches adjoint() for VolToSurfMorphism", {
  grid <- grid_spec(dims = c(2L, 2L, 2L), affine = diag(4))

  coords <- matrix(c(
    0.25, 0.25, 0.25,
    0.75, 0.75, 0.75
  ), ncol = 3, byrow = TRUE)

  v2s <- VolToSurfMorphism("vol", "surf", method = "trilinear", mid_coords = coords)
  y <- c(1, 2)

  bp1 <- backproject_surface_to_volume(y, v2s, grid, surface_coords = coords, method = "linear")

  adj <- adjoint(v2s, grid = grid, surface_coords = coords, method = "linear")
  bp2 <- adj(y)

  expect_equal(bp2, bp1, tolerance = 1e-12)
})

test_that("VolToSurf adjoint satisfies inner-product identity for linear sampling", {
  grid <- grid_spec(dims = c(2L, 2L, 2L), affine = diag(4))
  coords <- matrix(c(
    0.25, 0.25, 0.25,
    0.75, 0.75, 0.75
  ), ncol = 3, byrow = TRUE)

  v2s <- VolToSurfMorphism("vol", "surf", method = "trilinear", mid_coords = coords)

  x <- array(runif(8), dim = c(2, 2, 2))
  y <- runif(2)

  fx <- sample_volume_on_surface(x, v2s, surface_coords = coords, method = "linear")
  lhs <- sum(fx * y)

  aty <- backproject_surface_to_volume(y, v2s, grid, surface_coords = coords, method = "linear")
  rhs <- sum(x * aty)

  expect_equal(lhs, rhs, tolerance = 1e-6)
})

