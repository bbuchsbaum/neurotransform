test_that("compose and transform agree with sequential application", {
  A <- diag(4); A[1, 4] <- 1  # translate +1 in x
  B <- diag(4); B[2, 4] <- 2  # translate +2 in y

  f <- Affine3DMorphism("a", "b", A)
  g <- Affine3DMorphism("b", "c", B)

  fg <- compose(f, g)

  pts <- matrix(c(0, 0, 0,
                  1, 1, 1), ncol = 3, byrow = TRUE)

  seq_val <- transform(f, transform(g, pts))
  comp_val <- transform(fg, pts)

  expect_equal(comp_val, seq_val)
})

test_that("transform_path applies pullback in correct order", {
  A <- diag(4); A[1, 4] <- 1
  B <- diag(4); B[3, 4] <- -1
  f <- Affine3DMorphism("a", "b", A)
  g <- Affine3DMorphism("b", "c", B)
  path <- list(f, g)

  pts <- matrix(c(0, 0, 0), ncol = 3)
  direct <- transform_path(path, pts)
  manual <- transform(f, transform(g, pts))

  expect_equal(direct, manual)
})

test_that("invert on affine morphism returns identity when composed", {
  A <- diag(4); A[1:3, 4] <- c(1, -2, 3)
  f <- Affine3DMorphism("a", "b", A)
  invf <- invert(f)

  pts <- matrix(c(2, 2, 2), ncol = 3)
  back <- transform(f, transform(invf, pts))
  expect_equal(back, pts)
})

test_that("grid_spec and grid_coords generate correct centres", {
  grid <- grid_spec(dims = c(2L, 2L, 2L), affine = diag(4))
  coords <- grid_coords(grid)
  expect_equal(nrow(coords), prod(grid@dims))
  expect_true(any(coords[,1] == 0))
  expect_true(any(coords[,1] == 1))
})

test_that("transform on VolToSurfMorphism errors with guidance", {
  v2s <- VolToSurfMorphism("vol", "surf", method = "trilinear")

  expect_error(
    transform(v2s, matrix(c(1, 2, 3), ncol = 3)),
    "sample_volume_on_surface\\(\\), resample\\(\\), or adjoint\\(\\)"
  )
})
