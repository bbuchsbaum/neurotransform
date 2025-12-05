test_that("triplet assembly builds sparse matrix and rejects out-of-bounds", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  mat <- neurotransform:::cpp_triplets_to_dgC(
    i = c(1L, 2L),
    j = c(1L, 2L),
    x = c(1, 3),
    nrow = 2L,
    ncol = 2L
  )
  expect_s4_class(mat, "dgCMatrix")
  expect_equal(as.matrix(mat), matrix(c(1, 0, 0, 3), nrow = 2, byrow = TRUE))
  oob <- neurotransform:::cpp_triplets_to_dgC(i = 3L, j = 1L, x = 1, nrow = 2L, ncol = 2L)
  expect_true(all(as.matrix(oob) == 0))
  expect_error(neurotransform:::cpp_triplets_to_dgC(
    i = c(1L, 1L), j = c(1L, 1L), x = c(1, 2), nrow = 2L, ncol = 2L
  ))
})

test_that("resample applies jacobian modulation factors", {
  vol <- array(1, dim = c(2L, 2L, 2L))
  sampler <- volume_sampler(vol, method = "nearest", domain = "a")
  morph <- Affine3DMorphism("a", "a", diag(c(0.5, 0.5, 0.5, 1)))
  grid <- grid_from_data(vol)

  none <- resample(sampler, morph, grid, modulate = "none")
  jac <- resample(sampler, morph, grid, modulate = "jacobian")
  sjac <- resample(sampler, morph, grid, modulate = "sqrt_jacobian")

  expect_true(all(none == 1))
  expect_equal(jac, rep(0.125, length(jac)), tolerance = 1e-8)
  expect_equal(sjac, rep(sqrt(0.125), length(sjac)), tolerance = 1e-8)
})

test_that("warp loader registry returns expected loaders", {
  default_loader <- get_loader()
  expect_true(is.function(default_loader))
  expect_identical(default_loader, neurotransform:::load_warp_neuroim2)
  expect_error(get_loader("nonexistent"))

  called <- FALSE
  custom_loader <- function(path) {
    called <<- TRUE
    list(array = numeric(), dim = integer(3), world_to_vox = diag(4))
  }
  register_loader("test_loader", custom_loader)
  loader <- get_loader("test_loader")
  res <- loader("ignored.nii")

  expect_true(called)
  expect_equal(res$dim, integer(3))
})
