test_that("ribbon sampling averages along the line segment", {
  vol <- array(c(0, 1, 2), dim = c(3L, 1L, 1L))
  inner <- matrix(c(0, 0, 0), ncol = 3)
  outer <- matrix(c(2, 0, 0), ncol = 3)

  val <- neurotransform:::cpp_ribbon_sample_volume(
    data = vol,
    inner = inner,
    outer = outer,
    world_to_vox = diag(4),
    n_steps = 2L,
    method = "linear"
  )

  expect_equal(val, 1) # mean of samples at x = 0,1,2
})

test_that("ribbon weights respect mask and normalize", {
  inner <- matrix(c(0, 0, 0), ncol = 3)
  outer <- matrix(c(2, 0, 0), ncol = 3)
  dims <- c(3L, 1L, 1L)
  mask <- c(TRUE, FALSE, TRUE)

  w <- neurotransform:::cpp_ribbon_weights(inner, outer, diag(4), dims, n_steps = 2L, mask = mask)

  expect_length(w$rows, 2)
  expect_equal(w$rows, c(1L, 1L))
  expect_setequal(w$cols, c(1L, 3L))
  expect_equal(sum(w$vals), 1, tolerance = 1e-8)
  expect_equal(sort(w$vals), c(0.5, 0.5))
})

test_that("trilinear weights cover 8 neighbors and sum to one", {
  dims <- c(2L, 2L, 2L)
  coords <- matrix(c(0.2, 0.4, 0.6), ncol = 3)
  w <- neurotransform:::cpp_trilinear_weights(coords, dims)

  expect_equal(length(w$vals), 8)
  expect_equal(sum(w$vals), 1, tolerance = 1e-8)

  expected_cols <- c(
    1L, 2L, 3L, 4L, # z = 0 layer
    5L, 6L, 7L, 8L  # z = 1 layer
  )
  expect_setequal(w$cols, expected_cols)
})
