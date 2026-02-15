test_that("surface_resampling_plan nearest supports differing vertex counts", {
  moving <- surface_mesh(
    matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 1, 0
    ), ncol = 3, byrow = TRUE),
    faces = matrix(c(1, 2, 3), ncol = 3, byrow = TRUE)
  )

  reference <- surface_mesh(
    matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0.9, 0.1, 0
    ), ncol = 3, byrow = TRUE)
  )

  plan <- surface_resampling_plan(reference, moving, method = "nearest", spherical = FALSE)
  y <- apply_surface_resampling(plan, c(10, 20, 30), normalize = "none")
  expect_equal(y, c(10, 20, 30, 20))
})

test_that("surface_resampling_plan barycentric interpolates moving vertex data", {
  moving <- surface_mesh(
    matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 1, 0
    ), ncol = 3, byrow = TRUE),
    faces = matrix(c(1, 2, 3), ncol = 3, byrow = TRUE)
  )

  reference <- surface_mesh(
    matrix(c(
      1/3, 1/3, 0,
      0, 0, 0
    ), ncol = 3, byrow = TRUE)
  )

  plan <- surface_resampling_plan(reference, moving, method = "barycentric", spherical = FALSE)
  y <- apply_surface_resampling(plan, c(1, 2, 3), normalize = "element")
  expect_equal(y[1], 2, tolerance = 1e-6)
  expect_equal(y[2], 1, tolerance = 1e-6)
})

test_that("apply_surface_resampling sum normalization preserves total mass", {
  moving <- surface_mesh(
    matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 1, 0
    ), ncol = 3, byrow = TRUE)
  )
  reference <- surface_mesh(
    matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0.9, 0.1, 0
    ), ncol = 3, byrow = TRUE)
  )

  plan <- surface_resampling_plan(reference, moving, method = "nearest", spherical = FALSE)
  x <- c(1, 1, 1)
  y_none <- apply_surface_resampling(plan, x, normalize = "none")
  y_sum <- apply_surface_resampling(plan, x, normalize = "sum")

  expect_equal(sum(y_none), 4, tolerance = 1e-8)
  expect_equal(sum(y_sum), sum(x), tolerance = 1e-8)
})
