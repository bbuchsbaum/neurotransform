test_that("volume sampling matches linear function for multiple methods", {
  dim <- c(3L, 3L, 3L)
  vol <- array(0, dim = dim)
  for (x in 0:(dim[1] - 1)) {
    for (y in 0:(dim[2] - 1)) {
      for (z in 0:(dim[3] - 1)) {
        vol[x + 1, y + 1, z + 1] <- x + y + z
      }
    }
  }

  coords <- matrix(c(0.2, 0.1, 0.4), ncol = 3)
  expected <- sum(coords)

  lin <- neurotransform:::cpp_sample_volume(vol, coords, diag(4), "linear", outside = -99)
  cub <- neurotransform:::cpp_sample_volume(vol, coords, diag(4), "cubic", outside = -99)
  near <- neurotransform:::cpp_sample_volume(vol, coords, diag(4), "nearest", outside = -99)

  expect_equal(lin, expected, tolerance = 1e-6)
  expect_equal(cub, expected, tolerance = 1e-6)
  expect_equal(near, 0) # nearest to (0,0,0)

  oob <- neurotransform:::cpp_sample_volume(vol, matrix(c(5, 0, 0), ncol = 3), diag(4), "linear", outside = -99)
  expect_identical(oob, -99)
})

test_that("nearest vertex picks the closest index", {
  verts <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  idx <- neurotransform:::cpp_nearest_vertex(matrix(c(0.9, 0.1, 0), ncol = 3), verts)
  expect_identical(as.integer(idx), 2L)
})

test_that("barycentric sampling matches expected interpolation", {
  verts <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  faces <- matrix(c(1L, 2L, 3L), nrow = 1L)
  storage.mode(faces) <- "integer"
  data <- c(10, 20, 30)

  centroid <- matrix(c(1 / 3, 1 / 3, 0), ncol = 3)
  at_vertex <- matrix(c(1, 0, 0), ncol = 3)

  cent_val <- neurotransform:::cpp_barycentric_sample(centroid, verts, faces, data)
  vert_val <- neurotransform:::cpp_barycentric_sample(at_vertex, verts, faces, data)

  expect_equal(cent_val, 20, tolerance = 1e-8)
  expect_equal(vert_val, 20, tolerance = 1e-8)
})
