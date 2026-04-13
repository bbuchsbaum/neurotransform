test_that("SurfToSurfMorphism nearest maps by vertex correspondence", {
  src_v <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)

  # Target is translated in space, but correspondence is by vertex index.
  tgt_v <- src_v + 10
  faces <- matrix(c(1L, 2L, 3L), nrow = 1L)

  m <- SurfToSurfMorphism(
    source = "src",
    target = "tgt",
    mapping = "nearest",
    source_vertices = src_v,
    target_vertices = tgt_v,
    faces = faces
  )

  q <- matrix(c(10.9, 10.1, 10), ncol = 3)
  mapped <- transform(m, q)
  # Should map to source vertex 2: (1,0,0)
  expect_equal(mapped, matrix(c(1, 0, 0), ncol = 3))
})

test_that("SurfToSurfMorphism barycentric maps within matching topology", {
  src_v <- matrix(c(
    0, 0, 0,
    2, 0, 0,
    0, 2, 0
  ), ncol = 3, byrow = TRUE)

  # Target triangle is rotated 90 degrees around Z, then shifted.
  rot <- matrix(c(
    0, -1, 0,
    1,  0, 0,
    0,  0, 1
  ), 3, 3, byrow = TRUE)
  tgt_v <- (src_v %*% t(rot)) + 5

  faces <- matrix(c(1L, 2L, 3L), nrow = 1L)

  m <- SurfToSurfMorphism(
    source = "src",
    target = "tgt",
    mapping = "barycentric",
    source_vertices = src_v,
    target_vertices = tgt_v,
    faces = faces
  )

  # Query at centroid of target triangle
  q <- matrix(colMeans(tgt_v), ncol = 3)
  mapped <- transform(m, q)

  # Centroid should map to centroid of source triangle
  expect_equal(mapped, matrix(colMeans(src_v), ncol = 3), tolerance = 1e-8)
})

test_that("SurfToSurfMorphism barycentric marks uncovered points as NA", {
  src_v <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  tgt_v <- src_v
  faces <- matrix(c(1L, 2L, 3L), nrow = 1L)

  m <- SurfToSurfMorphism(
    source = "src",
    target = "tgt",
    mapping = "barycentric",
    source_vertices = src_v,
    target_vertices = tgt_v,
    faces = faces
  )

  q <- matrix(c(5, 5, 5), ncol = 3)
  mapped <- transform(m, q)

  expect_true(all(is.na(mapped)))
})
