octahedron <- function() {
  list(vertices = rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 1, 0),
                        c(0, -1, 0), c(0, 0, 1), c(0, 0, -1)),
       faces = rbind(c(1, 3, 5), c(3, 2, 5), c(2, 4, 5), c(4, 1, 5),
                     c(3, 1, 6), c(2, 3, 6), c(4, 2, 6), c(1, 4, 6)))
}

bary_matrix <- function(query, vertices, faces, closest = TRUE) {
  w <- neurotransform:::cpp_barycentric_weights(query, vertices, faces - 1L,
                                                closest = closest)
  out <- matrix(0, nrow(query), nrow(vertices))
  out[cbind(w$rows, w$cols)] <- w$vals
  out
}

test_that("octahedral face selection excludes the opposite side", {
  mesh <- octahedron()
  query <- as.matrix(expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1))) / sqrt(3)
  expected <- cbind(query[, 1] > 0, query[, 1] < 0,
                    query[, 2] > 0, query[, 2] < 0,
                    query[, 3] > 0, query[, 3] < 0) / 3
  for (order in list(1:8, 8:1)) {
    # Exercise the original three-argument compiled API as well as the public
    # spherical plan. Inspect every weight, not only a constant-field output.
    expect_equal(bary_matrix(query, mesh$vertices, mesh$faces[order, ]),
                 expected, tolerance = 1e-12, ignore_attr = TRUE)
    plan <- surface_resampling_plan(surface_mesh(query),
      surface_mesh(mesh$vertices, mesh$faces[order, ]))
    expect_equal(apply_surface_resampling(plan, mesh$vertices[, 3], normalize = "none"),
                 sign(query[, 3]) / 3, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(apply_surface_resampling(plan, diag(6), normalize = "none"),
                 expected, tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("closest spherical weights agree with an independent convex oracle", {
  mesh <- octahedron()
  set.seed(20926)
  query <- rbind(mesh$vertices, c(1, 1, 0), c(1, 0, 1), c(0, -1, 1),
                 c(0.97, 0.2, 0.05), c(1, 1, 1e-11), c(1, 1, -1e-11),
                 matrix(rnorm(600), ncol = 3))
  query <- query / sqrt(rowSums(query^2))
  # The solid octahedron is the unit L1 ball. Its Euclidean projection is
  # soft-thresholding with a threshold from sorted coordinates. This oracle
  # uses neither triangle planes nor the production edge-projection algorithm.
  expected <- t(apply(query, 1, function(q) {
    a <- sort(abs(q), decreasing = TRUE)
    thresholds <- (cumsum(a) - 1) / seq_along(a)
    rho <- max(which(a > thresholds))
    p <- pmax(abs(q) - thresholds[rho], 0)
    c(p[1] * (q[1] > 0), p[1] * (q[1] < 0),
      p[2] * (q[2] > 0), p[2] * (q[2] < 0),
      p[3] * (q[3] > 0), p[3] * (q[3] < 0))
  }))
  # The compiled API defaults to closest mesh projection, including boundaries.
  w <- neurotransform:::cpp_barycentric_weights(query, mesh$vertices, mesh$faces - 1L)
  actual <- matrix(0, nrow(query), 6)
  actual[cbind(w$rows, w$cols)] <- w$vals
  expect_equal(actual, expected, tolerance = 1e-12)
  windings <- rbind(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
                    c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  for (i in 1:18) {
    faces <- mesh$faces[sample.int(8), windings[(i - 1) %% 6 + 1, ]]
    scale <- c(1e-6, 1, 100, 1e6)[(i - 1) %% 4 + 1]
    got <- bary_matrix(query * scale, mesh$vertices * scale, faces, closest = TRUE)
    expect_equal(got, expected, tolerance = 1e-12)
    expect_equal(rowSums(got), rep(1, nrow(query)), tolerance = 1e-12)
    expect_true(all(got >= 0))
  }
  # Vertex numbering and query ordering are independent of the geometry.
  perm <- sample.int(6)
  faces <- matrix(match(mesh$faces, perm), ncol = 3)
  got <- bary_matrix(query, mesh$vertices[perm, ], faces, closest = TRUE)
  expect_equal(got[, order(perm)], expected, tolerance = 1e-12)
  plan <- surface_resampling_plan(surface_mesh(query), surface_mesh(mesh$vertices, mesh$faces))
  expect_equal(apply_surface_resampling(plan, diag(6), normalize = "none"),
               expected, tolerance = 1e-12)
})

test_that("generic projection uses distance and retains unsupported queries", {
  vertices <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0),
                    c(0, 0, 10), c(1, 0, 10), c(0, 1, 10))
  faces <- rbind(1:3, 4:6)
  query <- rbind(c(0.25, 0.25, 0.1), c(0.25, 0.25, 9.9), c(5, 5, 5))
  expected <- rbind(c(0.5, 0.25, 0.25, 0, 0, 0),
                    c(0, 0, 0, 0.5, 0.25, 0.25), rep(0, 6))
  expect_equal(bary_matrix(query, vertices, faces, closest = FALSE), expected, tolerance = 1e-12)
  expect_equal(bary_matrix(query, vertices, faces[2:1, ], closest = FALSE), expected, tolerance = 1e-12)
  # Exact ties between disjoint faces have a stable vertex-ID tie-break.
  tie <- matrix(c(0.25, 0.25, 5), nrow = 1)
  expect_equal(bary_matrix(tie, vertices, faces),
               bary_matrix(tie, vertices, faces[2:1, ]))
})

test_that("degenerate faces cannot fabricate spherical support", {
  mesh <- octahedron()
  query <- matrix(c(1, 1, 1) / sqrt(3), nrow = 1)
  degenerate <- rbind(c(1, 1, 1), c(1, 2, 1), c(1, 2, 2))
  expect_equal(bary_matrix(query, mesh$vertices, degenerate, closest = TRUE),
               matrix(0, 1, 6))
  expect_equal(bary_matrix(query, mesh$vertices, rbind(degenerate, mesh$faces), TRUE),
               bary_matrix(query, mesh$vertices, mesh$faces, TRUE))
  expect_error(surface_resampling_plan(surface_mesh(query),
    surface_mesh(mesh$vertices, degenerate)), "without valid triangle support")
  # Small but nonzero support must survive triplet emission.
  flat <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
  expect_equal(bary_matrix(matrix(c(1e-11, 0.5, 0), 1), flat, matrix(1:3, 1)),
               matrix(c(0.5 - 1e-11, 1e-11, 0.5), 1), tolerance = 1e-14)
  expect_error(neurotransform:::cpp_barycentric_weights(query, mesh$vertices,
    matrix(c(-1L, 0L, 1L), 1)), "zero-based")
})
