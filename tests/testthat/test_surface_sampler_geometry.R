test_that("surface samplers rank geometric distance rather than centroids", {
  vertices <- rbind(c(0, 0, 0), c(100, 0, 0), c(0, 100, 0),
                    c(0, 0, 1), c(1, 0, 1), c(0, 1, 1))
  faces <- rbind(1:3, 4:6)
  query <- matrix(c(0.2, 0.2, 0.1), 1)
  data <- c(0, 0, 0, 1, 1, 1)
  expect_equal(neurotransform:::cpp_barycentric_sample(query, vertices, faces, data), 0)
  sampler <- surface_sampler(vertices, data, faces, method = "barycentric")
  expect_equal(sampler@evaluate(query), 0)
})

test_that("samplers preserve data shape, mesh indices and projection policy", {
  vertices <- rbind(c(-10, 0, 0), c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
  faces <- matrix(2:4, 1)
  # The stored zero-based face does not contain zero: its known index base must
  # not be guessed again when the object is passed into a sampler.
  mesh <- surface_mesh(vertices, faces)
  q <- rbind(c(.25, .25, .1), c(5, 0, 0))
  for (data in list(c(99, 1, 2, 3), matrix(c(99, 1, 2, 3), 4),
                   cbind(c(99, 1, 2, 3), c(77, 2, 4, 6)))) {
    sampler <- surface_sampler(mesh, data, method = "barycentric")
    got <- sampler@evaluate(q)
    expected <- if (is.matrix(data)) rbind(c(1.75, 3.5)[seq_len(ncol(data))], rep(NA_real_, ncol(data))) else c(1.75, NA_real_)
    expect_equal(got, expected)
    closest <- surface_sampler(mesh, data, method = "barycentric", projection = "closest")
    if (is.matrix(data)) expect_equal(closest@evaluate(q)[2, ], data[3, ])
    else expect_equal(closest@evaluate(q)[2], data[3])
  }
  expect_error(neurotransform:::cpp_barycentric_sample(q, vertices,
    matrix(c(0L, 1L, 2L), 1), 1:4), "one-based")
  expect_error(neurotransform:::cpp_barycentric_sample(q, vertices, faces, 1:3), "per vertex")
})

test_that("closest samplers pass the independent spherical impulse oracle", {
  skip_if_not_installed("jsonlite")
  fixture <- jsonlite::read_json(system.file("extdata", "barycentric_oracle", "oracle.json",
                                            package = "neurotransform"), simplifyVector = TRUE)
  for (i in seq_len(nrow(fixture$cases))) {
    case <- fixture$cases[i, ]
    mesh <- surface_mesh(case$vertices[[1]], case$faces[[1]] + 1L)
    sampler <- surface_sampler(mesh, diag(nrow(mesh@coords)), method = "barycentric",
                                projection = "closest")
    got <- sampler@evaluate(case$query[[1]])
    expect_true(max(abs(got - case$weights[[1]])) < fixture$tolerance)
  }
})
