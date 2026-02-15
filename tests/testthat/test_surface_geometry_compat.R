# Compatibility tests for neurosurf::SurfaceGeometry-like objects

mock_surface_geometry <- function(vertices, faces0 = NULL, surf_to_world = diag(4),
                                  hemi = "lh", label = "pial") {
  stopifnot(is.matrix(vertices), ncol(vertices) == 3)
  if (is.null(faces0)) {
    it <- matrix(integer(0), nrow = 3, ncol = 0)
  } else {
    stopifnot(is.matrix(faces0), ncol(faces0) == 3)
    it <- t(faces0 + 1L) # mesh$it is 1-based in rgl mesh objects
  }
  vb <- rbind(t(vertices), rep(1, nrow(vertices)))
  structure(
    list(
      mesh = list(vb = vb, it = it),
      surf_to_world = surf_to_world,
      hemi = hemi,
      label = label
    ),
    class = "SurfaceGeometry"
  )
}

test_that("surface_sampler accepts SurfaceGeometry-like objects and applies surf_to_world", {
  vertices <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  faces0 <- matrix(c(0, 1, 2), ncol = 3, byrow = TRUE)
  stw <- diag(4)
  stw[1:3, 4] <- c(10, 20, 30)
  geom <- mock_surface_geometry(vertices, faces0, surf_to_world = stw)

  sampler <- surface_sampler(geom, data = c(10, 20, 30), method = "nearest")
  expect_s4_class(sampler, "Sampler")
  expect_equal(sampler@bounds$min, c(10, 20, 30))
  expect_equal(sampler@bounds$max, c(11, 21, 30))
})

test_that("sample_volume_on_surface accepts SurfaceGeometry-like coordinates", {
  vol <- array(0, dim = c(3, 3, 3))
  vol[1, 1, 1] <- 1
  vol[2, 1, 1] <- 2
  vol[1, 2, 1] <- 3

  vertices <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  geom <- mock_surface_geometry(vertices, faces0 = matrix(c(0, 1, 2), ncol = 3))
  m <- VolToSurfMorphism("vol", "surf", method = "trilinear")

  vals <- sample_volume_on_surface(vol, m, surface_coords = geom, method = "nearest")
  expect_length(vals, 3)
})

test_that("SurfToSurfMorphism accepts SurfaceGeometry-like inputs and infers faces", {
  vertices <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  faces0 <- matrix(c(0, 1, 2), ncol = 3, byrow = TRUE)

  src_geom <- mock_surface_geometry(vertices, faces0)
  tgt_geom <- mock_surface_geometry(vertices, faces0)

  s2s <- SurfToSurfMorphism(
    source = "src",
    target = "tgt",
    mapping = "barycentric",
    source_vertices = src_geom,
    target_vertices = tgt_geom
  )

  query <- matrix(c(0.25, 0.25, 0), ncol = 3)
  out <- transform(s2s, query)
  expect_equal(out, query, tolerance = 1e-6)
})
