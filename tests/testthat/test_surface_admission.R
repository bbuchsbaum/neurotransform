admission_octahedron <- function() {
  surface_mesh(rbind(c(1,0,0),c(-1,0,0),c(0,1,0),c(0,-1,0),c(0,0,1),c(0,0,-1)),
    rbind(c(1,3,5),c(3,2,5),c(2,4,5),c(4,1,5),c(3,1,6),c(2,3,6),c(4,2,6),c(1,4,6)))
}

test_that("sphere admission accepts winding, order and radius changes", {
  mesh <- admission_octahedron()
  set.seed(50926)
  for (scale in c(1e-100, 1e-6, 1, 100, 1e6, 1e100)) {
    m <- mesh
    m@coords <- mesh@coords * scale
    m@faces <- m@faces[sample.int(8), ]
    m@faces[c(1,4,6), ] <- m@faces[c(1,4,6), 3:1]
    d <- validate_surface_mesh(m)
    expect_true(d$valid)
    expect_equal(d$radial_degree, 1L)
    expect_equal(d$spherical_area, 4*pi, tolerance=1e-12)
    expect_equal(d$euler_characteristic, 2L)
    expect_equal(mesh_set_radius(m, 100)@coords, mesh@coords * 100)
  }
  # Face-free reference is a set of queries, not an invalid closed mesh.
  expect_s3_class(surface_resampling_plan(surface_mesh(matrix(c(1,1,1)/sqrt(3),1)), mesh),
                   "SurfaceResamplingPlan")
})

test_that("underflow cannot certify a radial determinant sign", {
  vertices <- rbind(c(.7661979452773544,-.6426046285646877,1e-323),
                    c(.07826187758764308,-.9969328355092216,-2e-323),
                    c(-.28286560123812593,.9591595548375638,1e-323), c(0,0,-1))
  mesh <- surface_mesh(vertices, rbind(c(1,3,2),c(1,2,4),c(1,4,3),c(2,3,4)))
  d <- validate_surface_mesh(mesh)
  expect_false(d$valid)
  expect_true("indeterminate_radial_faces" %in% names(d$issues))
})

test_that("spherical admission diagnoses malformed topology and geometry", {
  base <- admission_octahedron()
  hole <- base; hole@faces <- hole@faces[-1, ]
  expect_true("boundary_edge_faces" %in% names(validate_surface_mesh(hole)$issues))
  duplicated <- base; duplicated@faces <- rbind(duplicated@faces, duplicated@faces[1, ])
  expect_true("duplicate_faces" %in% names(validate_surface_mesh(duplicated)$issues))
  unused <- base; unused@coords <- rbind(unused@coords, c(1,0,0))
  expect_true("unused_vertices" %in% names(validate_surface_mesh(unused)$issues))
  split <- base; split@coords <- rbind(base@coords,base@coords)
  split@faces <- rbind(base@faces,base@faces+6L)
  expect_true("disconnected_faces" %in% names(validate_surface_mesh(split)$issues))
  # A pinched vertex joins two closed spheres, so its link has two cycles.
  pinch <- split; pinch@faces[pinch@faces==6L] <- 0L
  expect_true("nonmanifold_vertices" %in% names(validate_surface_mesh(pinch)$issues))
  folded <- base; folded@coords[5,] <- c(.1,.1,-1)/sqrt(1.02)
  expect_true("folded_faces" %in% names(validate_surface_mesh(folded)$issues))
  translated <- base; translated@coords <- translated@coords + .1
  expect_true("radius_spread" %in% names(validate_surface_mesh(translated)$issues))
  zero <- base; zero@coords[1,] <- 0
  expect_true("invalid_radius_vertices" %in% names(validate_surface_mesh(zero)$issues))
  expect_error(surface_resampling_plan(base,base,radius=0), "finite and positive")
  expect_error(neurotransform:::cpp_validate_surface(base@coords,base@faces,TRUE,Inf),
               "finite and greater")
  flat <- surface_mesh(rbind(c(0,0,0),c(1,0,0),c(0,1,0)), matrix(1:3,1))
  expect_true(validate_surface_mesh(flat,spherical=FALSE)$valid)
  expect_false(validate_surface_mesh(flat)$valid)
})

test_that("positive orientation and sphere topology do not hide degree-two overlap", {
  angle <- 4*pi*(0:4)/5
  vertices <- rbind(c(0,0,1),c(0,0,-1),cbind(cos(angle),sin(angle),0))
  faces <- do.call(rbind,lapply(0:4,function(i) {
    a <- i+3L; b <- (i+1L)%%5L+3L
    rbind(c(1L,a,b),c(2L,b,a))
  }))
  d <- validate_surface_mesh(surface_mesh(vertices,faces))
  expect_false(d$valid)
  expect_equal(d$euler_characteristic,2L)
  expect_equal(d$radial_degree,2L)
  expect_equal(d$spherical_area,8*pi,tolerance=1e-12)
  expect_true("radial_degree_not_one" %in% names(d$issues))
  expect_false("folded_faces" %in% names(d$issues))
})
