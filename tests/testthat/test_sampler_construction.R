# Test Sampler and Grid construction and properties

# =============================================================================
# IRREGULAR REFERENCE TESTS
# =============================================================================

test_that("sampled_points constructs lightweight irregular references", {
  coords <- matrix(c(
    0, 0, 0,
    1, 2, 3,
    -1, 4, 2
  ), ncol = 3, byrow = TRUE)

  pts <- sampled_points(coords)
  expect_s4_class(pts, "SampledPoints")
  expect_equal(pts@coords, coords)
})

test_that("surface_mesh accepts 1-based faces and stores 0-based faces", {
  verts <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  faces <- matrix(c(1, 2, 3), ncol = 3, byrow = TRUE)

  mesh <- surface_mesh(verts, faces)
  expect_s4_class(mesh, "SurfaceMesh")
  expect_equal(mesh@faces, matrix(c(0L, 1L, 2L), ncol = 3, byrow = TRUE))
})

test_that("mesh_is_sphere and mesh_set_radius work for spherical meshes", {
  verts <- matrix(c(
    1, 0, 0,
    -1, 0, 0,
    0, 1, 0,
    0, -1, 0,
    0, 0, 1,
    0, 0, -1
  ), ncol = 3, byrow = TRUE)
  mesh <- surface_mesh(verts)

  expect_true(mesh_is_sphere(mesh))
  mesh2 <- mesh_set_radius(mesh, radius = 100)
  radii <- sqrt(rowSums(mesh2@coords^2))
  expect_equal(radii, rep(100, length(radii)), tolerance = 1e-6)
})

# =============================================================================
# VOLUME SAMPLER TESTS
# =============================================================================

test_that("volume_sampler creates sampler for each interpolation method", {
  vol <- array(runif(27), dim = c(3, 3, 3))

  for (method in c("linear", "nearest", "cubic")) {
    sampler <- volume_sampler(vol, affine = diag(4), method = method)
    expect_s4_class(sampler, "Sampler")
    expect_equal(sampler@method, method)
    expect_equal(sampler@dims, c(3L, 3L, 3L))
    expect_equal(sampler@vdim, 1L)
  }
})

test_that("volume_sampler computes correct bounds for identity affine", {
  vol <- array(1, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4))

  # With 3 voxels (indices 0,1,2), 1mm spacing, origin at 0
  # Corners are at 0 and 2 in each dimension
  expect_equal(sampler@bounds$min, c(0, 0, 0))
  expect_equal(sampler@bounds$max, c(2, 2, 2))
})

test_that("volume_sampler computes correct bounds with offset origin", {
  vol <- array(1, dim = c(2, 2, 2))
  aff <- diag(4)
  aff[1:3, 4] <- c(10, 20, 30)  # Offset origin

  sampler <- volume_sampler(vol, affine = aff)

  # 2 voxels at 1mm spacing with offset
  expect_equal(sampler@bounds$min, c(10, 20, 30))
  expect_equal(sampler@bounds$max, c(11, 21, 31))
})

test_that("volume_sampler computes correct bounds with scaled affine", {
  vol <- array(1, dim = c(2, 2, 2))
  aff <- diag(c(2, 2, 2, 1))  # 2mm voxels

  sampler <- volume_sampler(vol, affine = aff)

  # 2 voxels at 2mm spacing
  expect_equal(sampler@bounds$min, c(0, 0, 0))
  expect_equal(sampler@bounds$max, c(2, 2, 2))  # (0-1)*2 = 0-2
})

test_that("volume_sampler handles 4D volumes correctly", {
  vol <- array(runif(3 * 3 * 3 * 5), dim = c(3, 3, 3, 5))
  sampler <- volume_sampler(vol, affine = diag(4))

  expect_equal(sampler@vdim, 5L)
  expect_equal(sampler@dims, c(3L, 3L, 3L))

  # Sampling should return 5 values per coordinate
  coords <- matrix(c(0, 0, 0), ncol = 3)
  result <- sampler@evaluate(coords)
  expect_equal(length(result), 5)
})

test_that("volume_sampler extracts affine from data if not provided", {
  vol <- array(1, dim = c(3, 3, 3))

  # Without explicit affine, should use identity
  sampler <- volume_sampler(vol)
  expect_equal(sampler@bounds$min, c(0, 0, 0))
  expect_equal(sampler@bounds$max, c(2, 2, 2))
})

test_that("volume_sampler respects outside parameter", {
  vol <- array(1, dim = c(3, 3, 3))

  sampler_na <- volume_sampler(vol, affine = diag(4), outside = NA_real_)
  sampler_zero <- volume_sampler(vol, affine = diag(4), outside = 0)
  sampler_neg <- volume_sampler(vol, affine = diag(4), outside = -999)

  expect_equal(sampler_na@outside, NA_real_)
  expect_equal(sampler_zero@outside, 0)
  expect_equal(sampler_neg@outside, -999)
})

test_that("volume_sampler evaluates correctly at voxel centers", {
  # Create volume with known values
  vol <- array(0, dim = c(3, 3, 3))
  vol[1, 1, 1] <- 1
  vol[2, 2, 2] <- 2
  vol[3, 3, 3] <- 3

  sampler <- volume_sampler(vol, affine = diag(4), method = "nearest")

  coords <- matrix(c(
    0, 0, 0,  # vol[1,1,1] = 1
    1, 1, 1,  # vol[2,2,2] = 2
    2, 2, 2   # vol[3,3,3] = 3
  ), ncol = 3, byrow = TRUE)

  result <- sampler@evaluate(coords)
  expect_equal(result, c(1, 2, 3))
})

test_that("volume_sampler returns outside value for out-of-bounds", {
  vol <- array(1, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4), outside = -99)

  coords <- matrix(c(
    -1, 0, 0,   # Out of bounds
    10, 10, 10  # Out of bounds
  ), ncol = 3, byrow = TRUE)

  result <- sampler@evaluate(coords)
  expect_equal(result, c(-99, -99))
})

# =============================================================================
# SURFACE SAMPLER TESTS
# =============================================================================

test_that("surface_sampler creates nearest sampler", {
  verts <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  data <- c(10, 20, 30)

  sampler <- surface_sampler(verts, data, method = "nearest")

  expect_s4_class(sampler, "Sampler")
  expect_equal(sampler@method, "nearest")
  expect_equal(sampler@vdim, 1L)
})

test_that("surface_sampler accepts SurfaceMesh and infers faces", {
  verts <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  faces <- matrix(c(0, 1, 2), ncol = 3, byrow = TRUE)
  mesh <- surface_mesh(verts, faces)
  data <- c(10, 20, 30)

  sampler <- surface_sampler(mesh, data, method = "nearest")
  expect_s4_class(sampler, "Sampler")
})

test_that("surface_sampler errors without faces for barycentric", {
  verts <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  data <- c(10, 20, 30)

  expect_error(
    surface_sampler(verts, data, method = "barycentric"),
    "Barycentric interpolation requires faces"
  )
})

test_that("surface_sampler computes bounds from vertices", {
  verts <- matrix(c(
    0, 0, 0,
    10, 20, 30,
    5, 5, 5
  ), ncol = 3, byrow = TRUE)
  data <- c(1, 2, 3)

  sampler <- surface_sampler(verts, data, method = "nearest")

  expect_equal(sampler@bounds$min, c(0, 0, 0))
  expect_equal(sampler@bounds$max, c(10, 20, 30))
})

test_that("surface_sampler handles multivariate data", {
  verts <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  data <- matrix(c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  ), ncol = 3, byrow = TRUE)

  sampler <- surface_sampler(verts, data, method = "nearest")

  expect_equal(sampler@vdim, 3L)
})

test_that("surface_sampler nearest returns closest vertex data", {
  verts <- matrix(c(
    0, 0, 0,
    10, 0, 0,
    0, 10, 0
  ), ncol = 3, byrow = TRUE)
  data <- c(100, 200, 300)

  sampler <- surface_sampler(verts, data, method = "nearest")

  # Query points closest to each vertex
  coords <- matrix(c(
    0.1, 0.1, 0,   # Closest to vertex 1 (0,0,0) -> 100
    9, 0, 0,       # Closest to vertex 2 (10,0,0) -> 200
    0, 9, 0        # Closest to vertex 3 (0,10,0) -> 300
  ), ncol = 3, byrow = TRUE)

  result <- sampler@evaluate(coords)
  expect_equal(result, c(100, 200, 300))
})

# =============================================================================
# GRID TESTS
# =============================================================================

test_that("grid_spec creates Grid with correct properties", {
  dims <- c(10L, 20L, 30L)
  aff <- diag(4)
  aff[1:3, 4] <- c(5, 10, 15)

  grid <- grid_spec(dims, aff)

  expect_s4_class(grid, "Grid")
  expect_equal(grid@dims, c(10L, 20L, 30L))
  expect_equal(grid@affine, aff)
})

test_that("grid_spec validates affine matrix", {
  expect_error(grid_spec(c(3, 3, 3), matrix(1:9, 3, 3)), "4x4")
})

test_that("grid_coords returns correct number of coordinates", {
  grid <- grid_spec(c(3L, 4L, 5L), diag(4))
  coords <- grid_coords(grid)

  expect_equal(nrow(coords), 3 * 4 * 5)
  expect_equal(ncol(coords), 3)
})

test_that("grid_coords with identity affine returns voxel indices", {
  grid <- grid_spec(c(2L, 2L, 2L), diag(4))
  coords <- grid_coords(grid)

  # Should contain all combinations of {0,1} x {0,1} x {0,1}
  expect_true(all(coords >= 0))
  expect_true(all(coords <= 1))
})

test_that("grid_coords with scaled affine returns world coordinates", {
  aff <- diag(c(2, 2, 2, 1))  # 2mm voxels
  grid <- grid_spec(c(2L, 2L, 2L), aff)
  coords <- grid_coords(grid)

  # With 2mm voxels, coordinates should be 0 and 2
  expect_true(all(coords %in% c(0, 2)))
})

test_that("grid_from_data extracts grid from array", {
  data <- array(1, dim = c(5, 6, 7))
  grid <- grid_from_data(data)

  expect_equal(grid@dims, c(5L, 6L, 7L))
  expect_equal(grid@affine, diag(4))  # Identity for plain array
})

# =============================================================================
# EXTRACT_AFFINE TESTS
# =============================================================================

test_that("extract_affine returns identity for plain array", {
  arr <- array(1, dim = c(3, 3, 3))
  aff <- extract_affine(arr)

  expect_equal(aff, diag(4))
})

test_that("extract_affine returns matrix itself for 4x4 matrix", {
  mat <- matrix(1:16, 4, 4)
  aff <- extract_affine(mat)

  expect_equal(aff, mat)
})

test_that("extract_affine returns identity for non-4x4 matrix", {
  mat <- matrix(1:9, 3, 3)
  aff <- extract_affine(mat)

  expect_equal(aff, diag(4))
})

# =============================================================================
# SHOW METHODS
# =============================================================================

test_that("Sampler show method works", {
  vol <- array(1, dim = c(3, 3, 3))
  sampler <- volume_sampler(vol, affine = diag(4))

  output <- capture.output(show(sampler))
  expect_match(output[1], "Sampler")
  expect_match(output[1], "3x3x3")
  expect_match(output[1], "linear")
})

test_that("Grid show method works", {
  grid <- grid_spec(c(10L, 20L, 30L), diag(4))

  output <- capture.output(show(grid))
  expect_match(output[1], "Grid")
  expect_match(output[1], "10x20x30")
})
