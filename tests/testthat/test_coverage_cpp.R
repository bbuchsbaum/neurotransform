# Coverage Tests for C++ Functions
#
# Targets:
# - rcpp_bary.cpp: 0.00% -> improve
# - rcpp_apply.cpp: 39.13% -> improve

# ==============================================================================
# CPP BARYCENTRIC WEIGHTS TESTS
# ==============================================================================

test_that("cpp_barycentric_weights computes weights for point on triangle", {
  # Simple right triangle in XY plane
  vertices <- matrix(c(0, 0, 0,
                       1, 0, 0,
                       0, 1, 0), ncol = 3, byrow = TRUE)
  # Face indices (0-based)
  faces <- matrix(c(0, 1, 2), ncol = 3)

  # Query point at centroid
  coords <- matrix(c(1/3, 1/3, 0), ncol = 3)

  result <- neurotransform:::cpp_barycentric_weights(coords, vertices, faces)

  expect_true(is.list(result))
  expect_true("rows" %in% names(result))
  expect_true("cols" %in% names(result))
  expect_true("vals" %in% names(result))

  # Should have 3 weights (one per vertex)
  expect_equal(length(result$rows), 3)

  # Weights should sum to 1
  expect_equal(sum(result$vals), 1, tolerance = 1e-6)
})

test_that("cpp_barycentric_weights handles point at vertex", {
  vertices <- matrix(c(0, 0, 0,
                       1, 0, 0,
                       0, 1, 0), ncol = 3, byrow = TRUE)
  faces <- matrix(c(0, 1, 2), ncol = 3)

  # Point exactly at first vertex
  coords <- matrix(c(0, 0, 0), ncol = 3)

  result <- neurotransform:::cpp_barycentric_weights(coords, vertices, faces)

  # Should have weight 1 at first vertex (index 1 in 1-based)
  expect_true(any(result$vals > 0.99))
})

test_that("cpp_barycentric_weights handles point on edge", {
  vertices <- matrix(c(0, 0, 0,
                       2, 0, 0,
                       0, 2, 0), ncol = 3, byrow = TRUE)
  faces <- matrix(c(0, 1, 2), ncol = 3)

  # Point on edge between v0 and v1
  coords <- matrix(c(1, 0, 0), ncol = 3)

  result <- neurotransform:::cpp_barycentric_weights(coords, vertices, faces)

  # Weights should exist and sum to 1
  expect_true(length(result$vals) > 0)
  expect_equal(sum(result$vals), 1, tolerance = 1e-6)
})

test_that("cpp_barycentric_weights handles multiple query points", {
  vertices <- matrix(c(0, 0, 0,
                       1, 0, 0,
                       0, 1, 0,
                       1, 1, 0), ncol = 3, byrow = TRUE)
  faces <- matrix(c(0, 1, 2,
                    1, 3, 2), ncol = 3, byrow = TRUE)

  # Multiple query points
  coords <- matrix(c(0.2, 0.2, 0,
                     0.8, 0.8, 0), ncol = 3, byrow = TRUE)

  result <- neurotransform:::cpp_barycentric_weights(coords, vertices, faces)

  # Should have results for both points
  expect_true(1 %in% result$rows)
  expect_true(2 %in% result$rows)
})

test_that("cpp_barycentric_weights handles point outside all triangles", {
  vertices <- matrix(c(0, 0, 0,
                       1, 0, 0,
                       0, 1, 0), ncol = 3, byrow = TRUE)
  faces <- matrix(c(0, 1, 2), ncol = 3)

  # Point far outside
  coords <- matrix(c(100, 100, 0), ncol = 3)

  result <- neurotransform:::cpp_barycentric_weights(coords, vertices, faces)

  # May have no weights or clipped weights
  expect_true(is.list(result))
})

# ==============================================================================
# CPP APPLY PROJECTOR TESTS
# ==============================================================================

test_that("cpp_apply_projector applies sparse matrix", {
  skip_if_not_installed("Matrix")

  # Create sparse identity-like projector
  proj <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = rep(1, 3),
                               dims = c(3, 3))

  data <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)

  result <- neurotransform:::cpp_apply_projector(proj, data, threads = 1)
  expected <- as.matrix(proj %*% data)

  expect_equal(dim(result), dim(data))
  expect_equal(as.matrix(result), expected)
})

test_that("cpp_apply_projector applies averaging projector", {
  skip_if_not_installed("Matrix")

  # Projector that averages first two rows
  proj <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2),
    j = c(1, 2, 3, 4),
    x = c(0.5, 0.5, 0.5, 0.5),
    dims = c(2, 4)
  )

  data <- matrix(c(1, 2, 3, 4,
                   5, 6, 7, 8), nrow = 4, ncol = 2)

  result <- neurotransform:::cpp_apply_projector(proj, data, threads = 1)
  expected <- as.matrix(proj %*% data)

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(as.matrix(result), expected, tolerance = 1e-12)
})

test_that("cpp_apply_projector errors when data rows do not match projector source", {
  skip_if_not_installed("Matrix")

  proj <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = c(1, 1), dims = c(2, 2))
  data <- matrix(1:3, nrow = 3, ncol = 1)

  expect_error(
    neurotransform:::cpp_apply_projector(proj, data, threads = 1),
    "projector expects 2 source elements"
  )
})

test_that("cpp_apply_projector returns the same result across thread counts", {
  skip_if_not_installed("Matrix")

  set.seed(1)
  proj <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 1, 3, 5, 6),
    j = c(1, 1, 2, 2, 3, 4, 5, 6),
    x = runif(8),
    dims = c(6, 6),
    repr = "C"
  )
  data <- matrix(runif(6 * 8), nrow = 6, ncol = 8)

  single <- neurotransform:::cpp_apply_projector(proj, data, threads = 1)
  multi <- neurotransform:::cpp_apply_projector(proj, data, threads = 4)

  expect_equal(multi, single, tolerance = 1e-12)
})

test_that("cpp_apply_projector handles realistic tall projector shapes", {
  skip_if_not_installed("Matrix")
  skip_on_cran()

  set.seed(42)
  n_rows <- 140000L
  n_cols <- 32000L
  nnz_per_row <- 2L

  proj <- Matrix::sparseMatrix(
    i = rep(seq_len(n_rows), each = nnz_per_row),
    j = sample.int(n_cols, n_rows * nnz_per_row, replace = TRUE),
    x = runif(n_rows * nnz_per_row),
    dims = c(n_rows, n_cols),
    repr = "C"
  )
  data <- matrix(runif(n_cols * 3L), nrow = n_cols, ncol = 3L)

  result <- neurotransform:::cpp_apply_projector(proj, data, threads = 1)
  expected <- as.matrix(proj %*% data)

  expect_equal(dim(result), c(n_rows, 3L))
  expect_equal(result, expected, tolerance = 1e-10)
})

# ==============================================================================
# CPP APPLY AFFINE CHAIN TESTS
# ==============================================================================

test_that("cpp_apply_affine_chain applies identity", {
  coords <- matrix(c(1, 2, 3,
                     4, 5, 6), ncol = 3, byrow = TRUE)

  result <- neurotransform:::cpp_apply_affine_chain(coords, list(diag(4)))

  expect_equal(dim(result), dim(coords))
  expect_equal(result, coords, tolerance = 1e-10)
})

test_that("cpp_apply_affine_chain applies scaling", {
  coords <- matrix(c(1, 1, 1), ncol = 3)
  mat <- diag(c(2, 3, 4, 1))

  result <- neurotransform:::cpp_apply_affine_chain(coords, list(mat))

  expect_equal(as.numeric(result), c(2, 3, 4))
})

test_that("cpp_apply_affine_chain applies translation", {
  coords <- matrix(c(0, 0, 0), ncol = 3)
  mat <- diag(4)
  mat[1:3, 4] <- c(10, 20, 30)

  result <- neurotransform:::cpp_apply_affine_chain(coords, list(mat))

  expect_equal(as.numeric(result), c(10, 20, 30))
})

test_that("cpp_apply_affine_chain composes multiple transforms", {
  coords <- matrix(c(0, 0, 0), ncol = 3)

  # First: translate by (1,0,0)
  mat1 <- diag(4)
  mat1[1, 4] <- 1

  # Second: scale by 2
  mat2 <- diag(c(2, 2, 2, 1))

  # Combined: 0 + 1 = 1, then 1 * 2 = 2
  result <- neurotransform:::cpp_apply_affine_chain(coords, list(mat1, mat2))

  expect_equal(as.numeric(result), c(2, 0, 0))
})

test_that("cpp_apply_affine_chain handles many matrices", {
  coords <- matrix(c(1, 1, 1), ncol = 3)

  # 5 identity matrices
  mats <- lapply(1:5, function(i) diag(4))

  result <- neurotransform:::cpp_apply_affine_chain(coords, mats)

  expect_equal(result, coords)
})

test_that("cpp_apply_affine_chain handles rotation", {
  coords <- matrix(c(1, 0, 0), ncol = 3)

  # 90 degree rotation around Z axis
  theta <- pi/2
  rot <- matrix(c(cos(theta), -sin(theta), 0, 0,
                  sin(theta), cos(theta), 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1), 4, 4, byrow = TRUE)

  result <- neurotransform:::cpp_apply_affine_chain(coords, list(rot))

  # (1,0,0) rotated 90 deg around Z -> (0,1,0)
  expect_equal(as.numeric(result), c(0, 1, 0), tolerance = 1e-10)
})

# ==============================================================================
# CPP NEAREST VERTEX TESTS
# ==============================================================================

test_that("cpp_nearest_vertex finds closest vertex", {
  vertices <- matrix(c(0, 0, 0,
                       10, 0, 0,
                       0, 10, 0,
                       0, 0, 10), ncol = 3, byrow = TRUE)

  # Query points near each vertex
  coords <- matrix(c(0.1, 0.1, 0.1,   # Near vertex 1
                     9.9, 0.1, 0.1,   # Near vertex 2
                     0.1, 9.9, 0.1,   # Near vertex 3
                     0.1, 0.1, 9.9),  # Near vertex 4
                   ncol = 3, byrow = TRUE)

  result <- neurotransform:::cpp_nearest_vertex(coords, vertices)

  expect_equal(length(result), 4)
  expect_equal(result, c(1, 2, 3, 4))
})

test_that("cpp_nearest_vertex handles single query", {
  vertices <- matrix(c(0, 0, 0,
                       1, 1, 1), ncol = 3, byrow = TRUE)

  coords <- matrix(c(0.1, 0.1, 0.1), ncol = 3)

  result <- neurotransform:::cpp_nearest_vertex(coords, vertices)

  expect_equal(result, 1)
})

# ==============================================================================
# CPP SAMPLE VOLUME TESTS
# ==============================================================================

test_that("cpp_sample_volume linear interpolation", {
  # 2x2x2 volume with known values
  vol <- array(c(0, 1, 2, 3, 4, 5, 6, 7), dim = c(2, 2, 2))

  # Identity world-to-vox
  w2v <- diag(4)

  # Query at center (0.5, 0.5, 0.5)
  coords <- matrix(c(0.5, 0.5, 0.5), ncol = 3)

  result <- neurotransform:::cpp_sample_volume(vol, coords, w2v,
                                                method = "linear", outside = NA)

  # Should be average of all 8 corners = 3.5
  expect_equal(result, 3.5, tolerance = 0.1)
})

test_that("cpp_sample_volume nearest interpolation", {
  vol <- array(c(0, 1, 2, 3, 4, 5, 6, 7), dim = c(2, 2, 2))
  w2v <- diag(4)

  # Query at (0.1, 0.1, 0.1) - nearest to (0,0,0)
  coords <- matrix(c(0.1, 0.1, 0.1), ncol = 3)

  result <- neurotransform:::cpp_sample_volume(vol, coords, w2v,
                                                method = "nearest", outside = NA)

  expect_equal(result, vol[1, 1, 1])
})

test_that("cpp_sample_volume returns outside for out-of-bounds", {
  vol <- array(1, dim = c(3, 3, 3))
  w2v <- diag(4)

  coords <- matrix(c(100, 100, 100), ncol = 3)

  result <- neurotransform:::cpp_sample_volume(vol, coords, w2v,
                                                method = "linear", outside = -999)

  expect_equal(result, -999)
})

test_that("cpp_sample_volume handles 4D volume", {
  vol <- array(1:16, dim = c(2, 2, 2, 2))
  w2v <- diag(4)

  coords <- matrix(c(0.5, 0.5, 0.5), ncol = 3)

  result <- neurotransform:::cpp_sample_volume(vol, coords, w2v,
                                                method = "linear", outside = NA)

  # Should return 2 values (one per 4th dimension)
  expect_equal(length(result), 2)
})

# ==============================================================================
# CPP VOLUME WORLD COORDS TESTS
# ==============================================================================

test_that("cpp_volume_world_coords generates correct grid", {
  dims <- c(2L, 2L, 2L)
  affine <- diag(4)

  coords <- neurotransform:::cpp_volume_world_coords(dims, affine)

  expect_equal(nrow(coords), 8)  # 2^3
  expect_equal(ncol(coords), 3)

  # First voxel (0,0,0) should be at origin
  expect_equal(coords[1, ], c(0, 0, 0))
})

test_that("cpp_volume_world_coords respects affine", {
  dims <- c(2L, 2L, 2L)
  affine <- diag(c(2, 2, 2, 1))  # 2mm spacing
  affine[1:3, 4] <- c(10, 20, 30)  # Origin offset

  coords <- neurotransform:::cpp_volume_world_coords(dims, affine)

  # First voxel should be at (10, 20, 30)
  expect_equal(coords[1, ], c(10, 20, 30))
})

# ==============================================================================
# CPP APPLY WARP FIELD TESTS
# ==============================================================================

test_that("cpp_apply_warp_field applies displacement field", {
  # Create simple 3x3x3 displacement field (zero displacement)
  field <- array(0, dim = c(3, 3, 3, 3))
  dims <- c(3L, 3L, 3L)

  # World to vox is identity
  w2v <- diag(4)

  coords <- matrix(c(1, 1, 1), ncol = 3)

  result <- neurotransform:::cpp_apply_warp_field(coords, field, dims, w2v)

  expect_equal(nrow(result), 1)
  # Zero displacement: output = input
  expect_equal(as.numeric(result), c(1, 1, 1), tolerance = 0.1)
})

test_that("cpp_apply_warp_field handles non-zero displacement", {
  # Create displacement field with constant (1,1,1) displacement
  field <- array(1, dim = c(3, 3, 3, 3))
  dims <- c(3L, 3L, 3L)
  w2v <- diag(4)

  coords <- matrix(c(1, 1, 1), ncol = 3)

  result <- neurotransform:::cpp_apply_warp_field(coords, field, dims, w2v)

  # Displacement of (1,1,1) added
  expect_equal(as.numeric(result), c(2, 2, 2), tolerance = 0.1)
})

test_that("cpp_apply_warp_field_cubic applies cubic interpolation", {
  field <- array(0, dim = c(5, 5, 5, 3))
  dims <- c(5L, 5L, 5L)
  w2v <- diag(4)

  coords <- matrix(c(2, 2, 2), ncol = 3)

  result <- neurotransform:::cpp_apply_warp_field_cubic(coords, field, dims, w2v)

  expect_equal(nrow(result), 1)
})

# ==============================================================================
# CPP TRILINEAR WEIGHTS TESTS
# ==============================================================================

test_that("cpp_trilinear_weights computes weights", {
  # Voxel coordinates
  vox_coords <- matrix(c(0.5, 0.5, 0.5), ncol = 3)
  dims <- c(3L, 3L, 3L)

  result <- neurotransform:::cpp_trilinear_weights(vox_coords, dims)

  expect_true(is.list(result))
  expect_true("rows" %in% names(result))
  expect_true("cols" %in% names(result))
  expect_true("vals" %in% names(result))

  # Should have 8 weights (one per corner of cube)
  expect_equal(length(result$vals), 8)

  # Weights should sum to 1
  expect_equal(sum(result$vals), 1, tolerance = 1e-6)
})

test_that("cpp_trilinear_weights handles integer coordinates", {
  # At exact voxel location
  vox_coords <- matrix(c(1, 1, 1), ncol = 3)
  dims <- c(3L, 3L, 3L)

  result <- neurotransform:::cpp_trilinear_weights(vox_coords, dims)

  # At integer location, only one weight should be 1
  expect_true(any(result$vals > 0.99))
})

# ==============================================================================
# CPP TRIPLETS TO DGC TESTS
# ==============================================================================

test_that("cpp_triplets_to_dgC creates sparse matrix", {
  i <- c(1L, 2L, 3L)  # Row indices (1-based)
  j <- c(1L, 2L, 3L)  # Col indices (1-based)
  x <- c(1.0, 2.0, 3.0)  # Values
  nrow <- 3L
  ncol <- 3L

  result <- neurotransform:::cpp_triplets_to_dgC(i, j, x, nrow, ncol)

  expect_s4_class(result, "dgCMatrix")
  expect_equal(dim(result), c(3, 3))
})

test_that("cpp_triplets_to_dgC handles off-diagonal", {
  i <- c(1L, 1L, 2L)
  j <- c(1L, 2L, 3L)
  x <- c(1.0, 2.0, 3.0)

  result <- neurotransform:::cpp_triplets_to_dgC(i, j, x, 3L, 3L)

  expect_s4_class(result, "dgCMatrix")
})

# ==============================================================================
# CPP COMPOSE WARP FIELDS TESTS
# ==============================================================================

test_that("cpp_compose_warp_fields combines two fields", {
  # Two identity-like warp fields (zero displacement)
  fieldA <- array(0, dim = c(3, 3, 3, 3))
  fieldB <- array(0, dim = c(3, 3, 3, 3))
  dims <- c(3L, 3L, 3L)
  w2v_A <- diag(4)
  v2w_B <- diag(4)

  result <- neurotransform:::cpp_compose_warp_fields(
    fieldA, dims, w2v_A,
    fieldB, dims, v2w_B
  )

  expect_true(is.list(result))
  expect_true("field" %in% names(result))
})

# ==============================================================================
# CPP ABSOLUTE TO DISPLACEMENT TESTS
# ==============================================================================

test_that("cpp_absolute_to_displacement converts coordinate field", {
  # Absolute coordinate field where each voxel points to itself
  dims <- c(3L, 3L, 3L)
  v2w <- diag(4)

  # Create field in correct format (flattened: X * Y * Z * 3)
  field <- numeric(3 * 3 * 3 * 3)

  # Fill with absolute coordinates (identity mapping)
  for (i in 0:2) {
    for (j in 0:2) {
      for (k in 0:2) {
        base <- ((k * 3 + j) * 3 + i) * 3 + 1  # 1-based index
        field[base]     <- i  # 0-based voxel x coord
        field[base + 1] <- j  # 0-based voxel y coord
        field[base + 2] <- k  # 0-based voxel z coord
      }
    }
  }

  result <- neurotransform:::cpp_absolute_to_displacement(field, dims, v2w)

  # Displacement should be approximately zero (identity mapping)
  expect_true(is.numeric(result))
  expect_equal(length(result), 3 * 3 * 3 * 3)
  # All displacements should be near zero
  expect_true(all(abs(result) < 1e-6))
})

# ==============================================================================
# CPP PATH APPLY FUNCTIONS TESTS
# ==============================================================================

test_that("cpp_path_apply_affine_warp applies affine then warp", {
  coords <- matrix(c(1, 1, 1), ncol = 3)
  affines <- list(diag(4))

  # Zero displacement warp field (flattened)
  field <- numeric(3 * 3 * 3 * 3)
  dims <- c(3L, 3L, 3L)
  w2v <- diag(4)

  result <- neurotransform:::cpp_path_apply_affine_warp(
    coords, affines, field, dims, w2v, cubic = FALSE
  )

  expect_equal(nrow(result), 1)
})

test_that("cpp_path_apply_affine_warp_affine handles full path", {
  coords <- matrix(c(1, 1, 1), ncol = 3)
  affines_pre <- list(diag(4))
  affines_post <- list(diag(4))

  field <- numeric(3 * 3 * 3 * 3)
  dims <- c(3L, 3L, 3L)
  w2v <- diag(4)

  result <- neurotransform:::cpp_path_apply_affine_warp_affine(
    coords, affines_pre, field, dims, w2v, cubic = FALSE, affines_post
  )

  expect_equal(nrow(result), 1)
})

# ==============================================================================
# CPP RESAMPLE PLAN TESTS
# ==============================================================================

test_that("cpp_make_resample_plan creates plan", {
  src_dim <- c(3L, 3L, 3L)
  src_coords_vox <- matrix(c(0, 0, 0,
                              1, 1, 1,
                              2, 2, 2), ncol = 3, byrow = TRUE)
  method <- "linear"

  plan <- neurotransform:::cpp_make_resample_plan(src_dim, src_coords_vox, method)

  # Plan should be an external pointer
  expect_true(inherits(plan, "externalptr"))
})

test_that("cpp_apply_resample_plan applies precomputed plan", {
  src_dim <- c(3L, 3L, 3L)
  src_coords_vox <- matrix(c(1, 1, 1), ncol = 3)

  plan <- neurotransform:::cpp_make_resample_plan(src_dim, src_coords_vox, "linear")

  # Create source data
  data <- array(1:27, dim = c(3, 3, 3))

  result <- neurotransform:::cpp_apply_resample_plan(plan, data, NA_real_)

  expect_true(is.numeric(result))
  expect_equal(length(result), 1)
})

test_that("cpp_apply_resample_plan matches direct sampling at boundary-sensitive coordinates", {
  src_dim <- c(5L, 5L, 5L)
  data <- array(0, dim = src_dim)
  for (x in 0:(src_dim[1] - 1)) {
    for (y in 0:(src_dim[2] - 1)) {
      for (z in 0:(src_dim[3] - 1)) {
        data[x + 1, y + 1, z + 1] <- x^2 + 2 * y^2 + 3 * z^2 + 0.5 * x * y
      }
    }
  }

  coords <- matrix(c(
    -0.51, 1.2, 2.0,
    -0.49, 1.2, 2.0,
    -0.25, 1.2, 2.0,
     0.49, 1.2, 2.0,
     3.49, 2.2, 1.4,
     3.51, 2.2, 1.4,
     4.49, 2.2, 1.4,
     4.51, 2.2, 1.4,
     2.3, -0.49, 1.7,
     2.3, 4.49, 1.7
  ), ncol = 3, byrow = TRUE)

  for (method in c("nearest", "linear", "cubic")) {
    plan <- neurotransform:::cpp_make_resample_plan(src_dim, coords, method)
    via_plan <- neurotransform:::cpp_apply_resample_plan(plan, data, outside = -999)
    direct <- neurotransform:::cpp_sample_volume(data, coords, diag(4), method, outside = -999)

    expect_equal(via_plan, direct, tolerance = 1e-10, info = method)
  }
})

# ==============================================================================
# CPP RIBBON FUNCTIONS TESTS
# ==============================================================================

test_that("cpp_ribbon_weights computes ribbon sampling weights", {
  inner <- matrix(c(1, 1, 1, 2, 2, 2), ncol = 3, byrow = TRUE)
  outer <- matrix(c(2, 2, 2, 3, 3, 3), ncol = 3, byrow = TRUE)
  dims <- c(5L, 5L, 5L)
  w2v <- diag(4)
  n_steps <- 3L
  mask <- rep(TRUE, 5 * 5 * 5)

  result <- neurotransform:::cpp_ribbon_weights(inner, outer, w2v, dims, n_steps, mask)

  expect_true(is.list(result))
  expect_true("rows" %in% names(result))
  expect_true("cols" %in% names(result))
  expect_true("vals" %in% names(result))
})

test_that("cpp_ribbon_sample_volume samples along ribbon", {
  data <- array(1:125, dim = c(5, 5, 5))
  inner <- matrix(c(2, 2, 2), ncol = 3)
  outer <- matrix(c(3, 3, 3), ncol = 3)
  w2v <- diag(4)
  n_steps <- 3L

  result <- neurotransform:::cpp_ribbon_sample_volume(data, inner, outer, w2v,
                                                       n_steps, "linear")

  expect_true(is.numeric(result))
})
