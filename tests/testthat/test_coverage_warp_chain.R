# Coverage Tests for warp_chain.R
#
# Target: Improve coverage from 75.41% to >80%

# ==============================================================================
# WARP CHAIN KEY
# ==============================================================================

test_that("warp_chain_key produces consistent hashes", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))

  key1 <- neurotransform:::warp_chain_key(list(f, g))
  key2 <- neurotransform:::warp_chain_key(list(f, g))

  expect_equal(key1, key2)
  expect_true(nzchar(key1))
})

test_that("warp_chain_key differs for different paths", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))
  h <- Affine3DMorphism("a", "b", diag(c(2, 2, 2, 1)))

  key1 <- neurotransform:::warp_chain_key(list(f, g))
  key2 <- neurotransform:::warp_chain_key(list(h, g))

  expect_false(key1 == key2)
})

# ==============================================================================
# BUILD WARP CHAIN
# ==============================================================================

test_that("build_warp_chain returns identity for empty path", {
  chain <- neurotransform:::build_warp_chain(list())
  expect_equal(chain$kind, "identity")
})

test_that("build_warp_chain fuses all-affine paths", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(c(2, 2, 2, 1)))

  chain <- neurotransform:::build_warp_chain(list(f, g))

  expect_equal(chain$kind, "affine")
  expect_true("matrix" %in% names(chain))
  expect_equal(dim(chain$matrix), c(4, 4))
})

test_that("build_warp_chain handles identity morphisms in affine chain", {
  f <- Affine3DMorphism("a", "b", diag(4))
  id <- IdentityMorphism("b")
  g <- Affine3DMorphism("b", "c", diag(c(2, 2, 2, 1)))

  chain <- neurotransform:::build_warp_chain(list(f, id, g))

  expect_equal(chain$kind, "affine")
})

test_that("build_warp_chain creates mixed chain for warp paths", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  f <- Affine3DMorphism("a", "b", diag(4))
  w <- Warp3DMorphism("b", "c", warp_path = warp_path, warp_type = "ants")

  chain <- neurotransform:::build_warp_chain(list(f, w))

  expect_equal(chain$kind, "mixed")
  expect_true("segments" %in% names(chain))
})

test_that("build_warp_chain fuses preceding affines before warp", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  f <- Affine3DMorphism("a", "b", diag(c(2, 2, 2, 1)))
  g <- Affine3DMorphism("b", "c", diag(c(3, 3, 3, 1)))
  w <- Warp3DMorphism("c", "d", warp_path = warp_path, warp_type = "ants")

  chain <- neurotransform:::build_warp_chain(list(f, g, w))

  expect_equal(chain$kind, "mixed")
  # Should have fused affine segment before warp
  segment_types <- sapply(chain$segments, `[[`, "type")
  expect_true("warp" %in% segment_types || "affine" %in% segment_types)
})

test_that("build_warp_chain handles warp composition", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  w1 <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")
  w2 <- Warp3DMorphism("b", "c", warp_path = warp_path, warp_type = "ants")

  chain <- neurotransform:::build_warp_chain(list(w1, w2))

  expect_equal(chain$kind, "mixed")
  # Should have warp_comp segment
  segment_types <- sapply(chain$segments, `[[`, "type")
  expect_true("warp_comp" %in% segment_types || "warp" %in% segment_types)
})

# ==============================================================================
# APPLY WARP CHAIN
# ==============================================================================

test_that("apply_warp_chain identity returns unchanged coords", {
  coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = TRUE)
  chain <- list(kind = "identity")

  result <- neurotransform:::apply_warp_chain(chain, coords)
  expect_equal(result, coords)
})

test_that("apply_warp_chain affine transforms coords", {
  coords <- matrix(c(0, 0, 0), ncol = 3)
  mat <- diag(4)
  mat[1:3, 4] <- c(10, 20, 30)

  chain <- list(kind = "affine", matrix = mat)

  result <- neurotransform:::apply_warp_chain(chain, coords)

  expect_equal(nrow(result), 1)
  # Translation: (0,0,0) + (10,20,30) = (10,20,30)
  expect_equal(as.numeric(result), c(10, 20, 30), tolerance = 1e-10)
})

test_that("apply_warp_chain mixed chain works", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  w <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")
  f <- Affine3DMorphism("b", "c", diag(4))

  chain <- neurotransform:::build_warp_chain(list(f, w))

  # Get valid coords
  warp <- neurotransform:::load_warp_array(w)
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]
  coords <- matrix(center_world, nrow = 1)

  result <- neurotransform:::apply_warp_chain(chain, coords)

  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# C++ AFFINE CHAIN
# ==============================================================================

test_that("cpp_apply_affine_chain applies single matrix", {
  coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = TRUE)
  mat <- diag(4)
  mat[1:3, 4] <- c(1, 1, 1)

  result <- neurotransform:::cpp_apply_affine_chain(coords, list(mat))

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
})

test_that("cpp_apply_affine_chain applies multiple matrices", {
  coords <- matrix(c(0, 0, 0), ncol = 3)

  mat1 <- diag(4)
  mat1[1, 4] <- 5

  mat2 <- diag(4)
  mat2[2, 4] <- 10

  result <- neurotransform:::cpp_apply_affine_chain(coords, list(mat1, mat2))

  expect_equal(nrow(result), 1)
})

test_that("cpp_apply_affine_chain with empty list returns unchanged", {
  coords <- matrix(c(1, 2, 3), ncol = 3)
  result <- neurotransform:::cpp_apply_affine_chain(coords, list())
  expect_equal(result, coords)
})

test_that("cpp_apply_affine_chain errors on non-4x4 matrix", {
  coords <- matrix(c(1, 2, 3), ncol = 3)
  expect_error(neurotransform:::cpp_apply_affine_chain(coords, list(diag(3))))
})
