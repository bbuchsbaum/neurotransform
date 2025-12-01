test_that("apply_affine works for vectors and matrices", {
  affine <- diag(4)
  affine[1:3, 4] <- c(1, -2, 3)

  vec <- apply_affine(c(0, 0, 0), affine)
  expect_equal(vec, c(1, -2, 3))

  pts <- matrix(c(1, 2, 3,
                  -1, 0, 5), ncol = 3, byrow = TRUE)
  res <- apply_affine(pts, affine)
  expect_equal(res[1, ], c(2, 0, 6))
  expect_equal(res[2, ], c(0, -2, 8))
})

test_that("compose_affines matches manual composition", {
  A <- diag(4); A[1, 4] <- 1
  B <- diag(4); B[2, 4] <- 2
  C <- compose_affines(A, B)
  pts <- matrix(c(0, 0, 0), ncol = 3)
  out <- apply_affine(apply_affine(pts, A), B)
  expect_equal(apply_affine(pts, C), out)
})

test_that("voxels_to_world and world_to_voxels are inverses", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 3; aff[3, 3] <- 4
  ijk <- matrix(c(0, 0, 0,
                  1, 1, 1), ncol = 3, byrow = TRUE)
  world <- voxels_to_world(ijk, aff)
  back <- world_to_voxels(world, aff)
  expect_equal(back, ijk)
})

test_that("world_to_voxels handles vector input", {
  aff <- diag(4)
  aff[1, 4] <- 1
  out <- world_to_voxels(c(3, 4, 5), aff)
  expect_equal(out, c(2, 4, 5))
})
