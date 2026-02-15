test_that("warp_from_field maps coordinates for displacement representation", {
  grid <- grid_spec(c(5, 5, 5), diag(c(2, 2, 2, 1)))
  field <- array(0, dim = c(5, 5, 5, 3))
  field[, , , 1] <- 1

  morph <- warp_from_field(
    source = "src",
    target = "tgt",
    field = field,
    grid = grid,
    representation = "displacements"
  )

  coords <- matrix(c(
    0, 0, 0,
    2, 4, 6,
    6, 2, 4
  ), ncol = 3, byrow = TRUE)

  out <- transform(morph, coords)
  expect_equal(out[, 1], coords[, 1] + 1, tolerance = 1e-6)
  expect_equal(out[, 2:3], coords[, 2:3], tolerance = 1e-6)
})

test_that("warp_from_field maps coordinates for deformation representation", {
  grid <- grid_spec(c(5, 5, 5), diag(c(2, 2, 2, 1)))

  dims <- grid@dims
  ijk <- as.matrix(expand.grid(0:(dims[1] - 1), 0:(dims[2] - 1), 0:(dims[3] - 1)))
  world <- t(grid@affine %*% t(cbind(ijk, 1)))[, 1:3, drop = FALSE]
  field <- array(world, dim = c(dims, 3))
  field[, , , 1] <- field[, , , 1] + 1

  morph <- warp_from_field(
    source = "src",
    target = "tgt",
    field = field,
    grid = grid,
    representation = "deformations"
  )

  coords <- matrix(c(
    0, 0, 0,
    2, 4, 6,
    6, 2, 4
  ), ncol = 3, byrow = TRUE)

  out <- transform(morph, coords)
  expect_equal(out[, 1], coords[, 1] + 1, tolerance = 1e-6)
  expect_equal(out[, 2:3], coords[, 2:3], tolerance = 1e-6)
})

test_that("write_warp_field round-trips displacement fields through NIfTI", {
  grid <- grid_spec(c(5, 5, 5), diag(c(2, 2, 2, 1)))
  field <- array(0, dim = c(5, 5, 5, 3))
  field[, , , 1] <- 1
  morph <- warp_from_field("src", "tgt", field, grid = grid, representation = "displacements")

  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))
  write_warp_field(morph, tmp, representation = "displacements")

  loaded <- read_transform(tmp, type = "ants", source = "src", target = "tgt")
  coords <- matrix(c(0, 0, 0, 2, 4, 6), ncol = 3, byrow = TRUE)
  expect_equal(transform(loaded, coords), transform(morph, coords), tolerance = 1e-5)
})

test_that("write_transform dispatches Warp3DMorphism exports", {
  skip_if_not_installed("neuroim2")

  grid <- grid_spec(c(4, 4, 4), diag(c(2, 2, 2, 1)))
  field <- array(0, dim = c(4, 4, 4, 3))
  field[, , , 1] <- 1
  morph <- warp_from_field("src", "tgt", field, grid = grid, representation = "displacements")

  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))
  expect_invisible(write_transform(morph, tmp, type = "deformations"))

  vec <- neuroim2::read_vec(tmp)
  arr <- as.array(vec)
  expect_equal(dim(arr), c(4, 4, 4, 3))
  # First voxel (0,0,0) deformation should be source coordinate (1,0,0).
  expect_equal(arr[1, 1, 1, 1], 1, tolerance = 1e-6)
  expect_equal(arr[1, 1, 1, 2], 0, tolerance = 1e-6)
  expect_equal(arr[1, 1, 1, 3], 0, tolerance = 1e-6)
})
