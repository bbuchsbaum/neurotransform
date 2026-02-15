test_that("get_linear_factory resolves aliases and returns handlers", {
  f1 <- get_linear_factory("itk", is_array = FALSE)
  expect_equal(f1$format, "itk")
  expect_true(is.function(f1$reader))
  expect_true(is.function(f1$writer))

  f2 <- get_linear_factory("ants", is_array = FALSE)
  expect_equal(f2$format, "itk")

  f3 <- get_linear_factory("fsl", is_array = TRUE)
  expect_equal(f3$format, "fsl")
  expect_identical(f3$reader, read_linear_transform_array)
  expect_identical(f3$writer, write_linear_transform_array)
})

test_that("get_linear_factory raises typed IO condition on bad format", {
  expect_error(
    get_linear_factory("fakepackage"),
    class = "TransformIOError"
  )
})

test_that("read/write_linear_transform_array supports FSL indexed files", {
  base <- tempfile(fileext = ".mat")
  on.exit(unlink(c(base, sprintf("%s.%03d", base, 0:4))))

  A <- diag(4)
  A[1:3, 4] <- c(1, 2, 3)
  B <- diag(4)
  B[1:3, 4] <- c(-2, 0.5, 1)

  write_linear_transform_array(
    list(A, B),
    path = base,
    format = "fsl",
    source_affine = diag(4),
    target_affine = diag(4)
  )

  expect_true(file.exists(sprintf("%s.%03d", base, 0)))
  expect_true(file.exists(sprintf("%s.%03d", base, 1)))

  arr <- read_linear_transform_array(
    base,
    format = "fsl",
    source = "src",
    target = "tgt",
    source_affine = diag(4),
    target_affine = diag(4)
  )
  expect_s3_class(arr, "LinearTransformArray")
  expect_length(arr$transforms, 2)
  expect_equal(arr$transforms[[1]]@matrix, A, tolerance = 1e-8)
  expect_equal(arr$transforms[[2]]@matrix, B, tolerance = 1e-8)
})

test_that("read_linear_transform_array raises typed file error on missing base", {
  expect_error(
    read_linear_transform_array(tempfile(fileext = ".mat"), format = "fsl"),
    class = "TransformFileError"
  )
})

test_that("read/write_linear_transform supports lta format", {
  tmp <- tempfile(fileext = ".lta")
  on.exit(unlink(tmp))

  A <- diag(4)
  A[1:3, 4] <- c(3, -2, 1)

  write_linear_transform(A, tmp, format = "lta")
  expect_equal(detect_transform_type(tmp), "lta_affine")
  morph <- read_linear_transform(tmp, format = "lta", source = "src", target = "tgt")

  expect_s4_class(morph, "Affine3DMorphism")
  expect_equal(morph@matrix, A, tolerance = 1e-8)
})

test_that("read/write_linear_transform_array supports lta format", {
  tmp <- tempfile(fileext = ".lta")
  on.exit(unlink(tmp))

  A <- diag(4)
  A[1:3, 4] <- c(1, 2, 3)
  B <- diag(4)
  B[1:3, 4] <- c(-4, 0.5, 2)

  write_linear_transform_array(list(A, B), path = tmp, format = "lta")
  arr <- read_linear_transform_array(tmp, format = "lta", source = "src", target = "tgt")

  expect_s3_class(arr, "LinearTransformArray")
  expect_length(arr$transforms, 2)
  expect_equal(arr$transforms[[1]]@matrix, A, tolerance = 1e-8)
  expect_equal(arr$transforms[[2]]@matrix, B, tolerance = 1e-8)
})

test_that("read/write_linear_transform_array supports itk h5 composite format", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp))

  A <- diag(4)
  A[1:3, 4] <- c(1, 2, 3)
  B <- diag(4)
  B[1:3, 4] <- c(-2, 0.5, 1)

  write_linear_transform_array(list(A, B), path = tmp, format = "itk")
  arr <- read_linear_transform_array(tmp, format = "itk", source = "src", target = "tgt")

  expect_s3_class(arr, "LinearTransformArray")
  expect_length(arr$transforms, 2)
  expect_equal(arr$transforms[[1]]@matrix, A, tolerance = 1e-8)
  expect_equal(arr$transforms[[2]]@matrix, B, tolerance = 1e-8)
})

test_that("read_linear_transform itk h5 errors with typed condition for multi-transform file", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp))

  A <- diag(4)
  B <- diag(4)
  B[1:3, 4] <- c(1, 0, 0)

  write_linear_transform_array(list(A, B), path = tmp, format = "itk")
  expect_error(
    read_linear_transform(tmp, format = "itk"),
    class = "TransformIOError"
  )
})

test_that("detect_transform_type distinguishes itk h5 and x5 h5", {
  skip_if_not_installed("hdf5r")

  tmp_itk <- tempfile(fileext = ".h5")
  tmp_x5 <- tempfile(fileext = ".h5")
  on.exit(unlink(c(tmp_itk, tmp_x5)))

  write_linear_transform_array(list(diag(4)), path = tmp_itk, format = "itk")
  write_x5(tmp_x5, list(
    x5_transform(type = "linear", transform = diag(4), dimension_kinds = c("space", "space"))
  ))

  expect_equal(detect_transform_type(tmp_itk), "itk_affine")
  expect_equal(detect_transform_type(tmp_x5), "x5")
})

test_that("read/write_linear_transform supports x5 linear format", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".x5")
  on.exit(unlink(tmp))

  A <- diag(4)
  A[1:3, 4] <- c(2, -1, 4)

  write_linear_transform(A, tmp, format = "x5")
  morph <- read_linear_transform(tmp, format = "x5", source = "src", target = "tgt")

  expect_s4_class(morph, "Affine3DMorphism")
  expect_equal(morph@matrix, A, tolerance = 1e-8)
})

test_that("read/write_linear_transform_array supports x5 linear arrays", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".x5")
  on.exit(unlink(tmp))

  A <- diag(4)
  A[1:3, 4] <- c(1, 2, 3)
  B <- diag(4)
  B[1:3, 4] <- c(-1, 0, 5)

  write_linear_transform_array(list(A, B), path = tmp, format = "x5")
  arr <- read_linear_transform_array(tmp, format = "x5", source = "src", target = "tgt")

  expect_s3_class(arr, "LinearTransformArray")
  expect_length(arr$transforms, 2)
  expect_equal(arr$transforms[[1]]@matrix, A, tolerance = 1e-8)
  expect_equal(arr$transforms[[2]]@matrix, B, tolerance = 1e-8)
})

test_that("read_linear_transform x5 errors with typed condition for multi-transform file", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".x5")
  on.exit(unlink(tmp))

  nodes <- list(
    x5_transform(type = "linear", transform = diag(4), dimension_kinds = c("space", "space")),
    x5_transform(type = "linear", transform = diag(4), dimension_kinds = c("space", "space"))
  )
  write_x5(tmp, nodes)

  expect_error(
    read_linear_transform(tmp, format = "x5"),
    class = "TransformIOError"
  )
})

test_that("read_transform detects and reads x5 linear paths", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".x5")
  on.exit(unlink(tmp))

  A <- diag(4)
  A[1:3, 4] <- c(1, 0, 0)
  B <- diag(4)
  B[2, 4] <- -2

  write_x5(tmp, list(
    x5_transform(type = "linear", transform = A, dimension_kinds = c("space", "space")),
    x5_transform(type = "linear", transform = B, dimension_kinds = c("space", "space"))
  ))

  expect_equal(detect_transform_type(tmp), "x5")

  morph <- read_transform(tmp, source = "src", target = "tgt")
  expect_s4_class(morph, "MorphismPath")
  expect_equal(length(morph@morphisms), 2L)
})

test_that("write_transform/read_transform round-trip x5 nonlinear warp", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".x5")
  on.exit(unlink(tmp))

  grid <- grid_spec(c(4, 4, 4), diag(4))
  field <- array(0, dim = c(4, 4, 4, 3))
  field[, , , 1] <- 1
  warp <- warp_from_field("src", "tgt", field, grid = grid, representation = "displacements")

  write_transform(warp, tmp, type = "x5")
  loaded <- read_transform(tmp, source = "src", target = "tgt")

  coords <- matrix(c(
    0, 0, 0,
    1, 2, 3
  ), byrow = TRUE, ncol = 3)
  expect_equal(transform(loaded, coords), transform(warp, coords), tolerance = 1e-6)
})
