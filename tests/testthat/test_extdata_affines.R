test_that("ANTs sample affine matches expected parameters", {
  path <- system.file("extdata/ants/sample_ANTs_0GenericAffine.mat", package = "neurotransform")
  lines <- readLines(path, warn = FALSE)
  params_line <- grep("^Parameters", lines, value = TRUE)
  stopifnot(length(params_line) == 1)
  vals <- scan(text = sub("Parameters:\\s*", "", params_line), quiet = TRUE)
  expect_equal(length(vals), 12)
  mat <- matrix(c(vals[1:3], 0,
                  vals[4:6], 0,
                  vals[7:9], 0,
                  vals[10:12], 1), nrow = 4, byrow = TRUE)
  expect_equal(mat[1, 1], 0.9961947, tolerance = 1e-6)
  expect_equal(mat[2, 1], 0.0871557, tolerance = 1e-6)
  expect_equal(mat[4, 1:3], c(1, -2, 0.5), tolerance = 1e-6)
})

test_that("FSL sample linear matrices match fixture values", {
  mat6_path <- system.file("extdata/fsl/S01_lin_6dof.mat", package = "neurotransform")
  mat12_path <- system.file("extdata/fsl/S01_lin_12dof.mat", package = "neurotransform")

  mat6 <- as.matrix(read.table(mat6_path, col.names = paste0("V", 1:4)))
  mat12 <- as.matrix(read.table(mat12_path, col.names = paste0("V", 1:4)))

  expect_equal(dim(mat6), c(4L, 4L))
  expect_equal(dim(mat12), c(4L, 4L))

  expect_equal(unname(mat6[1, 4]), 2, tolerance = 1e-8)
  expect_equal(unname(mat6[2, 4]), -3, tolerance = 1e-8)
  expect_equal(unname(mat6[3, 4]), 1.5, tolerance = 1e-8)

  expect_equal(unname(mat12[1, 1]), 1.02, tolerance = 1e-8)
  expect_equal(unname(mat12[2, 2]), 0.99, tolerance = 1e-8)
})

test_that("load_warp_array works on sample ANTs and FSL warps", {
  ants_path <- system.file("extdata/ants/sample_ANTs_1Warp.nii.gz", package = "neurotransform")
  fsl_path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")

  ants <- load_warp_array(Warp3DMorphism("src", "tgt", ants_path, warp_type = "ants"))
  fsl <- load_warp_array(Warp3DMorphism("src", "tgt", fsl_path, warp_type = "fsl"))

  expect_equal(ants$dim, c(3L, 3L, 3L))
  expect_equal(fsl$dim, c(3L, 3L, 3L))

  arr_ants <- array(ants$array, dim = c(3, ants$dim))
  arr_fsl <- array(fsl$array, dim = c(3, fsl$dim))

  expect_true(any(arr_ants != 0))
  expect_true(any(arr_fsl != 0))
})
