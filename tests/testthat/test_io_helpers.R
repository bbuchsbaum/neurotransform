test_that("read_image/write_image round trip 3D", {
  skip_if_not_installed("neuroim2")
  arr <- array(1:8, dim = c(2, 2, 2))
  space <- neuroim2::NeuroSpace(dim(arr), trans = diag(4))
  vol <- neuroim2::DenseNeuroVol(arr, space)
  tmp <- tempfile(fileext = ".nii.gz")
  write_image(vol, tmp)
  vol2 <- read_image(tmp)
  expect_equal(dim(vol2)[1:3], dim(vol))
  expect_equal(as.vector(as.array(vol2)), as.vector(arr))
})

test_that("resample_to identity matches input", {
  src <- array(0, dim = c(3, 3, 3))
  src[2, 2, 2] <- 1
  tgt <- array(0, dim = c(3, 3, 3))
  morph <- Affine3DMorphism("src", "tgt", diag(4))
  out <- resample_to(src, target = tgt, transform = morph, method = "nearest")
  vals <- if (requireNamespace("neuroim2", quietly = TRUE) &&
              (inherits(out, "DenseNeuroVol") || inherits(out, "DenseNeuroVec"))) {
    as.array(out)
  } else out
  expect_equal(dim(vals)[1:3], dim(src))
  expect_equal(sum(vals), sum(src))
  expect_equal(which.max(vals), which.max(src))
})
