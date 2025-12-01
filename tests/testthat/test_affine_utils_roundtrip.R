test_that("build_affine_matrix invertibility and compose round-trip", {
  mat <- build_affine_matrix(
    translation = c(1, -2, 3),
    scales = c(1.1, 0.9, 1.05),
    skews = c(0.1, -0.05, 0.02),
    angles = c(0.1, -0.05, 0.2)
  )
  inv <- invert_affine(mat)
  pts <- matrix(c(1, 2, 3,
                  -1, 0, 0), ncol = 3, byrow = TRUE)
  back <- apply_affine(apply_affine(pts, mat), inv)
  expect_equal(back, pts, tolerance = 1e-8)

  comp <- compose_affines(mat, inv)
  expect_true(is_affine_matrix(comp))
  expect_equal(comp, diag(4), tolerance = 1e-8)
})

test_that("decompose_affine_matrix recovers components approximately", {
  mat <- build_affine_matrix(
    translation = c(2, 3, -1),
    scales = c(1.2, 0.8, 1.1),
    skews = c(0.05, 0.02, -0.01),
    angles = c(0.05, 0.1, -0.02)
  )
  parts <- decompose_affine_matrix(mat)
  expect_equal(parts$translation, c(2, 3, -1), tolerance = 1e-6)
  expect_equal(parts$scales, c(1.2, 0.8, 1.1), tolerance = 1e-2)
  expect_equal(parts$skews, c(0.05, 0.02, -0.01), tolerance = 3e-2)
})

test_that("read/write affine matrix txt round-trip", {
  mat <- build_affine_matrix(translation = c(1, 2, 3))
  tmp <- tempfile(fileext = ".txt")
  write_affine_matrix_txt(mat, tmp, type = "generic")
  res <- read_affine_matrix_txt(tmp)
  expect_equal(res$matrix, mat, tolerance = 1e-8)
  expect_equal(res$type, "generic")
})

test_that("convert_affine_convention is invertible for identity voxel affines", {
  C <- diag(4)
  src_aff <- diag(4)
  tgt_aff <- diag(4)
  fsl <- convert_affine_convention(C, src_aff, tgt_aff, from = "generic", to = "fsl")
  back <- convert_affine_convention(fsl, src_aff, tgt_aff, from = "fsl", to = "generic")
  expect_equal(back, C, tolerance = 1e-8)
})
