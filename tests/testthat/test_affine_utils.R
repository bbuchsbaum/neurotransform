test_that("build and decompose affine round-trip basic components", {
  mat <- build_affine_matrix(
    translation = c(2, -3, 1.5),
    scales = c(1.1, 0.9, 1.05),
    skews = c(0.1, -0.05, 0.02),
    angles = c(0.1, -0.05, 0.2)
  )
  parts <- decompose_affine_matrix(mat)

  expect_equal(parts$translation, c(2, -3, 1.5), tolerance = 1e-6)
  expect_length(parts$scales, 3)
  expect_length(parts$skews, 3)
  expect_length(parts$angles, 3)
})

test_that("affine matrix I/O round trips", {
  mat <- diag(4)
  mat[1, 4] <- 1.23
  tmp <- tempfile(fileext = ".txt")
  write_affine_matrix_txt(mat, tmp, type = "generic")
  res <- read_affine_matrix_txt(tmp)
  expect_equal(res$matrix, mat)
  expect_true(is_affine_matrix(res$matrix))
})

test_that("convert_affine_convention flips FSL vox<->canonical with identity affines", {
  fsl <- diag(4)
  fsl[1, 4] <- 1
  src_aff <- diag(4)
  tgt_aff <- diag(4)

  canon <- convert_affine_convention(fsl, src_aff, tgt_aff, from = "fsl", to = "generic")
  back <- convert_affine_convention(canon, src_aff, tgt_aff, from = "generic", to = "fsl")

  expect_equal(back, fsl, tolerance = 1e-8)
})
