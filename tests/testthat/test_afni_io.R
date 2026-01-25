test_that("afni_read_aff12 loads AFNI affine and converts to RAS", {
  path <- system.file("extdata/afni/anatQQ.FT.aff12.1D", package = "neurotransform")
  skip_if_not(file.exists(path))
  mat <- afni_read_aff12(path)
  expect_equal(dim(mat), c(4, 4))

  ras <- afni_aff12_to_ras(mat)
  expect_equal(dim(ras), c(4, 4))
  expect_true(all(is.finite(ras)))
})
