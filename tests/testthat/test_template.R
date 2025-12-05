test_that("MNI template exists and is readable", {
  path <- system.file("extdata/template/MNI152NLin2009cAsym_T1w_2mm_brain.nii.gz",
                      package = "neurotransform")
  expect_true(file.exists(path))
  skip_if_not_installed("neuroim2")
  img <- neuroim2::read_vol(path)
  expect_equal(dim(img), c(97, 115, 97))
})
