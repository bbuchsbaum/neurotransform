test_that("chris_t1 can be loaded (if neuroim2 available)", {
  path <- system.file("extdata/chris/chris_t1.nii.gz", package = "neurotransform")
  skip_if_not(file.exists(path))
  skip_if_not_installed("neuroim2")
  img <- neuroim2::read_vol(path)
  expect_true(all(dim(img) > 0))
})
