test_that("extract_affine works for neuroim2 volumes", {
  skip_if_not_installed("neuroim2")

  vol_path <- testthat::test_path("..", "..", "inst", "extdata", "template",
                                  "MNI152NLin2009cAsym_T1w_2mm_brain.nii.gz")
  skip_if_not(file.exists(vol_path), "Template file missing")

  vol <- neuroim2::read_vol(vol_path)
  expect_s4_class(vol, "DenseNeuroVol")

  aff_expected <- neuroim2::trans(vol)
  aff_got <- extract_affine(vol)

  expect_equal(aff_got, aff_expected)
})
