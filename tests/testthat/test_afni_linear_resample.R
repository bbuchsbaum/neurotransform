test_that("AFNI affine resample matches AFNI output (high correlation)", {
  src_path <- system.file("extdata/afni/ss_sub-1001_T1w.nii.gz", package = "neurotransform")
  tgt_path <- system.file("extdata/afni/mni.nii", package = "neurotransform")
  aff_path <- system.file("extdata/afni/sub-1001_T1w_in_mni_Allin.aff12.1D", package = "neurotransform")
  ref_path <- system.file("extdata/afni/sub-1001_T1w_in_mni_Allin.nii", package = "neurotransform")

  skip_if_not(file.exists(src_path))
  skip_if_not(file.exists(tgt_path))
  skip_if_not(file.exists(aff_path))
  skip_if_not(file.exists(ref_path))

  src <- read_image(src_path)
  tgt <- read_image(tgt_path)
  ref <- read_image(ref_path)

  # 3dAllineate -1Dmatrix_save is base/target -> source/input in DICOM coords.
  morph <- afni_load_affine_morphism("src", "tgt", aff_path)

  out <- resample_to(src, target = tgt, transform = morph, method = "linear")

  out_arr <- as.array(out)
  ref_arr <- as.array(ref)
  if (length(dim(out_arr)) == 4) out_arr <- out_arr[, , , 1, drop = TRUE]
  if (length(dim(ref_arr)) == 4) ref_arr <- ref_arr[, , , 1, drop = TRUE]

  mask <- is.finite(out_arr) & is.finite(ref_arr)
  expect_gt(mean(mask), 0.97)
  expect_gt(cor(out_arr[mask], ref_arr[mask]), 0.995)
  expect_lt(mean(abs(out_arr[mask] - ref_arr[mask])), 35)
})
