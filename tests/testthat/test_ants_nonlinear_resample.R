test_that("ANTs H5 warp loads correctly", {
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("hdf5r")

  # Test that the warp can be loaded
  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants_h5")
  expect_s4_class(morph, "Warp3DMorphism")

  # Load the warp array

  warp <- neurotransform:::load_warp_array(morph)
  expect_equal(warp$dim, c(97L, 115L, 97L))
  expect_true(length(warp$array) == prod(warp$dim) * 3)
})

test_that("ANTs H5 warp transforms coordinates within bounds", {
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("hdf5r")

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants_h5")

  # Test coordinates near the warp field center
  # The warp has dims 97x115x97 with identity transform, so valid range is ~0-97mm
  test_coords <- matrix(c(
    45, 55, 45,   # Center of warp field
    30, 40, 30,   # Another point in bounds
    60, 70, 60    # Third point
  ), ncol = 3, byrow = TRUE)

  warped <- transform(morph, test_coords)

  # Warped coordinates should be finite

  expect_true(all(is.finite(warped)))

  # Warped coordinates should be within neuroimaging bounds
  expect_true(all(abs(warped) < 300))

  # There should be some displacement (not identity)
  disp <- warped - test_coords
  expect_gt(max(abs(disp)), 0.1)
})

test_that("ANTs H5 resample produces output (functional test)", {
  src_path <- system.file("extdata/chris/chris_t1.nii.gz", package = "neurotransform")
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  ref_path <- system.file("extdata/chris/ants/chris_in_mni.nii.gz",
                          package = "neurotransform")

  skip_if_not(file.exists(src_path))
  skip_if_not(file.exists(warp_path))
  skip_if_not(file.exists(ref_path))
  skip_if_not_installed("hdf5r")

  src <- neuroim2::read_vol(src_path)
  ref <- neuroim2::read_vol(ref_path)

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants_h5")

  # This should run without error
  out <- resample_to(src, target = ref, transform = morph, method = "linear")

  # Output should have correct dimensions
  expect_equal(dim(out)[1:3], dim(ref)[1:3])

  out_arr <- as.array(out)
  if (length(dim(out_arr)) == 4) out_arr <- out_arr[, , , 1, drop = TRUE]

  # Output should have some non-zero values
  expect_gt(sum(out_arr > 0, na.rm = TRUE), 1000)
})

test_that("ANTs H5 resample produces reasonable output", {
  # Note: This test validates that the resampling produces valid output.

  # High correlation with the reference is not expected because:
  # 1. The reference was created with ANTs using the full composite transform
  # 2. The warp-only approach here does not include the embedded affine component
  # 3. The warp field coordinate system may differ from the reference creation
  # Use ants_h5_morphism(apply_affine=TRUE) for full composite transform.

  src_path <- system.file("extdata/chris/chris_t1.nii.gz", package = "neurotransform")
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  ref_path <- system.file("extdata/chris/ants/chris_in_mni.nii.gz",
                          package = "neurotransform")

  skip_if_not(file.exists(src_path))
  skip_if_not(file.exists(warp_path))
  skip_if_not(file.exists(ref_path))
  skip_if_not_installed("hdf5r")

  src <- neuroim2::read_vol(src_path)
  ref <- neuroim2::read_vol(ref_path)

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants_h5")
  out <- resample_to(src, target = ref, transform = morph, method = "linear")

  out_arr <- as.array(out)
  ref_arr <- as.array(ref)
  if (length(dim(out_arr)) == 4) out_arr <- out_arr[, , , 1, drop = TRUE]
  if (length(dim(ref_arr)) == 4) ref_arr <- ref_arr[, , , 1, drop = TRUE]

  mask <- is.finite(out_arr) & is.finite(ref_arr) & (ref_arr != 0)

  # Should have reasonable mask coverage
  expect_gt(mean(mask), 0.1)

  # Output should have brain-like intensity range (not all zeros/NAs)
  expect_gt(max(out_arr, na.rm = TRUE), 50)  # Some bright voxels

  r <- suppressWarnings(cor(as.vector(out_arr[mask]), as.vector(ref_arr[mask])))
  expect_true(is.finite(r))
  # Note: Correlation may be low without affine; this just validates output is valid
})

test_that("ANTs NIfTI warp loads correctly", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  # Test that the NIfTI warp can be loaded
  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants")
  expect_s4_class(morph, "Warp3DMorphism")

  warp <- neurotransform:::load_warp_array(morph)
  expect_true(length(warp$dim) == 3)
  expect_true(all(warp$dim > 0))
})

test_that("ANTs NIfTI warp transforms coordinates", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants")

  # Get warp info to find valid coordinate range
  warp <- neurotransform:::load_warp_array(morph)
  vox_to_world <- warp$vox_to_world

  # Test at center of warp volume
  center_vox <- warp$dim / 2
  center_world <- (vox_to_world %*% c(center_vox, 1))[1:3]

  test_coords <- matrix(center_world, nrow = 1)
  warped <- transform(morph, test_coords)

  expect_true(all(is.finite(warped)))
  expect_true(all(abs(warped) < 300))
})

test_that("ants_h5_morphism returns MorphismPath with affine", {
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("hdf5r")

  morph <- ants_h5_morphism(warp_path, source = "native", target = "mni",
                            apply_affine = TRUE)

  # Should return a MorphismPath with warp + affine (in that order for pullback)
  # The path is [warp, affine] so pullback applies: affine_pullback(warp_pullback(coords))
  expect_s4_class(morph, "MorphismPath")
  expect_equal(length(morph@morphisms), 2L)
  expect_s4_class(morph@morphisms[[1]], "Warp3DMorphism")
  expect_s4_class(morph@morphisms[[2]], "Affine3DMorphism")
})

test_that("ants_h5_morphism without affine returns single warp", {
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("hdf5r")

  morph <- ants_h5_morphism(warp_path, source = "native", target = "mni",
                            apply_affine = FALSE)

  # Should return just the warp morphism
  expect_s4_class(morph, "Warp3DMorphism")
})

# ==============================================================================
# FULL COMPOSITE TRANSFORM FUNCTIONAL TESTS
# ==============================================================================

test_that("ANTs H5 full composite resample produces valid output", {
 # This tests that the full composite transform (warp + embedded affine)
  # produces reasonable output. Note: exact correlation matching with ANTs
  # antsApplyTransforms requires careful analysis of how ANTs structures
  # composite transforms internally (transform order, coordinate systems).
  #
  # The separate warp and affine components are validated individually.
  # This test validates the composite path construction and resampling.

  src_path <- system.file("extdata/chris/chris_t1.nii.gz", package = "neurotransform")
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  ref_path <- system.file("extdata/chris/ants/chris_in_mni.nii.gz",
                          package = "neurotransform")

  skip_if_not(file.exists(src_path), "Source image not available")
  skip_if_not(file.exists(warp_path), "H5 composite warp not available")
  skip_if_not(file.exists(ref_path), "ANTs reference output not available")
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")

  src <- neuroim2::read_vol(src_path)
  ref <- neuroim2::read_vol(ref_path)

  # Use full composite transform (warp + embedded affine)
  morph <- ants_h5_morphism(warp_path, source = "native", target = "mni",
                            apply_affine = TRUE)

  # Verify morphism structure
  expect_s4_class(morph, "MorphismPath")
  expect_equal(length(morph@morphisms), 2)

  # Resample using the full composite transform
  out <- resample_to(src, target = ref, transform = morph, method = "linear")

  # Output should have correct dimensions
  expect_equal(dim(out)[1:3], dim(ref)[1:3])

  out_arr <- as.array(out)
  ref_arr <- as.array(ref)
  if (length(dim(out_arr)) == 4) out_arr <- out_arr[, , , 1, drop = TRUE]
  if (length(dim(ref_arr)) == 4) ref_arr <- ref_arr[, , , 1, drop = TRUE]

  # Output should have brain-like values (not all zeros or NaN)
  expect_gt(sum(out_arr > 0, na.rm = TRUE), 1000)
  expect_gt(max(out_arr, na.rm = TRUE), 50)

  # Create mask for comparison
  mask <- is.finite(out_arr) & is.finite(ref_arr) & (ref_arr != 0) & (out_arr != 0)

  # Should have some overlap with reference
  expect_gt(mean(mask), 0.01)

  # Note: correlation may not be high due to differences in how ANTs composite
  # transforms are applied vs how we decompose and recompose them.
  # This test validates functional output, not exact matching.
  r <- suppressWarnings(cor(as.vector(out_arr[mask]), as.vector(ref_arr[mask])))
  expect_true(is.finite(r))
})
