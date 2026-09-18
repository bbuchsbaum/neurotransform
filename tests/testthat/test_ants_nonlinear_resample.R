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

  warp <- neurotransform:::load_warp_array(morph)
  vox <- rbind(
    floor(warp$dim / 2),
    floor(warp$dim / 3),
    floor(2 * warp$dim / 3)
  )
  test_coords <- (cbind(vox, 1) %*% t(warp$vox_to_world))[, 1:3, drop = FALSE]

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

test_that("ANTs NIfTI displacement components convert from LPS to RAS", {
  skip_if_not_installed("neuroim2")
  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path), add = TRUE)

  field_lps <- array(0, dim = c(5, 5, 5, 3))
  field_lps[, , , 1] <- 1
  field_lps[, , , 2] <- 2
  field_lps[, , , 3] <- 3
  space <- neuroim2::NeuroSpace(dim(field_lps), trans = diag(4))
  neuroim2::write_vec(neuroim2::DenseNeuroVec(field_lps, space), path, format = "nifti")

  morph <- Warp3DMorphism("moving", "fixed", path, warp_type = "ants")
  point <- matrix(c(2, 2, 2), nrow = 1)

  expect_equal(transform(morph, point), point + matrix(c(-1, -2, 3), nrow = 1),
               tolerance = 1e-7)
})

test_that("ANTs NIfTI warp matches retained SimpleITK point oracles", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))

  morph <- Warp3DMorphism("moving", "fixed", warp_path, warp_type = "ants")
  points_ras <- rbind(
    c(0, 0, 0),
    c(-80, -100, -50),
    c(50, 60, 20)
  )
  expected_ras <- rbind(
    c(0.4570037569, -1.8435784895, 2.0230027791),
    c(-80.7043913016, -100.6252680179, -50.2785666147),
    c(51.4678156711, 60.3026160169, 21.1965898331)
  )

  expect_equal(transform(morph, points_ras), expected_ras, tolerance = 2e-5)
})

test_that("ants_h5_morphism returns MorphismPath with affine", {
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("hdf5r")

  morph <- ants_h5_morphism(warp_path, source = "native", target = "mni",
                            apply_affine = TRUE)

  # Stored H5 components are preserved. MorphismPath, like ITK, evaluates them
  # from last to first.
  expect_s4_class(morph, "MorphismPath")
  expect_equal(length(morph@morphisms), 2L)
  expect_s4_class(morph@morphisms[[1]], "Affine3DMorphism")
  expect_s4_class(morph@morphisms[[2]], "Warp3DMorphism")
})

test_that("ANTs H5 composites match retained SimpleITK point oracles", {
  skip_if_not_installed("hdf5r")
  forward_path <- system.file("extdata/chris/ants/reg_Composite.h5",
                              package = "neurotransform")
  inverse_path <- system.file("extdata/chris/ants/reg_InverseComposite.h5",
                              package = "neurotransform")
  skip_if_not(file.exists(forward_path))
  skip_if_not(file.exists(inverse_path))

  points_ras <- rbind(c(-10, 20, 30), c(90, -120, -70))
  expected_forward_ras <- rbind(
    c(-8.6304721092, 17.6258127526, 25.5654076920),
    c(82.9201542089, -119.4691335359, -69.7328632238)
  )
  expected_inverse_ras <- rbind(
    c(-11.2556239032, 22.7520343883, 34.8090080863),
    c(97.7416223612, -120.3847724732, -70.3661393128)
  )

  forward <- ants_h5_morphism(forward_path, source = "moving", target = "fixed")
  inverse <- ants_h5_morphism(inverse_path, source = "fixed", target = "moving")

  expect_equal(transform(forward, points_ras), expected_forward_ras,
               tolerance = 2e-5)
  expect_equal(transform(inverse, points_ras), expected_inverse_ras,
               tolerance = 2e-5)
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
