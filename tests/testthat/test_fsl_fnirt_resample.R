# FSL FNIRT Nonlinear Warp Resampling Tests
#
# These tests validate FSL FNIRT warp handling against FSL reference outputs.
# Similar to test_afni_nonlinear_resample.R, these tests ensure that:
# 1. FNIRT warps can be loaded and applied correctly
# 2. Coordinate transformations produce reasonable results
# 3. Resampled volumes correlate with FSL's applywarp output
#
# Test data generation:
#   Run inst/extdata/fsl/register_to_mni.sh to generate real FNIRT warps
#   (requires FSL to be installed)

# ==============================================================================
# FNIRT WARP LOADING TESTS
# ==============================================================================

test_that("FSL FNIRT warp loads correctly", {
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available (run register_to_mni.sh)")

  # Test that the warp can be loaded

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "fsl")
  expect_s4_class(morph, "Warp3DMorphism")
  expect_equal(morph@warp_type, "fsl")

  # Load the warp array and verify structure
  warp <- neurotransform:::load_warp_array(morph)
  expect_true(length(warp$dim) == 3)
  expect_true(all(warp$dim > 0))
  expect_true(length(warp$array) == prod(warp$dim) * 3)
})

test_that("FSL FNIRT inverse warp loads correctly", {
  inv_warp_path <- system.file("extdata/fsl/standard2highres_warp.nii.gz",
                               package = "neurotransform")
  skip_if_not(file.exists(inv_warp_path), "Inverse FNIRT warp not available")

  morph <- Warp3DMorphism("mni", "native", warp_path = inv_warp_path, warp_type = "fsl")
  expect_s4_class(morph, "Warp3DMorphism")

  warp <- neurotransform:::load_warp_array(morph)
  expect_true(all(warp$dim > 0))
})

# ==============================================================================
# COORDINATE TRANSFORMATION TESTS
# ==============================================================================

test_that("FSL FNIRT warp transforms coordinates within expected range", {
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  inv_warp_path <- system.file("extdata/fsl/standard2highres_warp.nii.gz",
                               package = "neurotransform")

  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available")
  skip_if_not(file.exists(inv_warp_path), "Inverse FNIRT warp not available")

  # Create morphisms
  fwd_morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "fsl")
  inv_morph <- Warp3DMorphism("mni", "native", warp_path = inv_warp_path, warp_type = "fsl")

  # Test coords in MNI space (where the forward warp is defined)
  test_coords <- matrix(c(
     0,  0,  0,   # MNI origin
    10, 20, 30,   # Arbitrary point
   -20, 40, 50    # Another point
  ), ncol = 3, byrow = TRUE)

  # Both warps should transform coordinates

warped_fwd <- transform(fwd_morph, test_coords)
  warped_inv <- transform(inv_morph, test_coords)

  # Warped coordinates should be finite
  expect_true(all(is.finite(warped_fwd)))
  expect_true(all(is.finite(warped_inv)))

  # Warped coordinates should be within reasonable neuroimaging bounds (±200mm)
  expect_true(all(abs(warped_fwd) < 200))
  expect_true(all(abs(warped_inv) < 200))

  # Displacements should be non-trivial (warps are not identity)
  disp_fwd <- warped_fwd - test_coords
  disp_inv <- warped_inv - test_coords
  expect_gt(max(abs(disp_fwd)), 0.5)  # At least 0.5mm displacement somewhere
  expect_gt(max(abs(disp_inv)), 0.5)
})

test_that("FSL FNIRT forward/inverse warp round-trip is approximately correct", {
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  inv_warp_path <- system.file("extdata/fsl/standard2highres_warp.nii.gz",
                               package = "neurotransform")

  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available")
  skip_if_not(file.exists(inv_warp_path), "Inverse FNIRT warp not available")

  fwd_morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "fsl")
  inv_morph <- Warp3DMorphism("mni", "native", warp_path = inv_warp_path, warp_type = "fsl")

  # Test coordinates near center of brain in MNI space
  test_coords <- matrix(c(
     0,  0,  20,  # Near AC
    20, -20, 40,  # Parietal region
   -30,  10, 10   # Left frontal
  ), ncol = 3, byrow = TRUE)

  # Forward then inverse should approximately recover original
  warped <- transform(fwd_morph, test_coords)
  recovered <- transform(inv_morph, warped)

  # Should be close to original (within a few mm due to interpolation)
  max_error <- max(abs(recovered - test_coords))
  expect_lt(max_error, 5)  # Within 5mm
})

# ==============================================================================
# DEFORMATION TYPE DETECTION TESTS
# ==============================================================================

test_that("detect_fnirt_def_type correctly identifies real FNIRT warp as relative", {
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available")
  skip_if_not_installed("neuroim2")

  # Real FNIRT warps typically store relative displacements
  def_type <- detect_fnirt_def_type(warp_path, sample_n = 200, threshold_mm = 50)

  # FNIRT --fout produces relative displacement fields
  expect_equal(def_type, "relative")
})

# ==============================================================================
# VOLUME RESAMPLING VALIDATION TESTS
# ==============================================================================

test_that("FSL FNIRT resample produces output (functional test)", {
  src_path <- system.file("extdata/afni/ss_sub-1001_T1w.nii.gz", package = "neurotransform")
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  ref_path <- system.file("extdata/fsl/highres_in_mni.nii.gz",
                          package = "neurotransform")

  skip_if_not(file.exists(src_path), "Source image not available")
  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available")
  skip_if_not(file.exists(ref_path), "Reference output not available")
  skip_if_not_installed("neuroim2")

  src <- read_image(src_path)
  ref <- read_image(ref_path)

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "fsl")

  # This should run without error
  out <- resample_to(src, target = ref, transform = morph, method = "linear")

  # Output should have correct dimensions
  expect_equal(dim(out)[1:3], dim(ref)[1:3])

  out_arr <- as.array(out)
  if (length(dim(out_arr)) == 4) out_arr <- out_arr[, , , 1, drop = TRUE]

  # Output should have some non-zero values
  expect_gt(sum(out_arr > 0, na.rm = TRUE), 1000)
})

test_that("FSL FNIRT resample correlates with FSL applywarp reference", {
  # This is the key validation test: our resampling should match FSL's output

  src_path <- system.file("extdata/afni/ss_sub-1001_T1w.nii.gz", package = "neurotransform")
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  ref_path <- system.file("extdata/fsl/highres_in_mni_applywarp.nii.gz",
                          package = "neurotransform")

  skip_if_not(file.exists(src_path), "Source image not available")
  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available")
  skip_if_not(file.exists(ref_path), "FSL applywarp reference not available")

  src <- read_image(src_path)
  ref <- read_image(ref_path)

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "fsl")
  out <- resample_to(src, target = ref, transform = morph, method = "linear")

  out_arr <- as.array(out)
  ref_arr <- as.array(ref)
  if (length(dim(out_arr)) == 4) out_arr <- out_arr[, , , 1, drop = TRUE]
  if (length(dim(ref_arr)) == 4) ref_arr <- ref_arr[, , , 1, drop = TRUE]

  # Create mask of valid voxels
  mask <- is.finite(out_arr) & is.finite(ref_arr) & (ref_arr != 0)

  # Should have reasonable mask coverage (at least 10% of brain)
  expect_gt(mean(mask), 0.1)

  # Correlation with FSL reference output should be high
  r <- suppressWarnings(cor(as.vector(out_arr[mask]), as.vector(ref_arr[mask])))

  # Expect high correlation (similar to AFNI test threshold)
  expect_gt(r, 0.8)
})

# ==============================================================================
# FLIRT + FNIRT COMBINED PATH TESTS
# ==============================================================================

test_that("FSL FLIRT affine loads and transforms correctly", {
  flirt_path <- system.file("extdata/fsl/highres2standard.mat",
                            package = "neurotransform")
  src_path <- system.file("extdata/afni/ss_sub-1001_T1w.nii.gz",
                          package = "neurotransform")

  skip_if_not(file.exists(flirt_path), "Real FLIRT matrix not available")
  skip_if_not(file.exists(src_path), "Source image not available")

  # Read the FLIRT matrix and source image affine
  flirt_mat <- as.matrix(read.table(flirt_path))
  expect_equal(dim(flirt_mat), c(4, 4))

  # FLIRT matrix should be non-identity
  expect_false(identical(flirt_mat, diag(4)))

  # Should be invertible (not singular)
  expect_gt(abs(det(flirt_mat)), 0.01)
})

test_that("FSL FLIRT + FNIRT path produces valid composite transform", {
  flirt_path <- system.file("extdata/fsl/highres2standard.mat",
                            package = "neurotransform")
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  src_path <- system.file("extdata/afni/ss_sub-1001_T1w.nii.gz",
                          package = "neurotransform")

  skip_if_not(file.exists(flirt_path), "Real FLIRT matrix not available")
  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available")
  skip_if_not(file.exists(src_path), "Source image not available")
  skip_if_not_installed("neuroim2")

  # Load source to get affine
  src <- neuroim2::read_vol(src_path)
  src_aff <- neuroim2::trans(src)

  # MNI 2mm template affine (standard)
  mni_aff <- diag(c(2, 2, 2, 1))
  mni_aff[1:3, 4] <- c(-90, -126, -72)  # Standard MNI origin

  # Create FLIRT morphism
  flirt_mat <- as.matrix(read.table(flirt_path))
  internal_aff <- fsl_flirt_to_internal_affine(flirt_mat, src_aff, mni_aff)
  aff_morph <- Affine3DMorphism("native", "mni_lin", internal_aff)

  # Create FNIRT morphism
  warp_morph <- Warp3DMorphism("mni_lin", "mni", warp_path = warp_path, warp_type = "fsl")

  # Compose into path (FLIRT then FNIRT)
  path <- compose(aff_morph, warp_morph)

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)

  # Test coordinate transform through the full path
  test_coords <- matrix(c(0, 0, 0), ncol = 3)
  warped <- transform(path, test_coords)

  expect_true(all(is.finite(warped)))
  expect_true(all(abs(warped) < 200))
})

# ==============================================================================
# JACOBIAN TESTS FOR FNIRT WARPS
# ==============================================================================

test_that("FSL FNIRT warp Jacobian is computed correctly", {
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available")

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "fsl")

  # Test coordinates in the middle of the warp field
  test_coords <- matrix(c(
     0,  0,  20,
    10, 20, 30
  ), ncol = 3, byrow = TRUE)

  # Compute Jacobian matrices
  jac <- jacobian(morph, test_coords, mode = "pullback")

  # Should return a JacobianField object
  expect_s4_class(jac, "JacobianField")

  # Should have correct dimensions
  expect_equal(length(jac), nrow(test_coords))

  # Jacobian determinants should be positive and reasonable
  dets <- det(jac)
  expect_true(all(is.finite(dets)))
  expect_true(all(dets > 0))  # Orientation-preserving
  expect_true(all(dets < 10))  # Not extreme expansion
  expect_true(all(dets > 0.1)) # Not extreme contraction
})

test_that("FSL FNIRT warp Jacobian determinant field matches jacobian_det", {
  warp_path <- system.file("extdata/fsl/highres2standard_warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "Real FNIRT warp not available")

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "fsl")

  test_coords <- matrix(c(0, 0, 20, 10, 20, 30), ncol = 3, byrow = TRUE)

  # Two ways to get Jacobian determinants
  jac <- jacobian(morph, test_coords)
  dets_from_field <- det(jac)

  dets_direct <- jacobian_det(morph, test_coords)

  # Should match
  expect_equal(dets_from_field, dets_direct, tolerance = 1e-6)
})
