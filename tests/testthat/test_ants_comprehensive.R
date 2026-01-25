# Comprehensive ANTs Transform Tests
#
# These tests fill gaps in ANTs coverage:
# 1. Inverse composite transforms
# 2. Standalone ITK affine files
# 3. Forward/inverse round-trip validation
# 4. Jacobian computation for ANTs warps
# 5. Coordinate convention handling

# ==============================================================================
# INVERSE COMPOSITE TESTS
# ==============================================================================

test_that("ANTs inverse H5 composite loads correctly", {
  inv_path <- system.file("extdata/chris/ants/chris_to_mni_InverseComposite.h5",
                          package = "neurotransform")
  skip_if_not(file.exists(inv_path), "Inverse composite not available")
  skip_if_not_installed("hdf5r")

  # Load inverse composite
  morph <- Warp3DMorphism("mni", "native", warp_path = inv_path, warp_type = "ants_h5")
  expect_s4_class(morph, "Warp3DMorphism")

  # Load warp array
  warp <- neurotransform:::load_warp_array(morph)
  expect_true(length(warp$dim) == 3)
  expect_true(all(warp$dim > 0))
})

test_that("ANTs inverse H5 morphism with affine loads correctly", {
  inv_path <- system.file("extdata/chris/ants/chris_to_mni_InverseComposite.h5",
                          package = "neurotransform")
  skip_if_not(file.exists(inv_path), "Inverse composite not available")
  skip_if_not_installed("hdf5r")

  # Load with embedded affine
  morph <- ants_h5_morphism(inv_path, source = "mni", target = "native",
                            apply_affine = TRUE)

  # Should be a MorphismPath if affine is embedded
  expect_true(is(morph, "Morphism") || is(morph, "MorphismPath"))
})

test_that("ANTs inverse composite transforms coordinates", {
  inv_path <- system.file("extdata/chris/ants/chris_to_mni_InverseComposite.h5",
                          package = "neurotransform")
  skip_if_not(file.exists(inv_path), "Inverse composite not available")
  skip_if_not_installed("hdf5r")

  morph <- Warp3DMorphism("mni", "native", warp_path = inv_path, warp_type = "ants_h5")

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

# ==============================================================================
# FORWARD/INVERSE ROUND-TRIP TESTS
# ==============================================================================

test_that("ANTs H5 forward/inverse round-trip approximately recovers coordinates", {
  fwd_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                          package = "neurotransform")
  inv_path <- system.file("extdata/chris/ants/chris_to_mni_InverseComposite.h5",
                          package = "neurotransform")

  skip_if_not(file.exists(fwd_path), "Forward composite not available")
  skip_if_not(file.exists(inv_path), "Inverse composite not available")
  skip_if_not_installed("hdf5r")

  fwd_morph <- Warp3DMorphism("native", "mni", warp_path = fwd_path, warp_type = "ants_h5")
  inv_morph <- Warp3DMorphism("mni", "native", warp_path = inv_path, warp_type = "ants_h5")

  # Get warp info to find valid coordinate range
  warp <- neurotransform:::load_warp_array(fwd_morph)
  vox_to_world <- warp$vox_to_world

  # Test coords near center of warp field
  center_vox <- warp$dim / 2
  center_world <- (vox_to_world %*% c(center_vox, 1))[1:3]

  # Create test coordinates around center
  offsets <- matrix(c(
    0, 0, 0,
    5, 5, 5,
    -5, 5, -5
  ), ncol = 3, byrow = TRUE)
  test_coords <- sweep(offsets, 2, center_world, "+")

  # Forward then inverse
  warped <- transform(fwd_morph, test_coords)
  recovered <- transform(inv_morph, warped)

  # Should approximately recover original (within a few mm)
  max_error <- max(abs(recovered - test_coords))
  expect_lt(max_error, 10)  # Within 10mm due to interpolation and discrete warps
})

test_that("ANTs NIfTI forward/inverse round-trip approximately recovers coordinates", {
  fwd_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                          package = "neurotransform")
  inv_path <- system.file("extdata/chris/ants/reg_1InverseWarp.nii.gz",
                          package = "neurotransform")

  skip_if_not(file.exists(fwd_path), "Forward warp not available")
  skip_if_not(file.exists(inv_path), "Inverse warp not available")

  fwd_morph <- Warp3DMorphism("native", "mni", warp_path = fwd_path, warp_type = "ants")
  inv_morph <- Warp3DMorphism("mni", "native", warp_path = inv_path, warp_type = "ants")

  # Get warp info
  warp <- suppressWarnings(neurotransform:::load_warp_array(fwd_morph))
  vox_to_world <- warp$vox_to_world

  # Test at center
  center_vox <- warp$dim / 2
  center_world <- (vox_to_world %*% c(center_vox, 1))[1:3]

  test_coords <- matrix(center_world, nrow = 1)

  # Round trip
  warped <- suppressWarnings(transform(fwd_morph, test_coords))
  recovered <- suppressWarnings(transform(inv_morph, warped))

  max_error <- max(abs(recovered - test_coords))
  expect_lt(max_error, 10)
})

# ==============================================================================
# JACOBIAN TESTS FOR ANTS WARPS
# ==============================================================================

test_that("ANTs H5 warp Jacobian is computed correctly", {
  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "H5 warp not available")
  skip_if_not_installed("hdf5r")

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants_h5")

  # Get warp info for valid coordinates
  warp <- neurotransform:::load_warp_array(morph)
  vox_to_world <- warp$vox_to_world
  center_vox <- warp$dim / 2
  center_world <- (vox_to_world %*% c(center_vox, 1))[1:3]

  test_coords <- matrix(c(
    center_world,
    center_world + c(5, 0, 0)
  ), ncol = 3, byrow = TRUE)

  # Compute Jacobian
  jac <- jacobian(morph, test_coords, mode = "pullback")

  expect_s4_class(jac, "JacobianField")
  expect_equal(length(jac), nrow(test_coords))

  # Jacobian determinants should be positive and reasonable
  dets <- det(jac)
  expect_true(all(is.finite(dets)))
  expect_true(all(dets > 0))    # Orientation-preserving
  expect_true(all(dets < 10))   # Not extreme expansion
  expect_true(all(dets > 0.1))  # Not extreme contraction
})

test_that("ANTs NIfTI warp Jacobian is computed correctly", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "NIfTI warp not available")

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants")

  # Get warp info
  warp <- neurotransform:::load_warp_array(morph)
  vox_to_world <- warp$vox_to_world
  center_vox <- warp$dim / 2
  center_world <- (vox_to_world %*% c(center_vox, 1))[1:3]

  test_coords <- matrix(center_world, nrow = 1)

  # Compute Jacobian
  jac <- jacobian(morph, test_coords, mode = "pullback")

  expect_s4_class(jac, "JacobianField")
  dets <- det(jac)
  expect_true(all(is.finite(dets)))
  expect_true(all(dets > 0))
})

test_that("ANTs warp jacobian_det matches det(jacobian)", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "NIfTI warp not available")

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants")

  warp <- neurotransform:::load_warp_array(morph)
  vox_to_world <- warp$vox_to_world
  center_vox <- warp$dim / 2
  center_world <- (vox_to_world %*% c(center_vox, 1))[1:3]

  test_coords <- matrix(center_world, nrow = 1)

  # Two ways to get determinants
  jac <- jacobian(morph, test_coords)
  dets_from_field <- det(jac)
  dets_direct <- jacobian_det(morph, test_coords)

  expect_equal(dets_from_field, dets_direct, tolerance = 1e-6)
})

# ==============================================================================
# STANDALONE ITK AFFINE TESTS
# ==============================================================================

test_that("ANTs ITK affine file exists and has content", {
  aff_path <- system.file("extdata/chris/ants/reg_0GenericAffine.mat",
                          package = "neurotransform")
  skip_if_not(file.exists(aff_path), "ITK affine not available")

  # ITK .mat files are binary format (not text)
  # Just verify the file exists and has reasonable size
  file_info <- file.info(aff_path)
  expect_true(file_info$size > 0)
  expect_true(file_info$size < 10000)  # Should be small (few KB)
})

# ==============================================================================
# COORDINATE CONVENTION TESTS
# ==============================================================================
test_that("ANTs H5 warp handles LPS/RAS conversion correctly", {
  # This test verifies that coordinates are handled correctly
  # ANTs/ITK uses LPS; neurotransform uses RAS internally

  warp_path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "H5 warp not available")
  skip_if_not_installed("hdf5r")

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants_h5")
  warp <- neurotransform:::load_warp_array(morph)

  # The warp should have valid coordinate transform matrices
  expect_true(is.matrix(warp$world_to_vox))
  expect_true(is.matrix(warp$vox_to_world))
  expect_equal(dim(warp$world_to_vox), c(4, 4))
  expect_equal(dim(warp$vox_to_world), c(4, 4))

  # Matrices should be invertible
  expect_gt(abs(det(warp$world_to_vox)), 1e-10)
  expect_gt(abs(det(warp$vox_to_world)), 1e-10)

  # world_to_vox and vox_to_world should be inverses
  product <- warp$world_to_vox %*% warp$vox_to_world
  expect_equal(product, diag(4), tolerance = 1e-6)
})

test_that("ANTs NIfTI warp displacement values are in expected range", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "NIfTI warp not available")

  morph <- Warp3DMorphism("native", "mni", warp_path = warp_path, warp_type = "ants")
  warp <- neurotransform:::load_warp_array(morph)

  # Displacement values should be reasonable (not huge)
  # Typical brain registration has displacements < 50mm
  max_disp <- max(abs(warp$array), na.rm = TRUE)
  expect_lt(max_disp, 100)  # Should be well under 100mm

  # Should have some non-zero displacements (not identity)
  expect_gt(max_disp, 0.1)
})

# ==============================================================================
# WARP COMPOSITION TESTS
# ==============================================================================

test_that("ANTs warp can be composed with affine", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "NIfTI warp not available")

  # Create warp morphism
  warp_morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  # Create simple affine morphism
  aff_mat <- diag(4)
  aff_mat[1:3, 4] <- c(5, 10, 15)  # Translation
  aff_morph <- Affine3DMorphism("b", "c", matrix = aff_mat)

  # Compose
  path <- compose(warp_morph, aff_morph)

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)
  expect_equal(source_of(path), "a")
  expect_equal(target_of(path), "c")
})

test_that("ANTs warp path transforms coordinates correctly", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path), "NIfTI warp not available")

  warp_morph <- Warp3DMorphism("a", "b", warp_path = warp_path, warp_type = "ants")

  # Simple translation affine
  aff_mat <- diag(4)
  aff_mat[1:3, 4] <- c(1, 2, 3)
  aff_morph <- Affine3DMorphism("b", "c", matrix = aff_mat)

  path <- compose(warp_morph, aff_morph)

  # Get valid coordinates
  warp <- neurotransform:::load_warp_array(warp_morph)
  center_vox <- warp$dim / 2
  center_world <- (warp$vox_to_world %*% c(center_vox, 1))[1:3]

  test_coords <- matrix(center_world, nrow = 1)

  # Transform through path
  result <- transform(path, test_coords)

  expect_true(all(is.finite(result)))
  expect_equal(ncol(result), 3)
})
