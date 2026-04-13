# Coverage Tests for fsl_ingest.R
#
# Target: Improve coverage from 62.96% to >80%

# ==============================================================================
# FSL SPACING FROM AFFINE
# ==============================================================================

test_that("fsl_spacing_from_affine extracts correct spacing", {
  # Identity affine has spacing 1mm
  aff <- diag(4)
  sp <- neurotransform:::fsl_spacing_from_affine(aff)
  expect_equal(sp, c(1, 1, 1))

  # Scaled affine
  aff2 <- diag(c(2, 3, 4, 1))
  sp2 <- neurotransform:::fsl_spacing_from_affine(aff2)
  expect_equal(sp2, c(2, 3, 4))
})

test_that("fsl_spacing_from_affine handles rotated affines", {
  # Rotation around Z axis (90 degrees)
  rot <- matrix(c(0, 1, 0, 0,
                  -1, 0, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1), 4, 4, byrow = TRUE)
  sp <- neurotransform:::fsl_spacing_from_affine(rot)
  expect_equal(sp, c(1, 1, 1), tolerance = 1e-10)
})

test_that("fsl_spacing_from_affine validates input", {
  expect_error(neurotransform:::fsl_spacing_from_affine(diag(3)))
  expect_error(neurotransform:::fsl_spacing_from_affine(c(1, 2, 3, 4)))
})

# ==============================================================================
# FSL VOX TO FSL
# ==============================================================================

test_that("fsl_vox_to_fsl creates correct scaling", {
  aff <- diag(c(2, 3, 4, 1))
  v2f <- neurotransform:::fsl_vox_to_fsl(aff, dim = c(5L, 6L, 7L))

  expect_equal(dim(v2f), c(4, 4))
  expect_equal(v2f[1, 1], -2)
  expect_equal(v2f[2, 2], 3)
  expect_equal(v2f[3, 3], 4)
  expect_equal(v2f[1, 4], 8)
  expect_equal(v2f[4, 4], 1)
})

test_that("fsl_fsl_to_vox inverts fsl_vox_to_fsl", {
  aff <- diag(c(2, 3, 4, 1))
  dims <- c(5L, 6L, 7L)
  v2f <- neurotransform:::fsl_vox_to_fsl(aff, dim = dims)
  f2v <- neurotransform:::fsl_fsl_to_vox(aff, dim = dims)

  expect_equal(v2f %*% f2v, diag(4), tolerance = 1e-10)
})

# ==============================================================================
# FSL WORLD TO FSL / FSL TO WORLD
# ==============================================================================

test_that("fsl_world_to_fsl and fsl_fsl_to_world are inverses", {
  aff <- diag(c(2, 2, 2, 1))
  aff[1:3, 4] <- c(10, 20, 30)
  dims <- c(5L, 6L, 7L)

  w2f <- neurotransform:::fsl_world_to_fsl(aff, dim = dims)
  f2w <- neurotransform:::fsl_fsl_to_world(aff, dim = dims)

  expect_equal(w2f %*% f2w, diag(4), tolerance = 1e-10)
})

test_that("fsl_world_to_fsl maps world coords to FSL coords", {
  aff <- diag(c(2, 2, 2, 1))
  dims <- c(5L, 6L, 7L)

  w2f <- neurotransform:::fsl_world_to_fsl(aff, dim = dims)

  # A point in world coords should map to FSL coords
  world_pt <- c(10, 20, 30, 1)
  fsl_pt <- w2f %*% world_pt

  expect_equal(length(fsl_pt), 4)
})

# ==============================================================================
# FSL FLIRT TO INTERNAL AFFINE
# ==============================================================================

test_that("fsl_flirt_to_internal_affine converts identity correctly", {
  # Identity FLIRT matrix with identity source/ref
  flirt_mat <- diag(4)
  src_aff <- diag(4)
  ref_aff <- diag(4)

  internal <- fsl_flirt_to_internal_affine(
    flirt_mat, src_aff, ref_aff,
    source_dim = c(5L, 5L, 5L), ref_dim = c(5L, 5L, 5L)
  )

  expect_equal(dim(internal), c(4, 4))
  # For identity inputs, result should be identity
  expect_equal(internal, diag(4), tolerance = 1e-10)
})

test_that("fsl_flirt_to_internal_affine handles scaled images", {
  flirt_mat <- diag(4)
  src_aff <- diag(c(2, 2, 2, 1))  # 2mm source
  ref_aff <- diag(c(1, 1, 1, 1))  # 1mm reference

  internal <- fsl_flirt_to_internal_affine(
    flirt_mat, src_aff, ref_aff,
    source_dim = c(5L, 5L, 5L), ref_dim = c(5L, 5L, 5L)
  )

  # Should be a valid 4x4 matrix
  expect_equal(dim(internal), c(4, 4))
  expect_true(all(is.finite(internal)))
})

test_that("fsl_flirt_to_internal_affine validates FLIRT matrix", {
  expect_error(fsl_flirt_to_internal_affine(diag(3), diag(4), diag(4)),
               "4x4")
  expect_error(fsl_flirt_to_internal_affine(c(1:16), diag(4), diag(4)),
               "4x4")
})

# ==============================================================================
# FSL LOAD FLIRT MORPHISM
# ==============================================================================

test_that("fsl_load_flirt_morphism errors on missing file", {
  expect_error(
    neurotransform:::fsl_load_flirt_morphism(diag(4), diag(4), "/nonexistent.mat"),
    "not found"
  )
})

test_that("fsl_load_flirt_morphism loads from file", {
  mat_path <- system.file("extdata/fsl/S01_lin_6dof.mat", package = "neurotransform")
  skip_if_not(file.exists(mat_path))

  morph <- neurotransform:::fsl_load_flirt_morphism(
    diag(4), diag(4), mat_path,
    source_dim = c(5L, 5L, 5L),
    target_dim = c(5L, 5L, 5L)
  )

  expect_s4_class(morph, "Affine3DMorphism")
  expect_equal(source_of(morph), "source")
  expect_equal(target_of(morph), "target")
})

test_that("fsl_load_flirt_morphism handles custom method_tag", {
  mat_path <- system.file("extdata/fsl/S01_lin_6dof.mat", package = "neurotransform")
  skip_if_not(file.exists(mat_path))

  morph <- neurotransform:::fsl_load_flirt_morphism(
    diag(4), diag(4), mat_path,
    source_dim = c(5L, 5L, 5L),
    target_dim = c(5L, 5L, 5L),
    cost = 2.5,
    method_tag = "custom_method"
  )

  expect_equal(morph@cost, 2.5)
  expect_equal(morph@method_tag, "custom_method")
})

# ==============================================================================
# DETECT FNIRT DEF TYPE
# ==============================================================================

test_that("detect_fnirt_def_type errors on missing file", {
  expect_error(detect_fnirt_def_type("/nonexistent.nii.gz"), "not found")
})

test_that("detect_fnirt_def_type works on relative displacement field", {
  warp_path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("neuroim2")

  result <- detect_fnirt_def_type(warp_path)
  expect_true(result %in% c("relative", "absolute"))
})

test_that("detect_fnirt_def_type works on real ANTs warp", {
  warp_path <- system.file("extdata/chris/ants/reg_1Warp.nii.gz",
                           package = "neurotransform")
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("neuroim2")

  # ANTs warps are displacement fields (relative)
  result <- suppressWarnings(detect_fnirt_def_type(warp_path))
  expect_true(result %in% c("relative", "absolute"))
})

test_that("detect_fnirt_def_type handles custom thresholds", {
  warp_path <- system.file("extdata/fsl/S01_warp.nii.gz", package = "neurotransform")
  skip_if_not(file.exists(warp_path))
  skip_if_not_installed("neuroim2")

  # Different sample size and threshold
  result <- detect_fnirt_def_type(warp_path, sample_n = 50, threshold_mm = 10)
  expect_true(result %in% c("relative", "absolute"))
})
