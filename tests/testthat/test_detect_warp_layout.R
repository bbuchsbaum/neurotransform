write_vector_header <- function(name, shape = c(4L, 5L, 6L, 3L), intent = 0L) {
  dir <- tempfile("layout_")
  dir.create(dir)
  path <- file.path(dir, name)
  image <- RNifti::asNifti(array(0, dim = shape), reference = list(intent_code = intent))
  RNifti::writeNifti(image, path, datatype = "float")
  path
}

test_that("4D vector fields are FSL layout, not ANTs", {
  skip_if_not_installed("RNifti")
  expect_equal(detect_transform_type(write_vector_header("warp.nii.gz")), "fsl")
  expect_equal(detect_transform_type(write_vector_header("warpcoef.nii.gz")), "fsl_coef")
})

test_that("definitive intent codes outrank filename hints", {
  skip_if_not_installed("RNifti")
  expect_equal(detect_transform_type(write_vector_header("field.nii.gz", intent = 2006L)), "fsl")
  expect_equal(detect_transform_type(write_vector_header("field.nii.gz", intent = 2007L)), "fsl_coef")
  itk_named_fsl <- write_vector_header("fnirt_field.nii.gz", c(4L, 5L, 6L, 1L, 3L), 1007L)
  expect_equal(detect_transform_type(itk_named_fsl), "ants")
})

test_that("5D fields keep the existing ANTs default and AFNI name hints", {
  skip_if_not_installed("RNifti")
  expect_equal(detect_transform_type(write_vector_header("warp.nii.gz", c(4L, 5L, 6L, 1L, 3L))), "ants")
  expect_equal(detect_transform_type(write_vector_header("anat_qwarp.nii.gz", c(4L, 5L, 6L, 1L, 3L))), "afni")
  for (name in c("itk_oracle/warp.nii.gz", "ants/sample_ANTs_1Warp.nii.gz",
                 "ants/sample_ANTs_1InverseWarp.nii.gz")) {
    path <- system.file("extdata", name, package = "neurotransform")
    skip_if_not(file.exists(path))
    expect_equal(detect_transform_type(path), "ants")
  }
})

test_that("an inferred FSL field without geometry fails when read", {
  skip_if_not_installed("RNifti")
  path <- write_vector_header("warp.nii.gz")
  expect_error(read_transform(path), "Detected a dense FSL warp field")
  # Same source and reference geometry with a zero field is the identity.
  lattice <- neuroim2::trans(neuroim2::read_header(path))
  m <- read_transform(path, source_affine = lattice, source_dim = c(4L, 5L, 6L),
                      def_type = "relative")
  expect_equal(m@warp_type, "fsl")
  expect_equal(unname(transform(m, matrix(c(1, 2, 3), 1))), matrix(c(1, 2, 3), 1))
})
