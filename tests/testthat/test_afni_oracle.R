# AFNI conversion checked against AFNI itself.
#
# Every expected value below comes out of AFNI_26.1.04, not out of this package.
# See inst/extdata/afni_oracle/README.md for the exact commands. A round trip
# through our own reader and writer could hide a shared sign or direction error;
# these fixtures cannot.

oracle_path <- function(...) {
  system.file("extdata", "afni_oracle", ..., package = "neurotransform")
}

skip_without_oracle <- function() {
  skip_if(!nzchar(oracle_path("oracle.aff12.1D")), "AFNI oracle fixtures not installed")
}

# AFNI RAI/DICOM -> RAS: negate X and Y.
dicom_to_ras <- function(xyz) cbind(-xyz[, 1], -xyz[, 2], xyz[, 3])

unclass_attrs <- function(x) {
  attr(x, "comments") <- NULL
  x
}

read_landmarks <- function(file, cols) {
  m <- as.matrix(utils::read.table(oracle_path(file)))
  dimnames(m) <- NULL
  colnames(m) <- cols
  m
}

test_that("AFNI's own matrix file reads back with the layout AFNI wrote", {
  skip_without_oracle()
  mat <- afni_read_aff12(oracle_path("oracle.aff12.1D"))
  expect_equal(dim(mat), c(4L, 4L))
  # cat_matvec 'MATRIX(1,0,0,7,0,0,-1,-3,0,1,0,5)' is row-major 3x4.
  expect_equal(
    mat,
    rbind(c(1, 0, 0, 7), c(0, 0, -1, -3), c(0, 1, 0, 5), c(0, 0, 0, 1))
  )
})

test_that("our inversion agrees with AFNI's cat_matvec -I", {
  skip_without_oracle()
  mat <- afni_read_aff12(oracle_path("oracle.aff12.1D"))
  afni_inverse <- afni_read_aff12(oracle_path("oracle_inverse.aff12.1D"))
  expect_equal(invert_affine(mat), afni_inverse, tolerance = 1e-9)
})

test_that("the RAS pullback reproduces AFNI's resampled landmark positions", {
  skip_without_oracle()
  # AFNI placed these landmarks in the source volume ...
  src <- read_landmarks("landmarks_source_dicom.txt", c("x", "y", "z", "value"))
  # ... and 3dAllineate -1Dmatrix_apply put them here on the base grid.
  base <- read_landmarks("landmarks_base_dicom.txt",
                         c("i", "j", "k", "x", "y", "z", "value"))
  src <- src[order(src[, "value"]), , drop = FALSE]
  base <- base[order(base[, "value"]), , drop = FALSE]

  src_ras <- dicom_to_ras(src[, c("x", "y", "z"), drop = FALSE])
  base_ras <- dicom_to_ras(base[, c("x", "y", "z"), drop = FALSE])

  m <- read_linear_transform(oracle_path("oracle.aff12.1D"), format = "afni",
                             source = "src", target = "base")
  expect_s4_class(m, "Affine3DMorphism")

  # A morphism's matrix is the target-to-source pullback, which is exactly what
  # AFNI stores. Applying it to base coordinates must land on the source points.
  expect_lt(max(abs(transform(m, base_ras) - src_ras)), 1e-9)

  # Inverting it gives the forward source-to-base mapping.
  expect_lt(max(abs(transform(invert(m), src_ras) - base_ras)), 1e-9)
})

test_that("the pullback agrees with AFNI's own Vecwarp point mapping", {
  skip_without_oracle()
  base <- read_landmarks("landmarks_base_dicom.txt",
                         c("i", "j", "k", "x", "y", "z", "value"))
  vec <- read_landmarks("vecwarp_source_dicom.1D", c("x", "y", "z"))

  m <- read_linear_transform(oracle_path("oracle.aff12.1D"), format = "afni",
                             source = "src", target = "base")
  ours <- transform(m, dicom_to_ras(base[, c("x", "y", "z"), drop = FALSE]))
  expect_lt(max(abs(ours - dicom_to_ras(vec))), 1e-9)
})

test_that("a Z flip would have been caught by the oracle", {
  skip_without_oracle()
  # Guards the defect this fixture set was built to find: conjugating with
  # diag(1, 1, -1, 1) instead of diag(-1, -1, 1, 1) leaves the linear part of
  # this rotation unchanged but negates the translation.
  mat <- afni_read_aff12(oracle_path("oracle.aff12.1D"))
  correct <- afni_aff12_to_ras(mat, oblique_correction = FALSE)
  z_flip <- diag(c(1, 1, -1, 1))
  wrong <- z_flip %*% mat %*% z_flip
  expect_equal(correct[1:3, 4], c(-7, 3, 5))
  expect_equal(wrong[1:3, 4], c(7, -3, -5))
  expect_gt(max(abs(correct - wrong)), 1)
})

test_that("a real 3dvolreg series reads as one matrix per sub-brick", {
  skip_without_oracle()
  path <- oracle_path("volreg_series.aff12.1D")
  mats <- afni_read_aff12_array(path)
  # The synthetic series had four volumes.
  expect_length(mats, 4L)
  expect_true(all(vapply(mats, function(m) identical(dim(m), c(4L, 4L)), logical(1))))
  # AFNI's header comment is preserved rather than silently dropped.
  expect_match(attr(mats, "comments"), "DICOM-to-DICOM", all = FALSE)
  # Sub-brick 0 is the base, so its matrix is the identity (AFNI writes -0).
  expect_equal(mats[[1]], diag(4), tolerance = 1e-12)

  # The volumes were displaced by known RAS amounts; AFNI stores the
  # base-to-source mapping in DICOM, so X and Y appear negated.
  expect_lt(max(abs(mats[[2]][1:3, 4] - c(-1.5, 0, 0))), 1e-3)
  expect_lt(max(abs(mats[[3]][1:3, 4] - c(0, 2, 0.5))), 1e-3)
  expect_lt(max(abs(mats[[4]][1:3, 4] - c(1, -1, -1.5))), 1e-3)

  arr <- read_linear_transform_array(path, format = "afni",
                                     source = "vol", target = "base")
  expect_s3_class(arr, "LinearTransformArray")
  expect_length(arr$transforms, 4L)
})

test_that("a single-matrix read refuses an ambiguous multi-row file", {
  skip_without_oracle()
  expect_error(
    afni_read_aff12(oracle_path("volreg_series.aff12.1D")),
    "holds 4 matrices"
  )
  expect_equal(
    afni_read_aff12(oracle_path("volreg_series.aff12.1D"), row = 1L),
    diag(4),
    tolerance = 1e-12
  )
})

test_that("what we write is what AFNI writes", {
  skip_without_oracle()
  mats <- afni_read_aff12_array(oracle_path("volreg_series.aff12.1D"))
  path <- tempfile(fileext = ".aff12.1D")
  afni_write_aff12(mats, path)

  written <- readLines(path)
  # One comment banner, then exactly one row of twelve numbers per matrix --
  # the layout 3dAllineate -1Dmatrix_apply expects.
  body <- written[!grepl("^\\s*#", written)]
  expect_length(body, 4L)
  expect_true(all(vapply(body, function(l) length(scan(text = l, quiet = TRUE)), integer(1)) == 12L))
  expect_equal(unclass_attrs(afni_read_aff12_array(path)), unclass_attrs(mats),
               tolerance = 1e-9)
})

test_that("malformed and degenerate matrix files fail with a clear message", {
  bad <- tempfile(fileext = ".1D")

  writeLines(c("# header only"), bad)
  expect_error(afni_read_aff12(bad), "no matrix rows")

  writeLines("1 0 0 0 0 1 0 0", bad)
  expect_error(afni_read_aff12(bad), "unsupported layout")

  writeLines("1 0 0 7 0 0 NaN -3 0 1 0 5", bad)
  expect_error(afni_read_aff12(bad), "non-finite")

  # A collapsed registration: the linear part has no inverse.
  writeLines("0 0 0 7 0 0 0 -3 0 0 0 5", bad)
  expect_error(
    read_linear_transform(bad, format = "afni", source = "s", target = "t"),
    "singular"
  )
})

test_that("an AFNI matrix series round-trips through the array writer", {
  skip_without_oracle()
  arr <- read_linear_transform_array(oracle_path("volreg_series.aff12.1D"),
                                     format = "afni", source = "vol", target = "base")
  path <- tempfile(fileext = ".aff12.1D")
  write_linear_transform_array(arr, path, format = "afni")
  back <- read_linear_transform_array(path, format = "afni",
                                      source = "vol", target = "base")
  expect_length(back$transforms, 4L)
  for (i in seq_along(arr$transforms)) {
    expect_equal(back$transforms[[i]]@matrix, arr$transforms[[i]]@matrix,
                 tolerance = 1e-9)
  }
})
