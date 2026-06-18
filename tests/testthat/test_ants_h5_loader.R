write_ants_h5_fixture <- function(path, include_warp = TRUE, include_affine = TRUE,
                                  misspelled = FALSE) {
  dataset_name <- function(name) {
    if (isTRUE(misspelled)) sub("^Transform", "Tranform", name) else name
  }

  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all())

  tg <- h5$create_group("TransformGroup")
  g0 <- tg$create_group("0")
  g0$create_dataset("TransformType", robj = "CompositeTransform_double_3_3")

  key <- 1L
  if (isTRUE(include_affine)) {
    g <- tg$create_group(as.character(key))
    key <- key + 1L
    g$create_dataset("TransformType", robj = "AffineTransform_double_3_3")
    g$create_dataset(
      dataset_name("TransformParameters"),
      robj = c(as.numeric(t(diag(3))), 0, 0, 0)
    )
    g$create_dataset(dataset_name("TransformFixedParameters"), robj = c(0, 0, 0))
  }

  if (isTRUE(include_warp)) {
    g <- tg$create_group(as.character(key))
    size <- c(2L, 2L, 2L)
    nvox <- prod(size)
    fixed <- c(
      as.numeric(size),
      0, 0, 0,
      1, 1, 1,
      as.numeric(t(diag(3)))
    )
    params <- c(seq_len(nvox), 10 + seq_len(nvox), 20 + seq_len(nvox))

    g$create_dataset("TransformType", robj = "DisplacementFieldTransform_double_3_3")
    g$create_dataset(dataset_name("TransformFixedParameters"), robj = fixed)
    g$create_dataset(dataset_name("TransformParameters"), robj = as.numeric(params))
  }

  invisible(path)
}

test_that("ANTS composite H5 warp loads via hdf5r loader", {
  path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5", package = "neurotransform")
  skip_if_not(file.exists(path))
  skip_if_not_installed("hdf5r")
  loader <- get_loader("ants_h5")
  w <- loader(path)
  expect_equal(w$dim, c(97L, 115L, 97L))
  expect_true(length(w$array) == prod(w$dim) * 3)
})

test_that("Warp3DMorphism can load ants_h5 via default loader", {
  path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5", package = "neurotransform")
  skip_if_not(file.exists(path))
  skip_if_not_installed("hdf5r")
  m <- Warp3DMorphism("src", "tgt", path, warp_type = "ants_h5")
  w <- load_warp_array(m)
  expect_equal(w$dim, c(97L, 115L, 97L))
})

test_that("ants_h5_morphism returns path with embedded affine when requested", {
  path <- system.file("extdata/chris/ants/chris_to_mni_Composite.h5", package = "neurotransform")
  skip_if_not(file.exists(path))
  skip_if_not_installed("hdf5r")
  m <- ants_h5_morphism(path, source = "s", target = "t", apply_affine = TRUE)
  expect_true(is(m, "Morphism") || is(m, "MorphismPath"))
  if (is(m, "MorphismPath")) {
    expect_equal(length(m@morphisms), 2L)
    # Path is [warp, affine] for pullback: affine_pullback(warp_pullback(coords))
    expect_s4_class(m@morphisms[[1]], "Warp3DMorphism")
    expect_s4_class(m@morphisms[[2]], "Affine3DMorphism")
  }
})

test_that("ANTS H5 warp loader accepts TemplateFlow Tranform dataset names", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  write_ants_h5_fixture(path, misspelled = TRUE)

  w <- load_warp_ants_h5(path)

  expect_equal(w$dim, c(2L, 2L, 2L))
  expect_equal(length(w$array), prod(w$dim) * 3L)
  expect_equal(w$array[1:6], c(1, 11, 21, 2, 12, 22))
  expect_equal(w$affine, diag(4))
})

test_that("ants_h5_morphism and read_transform accept TemplateFlow Tranform dataset names", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  write_ants_h5_fixture(path, misspelled = TRUE)

  direct <- ants_h5_morphism(path, source = "native", target = "mni", apply_affine = TRUE)
  through_io <- read_transform(path, type = "ants_h5", source = "native", target = "mni")

  expect_s4_class(direct, "MorphismPath")
  expect_s4_class(through_io, "MorphismPath")
  expect_s4_class(direct@morphisms[[1]], "Warp3DMorphism")
  expect_s4_class(direct@morphisms[[2]], "Affine3DMorphism")
})

test_that("ITK H5 affine reader accepts TemplateFlow Tranform dataset names", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  write_ants_h5_fixture(path, include_warp = FALSE, misspelled = TRUE)

  morph <- read_linear_transform(path, format = "itk", source = "native", target = "mni")

  expect_s4_class(morph, "Affine3DMorphism")
  expect_equal(morph@matrix, diag(4))
})
