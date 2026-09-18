write_ants_h5_fixture <- function(path, include_warp = TRUE, include_affine = TRUE,
                                  misspelled = FALSE,
                                  component_order = c("affine", "warp"),
                                  affine_matrix = diag(3),
                                  affine_translation = c(0, 0, 0),
                                  affine_center = c(0, 0, 0),
                                  displacement = NULL,
                                  size = c(2L, 2L, 2L),
                                  origin = c(0, 0, 0),
                                  spacing = c(1, 1, 1),
                                  direction = diag(3)) {
  dataset_name <- function(name) {
    if (isTRUE(misspelled)) sub("^Transform", "Tranform", name) else name
  }

  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all())

  tg <- h5$create_group("TransformGroup")
  g0 <- tg$create_group("0")
  g0$create_dataset("TransformType", robj = "CompositeTransform_double_3_3")

  enabled <- component_order[
    (component_order != "affine" | isTRUE(include_affine)) &
      (component_order != "warp" | isTRUE(include_warp))
  ]
  for (i in seq_along(enabled)) {
    g <- tg$create_group(as.character(i))
    if (identical(enabled[[i]], "affine")) {
      g$create_dataset("TransformType", robj = "AffineTransform_double_3_3")
      g$create_dataset(
        dataset_name("TransformParameters"),
        robj = c(as.numeric(t(affine_matrix)), affine_translation)
      )
      g$create_dataset(
        dataset_name("TransformFixedParameters"), robj = affine_center
      )
    } else {
      nvox <- prod(size)
      fixed <- c(
        as.numeric(size), origin, spacing, as.numeric(t(direction))
      )
      params <- if (is.null(displacement)) {
        as.numeric(t(cbind(seq_len(nvox), 10 + seq_len(nvox), 20 + seq_len(nvox))))
      } else {
        rep(displacement, times = nvox)
      }
      g$create_dataset("TransformType", robj = "DisplacementFieldTransform_double_3_3")
      g$create_dataset(dataset_name("TransformFixedParameters"), robj = fixed)
      g$create_dataset(dataset_name("TransformParameters"), robj = as.numeric(params))
    }
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
    expect_s4_class(m@morphisms[[1]], "Affine3DMorphism")
    expect_s4_class(m@morphisms[[2]], "Warp3DMorphism")
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
  expect_equal(w$array[1:6], c(-1, -11, 21, -2, -12, 22))
  expect_equal(w$affine, diag(4))
  expect_equal(w$transform_order, c("affine", "warp"))
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
  expect_s4_class(direct@morphisms[[1]], "Affine3DMorphism")
  expect_s4_class(direct@morphisms[[2]], "Warp3DMorphism")
})

test_that("ANTs H5 composite matches analytic LPS point semantics", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)

  A_lps <- matrix(c(
    1.05, 0.10, 0.00,
    0.00, 0.95, 0.15,
    0.00, 0.00, 1.10
  ), 3, byrow = TRUE)
  translation <- c(1.5, -2, 0.75)
  center <- c(0.5, 1, -0.5)
  displacement_lps <- c(0.25, -0.5, 0.75)
  write_ants_h5_fixture(
    path,
    affine_matrix = A_lps,
    affine_translation = translation,
    affine_center = center,
    displacement = displacement_lps,
    size = c(4L, 4L, 4L)
  )

  morph <- ants_h5_morphism(path, source = "moving", target = "fixed")
  points_lps <- rbind(c(1, 1, 1), c(2, 1, 2))
  flip3 <- diag(c(-1, -1, 1))
  points_ras <- points_lps %*% flip3
  expected_lps <- t(vapply(seq_len(nrow(points_lps)), function(i) {
    displaced <- points_lps[i, ] + displacement_lps
    as.numeric(A_lps %*% (displaced - center) + center + translation)
  }, numeric(3)))

  expect_s4_class(morph, "MorphismPath")
  expect_equal(transform(morph, points_ras), expected_lps %*% flip3,
               tolerance = 1e-7)
})

test_that("ANTs H5 loader fails closed on ambiguous composites", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  write_ants_h5_fixture(path, component_order = c("affine", "affine", "warp"))

  expect_error(
    load_warp_ants_h5(path),
    "Multiple embedded affine transforms"
  )
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

test_that("H5 loading refuses unsupported components rather than dropping them", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path))
  write_ants_h5_fixture(path)
  h5 <- hdf5r::H5File$new(path, "r+")
  g <- h5[["TransformGroup"]]$create_group("3")
  g[["TransformType"]] <- "TranslationTransform_double_3_3"
  g[["TransformParameters"]] <- c(1, 2, 3)
  g[["TransformFixedParameters"]] <- numeric(0)
  h5$close_all()
  expect_error(read_transform(path), "Unsupported.*TranslationTransform")
  expect_error(read_linear_transform_array(path, format = "itk"), "Unsupported")
})

test_that("H5 loading refuses malformed affine and displacement parameters", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path))
  cases <- list(
    list("1", "TransformParameters", rep(0, 11), "affine parameters"),
    list("1", "TransformParameters", c(rep(0, 11), NaN), "affine parameters"),
    list("1", "TransformFixedParameters", c(0, 0), "affine fixed parameters"),
    list("2", "TransformParameters", rep(0, 25), "Displacement parameters"),
    list("2", "TransformParameters", rep(0, 23), "Displacement parameters"),
    list("2", "TransformFixedParameters", c(2.5, 2, 2, rep(0, 15)), "grid size"),
    list("2", "TransformFixedParameters", rep(0, 17), "fixed parameters")
  )
  for (case in cases) {
    write_ants_h5_fixture(path)
    h5 <- hdf5r::H5File$new(path, "r+")
    g <- h5[["TransformGroup"]][[case[[1]]]]
    g$link_delete(case[[2]])
    g[[case[[2]]]] <- case[[3]]
    h5$close_all()
    expect_error(read_transform(path), case[[4]])
  }
})
