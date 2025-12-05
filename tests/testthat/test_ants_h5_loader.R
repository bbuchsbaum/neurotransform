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
