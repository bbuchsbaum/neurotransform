# Test validate_path() and is_valid_path() functions

test_that("validate_path accepts empty path", {
  expect_silent(validate_path(list()))
})

test_that("validate_path accepts single morphism", {
  aff <- Affine3DMorphism("a", "b", diag(4))
  expect_silent(validate_path(list(aff)))
})

test_that("validate_path accepts composable two-element path", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))

  expect_silent(validate_path(list(f, g)))
})

test_that("validate_path accepts composable three-element path", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))
  h <- Affine3DMorphism("c", "d", diag(4))

  expect_silent(validate_path(list(f, g, h)))
})

test_that("validate_path accepts mixed morphism types", {
  aff <- Affine3DMorphism("a", "b", diag(4))
  warp <- Warp3DMorphism("b", "c", "dummy.nii", warp_type = "ants")
  v2s <- VolToSurfMorphism("c", "surf", method = "trilinear")

  expect_silent(validate_path(list(aff, warp, v2s)))
})

test_that("validate_path errors on incompatible domains", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("c", "d", diag(4))  # b != c

  expect_error(validate_path(list(f, g)), "not composable")
})

test_that("validate_path error message includes position", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))
  h <- Affine3DMorphism("x", "y", diag(4))  # c != x at position 2

  expect_error(validate_path(list(f, g, h)), "position 2")
})

test_that("validate_path error message includes domain names", {
  f <- Affine3DMorphism("domain_a", "domain_b", diag(4))
  g <- Affine3DMorphism("domain_x", "domain_y", diag(4))

  expect_error(validate_path(list(f, g)), "domain_b")
  expect_error(validate_path(list(f, g)), "domain_x")
})

test_that("validate_path errors on non-Morphism elements", {
  f <- Affine3DMorphism("a", "b", diag(4))

  expect_error(validate_path(list(f, "not a morphism")), "must be Morphism")
  expect_error(validate_path(list(f, 42)), "must be Morphism")
  expect_error(validate_path(list(f, list())), "must be Morphism")
  expect_error(validate_path(list(f, NULL)), "must be Morphism")
})

test_that("validate_path errors on non-Morphism first element", {
  expect_error(validate_path(list("not a morphism")), "must be Morphism")
})

test_that("validate_path returns path invisibly when valid", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))

  result <- validate_path(list(f, g))

  expect_identical(result, list(f, g))
})

test_that("is_valid_path returns TRUE for valid paths", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))
  h <- Affine3DMorphism("c", "d", diag(4))

  expect_true(is_valid_path(list()))
  expect_true(is_valid_path(list(f)))
  expect_true(is_valid_path(list(f, g)))
  expect_true(is_valid_path(list(f, g, h)))
})

test_that("is_valid_path returns FALSE for invalid paths", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("c", "d", diag(4))  # b != c

  expect_false(is_valid_path(list(f, g)))
})

test_that("is_valid_path returns FALSE for non-Morphism elements", {
  f <- Affine3DMorphism("a", "b", diag(4))

  expect_false(is_valid_path(list(f, "not a morphism")))
  expect_false(is_valid_path(list(f, 42)))
  expect_false(is_valid_path(list("string")))
})

test_that("is_valid_path does not throw errors (wraps validate_path)", {
  # Even for invalid inputs, is_valid_path should return FALSE, not error
  expect_false(is_valid_path(list("a", "b")))
  expect_false(is_valid_path(list(NULL)))
})

test_that("validate_path handles identity morphisms", {
  id1 <- IdentityMorphism("a")
  id2 <- IdentityMorphism("a")  # Same domain
  aff <- Affine3DMorphism("a", "b", diag(4))

  expect_silent(validate_path(list(id1, aff)))
  expect_silent(validate_path(list(id1, id2)))  # Both have domain "a"
})

test_that("validate_path handles path starting with warp", {
  warp <- Warp3DMorphism("a", "b", "dummy.nii", warp_type = "ants")
  aff <- Affine3DMorphism("b", "c", diag(4))

  expect_silent(validate_path(list(warp, aff)))
})

test_that("validate_path handles surface morphisms", {
  v2s <- VolToSurfMorphism("vol", "surf1", method = "trilinear")
  s2s <- SurfToSurfMorphism("surf1", "surf2", method = "sphere")

  expect_silent(validate_path(list(v2s, s2s)))
})

test_that("is_valid_path handles long paths", {
  # Create a chain of 10 affines
  morphisms <- list()
  domains <- paste0("d", 0:10)

  for (i in 1:10) {
    morphisms[[i]] <- Affine3DMorphism(domains[i], domains[i + 1], diag(4))
  }

  expect_true(is_valid_path(morphisms))
})

test_that("is_valid_path detects break in middle of long path", {
  morphisms <- list(
    Affine3DMorphism("a", "b", diag(4)),
    Affine3DMorphism("b", "c", diag(4)),
    Affine3DMorphism("c", "d", diag(4)),
    Affine3DMorphism("x", "y", diag(4)),  # Break: d != x
    Affine3DMorphism("y", "z", diag(4))
  )

  expect_false(is_valid_path(morphisms))
})
