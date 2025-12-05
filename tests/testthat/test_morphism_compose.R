# Test morphism compose() methods comprehensively

test_that("compose with identity returns other morphism (identity first)", {
  id <- IdentityMorphism("a")
  aff <- Affine3DMorphism("a", "b", diag(4))

  result <- compose(id, aff)
  expect_identical(result, aff)
})

test_that("compose with identity returns other morphism (identity second)", {
  aff <- Affine3DMorphism("a", "b", diag(4))
  id <- IdentityMorphism("b")

  result <- compose(aff, id)
  expect_identical(result, aff)
})

test_that("compose identity with identity returns second identity", {
  id1 <- IdentityMorphism("a")
  id2 <- IdentityMorphism("a")

  result <- compose(id1, id2)
  expect_identical(result, id2)
})

test_that("compose checks domain compatibility", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("c", "d", diag(4))  # b != c

  expect_error(compose(f, g), "target of first must equal source of second")
})

test_that("compose fuses consecutive affines into single affine", {
  A <- diag(4)
  A[1, 4] <- 10
  B <- diag(4)
  B[2, 4] <- 20

  f <- Affine3DMorphism("a", "b", A)
  g <- Affine3DMorphism("b", "c", B)

  fg <- compose(f, g)

  # Should be single affine, not path

  expect_s4_class(fg, "Affine3DMorphism")
  expect_equal(fg@source, "a")
  expect_equal(fg@target, "c")

  # Matrix should be A %*% B (pullback semantics)
  expect_equal(fg@matrix, A %*% B)
})

test_that("compose accumulates costs for fused affines", {
  f <- Affine3DMorphism("a", "b", diag(4), cost = 1.5)
  g <- Affine3DMorphism("b", "c", diag(4), cost = 2.5)

  fg <- compose(f, g)

  expect_equal(fg@cost, 4.0)  # 1.5 + 2.5
})

test_that("compose mixed types (affine + warp) creates MorphismPath", {
  aff <- Affine3DMorphism("a", "b", diag(4))
  warp <- Warp3DMorphism("b", "c", "dummy.nii", warp_type = "ants")

  path <- compose(aff, warp)

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)
  expect_equal(path@source, "a")
  expect_equal(path@target, "c")
})

test_that("compose warp + affine creates MorphismPath", {
  warp <- Warp3DMorphism("a", "b", "dummy.nii", warp_type = "ants")
  aff <- Affine3DMorphism("b", "c", diag(4))

  path <- compose(warp, aff)

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)
  expect_equal(path@source, "a")
  expect_equal(path@target, "c")
})

test_that("compose two warps creates MorphismPath", {
  warp1 <- Warp3DMorphism("a", "b", "dummy1.nii", warp_type = "ants")
  warp2 <- Warp3DMorphism("b", "c", "dummy2.nii", warp_type = "ants")

  path <- compose(warp1, warp2)

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)
})

test_that("compose MorphismPath with Morphism appends correctly", {
  f <- Affine3DMorphism("a", "b", diag(4))
  warp <- Warp3DMorphism("b", "c", "dummy.nii", warp_type = "ants")
  h <- Affine3DMorphism("c", "d", diag(4))

  # First compose f + warp to get a path
  path_fw <- compose(f, warp)
  expect_s4_class(path_fw, "MorphismPath")

  # Then append h
  path_fwh <- compose(path_fw, h)

  expect_s4_class(path_fwh, "MorphismPath")
  expect_equal(length(path_fwh@morphisms), 3)
  expect_equal(path_fwh@source, "a")
  expect_equal(path_fwh@target, "d")
})

test_that("compose MorphismPath with Morphism checks domain compatibility", {
  f <- Affine3DMorphism("a", "b", diag(4))
  warp <- Warp3DMorphism("b", "c", "dummy.nii", warp_type = "ants")
  h <- Affine3DMorphism("x", "y", diag(4))  # c != x

  path_fw <- compose(f, warp)

  expect_error(compose(path_fw, h), "target of path must equal source of morphism")
})

test_that("compose Morphism with MorphismPath prepends correctly", {
  f <- Affine3DMorphism("a", "b", diag(4))
  warp <- Warp3DMorphism("b", "c", "dummy.nii", warp_type = "ants")
  h <- Affine3DMorphism("c", "d", diag(4))

  # Create path warp + h
  path_wh <- compose(warp, h)
  expect_s4_class(path_wh, "MorphismPath")

  # Prepend f
  path_fwh <- compose(f, path_wh)

  expect_s4_class(path_fwh, "MorphismPath")
  expect_equal(length(path_fwh@morphisms), 3)
  expect_equal(path_fwh@source, "a")
  expect_equal(path_fwh@target, "d")
})

test_that("compose Morphism with MorphismPath checks domain compatibility", {
  warp <- Warp3DMorphism("b", "c", "dummy.nii", warp_type = "ants")
  h <- Affine3DMorphism("c", "d", diag(4))
  f <- Affine3DMorphism("x", "y", diag(4))  # y != b

  path_wh <- compose(warp, h)

  expect_error(compose(f, path_wh), "target of morphism must equal source of path")
})

test_that("compose two MorphismPaths concatenates correctly", {
  f <- Affine3DMorphism("a", "b", diag(4))
  warp1 <- Warp3DMorphism("b", "c", "dummy1.nii", warp_type = "ants")
  warp2 <- Warp3DMorphism("c", "d", "dummy2.nii", warp_type = "ants")
  g <- Affine3DMorphism("d", "e", diag(4))

  path1 <- compose(f, warp1)
  path2 <- compose(warp2, g)

  combined <- compose(path1, path2)

  expect_s4_class(combined, "MorphismPath")
  expect_equal(length(combined@morphisms), 4)
  expect_equal(combined@source, "a")
  expect_equal(combined@target, "e")
})

test_that("compose two MorphismPaths checks domain compatibility", {
  f <- Affine3DMorphism("a", "b", diag(4))
  warp1 <- Warp3DMorphism("b", "c", "dummy1.nii", warp_type = "ants")
  warp2 <- Warp3DMorphism("x", "y", "dummy2.nii", warp_type = "ants")  # x != c
  g <- Affine3DMorphism("y", "z", diag(4))

  path1 <- compose(f, warp1)
  path2 <- compose(warp2, g)

  expect_error(compose(path1, path2), "target of first path must equal source of second path")
})

test_that("compose three or more morphisms works correctly", {
  f <- Affine3DMorphism("a", "b", diag(4))
  g <- Affine3DMorphism("b", "c", diag(4))
  h <- Affine3DMorphism("c", "d", diag(4))

  # Compose all three (affines fuse)
  fg <- compose(f, g)
  fgh <- compose(fg, h)

  # Three fused affines should still be single affine
  expect_s4_class(fgh, "Affine3DMorphism")
  expect_equal(fgh@source, "a")
  expect_equal(fgh@target, "d")
})

test_that("compose preserves transformation semantics", {
  # Create translation affines
  A <- diag(4)
  A[1, 4] <- 5  # +5 in X
  B <- diag(4)
  B[2, 4] <- 10  # +10 in Y

  f <- Affine3DMorphism("a", "b", A)
  g <- Affine3DMorphism("b", "c", B)

  fg <- compose(f, g)

  # Transform test point
  coords <- matrix(c(0, 0, 0), ncol = 3)
  result <- transform(fg, coords)

  # Should get (5, 10, 0)
  expect_equal(result[1, 1], 5, tolerance = 1e-10)
  expect_equal(result[1, 2], 10, tolerance = 1e-10)
  expect_equal(result[1, 3], 0, tolerance = 1e-10)
})

test_that("composed affines agree with fused matrix application", {
  # With pullback semantics: compose(f, g) means "apply f's matrix then g's matrix"
  # The composed matrix is A %*% B (f first, g second in pullback order)
  A <- matrix(c(
    2, 0, 0, 1,
    0, 1, 0, 2,
    0, 0, 1, 3,
    0, 0, 0, 1
  ), 4, 4, byrow = TRUE)

  B <- matrix(c(
    1, 0, 0, -1,
    0, 2, 0, 0,
    0, 0, 3, 0,
    0, 0, 0, 1
  ), 4, 4, byrow = TRUE)

  f <- Affine3DMorphism("a", "b", A)
  g <- Affine3DMorphism("b", "c", B)
  fg <- compose(f, g)

  coords <- matrix(c(1, 2, 3, 10, 20, 30, -5, -5, -5), ncol = 3, byrow = TRUE)

  # Direct matrix application
  direct_result <- apply_affine(coords, A %*% B)

  # Composed
  comp_result <- transform(fg, coords)

  expect_equal(comp_result, direct_result, tolerance = 1e-10)
})

test_that("compose with VolToSurfMorphism creates path", {
  aff <- Affine3DMorphism("a", "b", diag(4))
  v2s <- VolToSurfMorphism("b", "surf", method = "trilinear")

  path <- compose(aff, v2s)

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)
  expect_equal(path@source, "a")
  expect_equal(path@target, "surf")
})

test_that("compose with SurfToSurfMorphism creates path", {
  v2s <- VolToSurfMorphism("vol", "surf1", method = "trilinear")
  s2s <- SurfToSurfMorphism("surf1", "surf2", method = "sphere")

  path <- compose(v2s, s2s)

  expect_s4_class(path, "MorphismPath")
  expect_equal(length(path@morphisms), 2)
})
