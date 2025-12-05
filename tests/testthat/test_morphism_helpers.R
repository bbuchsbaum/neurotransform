# Test morphism helper functions: morphism_hash, morphism_kind, is_linear_morphism, etc.

test_that("morphism_hash is stable and unique", {
  aff1 <- Affine3DMorphism("a", "b", diag(4))
  aff2 <- Affine3DMorphism("a", "b", diag(4))
  aff3 <- Affine3DMorphism("a", "c", diag(4))

  # Stability: same morphism config -> same hash components
  # (Note: hashes include id which uses digest, so we test via morphism_hash)
  hash1 <- morphism_hash(aff1)
  hash2 <- morphism_hash(aff1)  # Same object
  expect_equal(hash1, hash2)

  # Uniqueness: different targets -> different hashes
  hash_a <- morphism_hash(aff1)
  hash_c <- morphism_hash(aff3)
  expect_false(hash_a == hash_c)
})

test_that("morphism_hash errors on non-Morphism", {
  expect_error(morphism_hash("not a morphism"), "must be a Morphism")
  expect_error(morphism_hash(list()), "must be a Morphism")
  expect_error(morphism_hash(NULL), "must be a Morphism")
})

test_that("morphism_kind returns correct values for each type", {
  expect_equal(morphism_kind(IdentityMorphism("a")), "identity")
  expect_equal(morphism_kind(Affine3DMorphism("a", "b", diag(4))), "affine3d")

  # Create a minimal warp morphism (file doesn't need to exist for kind check)
  warp <- Warp3DMorphism("a", "b", "nonexistent.nii", warp_type = "ants")
  expect_equal(morphism_kind(warp), "warp3d")

  v2s <- VolToSurfMorphism("a", "b", method = "trilinear")
  expect_equal(morphism_kind(v2s), "vol2surf")

  s2s <- SurfToSurfMorphism("a", "b", method = "sphere")
  expect_equal(morphism_kind(s2s), "surf2surf")
})

test_that("morphism_kind errors on non-Morphism", {
  expect_error(morphism_kind("not a morphism"), "must be a Morphism")
  expect_error(morphism_kind(42), "must be a Morphism")
})

test_that("is_linear_morphism correctly identifies linear morphisms", {
  expect_true(is_linear_morphism(IdentityMorphism("a")))
  expect_true(is_linear_morphism(Affine3DMorphism("a", "b", diag(4))))
  expect_false(is_linear_morphism(Warp3DMorphism("a", "b", "x.nii", warp_type = "ants")))
  expect_false(is_linear_morphism(VolToSurfMorphism("a", "b", method = "trilinear")))
  expect_false(is_linear_morphism(SurfToSurfMorphism("a", "b", method = "sphere")))
})

test_that("is_warp_morphism correctly identifies warp morphisms", {
  expect_true(is_warp_morphism(Warp3DMorphism("a", "b", "x.nii", warp_type = "ants")))
  expect_false(is_warp_morphism(IdentityMorphism("a")))
  expect_false(is_warp_morphism(Affine3DMorphism("a", "b", diag(4))))
  expect_false(is_warp_morphism(VolToSurfMorphism("a", "b", method = "trilinear")))
  expect_false(is_warp_morphism(SurfToSurfMorphism("a", "b", method = "sphere")))
})

test_that("is_invertible checks inverse_type correctly", {
  # Identity has exact inverse
  expect_true(is_invertible(IdentityMorphism("a")))

  # Affine has exact inverse
  expect_true(is_invertible(Affine3DMorphism("a", "b", diag(4))))

  # Warp without inverse_path has none
  w_no_inv <- Warp3DMorphism("a", "b", "fwd.nii", warp_type = "ants")
  expect_false(is_invertible(w_no_inv))

  # Warp with inverse_path has "provided" not "exact"
  w_inv <- Warp3DMorphism("a", "b", "fwd.nii", warp_type = "ants", inverse_path = "inv.nii")
  expect_false(is_invertible(w_inv))  # "provided" != "exact"

  # VolToSurfMorphism has adjoint
  v2s <- VolToSurfMorphism("a", "b", method = "trilinear")
  expect_false(is_invertible(v2s))
})

test_that("is_invertible errors on non-Morphism", {
  expect_error(is_invertible("not a morphism"), "must be a Morphism")
  expect_error(is_invertible(NULL), "must be a Morphism")
})

test_that("has_adjoint checks inverse_type correctly", {
  # VolToSurfMorphism has adjoint
  v2s <- VolToSurfMorphism("a", "b", method = "trilinear")
  expect_true(has_adjoint(v2s))

  # Other types don't have adjoint
  expect_false(has_adjoint(IdentityMorphism("a")))
  expect_false(has_adjoint(Affine3DMorphism("a", "b", diag(4))))
  expect_false(has_adjoint(Warp3DMorphism("a", "b", "x.nii", warp_type = "ants")))
  expect_false(has_adjoint(SurfToSurfMorphism("a", "b", method = "sphere")))
})

test_that("has_adjoint errors on non-Morphism", {
  expect_error(has_adjoint("not a morphism"), "must be a Morphism")
  expect_error(has_adjoint(list(a = 1)), "must be a Morphism")
})

test_that("morphisms store hash in @hash slot", {
  id <- IdentityMorphism("test")
  expect_true(nzchar(id@hash))

  aff <- Affine3DMorphism("a", "b", diag(4))
  expect_true(nzchar(aff@hash))

  warp <- Warp3DMorphism("a", "b", "dummy.nii", warp_type = "ants")
  expect_true(nzchar(warp@hash))
})

test_that("morphism_hash differs by domain", {
  aff1 <- Affine3DMorphism("domain_a", "domain_b", diag(4))
  aff2 <- Affine3DMorphism("domain_x", "domain_y", diag(4))

  expect_false(morphism_hash(aff1) == morphism_hash(aff2))
})

test_that("morphism_hash differs by matrix", {
  mat1 <- diag(4)
  mat2 <- diag(4)
  mat2[1, 4] <- 10  # Translation

  aff1 <- Affine3DMorphism("a", "b", mat1)
  aff2 <- Affine3DMorphism("a", "b", mat2)

  expect_false(morphism_hash(aff1) == morphism_hash(aff2))
})

test_that("morphism_kind errors on MorphismPath", {
  # MorphismPath is not a base Morphism - it's a container
  # morphism_kind should error since it checks is(object, "Morphism")
  aff1 <- Affine3DMorphism("a", "b", diag(4))
  warp <- Warp3DMorphism("b", "c", "dummy.nii", warp_type = "ants")
  path <- compose(aff1, warp)

  # MorphismPath doesn't pass the is(object, "Morphism") check
  # so morphism_kind will error
  expect_error(morphism_kind(path), "must be a Morphism")
})
