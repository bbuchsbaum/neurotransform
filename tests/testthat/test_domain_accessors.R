test_that("domain accessors work on morphisms", {
  m <- Affine3DMorphism("a", "b", diag(4))
  expect_equal(source_domain(m), "a")
  expect_equal(target_domain(m), "b")
  expect_equal(source_of(m), "a")
  expect_equal(target_of(m), "b")
})
