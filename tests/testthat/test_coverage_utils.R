# Coverage Tests for aaa_utils.R
#
# Target: Improve coverage from 50% to >80%

# ==============================================================================
# NULL COALESCING OPERATOR
# ==============================================================================

test_that("%||% returns left if not NULL", {
  result <- "value" %||% "default"
  expect_equal(result, "value")

  result <- 42 %||% 0
  expect_equal(result, 42)

  result <- FALSE %||% TRUE
  expect_equal(result, FALSE)
})
test_that("%||% returns right if left is NULL", {
  result <- NULL %||% "default"
  expect_equal(result, "default")

  result <- NULL %||% 42
  expect_equal(result, 42)
})

# ==============================================================================
# COMPUTE HASH
# ==============================================================================

test_that("compute_hash produces consistent hashes", {
  h1 <- compute_hash("test")
  h2 <- compute_hash("test")
  expect_equal(h1, h2)
})

test_that("compute_hash produces different hashes for different inputs", {
  h1 <- compute_hash("test1")
  h2 <- compute_hash("test2")
  expect_false(h1 == h2)
})

test_that("compute_hash handles multiple arguments", {
  h1 <- compute_hash("a", "b", "c")
  h2 <- compute_hash("a", "b", "c")
  expect_equal(h1, h2)

  h3 <- compute_hash("a", "b")
  expect_false(h1 == h3)
})

test_that("compute_hash handles complex objects", {
  obj <- list(a = 1, b = c(2, 3), c = matrix(1:4, 2))
  h1 <- compute_hash(obj)
  h2 <- compute_hash(obj)
  expect_equal(h1, h2)
})

# ==============================================================================
# CHECK SLOTS
# ==============================================================================

test_that("check_slots returns TRUE for valid objects", {
  # Create a simple S4 object
  aff <- Affine3DMorphism("a", "b", diag(4))

  result <- neurotransform:::check_slots(aff, c("source", "target", "matrix"))
  expect_true(result)
})

test_that("check_slots returns error message for missing slots", {
  aff <- Affine3DMorphism("a", "b", diag(4))

  result <- neurotransform:::check_slots(aff, c("source", "nonexistent_slot"))
  expect_true(is.character(result))
  expect_true(grepl("Missing", result))
})

# ==============================================================================
# NEW CACHE ENV
# ==============================================================================

test_that("new_cache_env creates environment", {
  env <- neurotransform:::new_cache_env()
  expect_true(is.environment(env))
})

test_that("new_cache_env creates isolated environment", {
  env <- neurotransform:::new_cache_env()
  expect_identical(parent.env(env), emptyenv())
})

test_that("new_cache_env environments can store values", {
  env <- neurotransform:::new_cache_env()
  assign("key", "value", envir = env)
  expect_equal(get("key", envir = env), "value")
})

# ==============================================================================
# VALIDATE 4X4 MATRIX
# ==============================================================================

test_that("validate_4x4_matrix accepts valid 4x4 matrix", {
  mat <- diag(4)
  result <- neurotransform:::validate_4x4_matrix(mat)
  expect_true(result)
})

test_that("validate_4x4_matrix rejects non-matrix", {
  expect_error(neurotransform:::validate_4x4_matrix(c(1, 2, 3, 4)), "must be a matrix")
})

test_that("validate_4x4_matrix rejects wrong dimensions", {
  mat <- diag(3)
  expect_error(neurotransform:::validate_4x4_matrix(mat), "must be 4x4")

  mat <- matrix(1:20, nrow = 4)
  expect_error(neurotransform:::validate_4x4_matrix(mat), "must be 4x4")
})

test_that("validate_4x4_matrix uses custom name in error", {
  mat <- diag(3)
  expect_error(neurotransform:::validate_4x4_matrix(mat, "my_affine"), "my_affine")
})

# ==============================================================================
# VALIDATE POSITIVE INTEGER
# ==============================================================================

test_that("validate_positive_integer accepts valid input", {
  result <- neurotransform:::validate_positive_integer(5)
  expect_true(result)

  result <- neurotransform:::validate_positive_integer(1L)
  expect_true(result)
})

test_that("validate_positive_integer rejects non-positive", {
  expect_error(neurotransform:::validate_positive_integer(0), "positive integer")
  expect_error(neurotransform:::validate_positive_integer(-5), "positive integer")
})

test_that("validate_positive_integer rejects non-integer", {
  expect_error(neurotransform:::validate_positive_integer(5.5), "positive integer")
})

test_that("validate_positive_integer rejects non-numeric", {
  expect_error(neurotransform:::validate_positive_integer("5"), "positive integer")
})

test_that("validate_positive_integer rejects vectors", {
  expect_error(neurotransform:::validate_positive_integer(c(1, 2)), "positive integer")
})

# ==============================================================================
# VALIDATE NUMERIC VECTOR
# ==============================================================================

test_that("validate_numeric_vector accepts valid input", {
  result <- neurotransform:::validate_numeric_vector(c(1, 2, 3), 3)
  expect_true(result)
})

test_that("validate_numeric_vector rejects wrong length", {
  expect_error(neurotransform:::validate_numeric_vector(c(1, 2), 3), "length 3")
})

test_that("validate_numeric_vector rejects non-numeric", {
  expect_error(neurotransform:::validate_numeric_vector(c("a", "b", "c"), 3), "numeric")
})

test_that("validate_numeric_vector uses custom name in error", {
  expect_error(neurotransform:::validate_numeric_vector(c(1, 2), 3, "coords"), "coords")
})
