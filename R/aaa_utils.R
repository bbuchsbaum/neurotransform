#' @title Internal Utility Functions
#' @name aaa_utils
#' @description
#' Internal utility functions for the neurotransform package.
#' This file is loaded first (via Collate) to ensure utilities are available.
#' @keywords internal
NULL

#' Null-coalescing operator
#'
#' @param x Left operand
#' @param y Right operand (used if x is NULL)
#' @return x if not NULL, otherwise y
#' @name null-coalesce
#' @rdname null-coalesce
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Compute a stable hash for objects
#'
#' Creates a reproducible hash using digest with serialize=TRUE.
#'
#' @param ... Objects to hash
#' @param algo Hash algorithm (default "xxhash64" for speed)
#' @return Character string hash
#' @export
#' @examples
#' compute_hash("hello", "world")
#' compute_hash(list(a = 1, b = 2))
compute_hash <- function(..., algo = "xxhash64") {
  digest::digest(list(...), algo = algo, serialize = TRUE)
}

#' Check if object has required slots
#'
#' @param object S4 object to check
#' @param slots Character vector of required slot names
#' @return TRUE if valid, character error message otherwise
#' @keywords internal
check_slots <- function(object, slots) {
  for (s in slots) {
    if (!methods::.hasSlot(object, s)) {
      return(sprintf("Missing required slot: %s", s))
    }
  }
  TRUE
}

#' Create a new environment for caching
#'
#' @return A new environment with hash=TRUE for efficient lookup
#' @keywords internal
new_cache_env <- function() {
  new.env(hash = TRUE, parent = emptyenv())
}

#' Validate that a matrix is 4x4 (for affine transforms)
#'
#' @param mat Matrix to validate
#' @param name Name for error messages
#' @return TRUE if valid, error otherwise
#' @keywords internal
validate_4x4_matrix <- function(mat, name = "matrix") {
 if (!is.matrix(mat)) {
    stop(sprintf("%s must be a matrix", name))
  }
  if (!identical(dim(mat), c(4L, 4L))) {
    stop(sprintf("%s must be 4x4, got %dx%d", name, nrow(mat), ncol(mat)))
  }
  TRUE
}

#' Validate positive integer
#'
#' @param x Value to validate
#' @param name Name for error messages
#' @return TRUE if valid, error otherwise
#' @keywords internal
validate_positive_integer <- function(x, name = "value") {
  if (!is.numeric(x) || length(x) != 1 || x <= 0 || x != floor(x)) {
    stop(sprintf("%s must be a positive integer", name))
  }
  TRUE
}

#' Validate numeric vector of specific length
#'
#' @param x Vector to validate
#' @param len Expected length
#' @param name Name for error messages
#' @return TRUE if valid, error otherwise
#' @keywords internal
validate_numeric_vector <- function(x, len, name = "vector") {
  if (!is.numeric(x) || length(x) != len) {
    stop(sprintf("%s must be numeric vector of length %d", name, len))
  }
  TRUE
}
