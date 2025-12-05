#' @title Jacobian Computations
#' @name jacobian
#' @description
#' Compute Jacobian matrices and determinants for morphisms.
#' The Jacobian is the derivative of the coordinate transformation.
#'
#' @section Semantics:
#' For a morphism with pullback semantics (transform maps target -> source),
#' the pullback Jacobian is d(source)/d(target), evaluated at target coordinates.
#'
#' @section Modes:
#' \itemize{
#'   \item "pullback": d(source)/d(target) - natural for resampling
#'   \item "pushforward": d(target)/d(source) - requires invertible morphism
#' }
NULL

# =============================================================================
# JACOBIAN FIELD CLASS
# =============================================================================

#' Jacobian field: collection of 3x3 Jacobian matrices
#'
#' @slot values Numeric array (n x 3 x 3) of Jacobian matrices
#' @slot coords Numeric matrix (n x 3) of coordinate locations
#' @slot mode Character: "pullback" or "pushforward"
#' @slot morphism_hash Hash of source morphism
#' @param object A JacobianField object (for show method)
#' @param x A JacobianField object (for other methods)
#' @param i Numeric index for subsetting
#' @param j Second index (optional)
#' @param ... Additional arguments
#' @param drop Logical; if TRUE drop dimensions
#' @param a A JacobianField object (for solve method)
#' @param b Missing (solve method signature)
#' @export
setClass("JacobianField",
  slots = c(
    values = "array",
    coords = "matrix",
    mode = "character",
    morphism_hash = "character"
  ),
  prototype = list(
    values = array(NA_real_, dim = c(0, 3, 3)),
    coords = matrix(NA_real_, nrow = 0, ncol = 3),
    mode = "pullback",
    morphism_hash = ""
  )
)

#' @rdname JacobianField-class
#' @export
setMethod("show", "JacobianField", function(object) {
  n <- nrow(object@coords)
  cat(sprintf("<JacobianField | n=%d | mode=%s>\n", n, object@mode))
  if (n > 0) {
    dets <- apply(object@values, 1, det)
    cat(sprintf("  det range: [%.4f, %.4f]\n", min(dets), max(dets)))
  }
})

#' Extract single Jacobian matrix or subset
#' @rdname JacobianField-class
#' @export
setMethod("[", signature("JacobianField", "numeric"), function(x, i, j, ..., drop = TRUE) {
  if (missing(j)) {
    if (length(i) == 1 && drop) {
      return(x@values[i, , , drop = TRUE])
    }
    return(new("JacobianField",
               values = x@values[i, , , drop = FALSE],
               coords = x@coords[i, , drop = FALSE],
               mode = x@mode,
               morphism_hash = x@morphism_hash))
  }
  x@values[i, j, ..., drop = drop]
})

#' Number of Jacobians
#' @rdname JacobianField-class
#' @export
setMethod("length", "JacobianField", function(x) nrow(x@coords))

#' Compute determinants
#' @rdname JacobianField-class
#' @export
setMethod("det", "JacobianField", function(x, ...) {
  apply(x@values, 1, det)
})

#' Invert Jacobians (for mode switching)
#' @rdname JacobianField-class
#' @export
setMethod("solve", signature("JacobianField", "missing"), function(a, b, ...) {
  n <- length(a)
  inv_values <- array(NA_real_, dim = dim(a@values))
  for (i in seq_len(n)) {
    inv_values[i, , ] <- solve(a@values[i, , ])
  }
  new("JacobianField",
      values = inv_values,
      coords = a@coords,
      mode = if (a@mode == "pullback") "pushforward" else "pullback",
      morphism_hash = a@morphism_hash)
})

# =============================================================================
# JACOBIAN GENERIC
# =============================================================================

#' Compute Jacobian matrices at coordinate locations
#'
#' @param morphism A Morphism or MorphismPath
#' @param coords Numeric matrix (n x 3) of coordinates
#' @param mode "pullback" (default) or "pushforward"
#' @return JacobianField object
#'
#' @export
#' @examples
#' aff <- Affine3DMorphism("a", "b", diag(4))
#' coords <- matrix(c(0, 0, 0, 10, 20, 30), ncol = 3, byrow = TRUE)
#' J <- jacobian(aff, coords)
setGeneric("jacobian", function(morphism, coords, mode = c("pullback", "pushforward")) {
  standardGeneric("jacobian")
})

#' Compute Jacobian determinants
#'
#' @param morphism A Morphism or MorphismPath
#' @param coords Numeric matrix (n x 3) of coordinates
#' @param log Return log(|det(J)|)?
#' @param mode "pullback" or "pushforward"
#' @return Numeric vector of determinants
#'
#' @export
#' @examples
#' aff <- Affine3DMorphism("a", "b", diag(4) * 2)
#' coords <- matrix(c(0, 0, 0), ncol = 3)
#' jacobian_det(aff, coords)  # 8 (2^3)
#' jacobian_det(aff, coords, log = TRUE)  # log(8)
setGeneric("jacobian_det", function(morphism, coords, log = FALSE,
                                     mode = c("pullback", "pushforward")) {
  standardGeneric("jacobian_det")
})

# =============================================================================
# IDENTITY MORPHISM
# =============================================================================

#' @rdname jacobian
#' @export
setMethod("jacobian", "IdentityMorphism", function(morphism, coords,
                                                    mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)
  n <- nrow(coords)
  values <- array(0, dim = c(n, 3, 3))
  for (i in seq_len(n)) {
    values[i, , ] <- diag(3)
  }
  new("JacobianField",
      values = values,
      coords = coords,
      mode = mode,
      morphism_hash = morphism@hash)
})

#' @rdname jacobian_det
#' @export
setMethod("jacobian_det", "IdentityMorphism", function(morphism, coords, log = FALSE,
                                                        mode = c("pullback", "pushforward")) {
  n <- nrow(coords)
  if (log) rep(0, n) else rep(1, n)
})

# =============================================================================
# AFFINE MORPHISM
# =============================================================================

#' @rdname jacobian
#' @export
setMethod("jacobian", "Affine3DMorphism", function(morphism, coords,
                                                    mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)
  n <- nrow(coords)

  # Linear part is the Jacobian (constant everywhere)
  linear <- morphism@matrix[1:3, 1:3]
  J <- if (mode == "pullback") linear else solve(linear)

  values <- array(0, dim = c(n, 3, 3))
  for (i in seq_len(n)) {
    values[i, , ] <- J
  }

  new("JacobianField",
      values = values,
      coords = coords,
      mode = mode,
      morphism_hash = morphism@hash)
})

#' @rdname jacobian_det
#' @export
setMethod("jacobian_det", "Affine3DMorphism", function(morphism, coords, log = FALSE,
                                                        mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)
  linear <- morphism@matrix[1:3, 1:3]
  d <- det(linear)
  if (mode == "pushforward") d <- 1 / d

  n <- nrow(coords)
  if (log) rep(log(abs(d)), n) else rep(d, n)
})

# =============================================================================
# WARP MORPHISM
# =============================================================================

#' @rdname jacobian
#' @export
setMethod("jacobian", "Warp3DMorphism", function(morphism, coords,
                                                  mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)

  coords_in <- coords

  # Load warp field
  cache_env <- morphism@cache %||% new_cache_env()
  warp <- load_warp_array(morphism, cache_env = cache_env)

  method <- morphism@params$warp_method %||% "linear"

  # Compute via C++ (finite differences on displacement field)
  values <- if (identical(method, "cubic")) {
    cpp_warp_jacobian_cubic(
      coords = coords_in,
      field = warp$array,
      dim = warp$dim,
      world_to_vox = warp$world_to_vox,
      vox_to_world = warp$vox_to_world
    )
  } else {
    cpp_warp_jacobian(
      coords = coords_in,
      field = warp$array,
      dim = warp$dim,
      world_to_vox = warp$world_to_vox,
      vox_to_world = warp$vox_to_world
    )
  }

  jf <- new("JacobianField",
            values = values,
            coords = coords,
            mode = "pullback",
            morphism_hash = morphism@hash)

  if (mode == "pushforward") {
    if (morphism@inverse_type %in% c("none", "adjoint")) {
      stop("Pushforward Jacobian requires invertible morphism")
    }
    jf <- solve(jf)
  }

  jf
})

#' @rdname jacobian_det
#' @export
setMethod("jacobian_det", "Warp3DMorphism", function(morphism, coords, log = FALSE,
                                                      mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)

  coords_in <- coords

  cache_env <- morphism@cache %||% new_cache_env()
  warp <- load_warp_array(morphism, cache_env = cache_env)

  method <- morphism@params$warp_method %||% "linear"

  dets <- if (identical(method, "cubic")) {
    cpp_warp_jacobian_det_cubic(
      coords = coords_in,
      field = warp$array,
      dim = warp$dim,
      world_to_vox = warp$world_to_vox,
      vox_to_world = warp$vox_to_world
    )
  } else {
    cpp_warp_jacobian_det(
      coords = coords_in,
      field = warp$array,
      dim = warp$dim,
      world_to_vox = warp$world_to_vox,
      vox_to_world = warp$vox_to_world
    )
  }

  if (mode == "pushforward") dets <- 1 / dets
  if (log) log(abs(dets)) else dets
})

# =============================================================================
# MORPHISM PATH (CHAIN RULE)
# =============================================================================

#' Multiply JacobianFields element-wise
#' @keywords internal
jacobian_multiply <- function(J1, J2) {
  stopifnot(length(J1) == length(J2))
  n <- length(J1)
  result <- array(NA_real_, dim = c(n, 3, 3))
  for (i in seq_len(n)) {
    result[i, , ] <- J1@values[i, , ] %*% J2@values[i, , ]
  }
  new("JacobianField",
      values = result,
      coords = J2@coords,
      mode = J1@mode,
      morphism_hash = compute_hash(J1@morphism_hash, J2@morphism_hash))
}

#' @rdname jacobian
#' @export
setMethod("jacobian", "MorphismPath", function(morphism, coords,
                                                mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)
  path <- morphism@morphisms

  if (length(path) == 0) {
    return(jacobian(IdentityMorphism(morphism@source), coords, mode))
  }

  # For pullback: J(g . f)(x) = J(f)(g(x)) %*% J(g)(x)
  # Apply path in reverse, collect intermediate coords

  current_coords <- coords
  intermediate <- vector("list", length(path))
  intermediate[[length(path)]] <- current_coords

  for (i in rev(seq_along(path))) {
    if (i > 1) {
      current_coords <- transform(path[[i]], current_coords)
      intermediate[[i - 1]] <- current_coords
    }
  }

  # Multiply Jacobians in reverse order
  result <- NULL
  for (i in rev(seq_along(path))) {
    coords_for_J <- intermediate[[i]]
    J_i <- jacobian(path[[i]], coords_for_J, mode = "pullback")

    if (is.null(result)) {
      result <- J_i
    } else {
      result <- jacobian_multiply(J_i, result)
    }
  }

  result@mode <- mode
  if (mode == "pushforward") {
    result <- solve(result)
  }

  result
})

#' @rdname jacobian_det
#' @export
setMethod("jacobian_det", "MorphismPath", function(morphism, coords, log = FALSE,
                                                    mode = c("pullback", "pushforward")) {
  J <- jacobian(morphism, coords, mode = mode)
  dets <- det(J)
  if (log) log(abs(dets)) else dets
})

# =============================================================================
# SURFACE MORPHISMS (NOT IMPLEMENTED)
# =============================================================================

#' @rdname jacobian
#' @export
setMethod("jacobian", "VolToSurfMorphism", function(morphism, coords, mode) {
  stop("Jacobian not defined for VolToSurfMorphism.\n",
       "Volume-to-surface projections don't preserve dimensionality.")
})

#' @rdname jacobian
#' @export
setMethod("jacobian", "SurfToSurfMorphism", function(morphism, coords, mode) {
  stop("Jacobian not implemented for SurfToSurfMorphism.\n",
       "Surface mappings require tangent space computations.")
})
