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

#' Coerce surface-like objects to vertex/face matrices
#'
#' Supports plain matrices, internal SurfaceMesh/SampledPoints, and optional
#' neurosurf::SurfaceGeometry objects without requiring neurosurf as a hard
#' dependency.
#'
#' @param x Surface-like object
#' @param require_faces Logical; whether faces must be present
#' @param apply_surf_to_world Logical; apply SurfaceGeometry surf_to_world
#' @return List with \code{vertices} and \code{faces}
#' @keywords internal
coerce_surface_reference <- function(x, require_faces = FALSE, apply_surf_to_world = TRUE) {
  if (is.matrix(x)) {
    if (!is.numeric(x) || ncol(x) != 3L || !all(is.finite(x))) {
      stop("Surface vertices must be a finite numeric matrix with 3 columns")
    }
    return(list(vertices = x, faces = NULL))
  }

  if (inherits(x, "SurfaceMesh")) {
    return(list(vertices = x@coords, faces = x@faces))
  }
  if (inherits(x, "SampledPoints")) {
    return(list(vertices = x@coords, faces = NULL))
  }

  if (inherits(x, "SurfaceGeometry")) {
    verts <- NULL
    fcs <- NULL

    # Fast path for standard S4 SurfaceGeometry internals.
    if (isS4(x) && methods::.hasSlot(x, "mesh")) {
      mesh <- methods::slot(x, "mesh")
      if (is.list(mesh) && !is.null(mesh$vb)) {
        verts <- t(mesh$vb[1:3, , drop = FALSE])
      }
      if (is.list(mesh) && !is.null(mesh$it)) {
        fcs <- t(mesh$it)
      }
      if (isTRUE(apply_surf_to_world) && methods::.hasSlot(x, "surf_to_world") && !is.null(verts)) {
        stw <- methods::slot(x, "surf_to_world")
        if (is.matrix(stw) && identical(dim(stw), c(4L, 4L))) {
          verts <- apply_affine(verts, stw)
        }
      }
    }

    # Allow light-weight list-like SurfaceGeometry objects in tests/integration glue.
    if (is.null(verts) && is.list(x) && !is.null(x$mesh)) {
      mesh <- x$mesh
      if (is.list(mesh) && !is.null(mesh$vb)) {
        verts <- t(mesh$vb[1:3, , drop = FALSE])
      }
      if (is.list(mesh) && !is.null(mesh$it)) {
        fcs <- t(mesh$it)
      }
      if (isTRUE(apply_surf_to_world) && !is.null(verts) &&
          is.matrix(x$surf_to_world) && identical(dim(x$surf_to_world), c(4L, 4L))) {
        verts <- apply_affine(verts, x$surf_to_world)
      }
    }

    # Fallback via neurosurf accessors if needed.
    if ((is.null(verts) || (require_faces && is.null(fcs))) &&
        requireNamespace("neurosurf", quietly = TRUE)) {
      ns <- asNamespace("neurosurf")
      if (is.null(verts) && exists("coords", envir = ns, inherits = FALSE)) {
        verts <- get("coords", envir = ns)(x)
      }
      if (is.null(fcs) && exists("faces", envir = ns, inherits = FALSE)) {
        fcs <- get("faces", envir = ns)(x)
      }
      if (isTRUE(apply_surf_to_world) && !is.null(verts) && exists("surf_to_world", envir = ns, inherits = FALSE)) {
        stw <- tryCatch(get("surf_to_world", envir = ns)(x), error = function(e) NULL)
        if (is.matrix(stw) && identical(dim(stw), c(4L, 4L))) {
          verts <- apply_affine(verts, stw)
        }
      }
    }

    if (is.null(verts) || !is.matrix(verts) || ncol(verts) != 3L || !is.numeric(verts) || !all(is.finite(verts))) {
      stop("Could not extract finite Nx3 vertex coordinates from SurfaceGeometry")
    }

    if (!is.null(fcs)) {
      if (!is.matrix(fcs) || ncol(fcs) != 3L || !is.numeric(fcs) || !all(is.finite(fcs))) {
        stop("Surface faces must be a finite numeric matrix with 3 columns")
      }
      f0 <- matrix(as.integer(round(fcs)), ncol = 3)
      if (any(abs(fcs - f0) > 1e-8)) stop("Surface faces must be integer indices")
      if (nrow(f0) > 0L && min(f0) >= 1L) f0 <- f0 - 1L
      if (nrow(f0) > 0L && (min(f0) < 0L || max(f0) >= nrow(verts))) {
        stop("Surface faces are out of bounds for extracted vertices")
      }
      fcs <- f0
    }

    if (isTRUE(require_faces) && (is.null(fcs) || nrow(fcs) == 0L)) {
      stop("SurfaceGeometry does not contain usable triangle faces")
    }

    return(list(vertices = verts, faces = fcs))
  }

  stop("Unsupported surface reference type: ", paste(class(x), collapse = "/"))
}
