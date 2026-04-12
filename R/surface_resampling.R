#' @title Surface Resampling Plans
#' @name surface_resampling
#' @description
#' Lightweight surface-to-surface resampling using precomputed sparse triplets.
#' This keeps the core API small while enabling reusable vertex-data transport.
NULL

#' Build a surface resampling plan
#'
#' Computes reusable interpolation weights mapping vertex data from a moving
#' mesh onto a reference mesh.
#'
#' @param reference Reference/output surface (`SurfaceMesh` or surface-like object)
#' @param moving Moving/input surface (`SurfaceMesh` or surface-like object)
#' @param method Interpolation method: `"barycentric"` or `"nearest"`
#' @param spherical Logical; if `TRUE`, require both meshes to be approximately spherical
#' @param radius Radius used when `spherical=TRUE` (both meshes are rescaled to this radius)
#' @return A `SurfaceResamplingPlan` object
#' @export
surface_resampling_plan <- function(reference, moving,
                                    method = c("barycentric", "nearest"),
                                    spherical = TRUE,
                                    radius = 100) {
  method <- match.arg(method)

  ref <- if (inherits(reference, "SurfaceMesh")) reference else surface_mesh(reference)
  mov <- if (inherits(moving, "SurfaceMesh")) moving else surface_mesh(moving)

  if (isTRUE(spherical)) {
    if (!mesh_is_sphere(ref)) stop("reference mesh is not approximately spherical")
    if (!mesh_is_sphere(mov)) stop("moving mesh is not approximately spherical")
    ref <- mesh_set_radius(ref, radius = radius)
    mov <- mesh_set_radius(mov, radius = radius)
  }

  n_ref <- nrow(ref@coords)
  n_mov <- nrow(mov@coords)

  if (identical(method, "nearest")) {
    rows <- seq_len(n_ref)
    cols <- cpp_nearest_vertex(ref@coords, mov@coords)
    vals <- rep(1, n_ref)
  } else {
    if (nrow(mov@faces) == 0L) {
      stop("barycentric method requires triangle faces on the moving mesh")
    }

    w <- cpp_barycentric_weights(ref@coords, mov@coords, mov@faces)
    rows <- as.integer(w$rows)
    cols <- as.integer(w$cols)
    vals <- as.numeric(w$vals)

    # Robust fallback: if a point yields no barycentric support, use nearest vertex.
    covered <- rep(FALSE, n_ref)
    if (length(rows) > 0L) covered[rows] <- TRUE
    missing <- which(!covered)
    if (length(missing) > 0L) {
      rows <- c(rows, missing)
      cols <- c(cols, cpp_nearest_vertex(ref@coords[missing, , drop = FALSE], mov@coords))
      vals <- c(vals, rep(1, length(missing)))
    }
  }

  structure(
    list(
      rows = rows,
      cols = cols,
      vals = vals,
      n_reference = as.integer(n_ref),
      n_moving = as.integer(n_mov),
      method = method,
      spherical = isTRUE(spherical),
      radius = as.numeric(radius)
    ),
    class = "SurfaceResamplingPlan"
  )
}

#' Apply a surface resampling plan to vertex data
#'
#' @param plan A `SurfaceResamplingPlan`
#' @param x Vertex data (`length == n_moving`) or matrix (`n_moving x k`)
#' @param inverse Logical; apply inverse direction (`reference -> moving`) by transposing weights
#' @param normalize One of `"element"`, `"sum"`, or `"none"`
#' @return Resampled vector or matrix
#' @export
apply_surface_resampling <- function(plan, x, inverse = FALSE,
                                     normalize = c("element", "sum", "none")) {
  stopifnot(inherits(plan, "SurfaceResamplingPlan"))
  normalize <- match.arg(normalize)

  rows <- if (isTRUE(inverse)) plan$cols else plan$rows
  cols <- if (isTRUE(inverse)) plan$rows else plan$cols
  n_out <- if (isTRUE(inverse)) plan$n_moving else plan$n_reference
  n_in <- if (isTRUE(inverse)) plan$n_reference else plan$n_moving

  x_is_vec <- is.vector(x)
  x_mat <- if (x_is_vec) {
    matrix(as.numeric(x), ncol = 1L)
  } else {
    as.matrix(x)
  }
  if (!is.numeric(x_mat) || nrow(x_mat) != n_in) {
    stop("x must be numeric with nrow equal to the plan input vertex count")
  }

  vals <- .normalize_surface_triplets(rows, cols, plan$vals, n_out, n_in, normalize)

  out <- matrix(0, nrow = n_out, ncol = ncol(x_mat))
  for (k in seq_len(ncol(x_mat))) {
    contrib <- x_mat[cols, k, drop = TRUE] * vals
    tmp <- rowsum(contrib, group = rows, reorder = FALSE)
    out[as.integer(rownames(tmp)), k] <- tmp[, 1]
  }

  if (x_is_vec) out[, 1] else out
}

#' @export
print.SurfaceResamplingPlan <- function(x, ...) {
  cat(
    "<SurfaceResamplingPlan | method=", x$method,
    " | moving=", x$n_moving,
    " -> reference=", x$n_reference,
    " | nnz=", length(x$vals),
    ">\n",
    sep = ""
  )
}

.normalize_surface_triplets <- function(rows, cols, vals, n_out, n_in, normalize) {
  if (identical(normalize, "none")) return(vals)

  if (identical(normalize, "element")) {
    s <- rowsum(vals, group = rows, reorder = FALSE)
    denom <- as.numeric(s[as.character(rows), 1])
    denom[denom == 0] <- 1
    return(vals / denom)
  }

  # normalize == "sum": preserve total mass by normalizing contribution of each input element.
  cs <- rowsum(vals, group = cols, reorder = FALSE)
  cdenom <- as.numeric(cs[as.character(cols), 1])
  cdenom[cdenom == 0] <- 1
  vals / cdenom
}
