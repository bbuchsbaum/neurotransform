#' @title Surface Resampling Plans
#' @name surface_resampling
#' @description
#' Lightweight surface-to-surface resampling using precomputed sparse triplets.
#' This keeps the core API small while enabling reusable vertex-data transport.
NULL

#' Validate a surface mesh
#'
#' Reports invalid faces and topology without repairing them. Spherical
#' admission additionally requires origin-based radius spread within tolerance,
#' a closed connected orientable manifold of sphere topology, unambiguous radial
#' face orientations and radial mapping degree one. Input face winding is free.
#' Generic open meshes do not require closed topology or spherical geometry.
#'
#' @param mesh SurfaceMesh or surface-like object with faces.
#' @param spherical Require an embedded spherical triangulation.
#' @param error Stop when validation fails instead of returning diagnostics.
#' @param radius_tolerance Maximum ratio of vertex radii about the origin.
#' @return Diagnostic list with valid, issues (one-based offending indices),
#'   topology counts, radius range, radial degree and numerical margins.
#' @export
validate_surface_mesh <- function(mesh, spherical = TRUE, error = FALSE,
                                  radius_tolerance = 1.001) {
  mesh <- if (inherits(mesh, "SurfaceMesh")) mesh else surface_mesh(mesh)
  if (length(spherical) != 1L || is.na(spherical) || !is.logical(spherical))
    stop("spherical must be one logical value")
  if (length(radius_tolerance) != 1L || !is.finite(radius_tolerance) || radius_tolerance <= 1)
    stop("radius_tolerance must be finite and greater than one")
  result <- cpp_validate_surface(mesh@coords, mesh@faces, spherical, radius_tolerance)
  if (isTRUE(error) && !result$valid) {
    stop("invalid ", if (spherical) "spherical " else "", "mesh: ",
         paste(names(result$issues), collapse = ", "), call. = FALSE)
  }
  result
}

#' Build a surface resampling plan
#'
#' Computes reusable interpolation weights mapping vertex data from a moving
#' mesh onto a reference mesh.
#'
#' Spherical barycentric resampling uses the closest point on the rescaled
#' triangular mesh, including edges and vertices (ordinary barycentric, without
#' adaptive area correction). It does not use radial ray intersections. For
#' nonspherical meshes, only interior orthogonal projections are considered,
#' with an explicit policy for uncovered queries. The default is nearest-vertex
#' fallback for nonspherical meshes and an error for spherical meshes.
#'
#' @param reference Reference/output surface (`SurfaceMesh` or surface-like object)
#' @param moving Moving/input surface (`SurfaceMesh` or surface-like object)
#' @param source_mask Optional logical or finite numeric input ROI; positive values
#'   include vertices. Exclusion is applied before row normalization.
#' @param target_mask Optional output ROI, applied after normalization.
#' @param outside Uncovered-query policy: `"error"`, `"missing"`, or `"nearest"`.
#'   Defaults to `"error"` for spherical and `"nearest"` for nonspherical meshes.
#' @param method Interpolation method: `"barycentric"` or `"nearest"`
#' @param spherical Logical; if `TRUE`, require both meshes to be approximately spherical
#' @param radius Radius used when `spherical=TRUE` (both meshes are rescaled to this radius)
#' @return A `SurfaceResamplingPlan` object
#' @export
surface_resampling_plan <- function(reference, moving,
                                    method = c("barycentric", "nearest"),
                                    spherical = TRUE,
                                    radius = 100, source_mask = NULL,
                                    target_mask = NULL, outside = NULL) {
  if (is.null(outside)) outside <- if (isTRUE(spherical)) "error" else "nearest"
  outside <- match.arg(outside, c("error", "missing", "nearest"))
  method <- match.arg(method)

  ref <- if (inherits(reference, "SurfaceMesh")) reference else surface_mesh(reference)
  mov <- if (inherits(moving, "SurfaceMesh")) moving else surface_mesh(moving)

  geometry <- list(reference = .surface_identity(ref), moving = .surface_identity(mov))
  source_mask <- .surface_mask(source_mask, nrow(mov@coords), "source_mask")
  target_mask <- .surface_mask(target_mask, nrow(ref@coords), "target_mask")

  if (isTRUE(spherical)) {
    if (length(radius) != 1L || !is.finite(radius) || radius <= 0)
      stop("radius must be finite and positive")
    if (!mesh_is_sphere(ref)) stop("reference mesh is not approximately spherical")
    if (!mesh_is_sphere(mov)) stop("moving mesh is not approximately spherical")
    if (method != "nearest" || nrow(mov@faces)) validate_surface_mesh(mov, spherical = TRUE, error = TRUE)
    if (nrow(ref@faces)) validate_surface_mesh(ref, spherical = TRUE, error = TRUE)
    ref <- mesh_set_radius(ref, radius = radius)
    mov <- mesh_set_radius(mov, radius = radius)
  }

  n_ref <- nrow(ref@coords)
  n_mov <- nrow(mov@coords)
  timing <- NULL
  support <- rep("none", n_ref)

  if (identical(method, "nearest")) {
    rows <- seq_len(n_ref)
    cols <- cpp_nearest_vertex(ref@coords, mov@coords)
    vals <- rep(1, n_ref)
    support[] <- "nearest"
  } else {
    if (nrow(mov@faces) == 0L) {
      stop("barycentric method requires triangle faces on the moving mesh")
    }

    w <- cpp_barycentric_weights(ref@coords, mov@coords, mov@faces,
                                 closest = isTRUE(spherical))
    rows <- as.integer(w$rows)
    cols <- as.integer(w$cols)
    vals <- as.numeric(w$vals)
    timing <- w$timing

    # Preserve triangle support separately from any explicitly allowed fallback.
    covered <- rep(FALSE, n_ref)
    if (length(rows) > 0L) covered[rows] <- TRUE
    support[covered] <- "triangle"
    missing <- which(!covered)
    if (length(missing) > 0L) {
      if (outside == "error") stop("resampling has queries without valid triangle support")
      if (outside == "nearest") {
        rows <- c(rows, missing)
        cols <- c(cols, cpp_nearest_vertex(ref@coords[missing, , drop = FALSE], mov@coords))
        vals <- c(vals, rep(1, length(missing)))
        support[missing] <- "nearest_fallback"
      }
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
      radius = as.numeric(radius),
      timing = timing,
      geometry = geometry,
      support = support,
      source_mask = source_mask,
      target_mask = target_mask,
      outside = outside
    ),
    class = "SurfaceResamplingPlan"
  )
}

#' Apply a surface resampling plan to vertex data
#'
#' Geometry, source ROI mass and finite-data mass are distinct. Unsupported or
#' masked outputs are unavailable (`NA`), while supported numerical zero remains
#' zero. Zero-weight missing inputs never contaminate an output. Omission returns
#' the retained raw weight mass as an attribute, or in the detailed result.
#'
#' Label keys are never averaged. Aggregate voting sums weights for each key;
#' largest voting selects one source vertex. Ties select the smallest label key
#' for aggregate voting and the smallest source vertex index for largest voting.
#'
#' @param plan A `SurfaceResamplingPlan`.
#' @param x Numeric vertex vector or matrix with `n_moving` rows.
#' @param inverse Legacy transpose followed by normalization. Deprecated; use
#'   `apply_surface_adjoint()` or construct a reverse plan explicitly.
#' @param normalize `"element"` normalizes each output row; `"sum"` normalizes
#'   represented input columns; `"none"` retains stored weights. Column
#'   normalization preserves represented discrete mass before target masking;
#'   it is not anatomical area correction.
#' @param na_policy `"propagate"`, `"omit"` with finite-weight renormalization,
#'   or `"error"`. Omission requires `normalize="element"`.
#' @param data_type `"continuous"` or `"label"` (integer keys).
#' @param label_method `"aggregate"` or `"largest"`.
#' @param label_table Optional data frame with a unique integer `key` column.
#'   Other columns, such as names and colors, are preserved unchanged.
#' @param unassigned Integer output key for unavailable labels, or `NA`.
#' @param details Return values and per-column coverage diagnostics.
#' @return Vector/matrix, or a list when `details=TRUE`. Diagnostics include raw
#'   source and finite-data weight mass before normalization and target masking.
#' @export
apply_surface_resampling <- function(plan, x, inverse = FALSE,
                                     normalize = c("element", "sum", "none"),
                                     na_policy = c("propagate", "omit", "error"),
                                     data_type = c("continuous", "label"),
                                     label_method = c("aggregate", "largest"),
                                     label_table = NULL, unassigned = NA_integer_,
                                     details = FALSE) {
  stopifnot(inherits(plan, "SurfaceResamplingPlan"))
  normalize <- match.arg(normalize)
  na_policy <- match.arg(na_policy)
  data_type <- match.arg(data_type)
  label_method <- match.arg(label_method)
  if (isTRUE(inverse)) {
    if (na_policy != "propagate" || data_type != "continuous" || details ||
        !is.null(label_table) || !all(is.na(unassigned)))
      stop("legacy inverse cannot be combined with new data policies")
    if ((!is.null(plan$source_mask) && !all(plan$source_mask)) ||
        (!is.null(plan$target_mask) && !all(plan$target_mask)))
      stop("legacy inverse cannot be combined with fixed masks; use apply_surface_adjoint")
    warning("inverse=TRUE is deprecated; use apply_surface_adjoint() or a reverse plan", call. = FALSE)
    return(.surface_legacy_inverse(plan, x, normalize))
  }
  if (na_policy == "omit" && normalize != "element")
    stop("na_policy='omit' requires normalize='element'")
  input <- .surface_data(x, plan$n_moving)
  x_mat <- input$matrix
  if (data_type == "label") {
    .surface_label_keys(x_mat[is.finite(x_mat)], "x")
    if (length(unassigned) != 1L || (!is.na(unassigned) && !is.finite(unassigned)))
      stop("unassigned must be one integer key or NA")
    .surface_label_keys(unassigned[!is.na(unassigned)], "unassigned")
    if (!is.null(label_table)) {
      if (!is.data.frame(label_table) || !"key" %in% names(label_table) ||
          anyNA(label_table$key) || anyDuplicated(label_table$key))
        stop("label_table must have unique nonmissing integer keys")
      .surface_label_keys(label_table$key, "label_table$key")
      if (!all(c(x_mat[is.finite(x_mat)], unassigned[!is.na(unassigned)]) %in% label_table$key))
        stop("label_table does not contain every input and unassigned key")
    }
  }
  op <- .surface_operator(plan, normalize)
  n <- plan$n_reference
  source_mass <- .surface_rowsum(op$raw, op$rows, n)
  out <- finite_mass <- matrix(0, n, ncol(x_mat))
  status <- matrix("available", n, ncol(x_mat))
  for (k in seq_len(ncol(x_mat))) {
    finite <- is.finite(x_mat[op$cols, k])
    finite_mass[, k] <- .surface_rowsum(op$raw[finite], op$rows[finite], n)
    missing_mass <- .surface_rowsum(op$vals[!finite], op$rows[!finite], n)
    if (na_policy == "error" && any(missing_mass > 0))
      stop("nonfinite input has positive weight in column ", k)
    rows <- op$rows[finite]
    cols <- op$cols[finite]
    weights <- op$vals[finite]
    if (na_policy == "omit")
      weights <- .normalize_surface_triplets(rows, cols, weights, n, plan$n_moving, "element")
    if (data_type == "continuous") {
      out[, k] <- .surface_rowsum(weights * x_mat[cols, k], rows, n)
    } else {
      out[, k] <- .surface_vote(rows, cols, weights, x_mat[cols, k], n, label_method)
    }
    status[source_mass == 0, k] <- "no_source_support"
    status[finite_mass[, k] == 0 & source_mass > 0, k] <- "missing_data"
    if (na_policy == "propagate") status[missing_mass > 0, k] <- "missing_data"
    status[op$support == "none", k] <- "no_geometry"
    status[!op$target_mask, k] <- "masked_target"
    out[status[, k] != "available", k] <- if (data_type == "label") unassigned else NA_real_
  }
  shape <- function(z) if (input$vector) z[, 1L] else z
  values <- shape(out)
  if (data_type == "label") {
    storage.mode(values) <- "integer"
    if (!is.null(label_table)) attr(values, "label_table") <- label_table
  }
  if (details) return(list(values = values, support = op$support,
    geometric_support = op$support == "triangle", source_weight_mass = source_mass,
    finite_weight_mass = shape(finite_mass), target_mask = op$target_mask,
    available = shape(status == "available"), status = shape(status), label_table = label_table))
  if (na_policy == "omit") attr(values, "retained_weight_mass") <- shape(finite_mass)
  values
}

.surface_identity <- function(mesh) digest::digest(list(mesh@coords, mesh@faces), algo = "sha256")

.surface_mask <- function(mask, n, name) {
  if (is.null(mask)) return(rep(TRUE, n))
  if (!(is.logical(mask) || is.numeric(mask)) || !is.null(dim(mask)) ||
      length(mask) != n || any(!is.finite(mask)))
    stop(name, " must be a finite logical or numeric vector of length ", n)
  mask > 0
}

.surface_data <- function(x, n) {
  if (!is.numeric(x) || (!is.null(dim(x)) && !is.matrix(x)))
    stop("x must be a numeric vector or matrix")
  vector <- is.null(dim(x))
  mat <- if (vector) matrix(x, ncol = 1L) else x
  if (nrow(mat) != n) stop("x must have the plan input vertex count")
  list(matrix = mat, vector = vector)
}

.surface_rowsum <- function(values, rows, n) {
  out <- numeric(n)
  if (length(rows)) {
    sums <- rowsum(values, rows, reorder = FALSE)
    out[as.integer(rownames(sums))] <- sums[, 1L]
  }
  out
}

.surface_operator <- function(plan, normalize) {
  source <- .surface_mask(plan$source_mask, plan$n_moving, "source_mask")
  target <- .surface_mask(plan$target_mask, plan$n_reference, "target_mask")
  keep <- plan$vals > 0 & source[plan$cols]
  rows <- plan$rows[keep]
  cols <- plan$cols[keep]
  raw <- plan$vals[keep]
  vals <- .normalize_surface_triplets(rows, cols, raw, plan$n_reference, plan$n_moving, normalize)
  vals[!target[rows]] <- 0
  support <- plan$support
  if (is.null(support)) {
    support <- rep("none", plan$n_reference)
    support[plan$rows[plan$vals > 0]] <- if (plan$method == "nearest") "nearest" else "triangle"
  }
  list(rows = rows, cols = cols, raw = raw, vals = vals, support = support, target_mask = target)
}

.surface_label_keys <- function(x, name) {
  if (!is.numeric(x) || any(!is.finite(x) | x != trunc(x) | abs(x) > .Machine$integer.max))
    stop(name, " must contain integer label keys")
}

.surface_vote <- function(rows, cols, weights, keys, n, method) {
  out <- rep(NA_real_, n)
  positive <- weights > 0
  rows <- rows[positive]; cols <- cols[positive]
  weights <- weights[positive]; keys <- keys[positive]
  if (!length(rows)) return(out)
  if (method == "aggregate") {
    ord <- order(rows, keys)
    rows <- rows[ord]; keys <- keys[ord]; weights <- weights[ord]
    first <- c(TRUE, diff(rows) != 0 | diff(keys) != 0)
    weights <- as.numeric(rowsum(weights, cumsum(first), reorder = FALSE))
    rows <- rows[first]; keys <- keys[first]
    ord <- order(rows, -weights, keys)
  } else ord <- order(rows, -weights, cols)
  selected <- ord[!duplicated(rows[ord])]
  out[rows[selected]] <- keys[selected]
  out
}

.surface_legacy_inverse <- function(plan, x, normalize) {
  input <- .surface_data(x, plan$n_reference)
  vals <- .normalize_surface_triplets(plan$cols, plan$rows, plan$vals,
                                     plan$n_moving, plan$n_reference, normalize)
  out <- matrix(0, plan$n_moving, ncol(input$matrix))
  for (k in seq_len(ncol(out)))
    out[, k] <- .surface_rowsum(input$matrix[plan$rows, k] * vals, plan$cols, plan$n_moving)
  if (input$vector) out[, 1L] else out
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
  if (!length(vals) || identical(normalize, "none")) return(vals)

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

#' Apply the Euclidean adjoint of a surface plan
#'
#' Applies the transpose of the exact forward operator, including fixed masks
#' and the selected normalization, without renormalizing its transpose. This is
#' an algebraic adjoint, not reverse geometric interpolation or reconstruction.
#' Unsupported source columns receive algebraic zero; `details=TRUE` identifies
#' them. Missing-data omission is not a fixed linear operator and is unsupported.
#'
#' @param plan A `SurfaceResamplingPlan`.
#' @param x Numeric vector or matrix with one row per reference vertex.
#' @param normalize Forward-operator normalization: `"element"`, `"sum"`, or `"none"`.
#' @param details Return values and supported-source-column diagnostics.
#' @return Vector/matrix or a list with `values` and `available`.
#' @export
apply_surface_adjoint <- function(plan, x, normalize = c("element", "sum", "none"),
                                  details = FALSE) {
  stopifnot(inherits(plan, "SurfaceResamplingPlan"))
  normalize <- match.arg(normalize)
  input <- .surface_data(x, plan$n_reference)
  op <- .surface_operator(plan, normalize)
  keep <- op$vals > 0
  rows <- op$rows[keep]; cols <- op$cols[keep]; vals <- op$vals[keep]
  if (any(!is.finite(input$matrix[rows, , drop = FALSE])))
    stop("adjoint input must be finite at positively weighted reference vertices")
  out <- matrix(0, plan$n_moving, ncol(input$matrix))
  for (k in seq_len(ncol(out)))
    out[, k] <- .surface_rowsum(vals * input$matrix[rows, k], cols, plan$n_moving)
  values <- if (input$vector) out[, 1L] else out
  if (details) list(values = values, available = .surface_rowsum(vals, cols, plan$n_moving) > 0)
  else values
}

#' Construct geometric resampling in the reverse direction
#'
#' Rebuilds correspondence from the original reference to the original moving
#' mesh, checking both ordered geometry identities. The new source must contain
#' faces for barycentric interpolation. This generally differs from a transpose
#' and does not recover information discarded by downsampling. Masks swap sides.
#'
#' @param plan Original `SurfaceResamplingPlan` with recorded geometry identities.
#' @param reference Original reference mesh.
#' @param moving Original moving mesh.
#' @return A new `SurfaceResamplingPlan` with geometries swapped.
#' @export
reverse_surface_resampling_plan <- function(plan, reference, moving) {
  stopifnot(inherits(plan, "SurfaceResamplingPlan"))
  ref <- if (inherits(reference, "SurfaceMesh")) reference else surface_mesh(reference)
  mov <- if (inherits(moving, "SurfaceMesh")) moving else surface_mesh(moving)
  if (is.null(plan$geometry) ||
      !identical(.surface_identity(ref), plan$geometry$reference) ||
      !identical(.surface_identity(mov), plan$geometry$moving))
    stop("reference and moving must match the plan's recorded ordered geometries")
  surface_resampling_plan(mov, ref, method = plan$method,
    spherical = plan$spherical, radius = plan$radius,
    source_mask = plan$target_mask, target_mask = plan$source_mask, outside = plan$outside)
}
