#' @title Sampler Classes and Resampling
#' @name sampler
#' @description
#' Samplers encapsulate data + interpolation strategy, enabling elegant
#' composition with morphisms for resampling operations.
#'
#' @section Design:
#' A Sampler wraps volumetric or surface data with an interpolation method.
#' The core operation is: given coordinates, return interpolated values.
#'
#' Resampling through a morphism is then simply:
#' \preformatted{
#' source_coords <- transform(morphism, target_coords)
#' values <- sampler(source_coords)
#' }
NULL

# =============================================================================
# SAMPLER CLASS
# =============================================================================

#' Sampler: data + interpolation encapsulation
#'
#' @slot evaluate Function: (coords) -> values
#' @slot domain Character: domain hash
#' @slot dims Integer: data dimensions
#' @slot vdim Integer: output dimensionality per point
#' @slot bounds List: min/max coordinate bounds
#' @slot method Character: interpolation method
#' @slot outside Value for out-of-bounds queries
#' @param object A Sampler object (for show method)
#' @export
setClass("Sampler",
  slots = c(
    evaluate = "function",
    domain = "character",
    dims = "integer",
    vdim = "integer",
    bounds = "list",
    method = "character",
    outside = "ANY"
  ),
  prototype = list(
    evaluate = function(coords) NA,
    domain = "",
    dims = integer(0),
    vdim = 1L,
    bounds = list(min = c(-Inf, -Inf, -Inf), max = c(Inf, Inf, Inf)),
    method = "linear",
    outside = NA_real_
  )
)

#' @rdname Sampler-class
#' @export
setMethod("show", "Sampler", function(object) {
  cat(sprintf("<Sampler | dims=%s | vdim=%d | method=%s>\n",
              paste(object@dims, collapse = "x"),
              object@vdim,
              object@method))
})

# =============================================================================
# GRID CLASS
# =============================================================================

#' Grid: regular coordinate grid specification
#'
#' @slot dims Integer: grid dimensions
#' @slot affine 4x4 voxel-to-world matrix
#' @slot domain Character: domain hash
#' @param object A Grid object (for show method)
#' @export
setClass("Grid",
  slots = c(
    dims = "integer",
    affine = "matrix",
    domain = "character"
  ),
  prototype = list(
    dims = c(1L, 1L, 1L),
    affine = diag(4),
    domain = ""
  )
)

#' @rdname Grid-class
#' @export
setMethod("show", "Grid", function(object) {
  cat(sprintf("<Grid | dims=%s | spacing=%.2f x %.2f x %.2f mm>\n",
              paste(object@dims, collapse = "x"),
              sqrt(sum(object@affine[1:3, 1]^2)),
              sqrt(sum(object@affine[1:3, 2]^2)),
              sqrt(sum(object@affine[1:3, 3]^2))))
})

# =============================================================================
# LIGHTWEIGHT IRREGULAR REFERENCES
# =============================================================================

#' Sampled points (irregular spatial reference)
#'
#' Lightweight container for point samples in world coordinates.
#'
#' @slot coords Numeric matrix (N x 3)
#' @slot domain Character domain hash
#' @export
setClass("SampledPoints",
  slots = c(
    coords = "matrix",
    domain = "character"
  ),
  prototype = list(
    coords = matrix(0, nrow = 0, ncol = 3),
    domain = ""
  )
)

#' Surface mesh (vertices + faces)
#'
#' Lightweight mesh reference for surface operations.
#'
#' @slot coords Numeric matrix (V x 3) of vertices
#' @slot faces Integer matrix (F x 3) of 0-based triangle indices
#' @slot domain Character domain hash
#' @export
setClass("SurfaceMesh",
  slots = c(
    coords = "matrix",
    faces = "matrix",
    domain = "character"
  ),
  prototype = list(
    coords = matrix(0, nrow = 0, ncol = 3),
    faces = matrix(integer(0), nrow = 0, ncol = 3),
    domain = ""
  )
)

#' @rdname SampledPoints-class
#' @export
setMethod("show", "SampledPoints", function(object) {
  cat(sprintf("<SampledPoints | n=%d | domain=%s>\n",
              nrow(object@coords),
              if (nzchar(object@domain)) object@domain else "<none>"))
})

#' @rdname SurfaceMesh-class
#' @export
setMethod("show", "SurfaceMesh", function(object) {
  cat(sprintf("<SurfaceMesh | vertices=%d | faces=%d | domain=%s>\n",
              nrow(object@coords),
              nrow(object@faces),
              if (nzchar(object@domain)) object@domain else "<none>"))
})

#' Construct sampled points
#'
#' @param coords Numeric matrix with 3 columns (N x 3)
#' @param domain Optional domain identifier
#' @return SampledPoints object
#' @export
sampled_points <- function(coords, domain = NULL) {
  if (!is.matrix(coords)) {
    coords <- coerce_surface_reference(coords, require_faces = FALSE)$vertices
  }
  if (!is.matrix(coords) || ncol(coords) != 3 || !is.numeric(coords)) {
    stop("coords must be a numeric matrix with 3 columns")
  }
  if (!all(is.finite(coords))) stop("coords must contain only finite values")
  if (is.null(domain)) domain <- compute_hash("sampled_points", coords)
  new("SampledPoints", coords = coords, domain = domain)
}

#' Construct a surface mesh
#'
#' @param vertices Numeric matrix with 3 columns (V x 3)
#' @param faces Optional matrix with 3 columns (F x 3), 0- or 1-based indices
#' @param domain Optional domain identifier
#' @return SurfaceMesh object
#' @export
surface_mesh <- function(vertices, faces = NULL, domain = NULL) {
  if (!is.matrix(vertices)) {
    ref <- coerce_surface_reference(vertices, require_faces = FALSE)
    if (is.null(faces) && !is.null(ref$faces) && nrow(ref$faces) > 0L) faces <- ref$faces
    if (is.null(domain) && nzchar(ref$domain)) domain <- ref$domain
    vertices <- ref$vertices
  }

  if (!is.matrix(vertices) || ncol(vertices) != 3 || !is.numeric(vertices)) {
    stop("vertices must be a numeric matrix with 3 columns")
  }
  if (!all(is.finite(vertices))) stop("vertices must contain only finite values")
  if (nrow(vertices) < 1L) stop("vertices must be non-empty")

  if (is.null(faces)) {
    faces0 <- matrix(integer(0), nrow = 0, ncol = 3)
  } else {
    if (!is.matrix(faces) || ncol(faces) != 3 || !is.numeric(faces)) {
      stop("faces must be a numeric matrix with 3 columns")
    }
    if (!all(is.finite(faces))) stop("faces must contain only finite values")
    faces0 <- matrix(as.integer(round(faces)), ncol = 3)
    if (any(abs(faces - faces0) > 1e-8)) stop("faces must contain integer indices")
    if (nrow(faces0) > 0) {
      if (min(faces0) >= 1L) faces0 <- faces0 - 1L
      if (min(faces0) < 0L || max(faces0) >= nrow(vertices)) {
        stop("faces indices are out of bounds for vertices")
      }
    }
  }

  if (is.null(domain)) domain <- compute_hash("surface_mesh", vertices, faces0)
  new("SurfaceMesh", coords = vertices, faces = faces0, domain = domain)
}

#' Check whether a surface mesh is approximately spherical
#'
#' @param mesh SurfaceMesh object
#' @param tolerance Ratio tolerance for min/max radius comparison
#' @return Logical
#' @export
mesh_is_sphere <- function(mesh, tolerance = 1.001) {
  if (!inherits(mesh, "SurfaceMesh")) stop("mesh must be a SurfaceMesh")
  dists <- sqrt(rowSums(mesh@coords^2))
  if (!length(dists)) return(FALSE)
  (min(dists) * tolerance) > max(dists)
}

#' Set mesh radius for spherical meshes
#'
#' @param mesh SurfaceMesh object
#' @param radius Target radius
#' @return SurfaceMesh object with rescaled coordinates
#' @export
mesh_set_radius <- function(mesh, radius = 100) {
  if (!inherits(mesh, "SurfaceMesh")) stop("mesh must be a SurfaceMesh")
  if (!mesh_is_sphere(mesh)) {
    stop("mesh_set_radius() requires an approximately spherical mesh")
  }
  dists <- sqrt(rowSums(mesh@coords^2))
  mesh@coords <- mesh@coords * (radius / dists)
  mesh
}

# =============================================================================
# GRID UTILITIES
# =============================================================================

#' Create a grid specification
#'
#' @param dims Integer vector of dimensions
#' @param affine 4x4 voxel-to-world matrix
#' @param domain Optional domain hash
#' @return Grid object
#' @export
grid_spec <- function(dims, affine, domain = NULL) {
  validate_4x4_matrix(affine, "affine")
  if (is.null(domain)) {
    domain <- compute_hash("grid", dims, affine)
  }
  new("Grid",
      dims = as.integer(dims),
      affine = affine,
      domain = domain)
}

#' Get world coordinates for a grid
#'
#' @param grid Grid object
#' @return Numeric matrix (prod(dims) x 3)
#' @export
grid_coords <- function(grid) {
  cpp_volume_world_coords(grid@dims, grid@affine)
}

#' Get grid specification from volume data
#'
#' @param data Volume with geometry
#' @return Grid object
#' @export
grid_from_data <- function(data) {
  affine <- extract_affine(data)
  dims <- dim(data)[1:3]
  grid_spec(dims, affine)
}

# =============================================================================
# AFFINE EXTRACTION (BACKEND-AGNOSTIC)
# =============================================================================

#' Extract affine matrix from data
#' @param x Data object
#' @return 4x4 affine matrix
#' @export
setGeneric("extract_affine", function(x) standardGeneric("extract_affine"))

#' @rdname extract_affine
#' @export
setMethod("extract_affine", "array", function(x) {
  # Default: identity (1mm isotropic)
  diag(4)
})

#' @rdname extract_affine
#' @export
setMethod("extract_affine", "matrix", function(x) {
  if (all(dim(x) == c(4, 4))) {
    return(x)  # Already an affine

  }
  diag(4)
})

# Fallback for RNifti objects without declaring the class
#' @rdname extract_affine
#' @export
setMethod("extract_affine", "ANY", function(x) {
  if (inherits(x, "niftiImage") && requireNamespace("RNifti", quietly = TRUE)) {
    return(RNifti::xform(x))
  }
  diag(4)
})

# neuroim2 support (DenseNeuroVol / DenseNeuroVec are S3 classes; register as old-style S4)
if (!methods::isClass("DenseNeuroVol")) methods::setOldClass("DenseNeuroVol")
if (!methods::isClass("DenseNeuroVec")) methods::setOldClass("DenseNeuroVec")

#' @rdname extract_affine
#' @export
setMethod("extract_affine", "DenseNeuroVol", function(x) {
  neuroim2::trans(x)
})

#' @rdname extract_affine
#' @export
setMethod("extract_affine", "DenseNeuroVec", function(x) {
  neuroim2::trans(x)
})

# =============================================================================
# SAMPLER CONSTRUCTORS
# =============================================================================

#' Create a volume sampler
#'
#' @param data Numeric array (3D/4D) or volume object
#' @param affine 4x4 voxel-to-world (extracted if NULL)
#' @param method Interpolation: "linear", "nearest", "cubic"
#' @param outside Value for out-of-bounds
#' @param domain Optional domain hash
#' @return Sampler object
#' @export
#' @examples
#' vol <- array(rnorm(27), dim = c(3, 3, 3))
#' sampler <- volume_sampler(vol, affine = diag(4))
#' coords <- matrix(c(1, 1, 1), ncol = 3)
#' sampler@evaluate(coords)
volume_sampler <- function(data, affine = NULL,
                           method = c("linear", "nearest", "cubic"),
                           outside = NA_real_, domain = NULL) {
  method <- match.arg(method)

  if (is.null(affine)) {
    affine <- extract_affine(data)
  }
  validate_4x4_matrix(affine, "affine")

  dims <- dim(data)[1:3]
  vdim <- if (length(dim(data)) > 3) dim(data)[4] else 1L

  world_to_vox <- solve(affine)

  # Compute bounds
  corners <- as.matrix(expand.grid(
    i = c(0, dims[1] - 1),
    j = c(0, dims[2] - 1),
    k = c(0, dims[3] - 1)
  ))
  corners_world <- t(apply(corners, 1, function(ijk) {
    (affine %*% c(ijk, 1))[1:3]
  }))
  bounds <- list(
    min = apply(corners_world, 2, min),
    max = apply(corners_world, 2, max)
  )

  if (is.null(domain)) {
    domain <- compute_hash("volume", dims, affine)
  }

  # Evaluate function using C++
  data_array <- as.array(data)
  evaluate_fn <- function(coords) {
    cpp_sample_volume(
      data = data_array,
      coords = coords,
      world_to_vox = world_to_vox,
      method = method,
      outside = outside
    )
  }

  new("Sampler",
      evaluate = evaluate_fn,
      domain = domain,
      dims = as.integer(dims),
      vdim = as.integer(vdim),
      bounds = bounds,
      method = method,
      outside = outside)
}

#' Create a surface sampler
#'
#' @param vertices Vertex coordinates (V x 3), SampledPoints, or SurfaceMesh
#' @param data Vertex data (length V or V x k)
#' @param faces Face indices (F x 3, optional)
#' @param method "nearest" or "barycentric"
#' @param domain Optional domain hash
#' @return Sampler object
#' @export
surface_sampler <- function(vertices, data, faces = NULL,
                            method = c("nearest", "barycentric"),
                            domain = NULL) {
  method <- match.arg(method)

  if (!is.matrix(vertices)) {
    ref <- coerce_surface_reference(vertices, require_faces = FALSE)
    if (is.null(faces) && !is.null(ref$faces) && nrow(ref$faces) > 0L) faces <- ref$faces
    if (is.null(domain) && nzchar(ref$domain)) domain <- ref$domain
    vertices <- ref$vertices
  }

  if (!is.matrix(vertices) || ncol(vertices) != 3 || !is.numeric(vertices)) {
    stop("vertices must be a numeric matrix with 3 columns")
  }
  if (!all(is.finite(vertices))) stop("vertices must contain only finite values")

  if (method == "barycentric" && is.null(faces)) {
    stop("Barycentric interpolation requires faces")
  }

  nv <- nrow(vertices)
  vdim <- if (is.matrix(data)) ncol(data) else 1L
  if ((is.vector(data) && length(data) != nv) || (is.matrix(data) && nrow(data) != nv)) {
    stop("data must have one value (or row) per vertex")
  }

  if (is.null(domain)) {
    domain <- compute_hash("surface", nv)
  }

  # Build evaluate function
  evaluate_fn <- if (method == "nearest") {
    function(coords) {
      indices <- cpp_nearest_vertex(coords, vertices)
      if (vdim == 1L) data[indices] else data[indices, , drop = FALSE]
    }
  } else {
    function(coords) {
      cpp_barycentric_sample(coords, vertices, faces, data)
    }
  }

  new("Sampler",
      evaluate = evaluate_fn,
      domain = domain,
      dims = as.integer(nv),
      vdim = as.integer(vdim),
      bounds = list(
        min = apply(vertices, 2, min),
        max = apply(vertices, 2, max)
      ),
      method = method,
      outside = NA_real_)
}

# =============================================================================
# CORE RESAMPLING
# =============================================================================

#' Resample data through a morphism
#'
#' The fundamental resampling operation: transform coordinates via morphism,
#' then evaluate sampler at transformed coordinates.
#'
#' @param sampler Sampler object
#' @param morphism Morphism or MorphismPath
#' @param coords Target coords (n x 3 matrix) or Grid
#' @param modulate Jacobian modulation: "none", "jacobian", "sqrt_jacobian"
#' @return Sampled values
#' @export
#' @examples
#' # Create volume and sampler
#' vol <- array(1:27, dim = c(3, 3, 3))
#' sampler <- volume_sampler(vol, affine = diag(4))
#'
#' # Create identity morphism
#' morph <- IdentityMorphism("test")
#'
#' # Resample at specific coordinates
#' coords <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE)
#' resample(sampler, morph, coords)
resample <- function(sampler, morphism, coords,
                     modulate = c("none", "jacobian", "sqrt_jacobian")) {
  modulate <- match.arg(modulate)

  # Handle Grid input
  if (is(coords, "Grid")) {
    coords <- grid_coords(coords)
  }

  # Domain compatibility warning
  if (nzchar(sampler@domain) && sampler@domain != source_of(morphism)) {
    warning("Sampler domain doesn't match morphism source")
  }

  # Transform coordinates (pullback)
  source_coords <- transform(morphism, coords)

  # Evaluate sampler
  values <- sampler@evaluate(source_coords)

  # Apply modulation
  if (modulate != "none") {
    jdet <- abs(jacobian_det(morphism, coords, log = FALSE, mode = "pullback"))

    if (modulate == "sqrt_jacobian") {
      jdet <- sqrt(jdet)
    }

    if (is.matrix(values)) {
      values <- values * jdet
    } else {
      values <- values * jdet
    }
  }

  values
}

#' Resample volume to new geometry
#'
#' Convenience wrapper for common volume-to-volume resampling.
#'
#' @param data Source volume
#' @param morphism Transform
#' @param target Target geometry (Grid or volume)
#' @param method Interpolation method
#' @param modulate Jacobian modulation
#' @return Array with target dimensions
#' @export
resample_volume <- function(data, morphism, target,
                            method = "linear", modulate = "none") {
  target_grid <- if (is(target, "Grid")) target else grid_from_data(target)
  source_grid <- grid_from_data(data)

  # Use cached plan path when modulation is not requested
  if (identical(modulate, "none")) {
    plan <- make_resampling_plan(morphism, source_grid, target_grid,
                                 interpolation = method, reuse_count = 2L, cache = TRUE)
    return(apply_resampling_plan(plan, data, modulate = modulate))
  }

  sampler <- volume_sampler(data, method = method)
  values <- resample(sampler, morphism, target_grid, modulate = modulate)
  array(values, dim = target_grid@dims)
}

#' Sample volume on surface
#'
#' Project volumetric data onto surface vertices.
#'
#' @param data Source volume
#' @param morphism VolToSurfMorphism
#' @param surface_coords Surface vertex coordinates (optional)
#' @param method Interpolation method
#' @return Values at surface vertices
#' @export
sample_volume_on_surface <- function(data, morphism, surface_coords = NULL,
                                     method = "linear") {
  if (!is(morphism, "VolToSurfMorphism")) {
    stop("morphism must be a VolToSurfMorphism")
  }

  sampler <- volume_sampler(data, method = method)

  if (is.null(surface_coords)) {
    surface_coords <- morphism@params$mid_coords
    if (is.null(surface_coords)) {
      stop("Surface coordinates required")
    }
  } else if (!is.matrix(surface_coords)) {
    surface_coords <- coerce_surface_reference(surface_coords, require_faces = FALSE)$vertices
  }

  # Handle ribbon sampling specially
  if (morphism@method == "ribbon") {
    inner <- morphism@params$ribbon_inner_coords
    outer <- morphism@params$ribbon_outer_coords
    n_steps <- morphism@n_ribbon_samples

    return(cpp_ribbon_sample_volume(
      data = as.array(data),
      inner = inner,
      outer = outer,
      world_to_vox = solve(extract_affine(data)),
      n_steps = n_steps,
      method = method
    ))
  }

  # Standard sampling
  source_coords <- transform(morphism, surface_coords)
  sampler@evaluate(source_coords)
}

#' Backproject surface values into a volume (adjoint of VolToSurf sampling)
#'
#' Implements a simple backprojection operator that distributes each surface
#' value into the source volume grid using the same interpolation weights that
#' forward sampling would use.
#'
#' This is the transpose/adjoint of the linear sampling operator for
#' \code{sample_volume_on_surface(..., method="linear")} when the same
#' \code{surface_coords} are used.
#'
#' @param surface_values Numeric vector (length N) or matrix (N x k) of values at surface points
#' @param morphism VolToSurfMorphism
#' @param grid Source Grid (volume geometry)
#' @param surface_coords Optional N x 3 world coordinates; defaults to morphism@params$mid_coords
#' @param method "linear", "nearest", or "ribbon"
#' @param outside Value used for out-of-bounds (default 0)
#' @return 3D array (or 4D if multivariate surface_values) in grid geometry
#' @export
backproject_surface_to_volume <- function(surface_values, morphism, grid,
                                         surface_coords = NULL,
                                         method = c("linear", "nearest", "ribbon"),
                                         outside = 0) {
  method <- match.arg(method)
  if (!is(morphism, "VolToSurfMorphism")) stop("morphism must be a VolToSurfMorphism")
  if (!inherits(grid, "Grid")) stop("grid must be a Grid")

  if (is.null(surface_coords)) {
    surface_coords <- morphism@params$mid_coords
  } else if (!is.matrix(surface_coords)) {
    surface_coords <- coerce_surface_reference(surface_coords, require_faces = FALSE)$vertices
  }
  if (is.null(surface_coords)) stop("surface_coords required (or provide morphism@params$mid_coords)")
  if (!is.matrix(surface_coords) || ncol(surface_coords) != 3) stop("surface_coords must be an N x 3 matrix")

  n <- nrow(surface_coords)
  if (is.vector(surface_values)) {
    if (length(surface_values) != n) stop("surface_values length must match nrow(surface_coords)")
    surface_mat <- matrix(surface_values, ncol = 1)
  } else if (is.matrix(surface_values)) {
    if (nrow(surface_values) != n) stop("surface_values nrow must match nrow(surface_coords)")
    surface_mat <- surface_values
  } else {
    stop("surface_values must be a numeric vector or matrix")
  }

  dims <- as.integer(grid@dims)
  nvox <- prod(dims)
  vdim <- ncol(surface_mat)

  out <- matrix(0, nrow = nvox, ncol = vdim)

  # Convert world -> voxel (continuous, 0-based)
  w2v <- solve(grid@affine)
  vox <- apply_affine(surface_coords, w2v)

  if (identical(method, "nearest")) {
    xi <- as.integer(floor(vox[, 1] + 0.5))
    yi <- as.integer(floor(vox[, 2] + 0.5))
    zi <- as.integer(floor(vox[, 3] + 0.5))
    keep <- xi >= 0 & yi >= 0 & zi >= 0 & xi < dims[1] & yi < dims[2] & zi < dims[3]
    if (any(keep)) {
      lin <- xi[keep] + yi[keep] * dims[1] + zi[keep] * dims[1] * dims[2] + 1L
      for (k in seq_len(vdim)) {
        tmp <- rowsum(surface_mat[keep, k, drop = TRUE], group = lin, reorder = FALSE)
        out[as.integer(rownames(tmp)), k] <- out[as.integer(rownames(tmp)), k] + tmp[, 1]
      }
    }
  } else if (identical(method, "linear")) {
    w <- cpp_trilinear_weights(vox, dims)
    if (length(w$rows) > 0) {
      rows <- w$rows
      cols <- w$cols
      vals <- w$vals
      for (k in seq_len(vdim)) {
        contrib <- surface_mat[rows, k, drop = TRUE] * vals
        tmp <- rowsum(contrib, group = cols, reorder = FALSE)
        out[as.integer(rownames(tmp)), k] <- out[as.integer(rownames(tmp)), k] + tmp[, 1]
      }
    }
  } else {
    inner <- morphism@params$ribbon_inner_coords
    outer <- morphism@params$ribbon_outer_coords
    if (is.null(inner) || is.null(outer)) stop("ribbon backprojection requires ribbon_inner_coords and ribbon_outer_coords in morphism@params")
    mask <- rep(TRUE, nvox)
    w <- cpp_ribbon_weights(inner, outer, w2v, dims, n_steps = morphism@n_ribbon_samples, mask = mask)
    if (length(w$rows) > 0) {
      rows <- w$rows
      cols <- w$cols
      vals <- w$vals
      for (k in seq_len(vdim)) {
        contrib <- surface_mat[rows, k, drop = TRUE] * vals
        tmp <- rowsum(contrib, group = cols, reorder = FALSE)
        out[as.integer(rownames(tmp)), k] <- out[as.integer(rownames(tmp)), k] + tmp[, 1]
      }
    }
  }

  # Outside fill (currently only affects OOB points by contributing nothing)
  if (!identical(outside, 0)) {
    out[!is.finite(out)] <- outside
  }

  if (vdim == 1) {
    return(array(out[, 1], dim = dims))
  }
  array(out, dim = c(dims, vdim))
}
