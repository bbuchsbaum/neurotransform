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

# Define niftiImage as an old-style class if RNifti is present, to silence
# method registration warnings when the class has not been loaded yet.
if (!methods::isClass("niftiImage") && requireNamespace("RNifti", quietly = TRUE)) {
  methods::setOldClass("niftiImage")
}

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

# RNifti niftiImage support (conditional registration at load time)
if (requireNamespace("RNifti", quietly = TRUE)) {
  setMethod("extract_affine", "niftiImage", function(x) {
    RNifti::xform(x)
  })
}

# neuroim2 support (DenseNeuroVol / DenseNeuroVec are S3 classes)
# Methods registered conditionally at load time
if (requireNamespace("neuroim2", quietly = TRUE)) {
  methods::setOldClass("DenseNeuroVol")
  methods::setOldClass("DenseNeuroVec")

  setMethod("extract_affine", "DenseNeuroVol", function(x) {
    neuroim2::trans(x)
  })

  setMethod("extract_affine", "DenseNeuroVec", function(x) {
    neuroim2::trans(x)
  })
}

# RNifti support (conditional)
# setMethod("extract_affine", "niftiImage", function(x) {
#   RNifti::xform(x)
# })

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
#' @param vertices Vertex coordinates (V x 3)
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

  if (method == "barycentric" && is.null(faces)) {
    stop("Barycentric interpolation requires faces")
  }

  nv <- nrow(vertices)
  vdim <- if (is.matrix(data)) ncol(data) else 1L

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
