# neurotransform vNext: Jacobians and Resampling

## Overview

This document details the design for two major extensions to neurotransform:
1. **Jacobians** - derivatives of morphisms, enabling local geometry analysis
2. **Resampling** - elegant data interpolation through morphisms

Both maintain the core philosophy: morphisms are about coordinate transformations,
everything else composes on top.

---

## Part 1: Jacobian API

### 1.1 Core Types

```r
#' Jacobian field: collection of 3x3 Jacobian matrices at coordinate locations
#'
#' @slot values Numeric array (n x 3 x 3) of Jacobian matrices
#' @slot coords Numeric matrix (n x 3) of coordinate locations
#' @slot mode Character: "pullback" or "pushforward"
#' @slot morphism_hash Hash of source morphism for provenance
setClass("JacobianField",
  slots = c(
    values = "array",
    coords = "matrix",
    mode = "character",
    morphism_hash = "character"
  )
)
```

### 1.2 Core Generics

```r
#' Compute Jacobian matrices at coordinate locations
#'
#' For a morphism f: A -> B with pullback semantics (transform maps B coords to A coords),
#' the Jacobian is the derivative of this mapping.
#'
#' @param morphism A Morphism or MorphismPath
#' @param coords Numeric matrix (n x 3) of coordinates
#' @param mode "pullback" (d(source)/d(target)) or "pushforward" (d(target)/d(source))
#' @return JacobianField object with n x 3 x 3 array of matrices
#'
#' @details
#' Mode semantics:
#' - "pullback": Returns d(source_coords)/d(target_coords). This is the natural
#'   derivative of transform(), evaluated at target coordinates.
#' - "pushforward": Returns d(target_coords)/d(source_coords). Requires invertible
#'   morphism; computed as inverse of pullback Jacobian.
#'
#' @export
setGeneric("jacobian", function(morphism, coords, mode = c("pullback", "pushforward")) {

  standardGeneric("jacobian")
})

#' Compute Jacobian determinants
#'
#' @param morphism A Morphism or MorphismPath
#' @param coords Numeric matrix (n x 3) of coordinates
#' @param log Logical: return log(|det(J)|) instead of det(J)?
#' @param mode "pullback" or "pushforward"
#' @return Numeric vector of length n
#'
#' @details
#' The determinant measures local volume change:
#' - det(J) > 1: local expansion

#' - det(J) < 1: local contraction
#' - det(J) < 0: local folding (non-diffeomorphic)
#'
#' Log determinants are preferred for:
#' - Numerical stability with extreme deformations
#' - Statistical analysis (log-Jacobian is approximately normal)
#' - Additive composition: log|det(AB)| = log|det(A)| + log|det(B)|
#'
#' @export
setGeneric("jacobian_det", function(morphism, coords, log = FALSE,
                                     mode = c("pullback", "pushforward")) {
  standardGeneric("jacobian_det")
})
```

### 1.3 JacobianField Methods

```r
# Accessors
setMethod("[", "JacobianField", function(x, i) {

  # Return single 3x3 matrix or subset JacobianField
  if (length(i) == 1) {
    return(x@values[i, , , drop = TRUE])
  }
  new("JacobianField",
      values = x@values[i, , , drop = FALSE],
      coords = x@coords[i, , drop = FALSE],
      mode = x@mode,
      morphism_hash = x@morphism_hash)
})

# Determinants
setMethod("det", "JacobianField", function(x, log = FALSE) {
  # Vectorized determinant computation
  dets <- apply(x@values, 1, function(J) det(J))
  if (log) log(abs(dets)) else dets
})

# Inversion (for pushforward from pullback)
setMethod("solve", "JacobianField", function(a) {
  inv_values <- array(NA_real_, dim = dim(a@values))
  for (i in seq_len(nrow(a@coords))) {
    inv_values[i, , ] <- solve(a@values[i, , ])
  }
  new("JacobianField",
      values = inv_values,
      coords = a@coords,
      mode = if (a@mode == "pullback") "pushforward" else "pullback",
      morphism_hash = a@morphism_hash)
})

# Matrix multiplication for chain rule
setMethod("%*%", signature("JacobianField", "JacobianField"), function(x, y) {
  # J_composed[i] = x[i] %*% y[i]
  stopifnot(nrow(x@coords) == nrow(y@coords))
  result <- array(NA_real_, dim = dim(x@values))
  for (i in seq_len(nrow(x@coords))) {
    result[i, , ] <- x@values[i, , ] %*% y@values[i, , ]
  }
  new("JacobianField",
      values = result,
      coords = y@coords,  # Coords from the "inner" function
      mode = x@mode,
      morphism_hash = compute_hash(x@morphism_hash, y@morphism_hash))
})
```

### 1.4 Morphism-Specific Implementations

#### 1.4.1 IdentityMorphism

```r
setMethod("jacobian", "IdentityMorphism", function(morphism, coords,
                                                    mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)
  n <- nrow(coords)
  # Identity Jacobian is I_3 everywhere
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

setMethod("jacobian_det", "IdentityMorphism", function(morphism, coords, log = FALSE,
                                                        mode = c("pullback", "pushforward")) {
  n <- nrow(coords)
  if (log) rep(0, n) else rep(1, n)
})
```

#### 1.4.2 Affine3DMorphism

```r
setMethod("jacobian", "Affine3DMorphism", function(morphism, coords,
                                                    mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)
  n <- nrow(coords)

  # Extract linear part (upper-left 3x3)
  linear <- morphism@matrix[1:3, 1:3]

  # For pullback, the Jacobian IS the linear part (constant everywhere)
  # For pushforward, it's the inverse
  J <- if (mode == "pullback") linear else solve(linear)

  # Replicate for all coordinates
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

setMethod("jacobian_det", "Affine3DMorphism", function(morphism, coords, log = FALSE,
                                                        mode = c("pullback", "pushforward")) {
  linear <- morphism@matrix[1:3, 1:3]
  d <- det(linear)
  if (mode == "pushforward") d <- 1 / d

  n <- nrow(coords)
  if (log) rep(log(abs(d)), n) else rep(d, n)
})
```

#### 1.4.3 Warp3DMorphism

This is where the real work happens. We need C++ for performance.

```r
setMethod("jacobian", "Warp3DMorphism", function(morphism, coords,
                                                  mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)

  # Load warp field
  cache_env <- morphism@cache %||% new_cache_env()
  warp <- load_warp_array(morphism, cache_env = cache_env)

  # Compute Jacobians via C++
  # This interpolates the gradient of the displacement field
  values <- cpp_warp_jacobian(
    coords = coords,
    field = warp$array,
    dim = warp$dim,
    world_to_vox = warp$world_to_vox,
    vox_to_world = warp$vox_to_world
  )

  jf <- new("JacobianField",
            values = values,
            coords = coords,
            mode = "pullback",
            morphism_hash = morphism@hash)

  # If pushforward requested, invert
  if (mode == "pushforward") {
    if (morphism@inverse_type %in% c("none", "adjoint")) {
      stop("Pushforward Jacobian requires invertible morphism")
    }
    jf <- solve(jf)
  }

  jf
})

setMethod("jacobian_det", "Warp3DMorphism", function(morphism, coords, log = FALSE,
                                                      mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)

  cache_env <- morphism@cache %||% new_cache_env()
  warp <- load_warp_array(morphism, cache_env = cache_env)

  # Direct determinant computation is more efficient than full Jacobian
  dets <- cpp_warp_jacobian_det(
    coords = coords,
    field = warp$array,
    dim = warp$dim,
    world_to_vox = warp$world_to_vox,
    vox_to_world = warp$vox_to_world
  )

  if (mode == "pushforward") dets <- 1 / dets
  if (log) log(abs(dets)) else dets
})
```

#### 1.4.4 MorphismPath (Chain Rule)

```r
setMethod("jacobian", "MorphismPath", function(morphism, coords,
                                                mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)
  path <- morphism@morphisms

  if (length(path) == 0) {
    return(jacobian(IdentityMorphism(morphism@source), coords, mode))
  }

  # For pullback: apply morphisms last-to-first, accumulate Jacobians
  # Chain rule: J(g ∘ f)(x) = J(f)(g(x)) %*% J(g)(x)

  # Start with coordinates in target space
  current_coords <- coords

  # Apply path in reverse (pullback order), collecting intermediate coords
  intermediate <- vector("list", length(path))
  intermediate[[length(path)]] <- current_coords

  for (i in rev(seq_along(path))) {
    if (i > 1) {
      current_coords <- transform(path[[i]], current_coords)
      intermediate[[i - 1]] <- current_coords
    }
  }

  # Now compute Jacobians and multiply
  # For path [f, g] (f: A->B, g: B->C), pullback applies g then f
  # J_pullback = J_f(coords_in_B) %*% J_g(coords_in_C)

  result <- NULL
  for (i in rev(seq_along(path))) {
    coords_for_J <- intermediate[[i]]
    J_i <- jacobian(path[[i]], coords_for_J, mode = "pullback")

    if (is.null(result)) {
      result <- J_i
    } else {
      result <- jacobian_multiply(J_i, result)  # J_i on left
    }
  }

  result@mode <- mode
  if (mode == "pushforward") {
    result <- solve(result)
  }

  result
})
```

#### 1.4.5 Surface Morphisms (Not Implemented)

```r
setMethod("jacobian", "VolToSurfMorphism", function(morphism, coords, mode) {
  stop("Jacobian not defined for VolToSurfMorphism.\n",
       "Volume-to-surface projections don't preserve dimensionality.\n",
       "For tangential Jacobians on the surface, see `surface_jacobian()`.")
})

setMethod("jacobian", "SurfToSurfMorphism", function(morphism, coords, mode) {
  stop("Jacobian not implemented for SurfToSurfMorphism.\n",
       "Surface-to-surface mappings require tangent space computations.\n",
       "This may be supported in a future version.")
})
```

### 1.5 C++ Implementation for Warp Jacobians

```cpp
// [[Rcpp::export]]
arma::cube cpp_warp_jacobian(
    const arma::mat& coords,           // n x 3
    const arma::vec& field,            // flattened displacement field
    const arma::ivec& dim,             // field dimensions
    const arma::mat& world_to_vox,     // 4x4
    const arma::mat& vox_to_world      // 4x4
) {
    int n = coords.n_rows;
    int nx = dim(0), ny = dim(1), nz = dim(2);

    arma::cube result(n, 3, 3);

    // Reshape field to 4D: 3 x nx x ny x nz
    // field[c + 3*(i + nx*(j + ny*k))] = displacement component c at voxel (i,j,k)

    for (int p = 0; p < n; p++) {
        // Convert world coord to voxel coord
        arma::vec world_h = {coords(p, 0), coords(p, 1), coords(p, 2), 1.0};
        arma::vec vox_h = world_to_vox * world_h;
        double vx = vox_h(0), vy = vox_h(1), vz = vox_h(2);

        // Jacobian of the deformation: J = I + grad(displacement)
        // We compute grad(displacement) via trilinear interpolation of gradients

        arma::mat grad_disp(3, 3, arma::fill::zeros);  // d(disp_i)/d(x_j)

        // For each displacement component
        for (int c = 0; c < 3; c++) {
            // Compute gradient via central differences in the field
            // Then interpolate to the query point
            arma::vec grad_c = interpolate_gradient(field, dim, c, vx, vy, vz);

            // grad_c is in voxel space; convert to world space
            // d(disp)/d(world) = d(disp)/d(vox) * d(vox)/d(world)
            arma::mat vox_to_world_linear = vox_to_world.submat(0, 0, 2, 2);
            arma::mat world_to_vox_linear = world_to_vox.submat(0, 0, 2, 2);

            grad_disp.row(c) = (world_to_vox_linear * grad_c).t();
        }

        // Jacobian of pullback = I + grad(displacement)
        // (since transform(x) = x + displacement(x))
        result.slice(p) = arma::eye(3, 3) + grad_disp;
    }

    return result;
}

// Helper: interpolate gradient of displacement component at fractional voxel coords
arma::vec interpolate_gradient(
    const arma::vec& field,
    const arma::ivec& dim,
    int component,
    double vx, double vy, double vz
) {
    int nx = dim(0), ny = dim(1), nz = dim(2);

    // Clamp to valid range
    vx = std::max(1.0, std::min(vx, (double)(nx - 2)));
    vy = std::max(1.0, std::min(vy, (double)(ny - 2)));
    vz = std::max(1.0, std::min(vz, (double)(nz - 2)));

    int i0 = (int)floor(vx), j0 = (int)floor(vy), k0 = (int)floor(vz);
    double fx = vx - i0, fy = vy - j0, fz = vz - k0;

    // Compute gradients at the 8 neighboring voxels via central differences
    // Then trilinearly interpolate

    auto idx = [&](int c, int i, int j, int k) {
        return c + 3 * (i + nx * (j + ny * k));
    };

    auto central_diff_x = [&](int c, int i, int j, int k) {
        if (i <= 0 || i >= nx - 1) return 0.0;
        return 0.5 * (field(idx(c, i+1, j, k)) - field(idx(c, i-1, j, k)));
    };

    auto central_diff_y = [&](int c, int i, int j, int k) {
        if (j <= 0 || j >= ny - 1) return 0.0;
        return 0.5 * (field(idx(c, i, j+1, k)) - field(idx(c, i, j-1, k)));
    };

    auto central_diff_z = [&](int c, int i, int j, int k) {
        if (k <= 0 || k >= nz - 1) return 0.0;
        return 0.5 * (field(idx(c, i, j, k+1)) - field(idx(c, i, j, k-1)));
    };

    // Trilinear interpolation of each gradient component
    arma::vec grad(3);

    for (int g = 0; g < 3; g++) {
        auto diff_func = (g == 0) ? central_diff_x : (g == 1) ? central_diff_y : central_diff_z;

        double v000 = diff_func(component, i0,   j0,   k0);
        double v100 = diff_func(component, i0+1, j0,   k0);
        double v010 = diff_func(component, i0,   j0+1, k0);
        double v110 = diff_func(component, i0+1, j0+1, k0);
        double v001 = diff_func(component, i0,   j0,   k0+1);
        double v101 = diff_func(component, i0+1, j0,   k0+1);
        double v011 = diff_func(component, i0,   j0+1, k0+1);
        double v111 = diff_func(component, i0+1, j0+1, k0+1);

        grad(g) = (1-fx)*(1-fy)*(1-fz)*v000 + fx*(1-fy)*(1-fz)*v100 +
                  (1-fx)*fy*(1-fz)*v010     + fx*fy*(1-fz)*v110 +
                  (1-fx)*(1-fy)*fz*v001     + fx*(1-fy)*fz*v101 +
                  (1-fx)*fy*fz*v011         + fx*fy*fz*v111;
    }

    return grad;
}

// [[Rcpp::export]]
arma::vec cpp_warp_jacobian_det(
    const arma::mat& coords,
    const arma::vec& field,
    const arma::ivec& dim,
    const arma::mat& world_to_vox,
    const arma::mat& vox_to_world
) {
    // More efficient: compute det directly without forming full matrix
    arma::cube J = cpp_warp_jacobian(coords, field, dim, world_to_vox, vox_to_world);

    int n = coords.n_rows;
    arma::vec dets(n);

    for (int i = 0; i < n; i++) {
        dets(i) = arma::det(J.slice(i));
    }

    return dets;
}
```

---

## Part 2: Resampling API

### 2.1 Design Philosophy

The resampling API follows these principles:

1. **Separation of concerns**: Samplers know about data; morphisms know about coordinates
2. **Composability**: Complex operations built from simple primitives
3. **Backend agnosticism**: Works with neuroim2, RNifti, oro.nifti, arrays, etc.
4. **Lazy evaluation where beneficial**: Build up operations, execute efficiently

### 2.2 Core Types

```r
#' Sampler: encapsulates data + interpolation strategy
#'
#' A Sampler is fundamentally a function: coords → values
#' But wrapped in a class to carry metadata and enable method dispatch.
#'
#' @slot evaluate Function: (coords matrix) → values
#' @slot domain Character: domain hash this sampler is defined on
#' @slot dims Integer: dimensions of the data (for volume: c(nx, ny, nz))
#' @slot vdim Integer: dimensionality of output per point (1 for scalar, 3 for vector, etc.)
#' @slot bounds List with min/max coordinates (for bounds checking)
#' @slot method Character: interpolation method used
#' @slot outside Value to return for out-of-bounds queries
setClass("Sampler",
  slots = c(
    evaluate = "function",
    domain = "character",
    dims = "integer",
    vdim = "integer",
    bounds = "list",
    method = "character",
    outside = "ANY"
  )
)

#' Grid: regular coordinate grid specification
#'
#' Lightweight specification of a regular grid, without materializing coordinates.
#'
#' @slot dims Integer vector of grid dimensions
#' @slot affine 4x4 voxel-to-world matrix
#' @slot domain Character: domain hash
setClass("Grid",
  slots = c(
    dims = "integer",
    affine = "matrix",
    domain = "character"
  )
)
```

### 2.3 Sampler Constructors

```r
#' Create a volume sampler
#'
#' Wraps volumetric data with an interpolation strategy.
#'
#' @param data Numeric array (3D or 4D) or object with indexing
#' @param affine 4x4 voxel-to-world matrix (extracted automatically for neuroim2/RNifti)
#' @param method Interpolation: "nearest", "linear" (default), "cubic"
#' @param outside Value for out-of-bounds queries (default NA)
#' @param domain Optional domain hash (generated if not provided)
#' @return Sampler object
#'
#' @export
volume_sampler <- function(data, affine = NULL, method = c("linear", "nearest", "cubic"),
                           outside = NA, domain = NULL) {
  method <- match.arg(method)

  # Extract affine from known types
  if (is.null(affine)) {
    affine <- extract_affine(data)  # dispatches on class
  }
  validate_4x4_matrix(affine, "affine")

  # Get dimensions
  dims <- dim(data)[1:3]
  vdim <- if (length(dim(data)) > 3) dim(data)[4] else 1L

  # Compute world-to-voxel for the evaluate function

  world_to_vox <- solve(affine)

  # Build bounds in world coordinates
  corners <- expand.grid(i = c(0, dims[1]), j = c(0, dims[2]), k = c(0, dims[3]))
  corners_world <- t(apply(corners, 1, function(ijk) {
    (affine %*% c(ijk, 1))[1:3]
  }))
  bounds <- list(
    min = apply(corners_world, 2, min),
    max = apply(corners_world, 2, max)
  )

  # Generate domain hash if not provided
  if (is.null(domain)) {
    domain <- compute_hash("volume", dims, affine)
  }

  # Build the evaluate function
  # This is where we choose C++ vs R implementation
  evaluate_fn <- function(coords) {
    cpp_sample_volume(
      data = as.array(data),
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
#' Wraps surface mesh data with interpolation.
#'
#' @param vertices Numeric matrix (V x 3) of vertex coordinates
#' @param data Numeric vector (length V) or matrix (V x k) of vertex data
#' @param faces Integer matrix (F x 3) of face indices (optional, for barycentric)
#' @param method "nearest" (default) or "barycentric" (requires faces)
#' @param domain Optional domain hash
#' @return Sampler object
#'
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
    domain <- compute_hash("surface", nv, vertices[1:min(10, nv), ])
  }

  # Build evaluate function
  evaluate_fn <- if (method == "nearest") {
    function(coords) {
      # Find nearest vertex for each query point
      indices <- cpp_nearest_vertex(coords, vertices)
      if (vdim == 1L) data[indices] else data[indices, , drop = FALSE]
    }
  } else {
    function(coords) {
      # Project to surface and interpolate
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
      outside = NA)
}

#' Generic affine extraction
#' @keywords internal
setGeneric("extract_affine", function(x) standardGeneric("extract_affine"))

setMethod("extract_affine", "array", function(x) {
  # Default: assume RAS with 1mm spacing
  dims <- dim(x)[1:3]
  diag(c(1, 1, 1, 1))
})

# Method for RNifti
setMethod("extract_affine", "niftiImage", function(x) {
  RNifti::xform(x)
})

# Method for neuroim2 (if available)
# setMethod("extract_affine", "NeuroVol", function(x) {
#   neuroim2::trans(space(x))
# })
```

### 2.4 Grid Utilities

```r
#' Create a grid specification
#'
#' @param dims Integer vector of dimensions
#' @param affine 4x4 voxel-to-world matrix
#' @param domain Optional domain hash
#' @return Grid object
#'
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
#' @return Numeric matrix (prod(dims) x 3) of world coordinates
#'
#' @export
grid_coords <- function(grid) {
  cpp_volume_world_coords(grid@dims, grid@affine)
}

#' Get grid specification from data
#'
#' @param data Volume data with geometry
#' @return Grid object
#'
#' @export
grid_from_data <- function(data) {
  affine <- extract_affine(data)
  dims <- dim(data)[1:3]
  grid_spec(dims, affine)
}
```

### 2.5 Core Resampling Functions

```r
#' Resample data through a morphism
#'
#' The fundamental resampling operation: given a sampler defined on domain A,
#' a morphism from A to B, and coordinates in B, return values at those coordinates.
#'
#' @param sampler A Sampler object
#' @param morphism A Morphism or MorphismPath
#' @param coords Target coordinates (matrix n x 3) or Grid object
#' @param modulate Modulation by Jacobian determinant
#' @return Sampled values
#'
#' @details
#' The operation is:
#' 1. Transform target coords to source coords via morphism pullback
#' 2. Evaluate sampler at source coords
#' 3. Optionally modulate by Jacobian determinant
#'
#' Modulation options:
#' - "none": no modulation (default for most uses)
#' - "jacobian": multiply by |det(J)|  (preserves integrals)
#' - "sqrt_jacobian": multiply by sqrt(|det(J)|) (for correlation maps)
#'
#' @export
resample <- function(sampler, morphism, coords, modulate = c("none", "jacobian", "sqrt_jacobian")) {
  modulate <- match.arg(modulate)


  # Handle Grid input
  if (is(coords, "Grid")) {
    coords <- grid_coords(coords)
  }

  # Domain compatibility check (warning, not error)
  if (nzchar(sampler@domain) && sampler@domain != source_of(morphism)) {
    warning("Sampler domain '", sampler@domain, "' doesn't match morphism source '",
            source_of(morphism), "'")
  }

  # Transform coordinates (pullback: target -> source)
  source_coords <- transform(morphism, coords)

  # Evaluate sampler
  values <- sampler@evaluate(source_coords)

  # Apply modulation if requested
  if (modulate != "none") {
    jdet <- jacobian_det(morphism, coords, log = FALSE, mode = "pullback")
    jdet <- abs(jdet)

    if (modulate == "sqrt_jacobian") {
      jdet <- sqrt(jdet)
    }

    # Handle vector/matrix outputs
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
#' Convenience function for the common case of warping a volume.
#'
#' @param data Source volume data
#' @param morphism Transform from source to target space
#' @param target Target geometry (Grid, or object with extractable geometry)
#' @param method Interpolation method
#' @param modulate Jacobian modulation
#' @return Array with target geometry, filled with resampled values
#'
#' @export
resample_volume <- function(data, morphism, target, method = "linear",
                            modulate = "none") {
  # Create sampler from source data
  sampler <- volume_sampler(data, method = method)

  # Get target grid
  target_grid <- if (is(target, "Grid")) target else grid_from_data(target)

  # Resample
  values <- resample(sampler, morphism, target_grid, modulate = modulate)

  # Reshape to target dimensions
  result <- array(values, dim = target_grid@dims)

  # Attach geometry if we know how
  # (would depend on target type)

  result
}

#' Sample volume data on surface
#'
#' Project volumetric data onto a surface using a VolToSurf morphism.
#'
#' @param data Source volume data
#' @param morphism VolToSurfMorphism specifying the projection
#' @param surface_coords Surface vertex coordinates (n x 3)
#' @param method Volume interpolation method
#' @return Vector or matrix of values at surface vertices
#'
#' @details
#' This is a convenience wrapper that:
#' 1. Creates a volume sampler
#' 2. Uses the VolToSurf morphism to determine sample locations
#' 3. Handles ribbon sampling, mid-thickness, etc. based on morphism params
#'
#' @export
sample_volume_on_surface <- function(data, morphism, surface_coords = NULL,
                                      method = "linear") {
  if (!is(morphism, "VolToSurfMorphism")) {
    stop("morphism must be a VolToSurfMorphism")
  }

  sampler <- volume_sampler(data, method = method)

  # Get surface coordinates if not provided
  if (is.null(surface_coords)) {
    # Use coordinates from morphism params
    surface_coords <- morphism@params$mid_coords
    if (is.null(surface_coords)) {
      stop("Surface coordinates must be provided or embedded in morphism")
    }
  }

  # For ribbon sampling, we need special handling
  if (morphism@method == "ribbon") {
    # Sample along ribbon and average
    inner <- morphism@params$ribbon_inner_coords
    outer <- morphism@params$ribbon_outer_coords
    n_steps <- morphism@n_ribbon_samples

    # Build sampling coords along ribbon
    values <- cpp_ribbon_sample_volume(
      data = as.array(data),
      inner = inner,
      outer = outer,
      world_to_vox = solve(extract_affine(data)),
      n_steps = n_steps,
      method = method
    )

    return(values)
  }

  # Standard case: just sample at the coords
  source_coords <- transform(morphism, surface_coords)
  sampler@evaluate(source_coords)
}
```

### 2.6 C++ Sampling Implementation

```cpp
// [[Rcpp::export]]
Rcpp::NumericVector cpp_sample_volume(
    const Rcpp::NumericVector& data,      // with dim attribute
    const Rcpp::NumericMatrix& coords,     // n x 3 world coords
    const Rcpp::NumericMatrix& world_to_vox,  // 4x4
    const std::string& method,
    double outside
) {
    Rcpp::IntegerVector dims = data.attr("dim");
    int nx = dims[0], ny = dims[1], nz = dims[2];
    int n = coords.nrow();
    bool is_4d = dims.size() > 3;
    int nv = is_4d ? dims[3] : 1;

    Rcpp::NumericVector result(n * nv);
    if (is_4d) {
        result.attr("dim") = Rcpp::IntegerVector::create(n, nv);
    }

    for (int p = 0; p < n; p++) {
        // World to voxel
        double wx = coords(p, 0), wy = coords(p, 1), wz = coords(p, 2);
        double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy +
                    world_to_vox(0,2)*wz + world_to_vox(0,3);
        double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy +
                    world_to_vox(1,2)*wz + world_to_vox(1,3);
        double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy +
                    world_to_vox(2,2)*wz + world_to_vox(2,3);

        // Check bounds
        if (vx < 0 || vx > nx-1 || vy < 0 || vy > ny-1 || vz < 0 || vz > nz-1) {
            for (int v = 0; v < nv; v++) {
                result[p + v*n] = outside;
            }
            continue;
        }

        if (method == "nearest") {
            int i = (int)round(vx), j = (int)round(vy), k = (int)round(vz);
            i = std::max(0, std::min(i, nx-1));
            j = std::max(0, std::min(j, ny-1));
            k = std::max(0, std::min(k, nz-1));

            for (int v = 0; v < nv; v++) {
                result[p + v*n] = data[i + nx*(j + ny*(k + nz*v))];
            }
        } else if (method == "linear") {
            // Trilinear interpolation
            int i0 = (int)floor(vx), j0 = (int)floor(vy), k0 = (int)floor(vz);
            double fx = vx - i0, fy = vy - j0, fz = vz - k0;

            i0 = std::max(0, std::min(i0, nx-2));
            j0 = std::max(0, std::min(j0, ny-2));
            k0 = std::max(0, std::min(k0, nz-2));

            for (int v = 0; v < nv; v++) {
                auto idx = [&](int i, int j, int k) {
                    return i + nx*(j + ny*(k + nz*v));
                };

                double c000 = data[idx(i0, j0, k0)];
                double c100 = data[idx(i0+1, j0, k0)];
                double c010 = data[idx(i0, j0+1, k0)];
                double c110 = data[idx(i0+1, j0+1, k0)];
                double c001 = data[idx(i0, j0, k0+1)];
                double c101 = data[idx(i0+1, j0, k0+1)];
                double c011 = data[idx(i0, j0+1, k0+1)];
                double c111 = data[idx(i0+1, j0+1, k0+1)];

                result[p + v*n] =
                    (1-fx)*(1-fy)*(1-fz)*c000 + fx*(1-fy)*(1-fz)*c100 +
                    (1-fx)*fy*(1-fz)*c010     + fx*fy*(1-fz)*c110 +
                    (1-fx)*(1-fy)*fz*c001     + fx*(1-fy)*fz*c101 +
                    (1-fx)*fy*fz*c011         + fx*fy*fz*c111;
            }
        }
        // cubic would go here...
    }

    return result;
}
```

---

## Part 3: Vector and Tensor Transforms

### 3.1 API Design

```r
#' Transform vector field through morphism
#'
#' Transforms vectors (e.g., gradients, displacements, fiber directions) accounting
#' for the local geometry change.
#'
#' @param morphism A Morphism or MorphismPath
#' @param coords Coordinates where vectors are defined (n x 3)
#' @param vectors Vector values (n x 3)
#' @param type "contravariant" (displacements) or "covariant" (gradients)
#' @param mode "pullback" or "pushforward"
#' @return Transformed vectors (n x 3)
#'
#' @details
#' Vector transformation depends on vector type:
#'
#' Contravariant vectors (tangent vectors, displacements):
#' - Transform with the Jacobian directly
#' - v_target = J * v_source  (pushforward)
#' - v_source = J * v_target  (pullback, since our J is d(source)/d(target))
#'
#' Covariant vectors (gradients, normals):
#' - Transform with the inverse transpose of Jacobian
#' - g_target = J^{-T} * g_source  (pushforward)
#'
#' @export
transform_vectors <- function(morphism, coords, vectors,
                              type = c("contravariant", "covariant"),
                              mode = c("pullback", "pushforward")) {
  type <- match.arg(type)
  mode <- match.arg(mode)

  # Get Jacobians at the coordinates
  J <- jacobian(morphism, coords, mode = "pullback")

  n <- nrow(coords)
  result <- matrix(NA_real_, n, 3)

  for (i in seq_len(n)) {
    Ji <- J[i]  # 3x3 matrix

    if (type == "contravariant") {
      # Contravariant: transform directly with J
      if (mode == "pullback") {
        result[i, ] <- Ji %*% vectors[i, ]
      } else {
        result[i, ] <- solve(Ji) %*% vectors[i, ]
      }
    } else {
      # Covariant: transform with J^{-T}
      if (mode == "pullback") {
        result[i, ] <- solve(t(Ji)) %*% vectors[i, ]
      } else {
        result[i, ] <- t(Ji) %*% vectors[i, ]
      }
    }
  }

  result
}

#' Transform tensor field through morphism
#'
#' Transforms symmetric 3x3 tensors (e.g., diffusion tensors) through a morphism.
#'
#' @param morphism A Morphism or MorphismPath
#' @param coords Coordinates where tensors are defined (n x 3)
#' @param tensors Tensor values (n x 6 for upper triangular, or n x 3 x 3)
#' @param mode "pullback" or "pushforward"
#' @return Transformed tensors in same format as input
#'
#' @details
#' Tensors transform as: T' = J * T * J^T (for contravariant-contravariant tensors)
#'
#' For diffusion tensors in neuroimaging, we often want to preserve the
#' tensor's principal directions while accounting for local rotation, which
#' uses a slightly different formula involving the rotation part of the Jacobian.
#'
#' @export
transform_tensors <- function(morphism, coords, tensors,
                              mode = c("pullback", "pushforward"),
                              preserve_shape = FALSE) {
  mode <- match.arg(mode)

  # Get Jacobians
  J <- jacobian(morphism, coords, mode = "pullback")

  n <- nrow(coords)
  is_compact <- is.matrix(tensors) && ncol(tensors) == 6

  result <- if (is_compact) {
    matrix(NA_real_, n, 6)
  } else {
    array(NA_real_, dim = c(n, 3, 3))
  }

  for (i in seq_len(n)) {
    Ji <- J[i]

    # Extract full tensor
    Ti <- if (is_compact) {
      matrix(c(tensors[i, 1], tensors[i, 4], tensors[i, 5],
               tensors[i, 4], tensors[i, 2], tensors[i, 6],
               tensors[i, 5], tensors[i, 6], tensors[i, 3]), 3, 3)
    } else {
      tensors[i, , ]
    }

    # Transform
    if (preserve_shape) {
      # Use only rotation part (Finite Strain reorientation)
      R <- polar_rotation(Ji)
      Ti_new <- R %*% Ti %*% t(R)
    } else {
      # Full tensor transformation
      if (mode == "pullback") {
        Ti_new <- Ji %*% Ti %*% t(Ji)
      } else {
        Ji_inv <- solve(Ji)
        Ti_new <- Ji_inv %*% Ti %*% t(Ji_inv)
      }
    }

    # Store result
    if (is_compact) {
      result[i, ] <- c(Ti_new[1,1], Ti_new[2,2], Ti_new[3,3],
                       Ti_new[1,2], Ti_new[1,3], Ti_new[2,3])
    } else {
      result[i, , ] <- Ti_new
    }
  }

  result
}

#' Extract rotation component from matrix via polar decomposition
#' @keywords internal
polar_rotation <- function(A) {
  svd_A <- svd(A)
  svd_A$u %*% t(svd_A$v)
}
```

---

## Part 4: QC and Analysis Helpers

### 4.1 Jacobian Maps

```r
#' Compute Jacobian determinant map over a grid
#'
#' @param morphism A Morphism (typically Warp3DMorphism)
#' @param grid Grid specification or volume with geometry
#' @param log Return log(|det(J)|)?
#' @param mode "pullback" or "pushforward"
#' @return Array of Jacobian determinants with grid dimensions
#'
#' @export
jacobian_map <- function(morphism, grid, log = TRUE,
                         mode = c("pullback", "pushforward")) {
  mode <- match.arg(mode)

  if (!is(grid, "Grid")) {
    grid <- grid_from_data(grid)
  }

  coords <- grid_coords(grid)
  dets <- jacobian_det(morphism, coords, log = log, mode = mode)

  array(dets, dim = grid@dims)
}

#' Check for diffeomorphic violations (folding)
#'
#' Identifies regions where the Jacobian determinant is non-positive,
#' indicating local folding (the map is not invertible there).
#'
#' @param morphism A Morphism
#' @param grid Grid specification
#' @param threshold Minimum acceptable det(J) (default 0)
#' @return List with:
#'   - is_diffeomorphic: logical, TRUE if no violations
#'   - n_violations: count of voxels with det(J) <= threshold
#'   - fraction_violations: fraction of voxels with violations
#'   - min_det: minimum determinant value
#'   - violation_mask: logical array marking violation locations
#'
#' @export
check_diffeomorphic <- function(morphism, grid, threshold = 0) {
  if (!is(grid, "Grid")) {
    grid <- grid_from_data(grid)
  }

  coords <- grid_coords(grid)
  dets <- jacobian_det(morphism, coords, log = FALSE, mode = "pullback")

  violations <- dets <= threshold
  n_total <- length(dets)
  n_violations <- sum(violations)

  list(
    is_diffeomorphic = n_violations == 0,
    n_violations = n_violations,
    fraction_violations = n_violations / n_total,
    min_det = min(dets),
    max_det = max(dets),
    median_det = median(dets),
    violation_mask = array(violations, dim = grid@dims)
  )
}

#' Summarize morphism geometry
#'
#' Compute summary statistics about a morphism's geometric properties.
#'
#' @param morphism A Morphism
#' @param grid Optional grid for spatial sampling
#' @return List of summary statistics
#'
#' @export
summarize_morphism <- function(morphism, grid = NULL) {
  summary <- list(
    kind = morphism_kind(morphism),
    source = source_of(morphism),
    target = target_of(morphism),
    cost = morphism@cost,
    inverse_type = morphism@inverse_type,
    inverse_quality = morphism@inverse_quality
  )

  # Add Jacobian stats for spatial morphisms
  if (!is.null(grid) && morphism_kind(morphism) %in% c("affine3d", "warp3d")) {
    if (!is(grid, "Grid")) {
      grid <- grid_from_data(grid)
    }

    coords <- grid_coords(grid)
    dets <- jacobian_det(morphism, coords, log = FALSE)
    log_dets <- log(abs(dets))

    summary$jacobian_stats <- list(
      det_min = min(dets),
      det_max = max(dets),
      det_median = median(dets),
      log_det_mean = mean(log_dets),
      log_det_sd = sd(log_dets),
      n_folding = sum(dets <= 0),
      is_diffeomorphic = all(dets > 0)
    )
  }

  class(summary) <- "morphism_summary"
  summary
}

#' @export
print.morphism_summary <- function(x, ...) {
  cat("Morphism Summary\n")
  cat("================\n")
  cat(sprintf("  Kind:           %s\n", x$kind))
  cat(sprintf("  Source:         %s\n", x$source))
  cat(sprintf("  Target:         %s\n", x$target))
  cat(sprintf("  Cost:           %.3f\n", x$cost))
  cat(sprintf("  Inverse type:   %s\n", x$inverse_type))
  cat(sprintf("  Inverse qual:   %.3f\n", x$inverse_quality))

  if (!is.null(x$jacobian_stats)) {
    cat("\nJacobian Statistics:\n")
    js <- x$jacobian_stats
    cat(sprintf("  det range:      [%.4f, %.4f]\n", js$det_min, js$det_max))
    cat(sprintf("  det median:     %.4f\n", js$det_median))
    cat(sprintf("  log(det) mean:  %.4f (SD: %.4f)\n", js$log_det_mean, js$log_det_sd))
    cat(sprintf("  Folding voxels: %d\n", js$n_folding))
    cat(sprintf("  Diffeomorphic:  %s\n", if (js$is_diffeomorphic) "YES" else "NO"))
  }

  invisible(x)
}
```

---

## Part 5: Implementation Roadmap

### Phase 1: Core Jacobians (v0.2)
- [ ] `JacobianField` class
- [ ] `jacobian()` for IdentityMorphism, Affine3DMorphism
- [ ] `jacobian_det()` for same
- [ ] C++ `cpp_warp_jacobian()` for Warp3DMorphism
- [ ] `jacobian()` for MorphismPath (chain rule)
- [ ] Unit tests for all Jacobian computations

### Phase 2: Basic Resampling (v0.3)
- [ ] `Sampler` class
- [ ] `volume_sampler()` constructor
- [ ] `Grid` class and utilities
- [ ] `resample()` core function
- [ ] `resample_volume()` convenience wrapper
- [ ] C++ `cpp_sample_volume()` with linear/nearest

### Phase 3: Surface and Modulation (v0.4)
- [ ] `surface_sampler()` constructor
- [ ] `sample_volume_on_surface()`
- [ ] Jacobian modulation in `resample()`
- [ ] C++ `cpp_ribbon_sample_volume()`

### Phase 4: Vectors and Tensors (v0.5)
- [ ] `transform_vectors()`
- [ ] `transform_tensors()`
- [ ] Polar decomposition for FSL-style tensor reorientation

### Phase 5: QC and Analysis (v0.6)
- [ ] `jacobian_map()`
- [ ] `check_diffeomorphic()`
- [ ] `summarize_morphism()`
- [ ] Visualization helpers (if scope appropriate)

---

## Appendix: Design Decisions

### Why Samplers are functions, not data containers

The Sampler abstraction wraps `data + interpolation` into a callable. This:
1. Hides backend differences (RNifti vs neuroim2 vs array)
2. Allows lazy/streaming evaluation for huge datasets
3. Makes composition natural: `resample(sampler, morphism, coords)`
4. Keeps neurotransform agnostic to data representation

### Why JacobianField is a separate class

Rather than returning raw arrays, JacobianField:
1. Carries provenance (morphism hash, mode)
2. Enables clean operator overloading (`%*%`, `solve`)
3. Provides indexed access `J[i]` for single matrices
4. Makes intent clear in function signatures

### Why mode = "pullback" is the default

neurotransform's convention is pullback semantics:
- `transform(m, coords)` maps target → source
- So the natural Jacobian is d(source)/d(target)
- This is what you need for most resampling operations
- Pushforward is available when needed (e.g., vector field visualization)

### Why modulation is opt-in

Most neuroimaging resampling does NOT want Jacobian modulation:
- Anatomical images: preserve intensities
- Statistical maps: preserve values

Modulation is needed for:
- Modulated VBM (preserving local tissue amount)
- Some registration metrics

Making it explicit prevents silent errors.
