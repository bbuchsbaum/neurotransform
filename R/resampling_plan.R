#' @title Resampling plans for repeated volume resampling
#' @name resampling_plan
#' @description
#' Plan once, apply many times: precompute target voxel centres and flatten the
#' morphism path so repeated resampling of multiple volumes is cheap while
#' keeping the public API small.
NULL

.resampling_plan_cache <- new.env(parent = emptyenv())

# Extract modification times for file-backed morphisms
morphism_file_mtimes <- function(m) {
  times <- numeric(0)
  if (is(m, "Warp3DMorphism")) {
    if (nzchar(m@warp_path) && file.exists(m@warp_path)) {
      times <- c(times, unclass(file.info(m@warp_path)$mtime))
    }
    if (nzchar(m@inverse_path) && file.exists(m@inverse_path)) {
      times <- c(times, unclass(file.info(m@inverse_path)$mtime))
    }
  }
  if (is(m, "MorphismPath")) {
    for (p in m@morphisms) {
      times <- c(times, morphism_file_mtimes(p))
    }
  }
  times
}

#' Build a resampling plan
#'
#' @param morphism Morphism or MorphismPath
#' @param source_grid Grid describing the source volume geometry
#' @param target_grid Grid describing the target geometry
#' @param interpolation "linear", "cubic", or "nearest"
#' @param reuse_count Hint for expected reuse; >1 enables caching
#' @param cache Logical; cache the plan by hash
#' @return A ResamplingPlan object
#' @export
make_resampling_plan <- function(morphism, source_grid, target_grid,
                                 interpolation = c("linear", "cubic", "nearest"),
                                 reuse_count = 1L, cache = TRUE) {
  interpolation <- match.arg(interpolation)
  if (!inherits(source_grid, "Grid") || !inherits(target_grid, "Grid")) {
    stop("source_grid and target_grid must be Grid objects")
  }

  mtime_sig <- morphism_file_mtimes(morphism)
  key <- compute_hash("resample_plan",
                      morphism_hash(morphism),
                      mtime_sig,
                      source_grid@dims, source_grid@affine,
                      target_grid@dims, target_grid@affine,
                      interpolation)

  precompute <- reuse_count > 1L

  if (cache && precompute && exists(key, envir = .resampling_plan_cache, inherits = FALSE)) {
    return(get(key, envir = .resampling_plan_cache, inherits = FALSE))
  }

  coords_tgt <- grid_coords(target_grid)
  # Pull back into source world, then convert to source voxel coords (0-based)
  coords_src_world <- coords_tgt

  path <- if (is(morphism, "MorphismPath")) morphism@morphisms else list(morphism)
  # Attempt full flattening: build steps list for affine/warp sequence
  can_flatten <- all(vapply(path, function(m) morphism_kind(m) %in% c("affine3d", "warp3d"), logical(1))) &&
    !any(vapply(path, function(m) is(m, "Warp3DMorphism") &&
                  (identical(m@warp_type, "fsl_coef") || identical(m@warp_type, "ants_h5")),
                logical(1)))
  if (can_flatten) {
    path_pullback <- rev(path)
    steps <- vector("list", length(path_pullback))
    for (i in seq_along(path_pullback)) {
      m <- path_pullback[[i]]
      if (morphism_kind(m) == "affine3d") {
        steps[[i]] <- list(kind = "affine", matrix = m@matrix)
      } else {
        w <- load_warp_array(m)
        steps[[i]] <- list(kind = "warp",
                           field = w$array,
                           dim = w$dim,
                           world_to_vox = w$world_to_vox,
                           method = m@params$warp_method %||% "linear")
      }
    }
    coords_src_world <- cpp_path_apply_steps(coords_tgt, steps)
  } else {
    coords_src_world <- transform(morphism, coords_tgt)
  }
  w2v <- solve(source_grid@affine)
  hom <- cbind(coords_src_world, 1)
  coords_src_vox <- (hom %*% t(w2v))[, 1:3, drop = FALSE]

  plan_handle <- if (precompute) {
    cpp_make_resample_plan(source_grid@dims,
                           coords_src_vox,
                           method = interpolation)
  } else {
    NULL
  }

  plan <- structure(
    list(
      morphism = morphism,
      source_grid = source_grid,
      target_grid = target_grid,
      target_coords = coords_tgt,
      interpolation = interpolation,
      reuse_count = as.integer(reuse_count),
      key = key,
      handle = plan_handle,
      coords_src_vox = coords_src_vox
    ),
    class = "ResamplingPlan"
  )

  if (cache && precompute) {
    assign(key, plan, envir = .resampling_plan_cache)
  }
  plan
}

#' Apply a resampling plan to a volume
#'
#' @param plan ResamplingPlan
#' @param source_data 3D or 4D array in source_grid geometry
#' @param outside Value for out-of-bounds sampling
#' @param modulate Jacobian modulation: "none", "jacobian", "sqrt_jacobian"
#' @return Array with target_grid dimensions (and trailing data dim if 4D input)
#' @export
apply_resampling_plan <- function(plan, source_data, outside = NA_real_,
                                  modulate = c("none", "jacobian", "sqrt_jacobian")) {
  stopifnot(inherits(plan, "ResamplingPlan"))
  modulate <- match.arg(modulate)

  src_dims <- dim(source_data)[1:3]
  if (!all(src_dims == plan$source_grid@dims)) {
    stop("source_data dims do not match plan$source_grid dims")
  }

  if (is.null(plan$handle) || !identical(modulate, "none")) {
    # Fall back to existing path when modulation is required
    sampler <- volume_sampler(source_data, affine = plan$source_grid@affine,
                              method = plan$interpolation, outside = outside)
    vals <- resample(sampler, morphism = plan$morphism,
                     coords = plan$target_coords, modulate = modulate)
    tgt_dims <- plan$target_grid@dims
    if (length(dim(source_data)) > 3) {
      vdim <- dim(source_data)[4]
      return(array(vals, dim = c(tgt_dims, vdim)))
    }
    return(array(vals, dim = tgt_dims))
  }

  tgt_dims <- plan$target_grid@dims
  if (length(dim(source_data)) > 3) {
    vdim <- dim(source_data)[4]
    out <- array(NA_real_, dim = c(tgt_dims, vdim))
    for (i in seq_len(vdim)) {
      slice <- source_data[,,, i, drop = TRUE]
      out[,,, i] <- cpp_apply_resample_plan(plan$handle, slice, outside)
    }
    return(out)
  }

  vals <- cpp_apply_resample_plan(plan$handle, source_data, outside)
  array(vals, dim = tgt_dims)
}

#' @export
print.ResamplingPlan <- function(x, ...) {
  cat("<ResamplingPlan | interp=", x$interpolation,
      " | source dims=", paste(x$source_grid@dims, collapse = "x"),
      " | target dims=", paste(x$target_grid@dims, collapse = "x"), ">\n", sep = "")
}
