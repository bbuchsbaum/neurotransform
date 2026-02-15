#' @title Warp Transform Implementation
#' @name warp_transform
#' @description
#' Internal functions for applying warp displacement fields to coordinates.
#' Uses C++ for trilinear interpolation of displacement fields.
#' @keywords internal
NULL

#' Apply warp transform to coordinates
#'
#' Internal function that handles warp field loading, caching, and application.
#'
#' @param morphism A Warp3DMorphism object
#' @param coords Numeric matrix (N x 3) of world coordinates
#' @return Transformed coordinates
#' @keywords internal
warp_transform_coords <- function(morphism, coords) {
  if (morphism_kind(morphism) != "warp3d") {
    return(transform(morphism, coords)) # fallback
  }
  if (!nzchar(morphism@warp_path)) {
    stop("Warp3DMorphism missing warp_path")
  }

  # AFNI warps need special handling for RAI convention
  if (identical(morphism@warp_type, "afni")) {
    return(afni_warp_transform_coords(morphism, coords))
  }

  cache_env <- morphism@cache %||% new_cache_env()
  warp <- load_warp_array(morphism, cache_env = cache_env)

  method <- morphism@params$warp_method %||% "linear"

  # FNIRT coefficient fields are not dense displacements; evaluate cubic
  # B-spline basis directly at query coordinates.
  if (identical(warp$mode %||% "", "bspline_coefficients")) {
    return(cpp_apply_bspline_coeff_field(coords, warp$array, warp$dim, warp$world_to_vox))
  }

  coords_in <- coords

  # Note: For ANTs composite H5 files, the embedded affine is typically meant to be
  # applied AFTER the displacement field in the pullback direction. However, when
  # using ants_h5_morphism() with apply_affine=TRUE, the affine is extracted as a
  # separate MorphismPath component. Here we don't apply the embedded affine -
  # that should be handled by MorphismPath if needed.
  embedded_affine <- warp$affine
  warp$affine <- NULL  # Don't apply here; handled separately

  # Handle absolute vs relative displacement fields
  def_type <- morphism@params$def_type %||% "relative"
  if (identical(def_type, "absolute")) {
    # Convert absolute coords to displacement once and cache
    cache_key <- paste0(morphism@warp_path, "::relative")
    if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
      disp <- get(cache_key, envir = cache_env, inherits = FALSE)
    } else {
      disp_vec <- cpp_absolute_to_displacement(warp$array, warp$dim, solve(warp$world_to_vox))
      disp <- list(array = disp_vec, dim = warp$dim, world_to_vox = warp$world_to_vox)
      assign(cache_key, disp, envir = cache_env)
    }
    warp <- disp
  }

  result <- if (identical(method, "cubic")) {
    cpp_apply_warp_field_cubic(coords_in, warp$array, warp$dim, warp$world_to_vox)
  } else {
    cpp_apply_warp_field(coords_in, warp$array, warp$dim, warp$world_to_vox)
  }

  # Outside-warp locations default to identity displacement, mirroring common
  # dense-field behavior where undefined samples imply zero displacement.
  invalid <- !is.finite(result[, 1]) | !is.finite(result[, 2]) | !is.finite(result[, 3])
  if (any(invalid)) {
    result[invalid, ] <- coords_in[invalid, , drop = FALSE]
  }

  result
}

#' Compose two warp morphisms
#'
#' Composes warp B then A into a single displacement field.
#'
#' @param warpB First warp (applied first)
#' @param warpA Second warp (applied second)
#' @return List with composed array, dim, world_to_vox
#' @keywords internal
compose_warps <- function(warpB, warpA) {
  cache_envB <- warpB@cache %||% new_cache_env()
  cache_envA <- warpA@cache %||% new_cache_env()
  wB <- load_warp_array(warpB, cache_env = cache_envB)
  wA <- load_warp_array(warpA, cache_env = cache_envA)
  res <- cpp_compose_warp_fields(
    wA$array, wA$dim, wA$world_to_vox,
    wB$array, wB$dim, wB$world_to_vox
  )
  list(array = res$field, dim = res$dim, world_to_vox = wB$world_to_vox)
}
