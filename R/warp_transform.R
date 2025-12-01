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

  if (identical(method, "cubic")) {
    return(cpp_apply_warp_field_cubic(coords, warp$array, warp$dim, warp$world_to_vox))
  }
  cpp_apply_warp_field(coords, warp$array, warp$dim, warp$world_to_vox)
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
