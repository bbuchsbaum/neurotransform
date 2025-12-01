#' @title Warp Field Loader Registry
#' @name warp_loader
#' @description
#' Pluggable loader system for warp displacement fields. Provides a registry
#' for different file formats (NIfTI via RNifti, AFNI, FSL, etc.) with
#' morphism-level caching.
#'
#' @section Loader Function Contract:
#' A loader function must accept a file path and return a list with:
#' \itemize{
#'   \item \code{array}: Numeric vector of displacement values (flattened 3 x X x Y x Z)
#'   \item \code{dim}: Integer vector c(X, Y, Z)
#'   \item \code{world_to_vox}: 4x4 matrix mapping world coords to voxel indices
#' }
NULL

# Internal registry environment
.warp_loader_registry <- new.env(parent = emptyenv())

#' Default RNifti warp loader
#'
#' Loads a NIfTI displacement field using RNifti package.
#'
#' @param path Path to NIfTI warp file
#' @return List with array, dim, world_to_vox
#' @keywords internal
load_warp_rnifti <- function(path) {
  if (!requireNamespace("RNifti", quietly = TRUE)) {
    stop("RNifti required for warp transforms; please install RNifti")
  }
  vol <- RNifti::readNifti(path)
  dim4 <- dim(vol)
  # Permute to get displacement components first: (3, X, Y, Z) -> flattened
  arr <- as.numeric(aperm(vol, c(4, 1, 2, 3)))
  vox_to_world <- RNifti::xform(vol)
  inv_aff <- solve(vox_to_world)
  list(
    array = arr,
    dim = as.integer(dim4[1:3]),
    world_to_vox = inv_aff,
    vox_to_world = vox_to_world
  )
}

#' Register a named warp loader
#'
#' Adds a loader function to the registry. Loaders are used to read warp
#' displacement fields from various file formats.
#'
#' @param name Character identifier for the loader
#' @param loader Function(path) -> list(array, dim, world_to_vox)
#' @return Invisibly returns the name
#' @export
#' @examples
#' # Register a custom loader
#' my_loader <- function(path) {
#'   # Load and return list(array, dim, world_to_vox)
#' }
#' register_loader("my_format", my_loader)
register_loader <- function(name, loader) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  if (!is.function(loader)) stop("loader must be a function")
  assign(name, loader, envir = .warp_loader_registry)
  invisible(name)
}

#' Retrieve a warp loader by name
#'
#' Gets a loader function from the registry.
#'
#' @param name Loader name (default "rnifti")
#' @return Loader function
#' @export
#' @examples
#' loader <- get_loader("rnifti")
get_loader <- function(name = NULL) {
  if (is.null(name)) name <- "rnifti"
  if (!exists(name, envir = .warp_loader_registry, inherits = FALSE)) {
    stop("No warp loader registered for name: ", name)
  }
  get(name, envir = .warp_loader_registry, inherits = FALSE)
}

#' List all registered loaders
#'
#' @return Character vector of registered loader names
#' @export
#' @examples
#' list_loaders()
list_loaders <- function() {
  ls(envir = .warp_loader_registry)
}

#' Load warp array with caching
#'
#' Loads a warp displacement field, caching the result in the morphism's
#' cache environment for efficiency.
#'
#' @param morphism A Warp3DMorphism object
#' @param loader Loader function or name (default: registry default)
#' @param cache_env Optional cache environment (default: morphism's cache)
#' @return List with array, dim, world_to_vox
#' @keywords internal
load_warp_array <- function(morphism, loader = NULL, cache_env = NULL) {
  if (missing(morphism) || morphism_kind(morphism) != "warp3d") {
    stop("morphism must be a Warp3DMorphism")
  }

  # Resolve loader
  if (is.null(loader)) {
    loader <- get_loader()
  } else if (is.character(loader)) {
    loader <- get_loader(loader)
  }
  if (!is.function(loader)) stop("loader must be a function")

  cache_env <- cache_env %||% morphism@cache %||% new_cache_env()
  key <- morphism@warp_path
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }
  value <- loader(morphism@warp_path)
  assign(key, value, envir = cache_env)
  value
}

# =============================================================================
# COMPATIBILITY ALIASES
# =============================================================================

#' @rdname register_loader
#' @export
register_warp_loader <- register_loader

#' @rdname get_loader
#' @export
get_warp_loader <- get_loader
