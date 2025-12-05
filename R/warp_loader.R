#' @title Warp Field Loader Registry
#' @name warp_loader
#' @description
#' Pluggable loader system for warp displacement fields. Provides a registry
#' for different file formats (NIfTI via neuroim2, AFNI, FSL, ANTs H5)
#' with morphism-level caching.
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

#' Default neuroim2 warp loader
#'
#' Loads a NIfTI displacement field using neuroim2 (DenseNeuroVec).
#'
#' @param path Path to NIfTI warp file
#' @return List with array, dim, world_to_vox, vox_to_world
#' @keywords internal
load_warp_neuroim2 <- function(path) {
  if (!file.exists(path)) stop("Warp file not found: ", path)

  vol <- tryCatch(neuroim2::read_vec(path), error = function(e) NULL)
  if (is.null(vol)) vol <- tryCatch(neuroim2::read_vol(path), error = function(e) NULL)
  if (is.null(vol)) stop("Failed to read warp: ", path)

  dim4 <- dim(vol)
  if (length(dim4) < 4 || dim4[4] < 3) {
    stop("Warp must be 4D with last dimension length >= 3")
  }

  # C++ code expects interleaved (X, Y, Z, 3) layout where each voxel's
  # 3 displacement values are contiguous: [dx0, dy0, dz0, dx1, dy1, dz1, ...]
  # R's as.numeric() gives planar layout: [dx0, dx1, ..., dy0, dy1, ..., dz0, dz1, ...]
  # Need to explicitly interleave.
  raw <- as.array(vol)
  nvox <- prod(dim4[1:3])
  arr <- numeric(3 * nvox)
  # Interleave: position i*3+0=dx_i, i*3+1=dy_i, i*3+2=dz_i
  idx <- seq_len(nvox)
  arr[3L * (idx - 1L) + 1L] <- as.numeric(raw[, , , 1])  # X components
  arr[3L * (idx - 1L) + 2L] <- as.numeric(raw[, , , 2])  # Y components
  arr[3L * (idx - 1L) + 3L] <- as.numeric(raw[, , , 3])  # Z components

  vox_to_world <- neuroim2::trans(vol)
  inv_aff <- solve(vox_to_world)
  list(
    array = arr,
    dim = as.integer(dim4[1:3]),
    world_to_vox = inv_aff,
    vox_to_world = vox_to_world
  )
}

#' Load ANTs composite H5 displacement field
#'
#' Reads the displacement field parameters from an ANTs composite H5 file
#' (TransformType == "DisplacementFieldTransform_*") and returns an internal
#' warp list. Requires hdf5r.
#'
#' @section Coordinate System Notes:
#' ANTs/ITK documentation states that H5 files use LPS (Left-Posterior-Superior)
#' coordinates. However, the actual coordinate handling depends on the warp's
#' `vox_to_world` matrix:
#'
#' \itemize{
#'   \item If `vox_to_world` has negative diagonal elements, coordinates are in
#'         LPS physical space and explicit RAS<->LPS conversion may be needed.
#'   \item If `vox_to_world` is identity (origin=0, spacing=1, direction=I),
#'         coordinates are just voxel indices - no conversion is needed.
#' }
#'
#' The `world_to_vox` matrix returned by this loader handles coordinate
#' conversion automatically. Do NOT add explicit LPS flipping in the transform
#' code - this will cause out-of-bounds voxel lookups when `vox_to_world` is
#' identity.
#'
#' @section Embedded Affine:
#' ANTs composite H5 files often contain both a displacement field and an
#' affine transform. The embedded affine is extracted, converted from LPS to
#' RAS, and inverted for pullback semantics. Use `ants_h5_morphism()` with
#' `apply_affine=TRUE` to get a MorphismPath that applies both transforms.
#'
#' @param path Path to ANTs H5 file
#' @return List with array, dim, world_to_vox, vox_to_world, and optionally affine
#' @keywords internal
load_warp_ants_h5 <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r required to read ANTs H5 warps")
  }
  if (!file.exists(path)) stop("H5 warp not found: ", path)

  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  tg <- h5[["TransformGroup"]]
  if (is.null(tg)) stop("TransformGroup not found in H5 file: ", path)

  # Find displacement transform entry
  disp_key <- NULL
  for (k in names(tg)) {
    tt <- tg[[k]][["TransformType"]]$read()
    if (grepl("DisplacementFieldTransform", tt)) {
      disp_key <- k
      break
    }
  }
  if (is.null(disp_key)) stop("No DisplacementFieldTransform in H5: ", path)

  grp <- tg[[disp_key]]
  fixed <- grp[["TransformFixedParameters"]]$read()
  params <- grp[["TransformParameters"]]$read()

  size <- as.integer(fixed[1:3])
  origin <- fixed[4:6]
  spacing <- fixed[7:9]
  direction <- matrix(fixed[10:18], nrow = 3, byrow = TRUE)

  # Build voxel-to-world in LPS (ITK/ANTs native coordinate system)
  vox_to_world <- diag(4)
  vox_to_world[1:3, 1:3] <- direction %*% diag(spacing)
  vox_to_world[1:3, 4] <- origin
  world_to_vox <- solve(vox_to_world)

  # Parameters are stored voxel-major (dx, dy, dz) per voxel in LPS.
  # C++ expects interleaved layout: [dx0, dy0, dz0, dx1, dy1, dz1, ...]
  # ANTs H5 stores params as: X varies fastest, then Y, then Z, then component (planar)
  # So we need to interleave.
  nvox <- prod(size)
  arr <- numeric(3 * nvox)
  idx <- seq_len(nvox)
  arr[3L * (idx - 1L) + 1L] <- params[idx]              # X components
  arr[3L * (idx - 1L) + 2L] <- params[nvox + idx]       # Y components
  arr[3L * (idx - 1L) + 3L] <- params[2L * nvox + idx]  # Z components

  # Optional affine in the composite
  aff_key <- NULL
  aff_mat <- NULL
  for (k in names(tg)) {
    tt <- tg[[k]][["TransformType"]]$read()
    if (grepl("AffineTransform", tt)) {
      aff_key <- k; break
    }
  }
  if (!is.null(aff_key)) {
    grp_aff <- tg[[aff_key]]
    p_aff <- grp_aff[["TransformParameters"]]$read()
    fixed_aff <- grp_aff[["TransformFixedParameters"]]$read()
    if (length(p_aff) >= 12) {
      A <- matrix(p_aff[1:9], nrow = 3, byrow = TRUE)
      tvec <- p_aff[10:12]
      aff_mat_lps <- diag(4)
      aff_mat_lps[1:3, 1:3] <- A
      if (length(fixed_aff) == 3) {
        aff_mat_lps[1:3, 4] <- tvec + fixed_aff - A %*% fixed_aff
      } else {
        aff_mat_lps[1:3, 4] <- tvec
      }
      # Convert from LPS (ITK convention) to RAS (neuroim2 convention)
      flip <- diag(c(-1, -1, 1, 1))
      aff_mat_ras <- flip %*% aff_mat_lps %*% flip
      # ANTs stores forward transform; for pullback we need inverse
      aff_mat <- solve(aff_mat_ras)
    }
  }

  list(
    array = arr,
    dim = size,
    world_to_vox = world_to_vox,
    vox_to_world = vox_to_world,
    affine = aff_mat
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
#' @param name Loader name (default "neuroim2")
#' @return Loader function
#' @export
#' @examples
#' loader <- get_loader("neuroim2")
get_loader <- function(name = NULL) {
  if (is.null(name)) name <- "neuroim2"
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

  # Resolve loader based on warp_type
  if (is.null(loader)) {
    default_loader_name <- switch(
      morphism@warp_type,
      "ants_h5" = "ants_h5",
      "ants" = "neuroim2",
      "fsl" = "neuroim2",
      "freesurfer" = "neuroim2",
      "afni" = "neuroim2",
      "neuroim2"
    )
    if (!exists(default_loader_name, envir = .warp_loader_registry, inherits = FALSE)) {
      default_loader_name <- "neuroim2"
    }
    loader <- get_loader(default_loader_name)
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

  # AFNI 3dQwarp stores displacements in DICOM/LPS convention.

  # LPS to RAS: negate X and Y displacement components.
  # Do this at load time so all code paths (transform(), resample_volume(), etc.)
  # get RAS-convention displacements without needing special handling.
  if (morphism@warp_type == "afni") {
    # Array is stored as (X, Y, Z, 3) flattened - each voxel has 3 contiguous values
    # Layout: [dx0, dy0, dz0, dx1, dy1, dz1, ...]
    # To negate X: indices 1, 4, 7, ... (seq from 1 by 3)
    # To negate Y: indices 2, 5, 8, ... (seq from 2 by 3)
    nvox <- prod(value$dim)
    x_idx <- seq(1, 3 * nvox, by = 3)  # Indices for X component
    y_idx <- seq(2, 3 * nvox, by = 3)  # Indices for Y component
    value$array[x_idx] <- -value$array[x_idx]
    value$array[y_idx] <- -value$array[y_idx]
  }

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
