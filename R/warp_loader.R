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

.looks_like_coef_warp <- function(path) {
  lower <- tolower(basename(path))
  grepl("coef|coeff|warpcoef|fieldcoef", lower)
}

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

#' Default FSL coefficient warp loader
#'
#' Loads an FNIRT coefficient field as B-spline coefficients. This is not a
#' dense displacement field and must be evaluated through the cubic B-spline
#' basis at query coordinates.
#'
#' @param path Path to coefficient NIfTI file
#' @return List with array, dim, world_to_vox, vox_to_world, mode
#' @keywords internal
load_warp_fsl_coef <- function(path) {
  if (!file.exists(path)) stop("Coefficient warp file not found: ", path)

  vol <- tryCatch(neuroim2::read_vec(path), error = function(e) NULL)
  if (is.null(vol)) vol <- tryCatch(neuroim2::read_vol(path), error = function(e) NULL)
  if (is.null(vol)) stop("Failed to read coefficient warp: ", path)

  dim4 <- dim(vol)
  if (length(dim4) < 4 || dim4[4] < 3) {
    stop("Coefficient warp must be 4D with last dimension length >= 3")
  }

  raw <- as.array(vol)
  nvox <- prod(dim4[1:3])
  arr <- numeric(3 * nvox)
  idx <- seq_len(nvox)
  arr[3L * (idx - 1L) + 1L] <- as.numeric(raw[, , , 1])  # X coefficients
  arr[3L * (idx - 1L) + 2L] <- as.numeric(raw[, , , 2])  # Y coefficients
  arr[3L * (idx - 1L) + 3L] <- as.numeric(raw[, , , 3])  # Z coefficients

  vox_to_world <- neuroim2::trans(vol)
  inv_aff <- solve(vox_to_world)
  list(
    array = arr,
    dim = as.integer(dim4[1:3]),
    world_to_vox = inv_aff,
    vox_to_world = vox_to_world,
    mode = "bspline_coefficients"
  )
}

#' Load ANTs composite H5 displacement field
#'
#' Reads the displacement field parameters from an ANTs composite H5 file
#' (TransformType == "DisplacementFieldTransform_*") and returns an internal
#' warp list. Requires hdf5r.
#'
#' @section Coordinate System Notes:
#' ANTs/ITK H5 domains and displacement components are stored in LPS physical
#' coordinates. This loader converts both the grid mapping and vector components
#' to the package's internal RAS convention before returning them.
#'
#' @section Embedded Affine:
#' ANTs composite H5 files often contain both a displacement field and an
#' affine transform. The embedded affine is extracted and converted from LPS to
#' RAS. ITK image-resampling transforms already use pullback semantics, so the
#' affine is not inverted. Use `ants_h5_morphism()` with
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

  transform_types <- .ants_h5_components(tg, c("AffineTransform", "DisplacementFieldTransform"))
  keys <- names(transform_types)
  disp_keys <- keys[grepl("DisplacementFieldTransform", transform_types)]
  aff_keys <- keys[grepl("AffineTransform", transform_types)]
  if (length(disp_keys) != 1L) {
    stop("Expected exactly one DisplacementFieldTransform in H5, found ", length(disp_keys), ": ", path)
  }
  if (length(aff_keys) > 1L) {
    stop("Multiple embedded affine transforms are not yet supported: ", path)
  }
  disp_key <- disp_keys[[1L]]

  grp <- tg[[disp_key]]
  fixed <- .ants_h5_read_dataset(grp, "TransformFixedParameters")
  params <- .ants_h5_read_dataset(grp, "TransformParameters")

  if (length(fixed) != 18L || any(!is.finite(fixed))) {
    .stop_transform_file("Displacement fixed parameters must contain exactly 18 finite values.")
  }
  if (any(fixed[1:3] < 1 | fixed[1:3] > .Machine$integer.max | fixed[1:3] != floor(fixed[1:3]))) {
    .stop_transform_file("Displacement grid size must contain positive integer dimensions.")
  }
  size <- as.integer(fixed[1:3])
  origin <- fixed[4:6]
  spacing <- fixed[7:9]
  direction <- matrix(fixed[10:18], nrow = 3, byrow = TRUE)
  if (any(spacing <= 0) || abs(det(direction)) < 1e-12) {
    .stop_transform_file("Displacement grid must have positive spacing and an invertible direction matrix.")
  }

  # Build the native LPS grid mapping, then convert its world coordinates to RAS.
  vox_to_lps <- diag(4)
  vox_to_lps[1:3, 1:3] <- direction %*% diag(spacing)
  vox_to_lps[1:3, 4] <- origin
  lps_to_ras <- diag(c(-1, -1, 1, 1))
  vox_to_world <- lps_to_ras %*% vox_to_lps
  world_to_vox <- solve(vox_to_world)

  # Parameters are stored voxel-major and interleaved in LPS:
  # [dx0, dy0, dz0, dx1, dy1, dz1, ...].
  nvox <- prod(size)
  if (length(params) != 3 * nvox || any(!is.finite(params))) {
    .stop_transform_file("Displacement parameters must contain exactly 3 finite values per voxel.")
  }
  arr <- as.numeric(params[seq_len(3L * nvox)])
  x_idx <- seq.int(1L, length(arr), by = 3L)
  y_idx <- seq.int(2L, length(arr), by = 3L)
  arr[x_idx] <- -arr[x_idx]
  arr[y_idx] <- -arr[y_idx]

  # Optional affine in the composite
  aff_mat <- NULL
  if (length(aff_keys) == 1L) {
    aff_key <- aff_keys[[1L]]
    grp_aff <- tg[[aff_key]]
    aff_mat_lps <- .ants_h5_affine(grp_aff)
    aff_mat <- lps_to_ras %*% aff_mat_lps %*% lps_to_ras
  }

  component_keys <- keys[grepl("AffineTransform|DisplacementFieldTransform", transform_types)]
  component_types <- transform_types[match(component_keys, keys)]
  transform_order <- ifelse(
    grepl("AffineTransform", component_types),
    "affine",
    "warp"
  )

  list(
    array = arr,
    dim = size,
    world_to_vox = world_to_vox,
    vox_to_world = vox_to_world,
    affine = aff_mat,
    transform_order = unname(transform_order)
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
    if (identical(morphism@warp_type, "fsl") && .looks_like_coef_warp(morphism@warp_path)) {
      stop("Warp looks like an FNIRT coefficient field. Use warp_type='fsl_coef' (or read_transform(..., type='fsl_coef')).")
    }
    default_loader_name <- switch(
      morphism@warp_type,
      "ants_h5" = "ants_h5",
      "ants" = "neuroim2",
      "fsl" = "neuroim2",
      "fsl_coef" = "fsl_coef",
      "freesurfer" = "neuroim2",
      "afni" = "neuroim2",
      "dense" = "neuroim2",
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
  # FSL source/reference geometry changes the decoded values for the same file.
  key <- if (identical(morphism@warp_type, "fsl")) {
    paste0(morphism@warp_path, "::", morphism_hash(morphism))
  } else morphism@warp_path
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }
  value <- loader(morphism@warp_path)
  if (identical(morphism@warp_type, "fsl")) {
    value <- .fsl_dense_to_ras_displacement(value, morphism@params)
  }

  # ANTs and AFNI NIfTI warps store vector components in LPS/DICOM
  # coordinates even though neuroim2 exposes the NIfTI grid in RAS.
  # LPS to RAS negates X and Y displacement components.
  # Do this at load time so all code paths (transform(), resample_volume(), etc.)
  # get RAS-convention displacements without needing special handling.
  if (morphism@warp_type %in% c("ants", "afni")) {
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
