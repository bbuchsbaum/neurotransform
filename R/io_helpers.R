#' Simple image I/O using neuroim2
#'
#' @param path File path
#' @return neuroim2 volume (DenseNeuroVol or DenseNeuroVec) or array fallback
#' @export
read_image <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  # Prefer read_vec to preserve 4D
  vol <- tryCatch(neuroim2::read_vec(path), error = function(e) NULL)
  if (is.null(vol)) vol <- neuroim2::read_vol(path)
  vol
}

#' Write image using neuroim2
#'
#' @param x neuroim2 volume or array
#' @param path Output path
#' @export
write_image <- function(x, path) {
  if (inherits(x, "DenseNeuroVec")) {
    return(neuroim2::write_vec(x, path, format = "nifti"))
  }
  if (inherits(x, "DenseNeuroVol")) {
    return(neuroim2::write_vol(x, path, format = "nifti"))
  }
  # If array, build a minimal space assuming identity affine
  dims <- dim(x)
  space <- neuroim2::NeuroSpace(dims, trans = diag(4))
  if (length(dims) > 3) {
    return(neuroim2::write_vec(neuroim2::DenseNeuroVec(x, space), path, format = "nifti"))
  }
  neuroim2::write_vol(neuroim2::DenseNeuroVol(x, space), path, format = "nifti")
}

#' Extract Grid from volume or grid-like input
#'
#' @param x neuroim2 volume, array, or Grid
#' @return Grid
#' @export
grid_of <- function(x) {
  if (inherits(x, "Grid")) return(x)
  if (requireNamespace("neuroim2", quietly = TRUE) &&
      (inherits(x, "DenseNeuroVol") || inherits(x, "DenseNeuroVec"))) {
    return(grid_spec(dim(x)[1:3], neuroim2::trans(x)))
  }
  if (is.array(x)) {
    return(grid_spec(dim(x)[1:3], diag(4)))
  }
  stop("Unsupported type for grid_of")
}

#' Read transform from file into a Morphism
#'
#' @param path Path or Morphism
#' @param type Optional type hint: ants_h5, ants, fsl, afni
#' @param source Optional source id
#' @param target Optional target id
#' @param apply_affine Logical; for ANTs H5 files, whether to include embedded affine
#' @param ... Passed to morphism constructors
#' @return A Morphism or MorphismPath object
#' @export
read_transform <- function(path, type = NULL, source = NULL, target = NULL, apply_affine = TRUE, ...) {
  if (is(path, "Morphism")) return(path)
  if (!is.character(path) || length(path) != 1) stop("path must be a file path or Morphism")
  if (!file.exists(path)) stop("Transform file not found: ", path)
  if (is.null(type)) {
    if (grepl("\\.h5$", path, ignore.case = TRUE)) {
      type <- "ants_h5"
    } else if (grepl("\\.(mat|txt)$", path, ignore.case = TRUE)) {
      type <- "fsl"
    } else if (grepl("\\.1D$", path, ignore.case = TRUE)) {
      type <- "afni"
    } else {
      type <- "ants"
    }
  }
  if (identical(type, "ants_h5")) {
    return(ants_h5_morphism(path, source = source %||% "source", target = target %||% "target",
                             apply_affine = apply_affine, ...))
  }
  if (type %in% c("ants", "fsl", "afni")) {
    return(Warp3DMorphism(source %||% "source", target %||% "target",
                          warp_path = path, warp_type = type))
  }
  stop("Unsupported transform type: ", type)
}

#' Coerce list/paths to MorphismPath
#'
#' @param ... Morphism objects or file paths to transforms
#' @return MorphismPath object
#' @export
as_morphism_path <- function(...) {
  args <- list(...)
  flat <- unlist(args, recursive = FALSE)
  morphs <- lapply(flat, function(x) {
    if (is.character(x)) read_transform(x) else x
  })
  methods::new("MorphismPath",
               morphisms = morphs,
               source = morphs[[1]]@source,
               target = morphs[[length(morphs)]]@target)
}

#' High-level resampling helper
#'
#' @param moving Source volume (neuroim2 or array)
#' @param target Target volume or Grid
#' @param transform Morphism, MorphismPath, or file path
#' @param method Interpolation method
#' @param modulate Jacobian modulation
#' @return neuroim2 volume (matching target geometry)
#' @export
resample_to <- function(moving, target, transform,
                        method = "linear",
                        modulate = "none") {
  tgt_grid <- grid_of(target)
  morph <- if (is(transform, "Morphism") || is(transform, "MorphismPath")) {
    transform
  } else {
    read_transform(transform)
  }

  data_in <- if (requireNamespace("neuroim2", quietly = TRUE) &&
                 (inherits(moving, "DenseNeuroVol") || inherits(moving, "DenseNeuroVec"))) {
    moving
  } else moving

  out_array <- resample_volume(data_in, morphism = morph, target = tgt_grid,
                               method = method, modulate = modulate)

  if (requireNamespace("neuroim2", quietly = TRUE)) {
    space <- neuroim2::NeuroSpace(tgt_grid@dims, trans = tgt_grid@affine)
    if (length(dim(out_array)) == 4) {
      if (dim(out_array)[4] == 1) {
        out_array <- array(out_array[, , , 1, drop = FALSE], dim = tgt_grid@dims)
        return(neuroim2::DenseNeuroVol(out_array, space))
      }
      return(neuroim2::DenseNeuroVec(out_array, space))
    }
    return(neuroim2::DenseNeuroVol(out_array, space))
  }
  out_array
}

#' Load ANTs composite H5 as morphism or path
#'
#' Creates a morphism (or MorphismPath) from an ANTs composite H5 transform file.
#' These files typically contain both a displacement field and an affine transform.
#'
#' @section Pullback Ordering:
#' When `apply_affine=TRUE`, returns a MorphismPath ordered as `[warp, affine]`.
#' This ordering is critical for correct pullback (resampling) semantics:
#'
#' \itemize{
#'   \item MorphismPath applies transforms right-to-left for pullback
#'   \item `[warp, affine]` means: `affine_pullback(warp_pullback(coords))`
#'   \item Target coords -> warp lookup -> affine transform -> source coords
#' }
#'
#' The ANTs forward transform is "affine then warp", so pullback (inverse) is
#' "inverse_warp then inverse_affine" - hence the `[warp, affine]` order.
#'
#' @param path Path to ANTs H5 file
#' @param source Source coordinate space identifier
#' @param target Target coordinate space identifier
#' @param apply_affine If TRUE and an embedded affine is present, return a
#'   MorphismPath combining warp and affine. If FALSE, return only the warp.
#' @return Warp3DMorphism (if no affine or apply_affine=FALSE) or MorphismPath
#' @export
ants_h5_morphism <- function(path, source = "source", target = "target", apply_affine = TRUE) {
  warp_m <- Warp3DMorphism(source, target, warp_path = path, warp_type = "ants_h5")

  if (!apply_affine) {
    return(warp_m)
  }

  aff_mat <- NULL
  # Peek at embedded affine
  if (exists("load_warp_ants_h5", envir = asNamespace("neurotransform"))) {
    info <- load_warp_ants_h5(path)
    aff_mat <- info$affine
  }

  if (is.null(aff_mat)) {
    return(warp_m)
  }

  # For pullback semantics (target->source), the path is applied from last to first.
  # ANTs composite: forward is affine then warp, so pullback is inverse_warp then inverse_affine.
  # Since transform_path applies [f, g] as g_pullback(f_pullback(coords)),
  # we order as [warp, affine] so: affine_pullback(warp_pullback(coords))
  # This gives: coords -> warp lookup -> affine transform -> source coords
  #
  # The warp goes from source to an intermediate "warp_space" (its output coordinates).
  # The affine then goes from warp_space to target.
  # For pullback, we need: target -> warp_space -> source
  warp_space <- paste0(source, "_warp_space")
  warp_m <- Warp3DMorphism(source, warp_space, warp_path = path, warp_type = "ants_h5")
  aff_m <- Affine3DMorphism(warp_space, target, matrix = aff_mat)
  methods::new("MorphismPath",
               morphisms = list(warp_m, aff_m),
               source = source,
               target = target)
}
