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
#' @param type Optional type hint: ants_h5, ants, fsl, fsl_coef, afni, x5, linear,
#'   fsl_affine, afni_affine, itk_affine, lta_affine
#' @param source Optional source id
#' @param target Optional target id
#' @param apply_affine Logical; for ANTs H5 files, whether to include embedded affine
#' @param ... Passed to morphism constructors
#' @return A Morphism or MorphismPath object
#' @export
detect_transform_type <- function(path, source_affine = NULL, target_affine = NULL) {
  if (!is.character(path) || length(path) != 1L) stop("path must be a single file path")
  if (!file.exists(path)) stop("Transform file not found: ", path)

  lower <- tolower(path)
  base <- tolower(basename(path))

  if (grepl("\\.x5$", lower)) {
    return("x5")
  }
  if (grepl("\\.(h5|hdf5)$", lower)) {
    return(.detect_h5_transform_type(path))
  }
  if (grepl("\\.lta$", lower)) {
    return("lta_affine")
  }
  if (grepl("\\.aff12\\.1d$", lower) || grepl("\\.1d$", lower)) {
    return("afni_affine")
  }

  # Try content-based linear format detection for .mat/.txt files.
  if (grepl("\\.(mat|txt)$", lower)) {
    lfmt <- .detect_linear_format(path, source_affine = source_affine, target_affine = target_affine)
    return(switch(
      lfmt,
      "afni" = "afni_affine",
      "fsl" = "fsl_affine",
      "itk" = "itk_affine",
      "linear"
    ))
  }

  # NIfTI warps: use lightweight filename heuristics for better defaulting.
  if (grepl("\\.(nii|nii\\.gz)$", lower)) {
    if (grepl("coef|coeff|warpcoef|fieldcoef", lower)) {
      return("fsl_coef")
    }
    if (grepl("fnirt|_fsl_|/fsl/", lower)) {
      return("fsl")
    }
    if (grepl("qwarp|3dqwarp|_afni_|/afni/", lower)) {
      return("afni")
    }
    return("ants")
  }

  # Fallback: nonlinear ANTs-style warp path.
  if (grepl("coef|coeff", base)) return("fsl_coef")
  if (grepl("fnirt", base)) return("fsl")
  if (grepl("qwarp|3dqwarp", base)) return("afni")
  "ants"
}

.detect_h5_transform_type <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) return("ants_h5")
  h5 <- tryCatch(hdf5r::H5File$new(path, mode = "r"), error = function(e) NULL)
  if (is.null(h5)) return("ants_h5")
  on.exit(h5$close_all(), add = TRUE)

  attrs <- tryCatch(hdf5r::h5attributes(h5), error = function(e) NULL)
  if (!is.null(attrs) && "Format" %in% names(attrs) && identical(as.character(attrs$Format), "X5")) {
    return("x5")
  }

  tg <- tryCatch(h5[["TransformGroup"]], error = function(e) NULL)
  if (is.null(tg)) return("ants_h5")
  ttypes <- character(0)
  for (k in names(tg)) {
    g <- tg[[k]]
    if (!("TransformType" %in% names(g))) next
    tt <- tryCatch(g[["TransformType"]]$read(), error = function(e) NULL)
    if (is.null(tt)) next
    ttypes <- c(ttypes, as.character(tt))
  }
  if (!length(ttypes)) return("ants_h5")
  if (any(grepl("DisplacementFieldTransform", ttypes))) return("ants_h5")
  if (any(grepl("AffineTransform", ttypes))) return("itk_affine")
  "ants_h5"
}

#' @rdname read_transform
#' @export
read_transform <- function(path, type = NULL, source = NULL, target = NULL, apply_affine = TRUE, ...) {
  if (is(path, "Morphism")) return(path)
  if (!is.character(path) || length(path) != 1) stop("path must be a file path or Morphism")
  if (!file.exists(path)) stop("Transform file not found: ", path)
  extra <- list(...)
  inferred_type <- is.null(type)

  if (inferred_type) {
    type <- detect_transform_type(
      path,
      source_affine = extra$source_affine %||% NULL,
      target_affine = extra$target_affine %||% NULL
    )
  }
  if (identical(type, "ants_h5")) {
    return(ants_h5_morphism(path, source = source %||% "source", target = target %||% "target",
                             apply_affine = apply_affine, ...))
  }
  if (identical(type, "x5")) {
    return(.transform_from_x5(
      read_x5(path),
      source = source %||% "source",
      target = target %||% "target"
    ))
  }
  if (type %in% c("linear", "affine", "generic_affine")) {
    return(read_linear_transform(
      path,
      format = "generic",
      source = source %||% "source",
      target = target %||% "target",
      ...
    ))
  }
  if (type %in% c("itk_affine", "itk")) {
    return(read_linear_transform(
      path,
      format = "itk",
      source = source %||% "source",
      target = target %||% "target",
      ...
    ))
  }
  if (type %in% c("fsl_affine")) {
    return(read_linear_transform(
      path,
      format = "fsl",
      source = source %||% "source",
      target = target %||% "target",
      ...
    ))
  }
  if (type %in% c("afni_affine")) {
    return(read_linear_transform(
      path,
      format = "afni",
      source = source %||% "source",
      target = target %||% "target",
      ...
    ))
  }
  if (type %in% c("lta_affine", "lta", "fs_affine", "fs")) {
    return(read_linear_transform(
      path,
      format = "lta",
      source = source %||% "source",
      target = target %||% "target",
      ...
    ))
  }
  if (type %in% c("ants", "fsl", "fsl_coef", "afni")) {
    def_type <- extra$def_type %||% NULL
    if (isTRUE(inferred_type) && identical(type, "fsl") && is.null(def_type)) {
      def_type <- tryCatch(
        detect_fnirt_def_type(path),
        error = function(e) "relative"
      )
    }
    return(Warp3DMorphism(
      source %||% "source",
      target %||% "target",
      warp_path = path,
      warp_type = type,
      def_type = def_type %||% "relative"
    ))
  }
  stop("Unsupported transform type: ", type)
}

#' Read an affine transform file into an Affine3DMorphism
#'
#' @param path Path to transform file
#' @param format One of "generic", "fsl", "afni", or "itk"
#' @param source Source domain id
#' @param target Target domain id
#' @param source_affine Source voxel-to-world affine (required for fsl)
#' @param target_affine Target voxel-to-world affine (required for fsl)
#' @param source_dim Optional source image dimensions for FSL handedness-aware conversion
#' @param target_dim Optional target image dimensions for FSL handedness-aware conversion
#' @param oblique_correction Logical; for AFNI format, apply cardinal/deoblique correction
#' @param invert Logical; if TRUE, invert the read matrix before wrapping
#' @return Affine3DMorphism
#' @export
read_linear_transform <- function(path,
                                  format = c("generic", "fsl", "afni", "itk", "lta", "x5"),
                                  source = "source",
                                  target = "target",
                                  source_affine = NULL,
                                  target_affine = NULL,
                                  source_dim = NULL,
                                  target_dim = NULL,
                                  oblique_correction = TRUE,
                                  invert = FALSE) {
  if (!is.character(path) || length(path) != 1L) stop("path must be a single file path")
  if (!file.exists(path)) stop("Linear transform file not found: ", path)
  format <- match.arg(format)

  mat <- switch(
    format,
    generic = read_affine_matrix_txt(path)$matrix,
    afni = afni_aff12_to_ras(
      afni_read_aff12(path),
      source_affine = source_affine,
      target_affine = target_affine,
      oblique_correction = oblique_correction
    ),
    itk = {
      itk_lps <- .read_itk_affine(path)
      flip <- diag(c(-1, -1, 1, 1))
      itk_ras <- flip %*% itk_lps %*% flip
      # ITK stores forward transforms. Internal representation is pullback.
      solve(itk_ras)
    },
    lta = {
      lta <- .read_lta(path)
      .lta_matrix_to_internal(
        lta$xforms[[1]]$m_L,
        type = lta$type,
        src_vg = lta$xforms[[1]]$src,
        dst_vg = lta$xforms[[1]]$dst
      )
    },
    x5 = {
      nodes <- read_x5(path)
      linear_nodes <- Filter(function(n) identical(n$type, "linear"), nodes)
      if (!length(linear_nodes)) {
        .stop_transform_file("No linear transform found in X5 file.")
      }
      if (length(linear_nodes) > 1L) {
        .stop_transform_io("X5 file contains multiple linear transforms; use read_linear_transform_array(format='x5').")
      }
      as.matrix(linear_nodes[[1]]$transform)
    },
    fsl = {
      if (is.null(source_affine) || is.null(target_affine)) {
        stop("source_affine and target_affine are required for format='fsl'")
      }
      flirt <- read_affine_matrix_txt(path)$matrix
      fsl_flirt_to_internal_affine(
        flirt,
        source_affine = source_affine,
        ref_affine = target_affine,
        source_dim = source_dim,
        ref_dim = target_dim
      )
    }
  )
  if (isTRUE(invert)) mat <- invert_affine(mat)
  Affine3DMorphism(source = source, target = target, matrix = mat)
}

#' Write an affine transform to disk
#'
#' @param x Affine3DMorphism or 4x4 matrix (internal pullback convention)
#' @param path Output file path
#' @param format One of "generic", "fsl", "afni", or "itk"
#' @param source_affine Source voxel-to-world affine (required for fsl)
#' @param target_affine Target voxel-to-world affine (required for fsl)
#' @param source_dim Optional source image dimensions for FSL handedness-aware conversion
#' @param target_dim Optional target image dimensions for FSL handedness-aware conversion
#' @param oblique_correction Logical; for AFNI format, apply cardinal/deoblique correction
#' @param invert Logical; if TRUE, invert before writing
#' @return Invisibly returns \code{path}
#' @export
write_linear_transform <- function(x,
                                   path,
                                   format = c("generic", "fsl", "afni", "itk", "lta", "x5"),
                                   source_affine = NULL,
                                   target_affine = NULL,
                                   source_dim = NULL,
                                   target_dim = NULL,
                                   oblique_correction = TRUE,
                                   invert = FALSE) {
  format <- match.arg(format)
  mat <- if (is(x, "Affine3DMorphism")) x@matrix else x
  validate_4x4_matrix(mat, "x")
  if (isTRUE(invert)) mat <- invert_affine(mat)

  if (identical(format, "generic")) {
    write_affine_matrix_txt(mat, path)
    return(invisible(path))
  }

  if (identical(format, "fsl")) {
    if (is.null(source_affine) || is.null(target_affine)) {
      stop("source_affine and target_affine are required for format='fsl'")
    }
    flirt <- convert_affine_convention(
      mat,
      source_affine = source_affine,
      target_affine = target_affine,
      source_dim = source_dim,
      target_dim = target_dim,
      from = "generic",
      to = "fsl"
    )
    utils::write.table(format(flirt, digits = 10), file = path,
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
    return(invisible(path))
  }

  if (identical(format, "itk")) {
    # Internal is pullback (target->source); ITK file expects forward (source->target).
    forward_ras <- solve(mat)
    flip <- diag(c(-1, -1, 1, 1))
    forward_lps <- flip %*% forward_ras %*% flip
    .write_itk_affine(forward_lps, path)
    return(invisible(path))
  }

  if (identical(format, "lta")) {
    .write_lta(path, mats_internal = list(mat), type = 1L)
    return(invisible(path))
  }

  if (identical(format, "x5")) {
    node <- x5_transform(
      type = "linear",
      transform = mat,
      dimension_kinds = c("space", "space")
    )
    write_x5(path, list(node))
    return(invisible(path))
  }

  # afni
  mat_rai <- afni_ras_to_aff12(
    mat,
    source_affine = source_affine,
    target_affine = target_affine,
    oblique_correction = oblique_correction
  )
  utils::write.table(format(mat_rai[1:3, , drop = FALSE], digits = 10), file = path,
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  invisible(path)
}

#' Create an in-memory dense-field warp morphism
#'
#' Builds a \code{Warp3DMorphism} directly from a 4D field array and a grid,
#' without requiring an on-disk warp file.
#'
#' @param source Source domain id
#' @param target Target domain id
#' @param field Numeric array with dims \code{c(X, Y, Z, 3)}
#' @param grid Grid describing voxel geometry. If \code{NULL}, uses
#'   \code{attr(field, "grid")} when available (e.g., from \code{deformation_field()}).
#' @param representation Either \code{"displacements"} (relative offsets) or
#'   \code{"deformations"} (absolute source coordinates)
#' @param warp_method Interpolation method for warp lookup
#' @param cost Path cost
#' @param method_tag Method tag
#' @param id Optional stable id for the in-memory warp payload
#' @return Warp3DMorphism
#' @export
warp_from_field <- function(source, target, field, grid = NULL,
                            representation = c("displacements", "deformations"),
                            warp_method = c("linear", "cubic"),
                            cost = 1.5, method_tag = "anatomical",
                            id = NULL) {
  representation <- match.arg(representation)
  warp_method <- match.arg(warp_method)

  if (!is.character(source) || length(source) != 1L || !nzchar(source)) {
    stop("source must be a non-empty character string")
  }
  if (!is.character(target) || length(target) != 1L || !nzchar(target)) {
    stop("target must be a non-empty character string")
  }

  if (!is.array(field) || !is.numeric(field) || length(dim(field)) != 4L || dim(field)[4] != 3L) {
    stop("field must be a numeric array with dims c(X, Y, Z, 3)")
  }
  if (!all(is.finite(field))) stop("field must contain finite values only")

  grid <- grid %||% attr(field, "grid")
  if (!inherits(grid, "Grid")) {
    stop("grid must be a Grid object (or present as attr(field, 'grid'))")
  }
  if (!identical(as.integer(dim(field)[1:3]), as.integer(grid@dims))) {
    stop("field spatial dimensions must match grid@dims")
  }

  dims <- as.integer(grid@dims)
  arr <- .flatten_warp_components(field)

  payload <- list(
    array = arr,
    dim = dims,
    world_to_vox = solve(grid@affine),
    vox_to_world = grid@affine
  )

  key <- id %||% paste0("inline_", compute_hash(source, target, dims, grid@affine, runif(1)))
  def_type <- if (identical(representation, "deformations")) "absolute" else "relative"

  m <- Warp3DMorphism(
    source = source,
    target = target,
    warp_path = key,
    warp_type = "dense",
    def_type = def_type,
    warp_method = warp_method,
    cost = cost,
    method_tag = method_tag
  )

  assign(key, payload, envir = m@cache)
  m@params$representation <- representation
  m@params$inline <- TRUE
  m@hash <- morphism_hash(m)
  m
}

#' Write a warp morphism as a NIfTI vector field
#'
#' @param x Warp3DMorphism
#' @param path Output NIfTI path
#' @param representation Output representation: \code{"displacements"},
#'   \code{"deformations"}, or \code{"auto"} (use morphism native representation)
#' @return Invisibly returns \code{path}
#' @export
write_warp_field <- function(x, path, representation = c("auto", "displacements", "deformations")) {
  if (!is(x, "Warp3DMorphism")) stop("x must be a Warp3DMorphism")
  representation <- match.arg(representation)

  warp <- load_warp_array(x)
  native <- if (identical(x@params$def_type %||% "relative", "absolute")) {
    "deformations"
  } else {
    "displacements"
  }
  target_repr <- if (identical(representation, "auto")) native else representation

  field <- .unflatten_warp_components(warp$array, warp$dim)
  if (!identical(native, target_repr)) {
    world <- .grid_world_coords_matrix(warp$dim, warp$vox_to_world)
    if (identical(native, "displacements") && identical(target_repr, "deformations")) {
      field <- field + array(world, dim = dim(field))
    } else if (identical(native, "deformations") && identical(target_repr, "displacements")) {
      field <- field - array(world, dim = dim(field))
    }
  }

  # Pass spacing/origin explicitly for robust qform/sform consistency.
  D <- min(length(warp$dim), 3L)
  spacing <- sqrt(colSums(warp$vox_to_world[1:3, 1:3, drop = FALSE]^2))
  origin <- warp$vox_to_world[1:3, 4]
  space <- neuroim2::NeuroSpace(
    c(warp$dim, 3L),
    spacing = spacing[seq_len(D)],
    origin = origin[seq_len(D)],
    trans = warp$vox_to_world
  )
  vec <- neuroim2::DenseNeuroVec(field, space)
  neuroim2::write_vec(vec, path, format = "nifti")
  invisible(path)
}

#' Write a transform to disk
#'
#' Convenience wrapper around \code{write_linear_transform()} for affine morphisms
#' and raw 4x4 matrices. For \code{Warp3DMorphism}, writes a NIfTI vector field.
#'
#' @param x Affine3DMorphism, Warp3DMorphism, or 4x4 numeric matrix
#' @param path Output file path
#' @param type Optional output type:
#'   \itemize{
#'     \item linear: \code{generic}, \code{fsl}, \code{afni}, \code{itk}
#'     \item warp: \code{nifti}, \code{displacements}, \code{deformations}
#'   }
#' @param ... Passed to \code{write_linear_transform()} or \code{write_warp_field()}
#' @return Invisibly returns \code{path}
#' @export
write_transform <- function(x, path, type = NULL, ...) {
  fmt_hint <- tolower(type %||% "")
  lower_path <- tolower(path)
  if (!nzchar(fmt_hint) && grepl("\\.x5$", lower_path)) fmt_hint <- "x5"

  if (identical(fmt_hint, "x5")) {
    nodes <- .x5_nodes_from_transform(x)
    write_x5(path, nodes)
    return(invisible(path))
  }

  if (is(x, "Warp3DMorphism")) {
    wtype <- tolower(type %||% "auto")
    if (wtype %in% c("nifti", "auto")) {
      return(write_warp_field(x, path, representation = "auto"))
    }
    if (wtype %in% c("displacements", "deformations")) {
      return(write_warp_field(x, path, representation = wtype))
    }
    stop("Unsupported warp write type: ", wtype)
  }

  format <- tolower(type %||% "")
  if (!nzchar(format)) {
    format <- if (grepl("\\.aff12\\.1d$", lower_path)) {
      "afni"
    } else if (grepl("\\.lta$", lower_path)) {
      "lta"
    } else if (grepl("\\.mat$", lower_path) || grepl("\\.h5$", lower_path)) {
      "itk"
    } else {
      "generic"
    }
  }
  if (!format %in% c("generic", "fsl", "afni", "itk", "lta", "x5")) {
    stop("Unsupported write type: ", format)
  }
  if (is(x, "Affine3DMorphism") || is_affine_matrix(x)) {
    return(write_linear_transform(x, path = path, format = format, ...))
  }
  stop("write_transform currently supports Affine3DMorphism, Warp3DMorphism, or 4x4 matrices.")
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

# Determine linear affine format from file content and light context.
.detect_linear_format <- function(path, source_affine = NULL, target_affine = NULL) {
  lower <- tolower(path)
  lines <- tryCatch(readLines(path, warn = FALSE, n = 50L), error = function(e) character())
  first <- if (length(lines) > 0L) tolower(lines[1]) else ""
  text <- tolower(paste(lines, collapse = "\n"))

  if (startsWith(first, "#insight transform file")) return("itk")
  if (grepl("affinetransform_(double|float)_3_3", text)) return("itk")
  if (grepl("^\\s*transform\\s*:\\s*affinetransform", text)) return("itk")
  if (grepl("affinetype\\s*:\\s*(itk|ants|insight)", text)) return("itk")
  if (grepl("affinetype\\s*:\\s*afni", text)) return("afni")
  if (grepl("affinetype\\s*:\\s*(fsl|flirt)", text)) return("fsl")

  # Binary MATLAB v4 ANTs affine files include this variable name.
  raw_head <- tryCatch(readBin(path, "raw", n = 256L), error = function(e) raw())
  raw_int <- as.integer(raw_head)
  raw_int[raw_int < 32L | raw_int > 126L] <- 32L
  raw_txt <- tolower(rawToChar(as.raw(raw_int)))
  if (grepl("affinetransform_(double|float)_3_3", raw_txt)) return("itk")

  if (grepl("\\.aff12\\.1d$", lower) || grepl("\\.1d$", lower)) return("afni")

  # FSL .mat and generic 4x4 text are ambiguous; only pick FSL when caller
  # provides image affines (required for correct convention conversion).
  if (grepl("\\.mat$", lower) && !is.null(source_affine) && !is.null(target_affine)) {
    return("fsl")
  }

  "generic"
}

# Read ITK affine from either Insight text file or MATLAB v4 binary .mat.
.read_itk_affine <- function(path) {
  if (!file.exists(path)) stop("ITK affine file not found: ", path)
  if (grepl("\\.h5$", tolower(path))) {
    mats <- .read_itk_affine_h5(path)
    if (length(mats) != 1L) {
      .stop_transform_io(
        sprintf(
          "ITK H5 file contains %d linear transforms; use read_linear_transform_array(format='itk').",
          length(mats)
        )
      )
    }
    return(mats[[1]])
  }
  txt <- tryCatch(readLines(path, warn = FALSE, n = 80L), error = function(e) character())
  if (length(txt) > 0L && startsWith(tolower(txt[1]), "#insight transform file")) {
    return(.read_itk_affine_text(txt, path))
  }
  .read_itk_affine_mat_v4(path)
}

# Read all affine transforms from an ITK Composite H5 file.
.read_itk_affine_h5 <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    .stop_transform_io("hdf5r is required to read ITK H5 transforms.")
  }
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  tg <- h5[["TransformGroup"]]
  if (is.null(tg)) {
    .stop_transform_file("TransformGroup not found in ITK H5 file.")
  }

  mats <- list()
  keys <- names(tg)
  keys <- keys[order(suppressWarnings(as.integer(keys)), na.last = TRUE)]
  for (k in keys) {
    g <- tg[[k]]
    if (!"TransformType" %in% names(g)) next
    tt <- g[["TransformType"]]$read()
    tt <- if (length(tt)) as.character(tt[[1]]) else ""
    if (!grepl("AffineTransform", tt)) next
    if (!all(c("TransformParameters", "TransformFixedParameters") %in% names(g))) next

    p <- as.numeric(g[["TransformParameters"]]$read())
    fixed <- as.numeric(g[["TransformFixedParameters"]]$read())
    if (length(p) < 12L) next
    if (length(fixed) < 3L) fixed <- c(fixed, rep(0, 3L - length(fixed)))
    mats[[length(mats) + 1L]] <- .itk_params_to_affine(p, fixed)
  }

  if (!length(mats)) {
    .stop_transform_file("No affine transforms found in ITK H5 file.")
  }
  mats
}

.read_itk_affine_text <- function(lines, path) {
  params_line <- grep("^\\s*Parameters\\s*:", lines, value = TRUE, ignore.case = TRUE)
  if (length(params_line) == 0L) stop("ITK text affine missing Parameters line: ", path)
  params <- scan(text = sub("^[^:]*:\\s*", "", params_line[1]), quiet = TRUE)
  if (length(params) < 12L) stop("ITK text affine needs at least 12 parameters: ", path)

  fixed_line <- grep("^\\s*FixedParameters\\s*:", lines, value = TRUE, ignore.case = TRUE)
  fixed <- if (length(fixed_line) > 0L) {
    scan(text = sub("^[^:]*:\\s*", "", fixed_line[1]), quiet = TRUE)
  } else {
    c(0, 0, 0)
  }
  if (length(fixed) < 3L) fixed <- c(fixed, rep(0, 3L - length(fixed)))
  .itk_params_to_affine(params, fixed)
}

.read_itk_affine_mat_v4 <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con))

  vars <- list()
  repeat {
    hdr <- readBin(con, integer(), n = 5L, size = 4L, endian = "little")
    if (length(hdr) < 5L) break
    mopt <- hdr[1]
    mrows <- hdr[2]
    ncols <- hdr[3]
    imagf <- hdr[4]
    namelen <- hdr[5]
    if (mrows <= 0L || ncols <= 0L || namelen <= 0L) break

    name_raw <- readBin(con, "raw", n = namelen)
    if (length(name_raw) < namelen) stop("Truncated MATLAB v4 variable name in ", path)
    name <- paste(rawToChar(name_raw, multiple = TRUE), collapse = "")
    name <- sub("[[:cntrl:]]+$", "", name)

    precision_code <- (mopt %/% 10L) %% 10L
    n <- as.integer(mrows * ncols)
    vals <- switch(
      as.character(precision_code),
      "0" = readBin(con, numeric(), n = n, size = 8L, endian = "little"),
      "1" = readBin(con, numeric(), n = n, size = 4L, endian = "little"),
      "2" = as.numeric(readBin(con, integer(), n = n, size = 4L, signed = TRUE, endian = "little")),
      stop("Unsupported MATLAB v4 precision code ", precision_code, " in ", path)
    )
    if (length(vals) < n) stop("Truncated MATLAB v4 matrix in ", path)

    if (!identical(imagf, 0L)) {
      # Skip imaginary part if present.
      bytes <- if (precision_code == 0L) 8L else 4L
      readBin(con, "raw", n = n * bytes)
    }

    vars[[name]] <- matrix(vals, nrow = mrows, ncol = ncols, byrow = FALSE)
  }

  if (length(vars) == 0L) {
    stop("No MATLAB v4 variables found in ITK affine file: ", path)
  }
  keys <- names(vars)
  aff_key <- keys[grepl("affinetransform_(double|float)_3_3", tolower(keys))]
  if (length(aff_key) == 0L) {
    stop("Could not locate ITK affine parameters in MATLAB v4 file: ", path)
  }
  params <- as.numeric(vars[[aff_key[1]]])
  if (length(params) < 12L) {
    stop("ITK MATLAB v4 affine has insufficient parameters: ", path)
  }

  fixed_key <- keys[grepl("^fixed$", tolower(keys))]
  fixed <- if (length(fixed_key) > 0L) as.numeric(vars[[fixed_key[1]]]) else c(0, 0, 0)
  if (length(fixed) < 3L) fixed <- c(fixed, rep(0, 3L - length(fixed)))

  .itk_params_to_affine(params, fixed)
}

.itk_params_to_affine <- function(params, fixed) {
  A <- matrix(params[1:9], nrow = 3L, byrow = TRUE)
  tvec <- params[10:12]
  center <- fixed[1:3]
  M <- diag(4)
  M[1:3, 1:3] <- A
  M[1:3, 4] <- as.numeric(tvec + center - A %*% center)
  M
}

.write_itk_affine <- function(matrix_lps, path) {
  if (grepl("\\.h5$", tolower(path))) {
    .write_itk_affine_h5(list(matrix_lps), path)
    return(invisible(path))
  }
  if (!is_affine_matrix(matrix_lps)) stop("matrix_lps must be a numeric 4x4 affine matrix")
  A <- matrix_lps[1:3, 1:3, drop = FALSE]
  tvec <- matrix_lps[1:3, 4]
  params <- c(as.numeric(t(A)), as.numeric(tvec))
  lines <- c(
    "#Insight Transform File V1.0",
    "Transform: AffineTransform_double_3_3",
    paste("Parameters:", paste(format(params, digits = 10), collapse = " ")),
    "FixedParameters: 0 0 0"
  )
  writeLines(lines, con = path)
  invisible(path)
}

# Write one or more ITK affine transforms to Composite H5 format.
.write_itk_affine_h5 <- function(mats_lps, path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    .stop_transform_io("hdf5r is required to write ITK H5 transforms.")
  }
  if (!is.list(mats_lps)) mats_lps <- list(mats_lps)
  if (!length(mats_lps)) .stop_transform_io("No matrices provided for ITK H5 writing.")
  for (m in mats_lps) validate_4x4_matrix(m, "itk_lps_matrix")

  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all())
  tg <- h5$create_group("TransformGroup")

  # Composite sentinel entry mirrors ITK's conventional layout.
  g0 <- tg$create_group("0")
  g0$create_dataset("TransformType", robj = "CompositeTransform_double_3_3")

  for (i in seq_along(mats_lps)) {
    g <- tg$create_group(as.character(i))
    m <- mats_lps[[i]]
    A <- m[1:3, 1:3, drop = FALSE]
    tvec <- as.numeric(m[1:3, 4])
    params <- c(as.numeric(t(A)), tvec)
    fixed <- c(0, 0, 0)

    g$create_dataset("TransformType", robj = "AffineTransform_double_3_3")
    g$create_dataset("TransformParameters", robj = params)
    g$create_dataset("TransformFixedParameters", robj = fixed)
  }
  invisible(path)
}

.flatten_warp_components <- function(field) {
  dims <- dim(field)[1:3]
  nvox <- prod(dims)
  idx <- seq_len(nvox)
  out <- numeric(3L * nvox)
  out[3L * (idx - 1L) + 1L] <- as.numeric(field[, , , 1, drop = TRUE])
  out[3L * (idx - 1L) + 2L] <- as.numeric(field[, , , 2, drop = TRUE])
  out[3L * (idx - 1L) + 3L] <- as.numeric(field[, , , 3, drop = TRUE])
  out
}

.unflatten_warp_components <- function(array_interleaved, dims) {
  dims <- as.integer(dims)
  nvox <- prod(dims)
  if (length(array_interleaved) != 3L * nvox) {
    stop("Expected interleaved array length 3 * prod(dims)")
  }
  idx <- seq_len(nvox)
  out <- array(0, dim = c(dims, 3L))
  out[, , , 1] <- array(array_interleaved[3L * (idx - 1L) + 1L], dim = dims)
  out[, , , 2] <- array(array_interleaved[3L * (idx - 1L) + 2L], dim = dims)
  out[, , , 3] <- array(array_interleaved[3L * (idx - 1L) + 3L], dim = dims)
  out
}

.grid_world_coords_matrix <- function(dims, vox_to_world) {
  dims <- as.integer(dims)
  ijk <- as.matrix(expand.grid(
    seq.int(0L, dims[1] - 1L),
    seq.int(0L, dims[2] - 1L),
    seq.int(0L, dims[3] - 1L)
  ))
  world <- t(vox_to_world %*% t(cbind(ijk, 1)))
  world[, 1:3, drop = FALSE]
}
