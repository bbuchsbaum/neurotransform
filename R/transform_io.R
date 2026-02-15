#' @title Transform IO Helpers
#' @name transform_io
#' @description
#' Lightweight IO helpers for linear transform backends, including format
#' factories, typed IO conditions, and array-style linear transform IO.
NULL

#' Construct a transform IO error condition
#'
#' @param message Error message
#' @return Condition object of class `TransformIOError`
#' @export
TransformIOError <- function(message) {
  structure(
    list(message = as.character(message), call = NULL),
    class = c("TransformIOError", "error", "condition")
  )
}

#' Construct a transform file-format error condition
#'
#' @param message Error message
#' @return Condition object of class `TransformFileError`
#' @export
TransformFileError <- function(message) {
  structure(
    list(message = as.character(message), call = NULL),
    class = c("TransformFileError", "TransformIOError", "error", "condition")
  )
}

.stop_transform_io <- function(message) stop(TransformIOError(message))
.stop_transform_file <- function(message) stop(TransformFileError(message))

.canonical_linear_format <- function(fmt) {
  key <- tolower(fmt %||% "")
  if (!nzchar(key)) .stop_transform_io("Format must be a non-empty string.")
  aliases <- c(
    itk = "itk",
    ants = "itk",
    elastix = "itk",
    fsl = "fsl",
    afni = "afni",
    lta = "lta",
    fs = "lta",
    x5 = "x5",
    generic = "generic"
  )
  out <- unname(aliases[key])
  if (is.na(out) || !nzchar(out)) {
    .stop_transform_io(sprintf("Unsupported transform format <%s>.", fmt))
  }
  out
}

#' Get linear transform IO backend factory
#'
#' Returns a small factory descriptor for linear transform IO, including
#' canonical format and reader/writer function handles.
#'
#' @param fmt Format name (e.g. `"itk"`, `"ants"`, `"fsl"`, `"afni"`, `"generic"`)
#' @param is_array Logical; if `TRUE`, returns array-IO handlers
#' @return List with `format`, `reader`, and `writer`
#' @export
get_linear_factory <- function(fmt, is_array = TRUE) {
  format <- .canonical_linear_format(fmt)
  if (isTRUE(is_array)) {
    return(list(
      format = format,
      reader = read_linear_transform_array,
      writer = write_linear_transform_array
    ))
  }

  list(
    format = format,
    reader = read_linear_transform,
    writer = write_linear_transform
  )
}

#' Read a series of linear transforms
#'
#' Supports single-file transforms and FSL/MCFLIRT-style indexed files
#' (`path.000`, `path.001`, ...).
#'
#' @param path Base file path
#' @param format One of `"generic"`, `"fsl"`, `"afni"`, or `"itk"`
#' @param source Source domain id
#' @param target Target domain id
#' @param ... Passed to `read_linear_transform()`
#' @return A `LinearTransformArray` object
#' @export
read_linear_transform_array <- function(path,
                                        format = c("generic", "fsl", "afni", "itk", "lta", "x5"),
                                        source = "source",
                                        target = "target",
                                        ...) {
  format <- .canonical_linear_format(match.arg(format))
  extra <- list(...)

  if (identical(format, "itk") && grepl("\\.h5$", tolower(path))) {
    mats_lps <- .read_itk_affine_h5(path)
    flip <- diag(c(-1, -1, 1, 1))
    tx <- lapply(mats_lps, function(itk_lps) {
      itk_ras <- flip %*% itk_lps %*% flip
      Affine3DMorphism(source = source, target = target, matrix = solve(itk_ras))
    })
    return(structure(list(format = format, transforms = tx, paths = path), class = "LinearTransformArray"))
  }

  if (identical(format, "lta")) {
    lta <- .read_lta(path)
    tx <- lapply(lta$xforms, function(xf) {
      mat <- .lta_matrix_to_internal(xf$m_L, lta$type, src_vg = xf$src, dst_vg = xf$dst)
      Affine3DMorphism(source = source, target = target, matrix = mat)
    })
    return(structure(list(format = format, transforms = tx, paths = path), class = "LinearTransformArray"))
  }

  if (identical(format, "x5")) {
    nodes <- read_x5(path)
    linear_nodes <- Filter(function(n) identical(n$type, "linear"), nodes)
    if (!length(linear_nodes)) .stop_transform_file("No linear transforms found in X5 file.")
    tx <- lapply(linear_nodes, function(node) {
      Affine3DMorphism(source = source, target = target, matrix = as.matrix(node$transform))
    })
    return(structure(list(format = format, transforms = tx, paths = path), class = "LinearTransformArray"))
  }

  paths <- character(0)
  if (file.exists(path)) {
    paths <- c(paths, path)
  } else if (identical(format, "fsl")) {
    idx <- 0L
    repeat {
      p <- sprintf("%s.%03d", path, idx)
      if (!file.exists(p)) break
      paths <- c(paths, p)
      idx <- idx + 1L
    }
  }

  if (!length(paths)) {
    .stop_transform_file(sprintf("Linear transform file(s) not found for base path: %s", path))
  }

  tx <- lapply(paths, function(p) {
    do.call(read_linear_transform, c(list(
      path = p,
      format = format,
      source = source,
      target = target
    ), extra))
  })

  structure(
    list(format = format, transforms = tx, paths = paths),
    class = "LinearTransformArray"
  )
}

#' Write a series of linear transforms
#'
#' @param x A `LinearTransformArray`, list of affine morphisms/matrices, or single transform
#' @param path Base output path
#' @param format One of `"generic"`, `"fsl"`, `"afni"`, or `"itk"`
#' @param ... Passed to `write_linear_transform()`
#' @return Invisibly returns written path(s)
#' @export
write_linear_transform_array <- function(x,
                                         path,
                                         format = c("generic", "fsl", "afni", "itk", "lta", "x5"),
                                         ...) {
  format <- .canonical_linear_format(match.arg(format))

  tx <- if (inherits(x, "LinearTransformArray")) {
    x$transforms
  } else if (is.list(x)) {
    x
  } else {
    list(x)
  }

  if (!length(tx)) .stop_transform_io("No transforms provided.")

  if (length(tx) == 1L) {
    write_linear_transform(tx[[1]], path = path, format = format, ...)
    return(invisible(path))
  }

  mats <- lapply(tx, function(txi) {
    if (is(txi, "Affine3DMorphism")) txi@matrix else txi
  })
  for (m in mats) validate_4x4_matrix(m, "linear_transform_array_matrix")

  if (identical(format, "itk")) {
    # Composite H5 writer for linear arrays.
    if (!grepl("\\.h5$", tolower(path))) {
      .stop_transform_io("ITK linear transform arrays must be written to a .h5 path.")
    }
    flip <- diag(c(-1, -1, 1, 1))
    mats_lps <- lapply(mats, function(mat) {
      forward_ras <- solve(mat)
      flip %*% forward_ras %*% flip
    })
    .write_itk_affine_h5(mats_lps, path)
    return(invisible(path))
  }

  if (identical(format, "lta")) {
    .write_lta(path, mats_internal = mats, type = 1L)
    return(invisible(path))
  }

  if (identical(format, "x5")) {
    nodes <- lapply(mats, function(mat) {
      x5_transform(type = "linear", transform = mat, dimension_kinds = c("space", "space"))
    })
    write_x5(path, nodes)
    return(invisible(path))
  }

  if (!identical(format, "fsl")) {
    .stop_transform_io("Multi-transform array writing is currently supported only for format='fsl'.")
  }

  out_paths <- character(length(tx))
  for (i in seq_along(tx)) {
    p <- sprintf("%s.%03d", path, i - 1L)
    write_linear_transform(tx[[i]], path = p, format = format, ...)
    out_paths[i] <- p
  }
  invisible(out_paths)
}

#' @export
print.LinearTransformArray <- function(x, ...) {
  cat(
    "<LinearTransformArray | format=", x$format,
    " | n=", length(x$transforms),
    ">\n",
    sep = ""
  )
}
