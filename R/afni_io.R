#' AFNI .aff12.1D matrix I/O
#'
#' @description
#' Readers and writers for the matrix files produced by
#' \verb{3dvolreg -1Dmatrix_save} and \verb{3dAllineate -1Dmatrix_save}.
#'
#' @section File layout:
#' AFNI stores one affine per sub-brick as a single row of twelve row-major
#' numbers, \verb{u11 u12 u13 v1 u21 u22 u23 v2 u31 u32 u33 v3}, optionally
#' preceded by comment lines such as
#' \verb{# 3dvolreg matrices (DICOM-to-DICOM, row-by-row):}. A registration of a
#' four-volume series therefore yields four rows. The three-lines-of-four layout
#' accepted by \code{cat_matvec} FORM 1 is also read, and always denotes a single
#' matrix.
#'
#' @section Coordinates and direction:
#' The numbers are in AFNI RAI/DICOM order, whose axis letters name the negative
#' end of each axis, so numerically \verb{+x = Left}, \verb{+y = Posterior},
#' \verb{+z = Superior}. The stored mapping is base-to-source: \verb{Xsource =
#' M Xbase}. Use \code{\link{afni_aff12_to_ras}} to reach the package's RAS
#' pullback convention.
#'
#' @name afni_io
#' @keywords internal
NULL

.afni_aff12_banner <- "# AFNI affine matrices (DICOM-to-DICOM, row-by-row):"

# Split a .1D file into comment lines and numeric rows, preserving order.
.afni_read_1d <- function(path) {
  if (!is.character(path) || length(path) != 1L) {
    .stop_transform_io("path must be a single file path")
  }
  if (!file.exists(path)) {
    .stop_transform_file(paste0("AFNI aff12 file not found: ", path))
  }
  lines <- readLines(path, warn = FALSE)
  is_comment <- grepl("^\\s*#", lines)
  body <- lines[!is_comment]
  body <- body[nzchar(trimws(body))]
  list(comments = lines[is_comment], body = body)
}

# Parse the numeric body into a list of 3x4 row-major blocks.
.afni_parse_aff12_body <- function(body, path) {
  if (!length(body)) {
    .stop_transform_file(paste0("AFNI aff12 file contains no matrix rows: ", path))
  }
  counts <- vapply(body, function(line) {
    length(suppressWarnings(scan(text = line, quiet = TRUE, what = numeric())))
  }, integer(1L), USE.NAMES = FALSE)

  if (all(counts == 12L)) {
    vals <- lapply(body, function(line) scan(text = line, quiet = TRUE))
  } else if (length(body) == 3L && all(counts == 4L)) {
    # cat_matvec FORM 1: one matrix written as three lines of four numbers.
    vals <- list(scan(text = paste(body, collapse = " "), quiet = TRUE))
  } else {
    .stop_transform_file(sprintf(
      paste0(
        "AFNI aff12 file has an unsupported layout: %s. Expected one row of 12 ",
        "numbers per matrix, or three rows of 4 numbers for a single matrix; ",
        "found %d row(s) with %s number(s)."
      ),
      path, length(body), paste(unique(counts), collapse = "/")
    ))
  }

  lapply(seq_along(vals), function(i) {
    v <- vals[[i]]
    if (!all(is.finite(v))) {
      .stop_transform_file(sprintf(
        "AFNI aff12 row %d of %s contains a non-finite value.", i, path
      ))
    }
    rbind(matrix(v, nrow = 3L, ncol = 4L, byrow = TRUE), c(0, 0, 0, 1))
  })
}

# AFNI emits degenerate rows when a registration collapses; they are not
# invertible mappings and must not reach a morphism.
.afni_require_nonsingular <- function(mat, index, path) {
  if (abs(det(mat[1:3, 1:3, drop = FALSE])) < .Machine$double.eps^0.5) {
    .stop_transform_file(sprintf(
      "AFNI aff12 row %d of %s is singular and cannot be inverted.", index, path
    ))
  }
  invisible(mat)
}

#' Read every affine row of an AFNI .aff12.1D file
#'
#' Reads all matrix rows, returning one 4x4 matrix per row in AFNI RAI/DICOM
#' coordinates. Comment lines are returned in the \code{"comments"} attribute.
#'
#' @param path Path to a \code{.aff12.1D} (or \code{.1D}) matrix file
#' @return A list of 4x4 matrices in RAI, with a \code{"comments"} attribute
#'   holding the file's comment lines
#' @seealso \code{\link{afni_read_aff12}} for a single-matrix file
#' @export
#' @examples
#' path <- tempfile(fileext = ".aff12.1D")
#' afni_write_aff12(list(diag(4), diag(4)), path)
#' length(afni_read_aff12_array(path))
afni_read_aff12_array <- function(path) {
  parsed <- .afni_read_1d(path)
  mats <- .afni_parse_aff12_body(parsed$body, path)
  attr(mats, "comments") <- parsed$comments
  mats
}

#' Read AFNI .aff12.1D affine matrix
#'
#' Reads a single affine from a \code{.aff12.1D} file and returns it as a 4x4
#' matrix in AFNI RAI/DICOM coordinates. AFNI axis labels describe the negative
#' end of each axis, so numeric RAI coordinates increase toward Left, Posterior,
#' and Superior (LPS).
#'
#' @param path Path to \code{.aff12.1D} file
#' @param row Optional 1-based row to select from a multi-row file. Required
#'   when the file holds more than one matrix.
#' @return 4x4 affine matrix (in RAI)
#' @seealso \code{\link{afni_read_aff12_array}} for per-volume matrix series
#' @export
#' @examples
#' \dontrun{
#' mat_rai <- afni_read_aff12("transform.aff12.1D")
#' }
afni_read_aff12 <- function(path, row = NULL) {
  mats <- afni_read_aff12_array(path)
  if (is.null(row)) {
    if (length(mats) != 1L) {
      .stop_transform_file(sprintf(
        paste0(
          "AFNI aff12 file %s holds %d matrices. Select one with `row=`, or ",
          "read the whole series with afni_read_aff12_array()."
        ),
        path, length(mats)
      ))
    }
    return(mats[[1L]])
  }
  if (!is.numeric(row) || length(row) != 1L || is.na(row) ||
      row < 1 || row > length(mats) || row != as.integer(row)) {
    .stop_transform_io(sprintf(
      "row must be a single index between 1 and %d.", length(mats)
    ))
  }
  mats[[as.integer(row)]]
}

#' Write AFNI .aff12.1D affine matrices
#'
#' Writes matrices in AFNI's canonical layout: one row of twelve row-major
#' numbers per matrix, optionally preceded by a banner comment. This is the
#' layout \verb{3dAllineate -1Dmatrix_apply} expects.
#'
#' @param x A 4x4 (or 3x4) matrix in RAI, or a list of them
#' @param path Output path
#' @param comment Logical, or a character vector of comment lines to write as a
#'   header. \code{TRUE} writes the package's banner; \code{FALSE} writes none.
#' @return Invisibly, \code{path}
#' @export
#' @examples
#' path <- tempfile(fileext = ".aff12.1D")
#' afni_write_aff12(diag(4), path)
afni_write_aff12 <- function(x, path, comment = TRUE) {
  mats <- if (is.matrix(x)) list(x) else x
  if (!is.list(mats) || !length(mats)) {
    .stop_transform_io("x must be a 4x4 matrix or a non-empty list of matrices")
  }
  rows <- vapply(seq_along(mats), function(i) {
    m <- mats[[i]]
    if (is(m, "Affine3DMorphism")) m <- m@matrix
    if (!is.matrix(m) || !is.numeric(m) || ncol(m) != 4L || !nrow(m) %in% c(3L, 4L)) {
      .stop_transform_io(sprintf(
        "x[[%d]] must be a numeric 3x4 or 4x4 matrix", i
      ))
    }
    if (!all(is.finite(m[1:3, ]))) {
      .stop_transform_io(sprintf("x[[%d]] contains a non-finite value", i))
    }
    paste(formatC(as.vector(t(m[1:3, , drop = FALSE])),
                  format = "g", digits = 12L, width = 18L),
          collapse = "")
  }, character(1L))

  header <- if (isTRUE(comment)) {
    .afni_aff12_banner
  } else if (is.character(comment)) {
    ifelse(grepl("^\\s*#", comment), comment, paste("#", comment))
  } else {
    character(0)
  }
  writeLines(c(header, rows), con = path)
  invisible(path)
}

# Read one AFNI matrix file as a RAS pullback matrix, rejecting degenerate rows.
.afni_read_linear_matrix <- function(path, source_affine = NULL,
                                     target_affine = NULL,
                                     oblique_correction = TRUE, row = NULL) {
  mat_rai <- afni_read_aff12(path, row = row)
  .afni_require_nonsingular(mat_rai, if (is.null(row)) 1L else as.integer(row), path)
  afni_aff12_to_ras(
    mat_rai,
    source_affine = source_affine,
    target_affine = target_affine,
    oblique_correction = oblique_correction
  )
}

# Keep only the arguments the RAI/RAS conversion itself understands, so callers
# can pass the wider read/write argument set through `...`.
.afni_conversion_args <- function(args) {
  keep <- c("source_affine", "target_affine", "oblique_correction")
  args[intersect(names(args), keep)]
}
