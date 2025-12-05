#' @title Affine utilities (matrix-level and morphism helpers)
#' @name affine_utils
#' @description
#' Lightweight helpers for building, decomposing, and testing affine transforms
#' without introducing external dependencies or tool-specific semantics.
NULL

#' Check if an object is an Affine3DMorphism
#' @param x Object to test
#' @return Logical
#' @export
is_affine_morphism <- function(x) inherits(x, "Affine3DMorphism")

#' Check if an object is a numeric 4x4 affine matrix
#' @param x Object to test
#' @return Logical
#' @export
is_affine_matrix <- function(x) {
  is.matrix(x) && all(dim(x) == c(4L, 4L)) && is.numeric(x) && all(is.finite(x))
}

#' Build a 4x4 affine matrix from components
#'
#' @param translation length-3 translation vector (applied last)
#' @param scales length-3 scaling
#' @param skews length-3 shear terms (xy, xz, yz)
#' @param angles length-3 rotation angles (radians) about x, y, z (applied in Z-Y-X order)
#' @param anchor Either "none"/"origin"/"centre"/"center" or a numeric length-3 point
#' @return 4x4 numeric matrix
#' @export
build_affine_matrix <- function(
    translation = c(0, 0, 0),
    scales = c(1, 1, 1),
    skews = c(0, 0, 0),
    angles = c(0, 0, 0),
    anchor = c("none", "origin", "centre", "center")) {
  stopifnot(length(translation) == 3, length(scales) == 3,
            length(skews) == 3, length(angles) == 3)
  if (is.character(anchor)) {
    anchor <- match.arg(anchor)
    if (anchor %in% c("none", "origin")) {
      anchor_vec <- c(0, 0, 0)
    } else {
      anchor_vec <- c(0, 0, 0)
    }
  } else {
    if (length(anchor) != 3) stop("anchor must be length 3 or a known keyword")
    anchor_vec <- anchor
  }

  # Rotations (Z-Y-X / yaw-pitch-roll)
  cx <- cos(angles[1]); sx <- sin(angles[1])
  cy <- cos(angles[2]); sy <- sin(angles[2])
  cz <- cos(angles[3]); sz <- sin(angles[3])

  Rx <- matrix(c(
    1, 0, 0,
    0, cx, -sx,
    0, sx, cx
  ), nrow = 3, byrow = TRUE)
  Ry <- matrix(c(
    cy, 0, sy,
    0, 1, 0,
    -sy, 0, cy
  ), nrow = 3, byrow = TRUE)
  Rz <- matrix(c(
    cz, -sz, 0,
    sz, cz, 0,
    0, 0, 1
  ), nrow = 3, byrow = TRUE)
  R <- Rz %*% Ry %*% Rx

  # Shear matrix (xy, xz, yz)
  Sh <- matrix(c(
    1, skews[1], skews[2],
    0, 1,        skews[3],
    0, 0,        1
  ), nrow = 3, byrow = TRUE)

  S <- diag(scales)
  linear <- R %*% Sh %*% S

  # Assemble anchored transform
  M_lin <- diag(4)
  M_lin[1:3, 1:3] <- linear

  T_anchor <- diag(4); T_anchor[1:3, 4] <- anchor_vec
  T_anchor_inv <- diag(4); T_anchor_inv[1:3, 4] <- -anchor_vec
  T_trans <- diag(4); T_trans[1:3, 4] <- translation

  anchored <- T_anchor %*% M_lin %*% T_anchor_inv
  T_trans %*% anchored
}

#' Decompose a 4x4 affine matrix into components
#'
#' @param matrix 4x4 numeric matrix
#' @return list with translation, rotation_matrix, scales, skews, angles
#' @export
decompose_affine_matrix <- function(matrix) {
  if (!is_affine_matrix(matrix)) stop("Input must be a numeric 4x4 affine matrix")
  A <- matrix[1:3, 1:3]
  tvec <- matrix[1:3, 4]

  # Polar decomposition to get closest rotation
  sv <- svd(A)
  R <- sv$u %*% t(sv$v)

  # Handle reflection
  if (det(R) < 0) {
    sv$v[, 3] <- -sv$v[, 3]
    sv$d[3] <- -sv$d[3]
    R <- sv$u %*% t(sv$v)
  }

  C <- t(R) %*% A
  scales <- diag(C)
  scales[abs(scales) < .Machine$double.eps] <- .Machine$double.eps

  # Shears relative to scaled axes
  skews <- c(
    C[1, 2] / scales[2],
    C[1, 3] / scales[3],
    C[2, 3] / scales[3]
  )

  # Extract Z-Y-X angles from R
  angles <- c(
    atan2(R[3, 2], R[3, 3]),           # x-axis rotation
    asin(-R[3, 1]),                    # y-axis rotation
    atan2(R[2, 1], R[1, 1])            # z-axis rotation
  )

  list(
    translation = as.numeric(tvec),
    rotation_matrix = R,
    scales = as.numeric(scales),
    skews = as.numeric(skews),
    angles = as.numeric(angles)
  )
}

#' Wrap an affine matrix built from components into an Affine3DMorphism
#'
#' @param source Source domain identifier
#' @param target Target domain identifier
#' @param ... Arguments passed to \code{\link{build_affine_matrix}}
#' @return Affine3DMorphism object
#' @export
affine_from_components <- function(source, target, ...) {
  mat <- build_affine_matrix(...)
  Affine3DMorphism(source = source, target = target, matrix = mat)
}

#' Read a 4x4 affine matrix from a text file
#'
#' Comment lines (starting with "#") are ignored. If a line contains
#' "affineType:" the value is returned as the detected type.
#'
#' @param file Path to text file containing affine matrix
#' @param type Optional type hint (overrides auto-detection)
#' @return list(matrix = 4x4, type = character|NULL)
#' @export
read_affine_matrix_txt <- function(file, type = NULL) {
  lines <- readLines(file, warn = FALSE)
  inferred <- NULL
  for (ln in lines) {
    if (grepl("affineType", ln, ignore.case = TRUE)) {
      val <- sub(".*affineType:\\s*", "", ln, ignore.case = TRUE)
      inferred <- trimws(val)
    }
  }
  dat <- read.table(text = paste(lines, collapse = "\n"),
                    comment.char = "#", colClasses = "numeric")
  vals <- as.numeric(t(dat))
  if (length(vals) < 16) stop("Expected at least 16 numeric values for a 4x4 matrix")
  mat <- matrix(vals[1:16], nrow = 4, byrow = TRUE)
  list(matrix = mat, type = type %||% inferred)
}

#' Write a 4x4 affine matrix to a text file
#'
#' @param matrix 4x4 numeric affine matrix
#' @param file Output file path
#' @param type Optional type label to write as comment
#' @param comment Logical; if TRUE, write type as comment header
#' @return Invisibly returns the file path
#' @export
write_affine_matrix_txt <- function(matrix, file, type = NULL, comment = TRUE) {
  if (!is_affine_matrix(matrix)) stop("matrix must be 4x4 numeric")
  con <- file(file, open = "w")
  on.exit(close(con))
  if (!is.null(type) && isTRUE(comment)) {
    writeLines(sprintf("# affineType: %s", type), con = con)
  }
  write.table(format(matrix, digits = 10), file = con,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  invisible(file)
}

#' Convert affine conventions between "generic" (target→source world) and "fsl" (source→target vox)
#'
#' @param matrix 4x4 matrix in the `from` convention
#' @param source_affine 4x4 voxel-to-world for source image
#' @param target_affine 4x4 voxel-to-world for target image
#' @param from Input convention ("generic", "fsl")
#' @param to Output convention ("generic", "fsl")
#' @export
convert_affine_convention <- function(matrix,
                                      source_affine,
                                      target_affine,
                                      from = c("generic", "fsl"),
                                      to = c("generic", "fsl")) {
  from <- match.arg(from)
  to <- match.arg(to)
  if (from == to) return(matrix)

  if (from == "fsl" && to == "generic") {
    # FSL FLIRT: target_vox = F * source_vox
    # Canonical: target_world -> source_world
    W_s <- source_affine
    W_t <- target_affine
    F <- matrix
    return(W_s %*% solve(F) %*% solve(W_t))
  }
  if (from == "generic" && to == "fsl") {
    # Invert the above relationship
    W_s <- source_affine
    W_t <- target_affine
    C <- matrix
    return(solve(W_s) %*% solve(C) %*% W_t)
  }
  # For now other conventions are passthrough; add as needed.
  matrix
}
