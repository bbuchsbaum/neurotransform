#' @title Coordinate Convention Utilities
#' @name coordinates
#' @description
#' Functions for converting between coordinate conventions used in neuroimaging.
#'
#' Internal convention: RAS (Right-Anterior-Superior) in millimeters.
#'
#' @section Supported conventions:
#' \itemize{
#'   \item \strong{RAS}: Right-Anterior-Superior (internal standard)
#'   \item \strong{LPS}: Left-Posterior-Superior (ANTs/ITK convention)
#'   \item \strong{tkRAS}: FreeSurfer tkregister RAS
#'   \item \strong{FSL}: FSL world coordinates (mm, may need scaling)
#' }
NULL

#' Convert LPS coordinates to RAS
#'
#' ANTs and ITK use LPS (Left-Posterior-Superior) convention. This function
#' converts to our internal RAS convention by flipping the first two axes.
#'
#' @param coords Numeric matrix (N x 3) or vector of length 3 in LPS
#' @return Coordinates in RAS convention
#' @export
#' @examples
#' lps <- c(-10, -20, 30)
#' ras <- lps_to_ras(lps)
lps_to_ras <- function(coords) {
  if (is.vector(coords)) {
    if (length(coords) != 3) {
      stop("Vector must have length 3")
    }
    coords[1:2] <- -coords[1:2]
    return(coords)
  }
  if (is.matrix(coords)) {
    if (ncol(coords) != 3) {
      stop("Matrix must have 3 columns")
    }
    coords[, 1:2] <- -coords[, 1:2]
    return(coords)
  }
  stop("coords must be a vector or matrix")
}

#' Convert RAS coordinates to LPS
#'
#' Converts from our internal RAS convention to ANTs/ITK LPS convention.
#'
#' @param coords Numeric matrix (N x 3) or vector of length 3 in RAS
#' @return Coordinates in LPS convention
#' @export
#' @examples
#' ras <- c(10, 20, 30)
#' lps <- ras_to_lps(ras)
ras_to_lps <- function(coords) {
  # Same operation - flip first two axes
  lps_to_ras(coords)
}

#' Convert FreeSurfer tkRAS to scanner RAS
#'
#' FreeSurfer's tkregister uses a shifted coordinate system (tkRAS) that differs
#' from scanner RAS by the c_ras offset stored in the surface header.
#'
#' @param tkras Numeric matrix (N x 3) or vector of length 3 in tkRAS
#' @param c_ras Numeric vector of length 3: the c_ras offset from FreeSurfer
#' @return Coordinates in scanner RAS
#' @export
#' @examples
#' tkras <- c(0, 0, 0)
#' c_ras <- c(0.5, -2.3, 1.1)
#' scanner_ras <- tkras_to_ras(tkras, c_ras)
tkras_to_ras <- function(tkras, c_ras) {
  validate_numeric_vector(c_ras, 3, "c_ras")

  if (is.vector(tkras)) {
    validate_numeric_vector(tkras, 3, "tkras")
    return(tkras + c_ras)
  }
  if (is.matrix(tkras)) {
    if (ncol(tkras) != 3) {
      stop("Matrix must have 3 columns")
    }
    return(sweep(tkras, 2, c_ras, "+"))
  }
  stop("tkras must be a vector or matrix")
}

#' Convert scanner RAS to FreeSurfer tkRAS
#'
#' Converts from scanner RAS to FreeSurfer's tkRAS coordinate system.
#'
#' @param ras Numeric matrix (N x 3) or vector of length 3 in scanner RAS
#' @param c_ras Numeric vector of length 3: the c_ras offset from FreeSurfer
#' @return Coordinates in tkRAS
#' @export
#' @examples
#' scanner_ras <- c(0.5, -2.3, 1.1)
#' c_ras <- c(0.5, -2.3, 1.1)
#' tkras <- ras_to_tkras(scanner_ras, c_ras)
ras_to_tkras <- function(ras, c_ras) {
  validate_numeric_vector(c_ras, 3, "c_ras")

  if (is.vector(ras)) {
    validate_numeric_vector(ras, 3, "ras")
    return(ras - c_ras)
  }
  if (is.matrix(ras)) {
    if (ncol(ras) != 3) {
      stop("Matrix must have 3 columns")
    }
    return(sweep(ras, 2, c_ras, "-"))
  }
  stop("ras must be a vector or matrix")
}

#' Apply affine transformation to coordinates
#'
#' Generic function to apply a 4x4 affine matrix to 3D coordinates.
#'
#' @param coords Numeric matrix (N x 3) or vector of length 3
#' @param affine 4x4 affine transformation matrix
#' @return Transformed coordinates
#' @export
#' @examples
#' coords <- matrix(c(10, 20, 30), ncol = 3)
#' affine <- diag(4)
#' affine[1:3, 4] <- c(5, 10, 15)  # translation
#' apply_affine(coords, affine)
apply_affine <- function(coords, affine) {
  validate_4x4_matrix(affine, "affine")

  if (is.vector(coords)) {
    validate_numeric_vector(coords, 3, "coords")
    hom <- c(coords, 1)
    result <- as.numeric(affine %*% hom)
    return(result[1:3])
  }
  if (is.matrix(coords)) {
    if (ncol(coords) != 3) {
      stop("Matrix must have 3 columns")
    }
    hom <- cbind(coords, 1)
    result <- t(affine %*% t(hom))
    return(result[, 1:3, drop = FALSE])
  }
  stop("coords must be a vector or matrix")
}

#' Compose two affine transformations
#'
#' Returns the composition B(A(x)) = BA * x.
#'
#' @param A First affine (applied first)
#' @param B Second affine (applied second)
#' @return Composed 4x4 affine matrix
#' @export
#' @examples
#' A <- diag(4); A[1:3, 4] <- c(1, 2, 3)
#' B <- diag(4); B[1:3, 4] <- c(4, 5, 6)
#' compose_affines(A, B)
compose_affines <- function(A, B) {
  validate_4x4_matrix(A, "A")
  validate_4x4_matrix(B, "B")
  B %*% A
}

#' Invert an affine transformation
#'
#' @param affine 4x4 affine matrix
#' @return Inverted 4x4 affine matrix
#' @export
#' @examples
#' affine <- diag(4)
#' affine[1:3, 4] <- c(5, 10, 15)
#' inv <- invert_affine(affine)
invert_affine <- function(affine) {
  validate_4x4_matrix(affine, "affine")
  solve(affine)
}

#' Voxel indices to world coordinates
#'
#' Convert voxel indices to world (RAS mm) coordinates using an affine.
#'
#' @param voxels Integer matrix (N x 3) of voxel indices (0-based)
#' @param affine 4x4 voxel-to-world affine matrix
#' @return Matrix of world coordinates (N x 3)
#' @keywords internal
voxels_to_world <- function(voxels, affine) {
  apply_affine(voxels, affine)
}

#' World coordinates to voxel indices
#'
#' Convert world (RAS mm) coordinates to voxel indices using an affine.
#'
#' @param world Numeric matrix (N x 3) of world coordinates
#' @param affine 4x4 voxel-to-world affine matrix
#' @return Matrix of voxel indices (0-based, continuous - round for discrete)
#' @keywords internal
world_to_voxels <- function(world, affine) {
  inv_affine <- invert_affine(affine)
  apply_affine(world, inv_affine)
}
