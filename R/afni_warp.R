#' @title AFNI Warp and Affine Handlers
#' @name afni_warp
#' @description
#' AFNI-specific coordinate convention handling. AFNI 3dQwarp stores
#' displacement fields in DICOM/LPS convention. These functions handle
#' the LPS <-> RAS conversion for warps and RAI <-> RAS for affines.
#'
#' @section Coordinate Systems:
#' \itemize{
#'   \item LPS (DICOM, used in AFNI warp displacements): +X=Left, +Y=Posterior, +Z=Superior
#'   \item RAI (AFNI affines): +X=Right, +Y=Anterior, +Z=Inferior
#'   \item RAS (NIfTI/neuroim2): +X=Right, +Y=Anterior, +Z=Superior
#' }
#'
#' LPS to RAS: negate X and Y (Z is the same).
#' RAI to RAS: negate Z only (X and Y are the same).
NULL

#' Apply AFNI warp displacement to coordinates
#'
#' AFNI 3dQwarp displacement fields store displacements in DICOM/LPS coordinate
#' convention (Left-Posterior-Superior), but neuroim2/NIfTI uses RAS
#' (Right-Anterior-Superior).
#'
#' The approach:
#' 1. Use RAS world_to_vox to get voxel indices (NIfTI geometry is RAS)
#' 2. Sample the displacement values (stored in DICOM/LPS convention)
#' 3. Convert displacement from LPS to RAS (negate X and Y components)
#' 4. Add RAS displacement to input RAS coords
#'
#' @param morphism A Warp3DMorphism with warp_type="afni"
#' @param coords Numeric matrix (N x 3) of RAS world coordinates
#' @return Transformed RAS coordinates
#' @keywords internal
afni_warp_transform_coords <- function(morphism, coords) {
  if (!nzchar(morphism@warp_path)) stop("Warp3DMorphism missing warp_path")

  cache_env <- morphism@cache %||% new_cache_env()
  warp <- load_warp_array(morphism, cache_env = cache_env)

  # LPS->RAS conversion is now done at load time in load_warp_array().
  # The warp array already contains RAS-convention displacements,
  # so we can apply directly.
  cpp_apply_warp_field(coords, warp$array, warp$dim, warp$world_to_vox)
}

#' Read AFNI .aff12.1D affine matrix
#'
#' Reads a 3x4 affine matrix in RAI coordinates and returns a 4x4 matrix.
#'
#' @param path Path to .aff12.1D file
#' @return 4x4 affine matrix (in RAI)
#' @export
#' @examples
#' \dontrun{
#' mat_rai <- afni_read_aff12("transform.aff12.1D")
#' }
afni_read_aff12 <- function(path) {
  if (!file.exists(path)) stop("AFNI aff12 file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  vals <- scan(text = lines, quiet = TRUE)
  if (length(vals) != 12) stop("AFNI aff12 file must contain 12 numbers")
  mat <- matrix(vals, nrow = 3, ncol = 4, byrow = TRUE)
  rbind(mat, c(0, 0, 0, 1))
}

#' Convert RAI affine to RAS
#'
#' Converts an AFNI RAI affine matrix to our internal RAS convention
#' by flipping the Z axis. Optionally applies AFNI-style deobliquing
#' compensation when source/target image affines are provided.
#'
#' @param mat_rai 4x4 affine in RAI coordinates
#' @param source_affine Optional 4x4 source (moving) voxel-to-world affine
#' @param target_affine Optional 4x4 target (reference) voxel-to-world affine
#' @param oblique_correction Logical; apply AFNI cardinal rotation correction
#' @return 4x4 affine in RAS coordinates
#' @export
#' @examples
#' mat_rai <- diag(4)
#' mat_ras <- afni_aff12_to_ras(mat_rai)
afni_aff12_to_ras <- function(mat_rai,
                              source_affine = NULL,
                              target_affine = NULL,
                              oblique_correction = TRUE) {
  stopifnot(is.matrix(mat_rai), all(dim(mat_rai) == 4))
  flip <- diag(c(1, 1, -1, 1))
  mat_ras <- flip %*% mat_rai %*% flip

  if (isTRUE(oblique_correction)) {
    if (!is.null(target_affine) && afni_is_oblique(target_affine)) {
      mat_ras <- mat_ras %*% afni_cardinal_rotation(target_affine, real_to_card = TRUE)
    }
    if (!is.null(source_affine) && afni_is_oblique(source_affine)) {
      mat_ras <- afni_cardinal_rotation(source_affine, real_to_card = FALSE) %*% mat_ras
    }
  }

  mat_ras
}

# Convert internal RAS pullback affine to AFNI RAI affine, with optional
# AFNI-style deobliquing compensation.
afni_ras_to_aff12 <- function(mat_ras,
                              source_affine = NULL,
                              target_affine = NULL,
                              oblique_correction = TRUE) {
  stopifnot(is.matrix(mat_ras), all(dim(mat_ras) == 4))
  out <- mat_ras

  if (isTRUE(oblique_correction)) {
    if (!is.null(target_affine) && afni_is_oblique(target_affine)) {
      out <- out %*% afni_cardinal_rotation(target_affine, real_to_card = FALSE)
    }
    if (!is.null(source_affine) && afni_is_oblique(source_affine)) {
      out <- afni_cardinal_rotation(source_affine, real_to_card = TRUE) %*% out
    }
  }

  flip <- diag(c(1, 1, -1, 1))
  flip %*% out %*% flip
}

#' Detect whether an affine is oblique (not cardinal axis-aligned)
#'
#' @param affine 4x4 voxel-to-world affine
#' @param threshold_deg Angular threshold in degrees
#' @return Logical
#' @export
afni_is_oblique <- function(affine, threshold_deg = 0.01) {
  stopifnot(is.matrix(affine), all(dim(affine) == 4))
  A <- affine[1:3, 1:3, drop = FALSE]
  norms <- sqrt(colSums(A^2))
  norms[norms == 0] <- 1
  dirs <- sweep(A, 2, norms, "/")
  max_abs <- pmin(1, pmax(0, apply(abs(dirs), 2, max)))
  ang <- acos(max_abs) * 180 / pi
  max(ang, na.rm = TRUE) > threshold_deg
}

# Compute AFNI "DICOM cardinal" matrix by dropping obliquity and preserving
# voxel sizes and origin.
afni_dicom_real_to_card <- function(oblique) {
  stopifnot(is.matrix(oblique), all(dim(oblique) == 4))

  retval <- diag(4)
  retval[1:3, 4] <- oblique[1:3, 4]

  A <- oblique[1:3, 1:3, drop = FALSE]
  vs <- sqrt(colSums(A^2))
  maxabs <- apply(abs(A), 2, max)
  maxabs[maxabs == 0] <- 1

  cosines <- sweep(A, 2, maxabs, "/")
  cosines[abs(cosines) < 1] <- 0
  retval[1:3, 1:3] <- sweep(cosines, 2, round(vs, 4), "*")
  retval
}

# Rotation matrix that maps between AFNI "real" and "cardinal" spaces.
afni_cardinal_rotation <- function(oblique, real_to_card = TRUE) {
  card <- afni_dicom_real_to_card(oblique)
  if (isTRUE(real_to_card)) card %*% solve(oblique) else oblique %*% solve(card)
}

#' Build Affine3DMorphism from AFNI .aff12.1D
#'
#' Creates a morphism from an AFNI affine file, handling RAI->RAS conversion.
#'
#' @param source Source domain hash
#' @param target Target domain hash
#' @param aff12_path Path to .aff12.1D file
#' @param direction Whether matrix maps source_to_target or target_to_source
#' @param cost Path cost
#' @param method_tag Method tag
#' @return Affine3DMorphism object
#' @keywords internal
afni_load_affine_morphism <- function(source, target, aff12_path,
                                      direction = c("source_to_target", "target_to_source"),
                                      cost = 1.0, method_tag = "anatomical") {
  direction <- match.arg(direction)
  m_rai <- afni_read_aff12(aff12_path)
  if (direction == "target_to_source") {
    m_rai <- invert_affine(m_rai)
  }
  m_ras <- afni_aff12_to_ras(m_rai)
  Affine3DMorphism(
    source = source,
    target = target,
    matrix = m_ras,
    cost = cost,
    method_tag = method_tag
  )
}
