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
#' by flipping the Z axis.
#'
#' @param mat_rai 4x4 affine in RAI coordinates
#' @return 4x4 affine in RAS coordinates
#' @export
#' @examples
#' mat_rai <- diag(4)
#' mat_ras <- afni_aff12_to_ras(mat_rai)
afni_aff12_to_ras <- function(mat_rai) {
  stopifnot(is.matrix(mat_rai), all(dim(mat_rai) == 4))
  flip <- diag(c(1, 1, -1, 1))
  flip %*% mat_rai %*% flip
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
