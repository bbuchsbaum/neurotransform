#' @title AFNI Warp and Affine Handlers
#' @name afni_warp
#' @description
#' AFNI-specific coordinate convention handling. AFNI uses RAI
#' (Right-Anterior-Inferior) convention internally. These functions
#' handle the RAI <-> RAS conversion.
NULL

#' Apply AFNI warp (RAI displacement) to coordinates
#'
#' AFNI 3dQwarp/Nwarp displacement fields are defined in RAI space.
#' This function handles the RAS <-> RAI conversion automatically.
#'
#' @param morphism A Warp3DMorphism with warp_type="afni"
#' @param coords Numeric matrix (N x 3) of RAS world coordinates
#' @return Transformed RAS coordinates
#' @keywords internal
afni_warp_transform_coords <- function(morphism, coords) {
  if (!nzchar(morphism@warp_path)) stop("Warp3DMorphism missing warp_path")

  cache_env <- morphism@cache %||% new_cache_env()
  warp <- load_warp_array(morphism, cache_env = cache_env)

  # AFNI warps are defined in RAI; if the header world_to_vox looks RAS (positive z scale), flip it
  if (warp$world_to_vox[3, 3] > 0) {
    rai_flip <- diag(c(1, 1, -1, 1))
    warp$world_to_vox <- rai_flip %*% warp$world_to_vox
  }

  # AFNI displacements are in RAI; convert coords (RAS) to RAI, apply, flip back
  ras_to_rai <- diag(c(1, 1, -1, 1))
  rai_to_ras <- ras_to_rai  # same matrix (self-inverse for Z flip)

  coords_h <- cbind(coords, 1)
  coords_rai <- (coords_h %*% t(ras_to_rai))[, 1:3, drop = FALSE]
  out_rai <- cpp_apply_warp_field(coords_rai, warp$array, warp$dim, warp$world_to_vox)
  out_h <- cbind(out_rai, 1)
  out_ras <- (out_h %*% t(rai_to_ras))[, 1:3, drop = FALSE]
  out_ras
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
