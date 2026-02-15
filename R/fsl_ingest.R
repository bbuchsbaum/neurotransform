#' @title FSL Coordinate Convention Handlers
#' @name fsl_ingest
#' @description
#' FSL-specific coordinate convention handling. FSL uses scaled voxel
#' coordinates for FLIRT matrices and has specific conventions for
#' FNIRT warp fields (relative vs absolute).
NULL

#' Compute voxel spacing from affine
#'
#' @param affine 4x4 voxel-to-world affine
#' @return Numeric vector of spacing in each dimension
#' @keywords internal
fsl_spacing_from_affine <- function(affine) {
  stopifnot(is.matrix(affine), all(dim(affine) == 4))
  sqrt(colSums(affine[1:3, 1:3, drop = FALSE]^2))
}

#' Build voxel-to-FSL scaling matrix
#'
#' @param affine 4x4 voxel-to-world affine
#' @param dim Optional image dimensions (required to apply FSL handedness swap)
#' @return 4x4 scaling matrix
#' @keywords internal
fsl_vox_to_fsl <- function(affine, dim = NULL) {
  stopifnot(is.matrix(affine), all(dim(affine) == 4))
  sp <- fsl_spacing_from_affine(affine)
  swp <- diag(4)
  # Preserve historical behavior unless dimensions are explicitly provided.
  # Handedness swap in FSL needs image extent to set the x-offset correctly.
  if (!is.null(dim) && det(affine[1:3, 1:3, drop = FALSE]) > 0) {
    swp[1, 1] <- -1
    swp[1, 4] <- (as.integer(dim)[1] - 1) * sp[1]
  }
  swp %*% diag(c(sp, 1))
}

#' Build FSL-to-voxel scaling matrix
#'
#' @param affine 4x4 voxel-to-world affine
#' @param dim Optional image dimensions (required to apply FSL handedness swap)
#' @return 4x4 inverse scaling matrix
#' @keywords internal
fsl_fsl_to_vox <- function(affine, dim = NULL) {
  solve(fsl_vox_to_fsl(affine, dim = dim))
}

#' Convert world coords to FSL coords
#'
#' @param affine 4x4 voxel-to-world affine
#' @param dim Optional image dimensions
#' @return 4x4 world-to-FSL transform
#' @keywords internal
fsl_world_to_fsl <- function(affine, dim = NULL) {
  fsl_vox_to_fsl(affine, dim = dim) %*% invert_affine(affine)
}

#' Convert FSL coords to world coords
#'
#' @param affine 4x4 voxel-to-world affine
#' @param dim Optional image dimensions
#' @return 4x4 FSL-to-world transform
#' @keywords internal
fsl_fsl_to_world <- function(affine, dim = NULL) {
  affine %*% fsl_fsl_to_vox(affine, dim = dim)
}

#' Convert FLIRT matrix to internal affine
#'
#' Given a FLIRT matrix (source_FSL -> ref_FSL), returns the internal
#' ref_world -> src_world affine for pullback semantics.
#'
#' @param flirt_mat 4x4 FLIRT matrix
#' @param source_affine 4x4 voxel-to-world for source image
#' @param ref_affine 4x4 voxel-to-world for reference image
#' @param source_dim Optional source image dimensions (needed for handedness swap)
#' @param ref_dim Optional reference image dimensions (needed for handedness swap)
#' @return 4x4 internal affine (ref_world -> src_world)
#' @export
#' @examples
#' \dontrun{
#' flirt_mat <- as.matrix(read.table("xform.mat"))
#' src_aff <- diag(4)  # from source image header
#' ref_aff <- diag(4)  # from reference image header
#' internal <- fsl_flirt_to_internal_affine(flirt_mat, src_aff, ref_aff)
#' }
fsl_flirt_to_internal_affine <- function(flirt_mat, source_affine, ref_affine,
                                         source_dim = NULL, ref_dim = NULL) {
  if (!is.matrix(flirt_mat) || any(dim(flirt_mat) != 4)) {
    stop("FLIRT matrix must be 4x4")
  }
  W_ref_to_fsl <- fsl_world_to_fsl(ref_affine, dim = ref_dim)
  W_fsl_to_src <- fsl_fsl_to_world(source_affine, dim = source_dim)
  # FLIRT matrix maps source_fsl -> ref_fsl. Pullback needs inverse.
  phi <- W_fsl_to_src %*% invert_affine(flirt_mat) %*% W_ref_to_fsl
  phi
}

#' Build Affine3DMorphism from FLIRT matrix
#'
#' Creates a morphism from an FSL FLIRT/MCFLIRT matrix file.
#'
#' @param source Source domain (with @geometry@affine) or affine matrix
#' @param target Target domain (with @geometry@affine) or affine matrix
#' @param mat_path Path to FLIRT .mat file
#' @param cost Path cost
#' @param method_tag Method tag
#' @return Affine3DMorphism object
#' @keywords internal
fsl_load_flirt_morphism <- function(source, target, mat_path, cost = 1.0, method_tag = "anatomical") {
  if (!file.exists(mat_path)) stop("FLIRT matrix not found: ", mat_path)
  mat <- as.matrix(read.table(mat_path))

  # Extract affines - caller should provide 4x4 matrices or objects with @geometry@affine
  src_aff <- if (is.matrix(source)) source else source@geometry@affine
  ref_aff <- if (is.matrix(target)) target else target@geometry@affine

  phi <- fsl_flirt_to_internal_affine(mat, src_aff, ref_aff)

  source_id <- if (is.matrix(source)) "source" else source@domain_hash

  target_id <- if (is.matrix(target)) "target" else target@domain_hash

  Affine3DMorphism(
    source = source_id,
    target = target_id,
    matrix = phi,
    cost = cost,
    method_tag = method_tag
  )
}

#' Detect FNIRT deformation type (relative vs absolute)
#'
#' Heuristically determines whether an FNIRT warp stores relative
#' displacements or absolute coordinates.
#'
#' @param warp_path Path to FNIRT warp file
#' @param sample_n Number of voxels to sample for heuristic
#' @param threshold_mm Threshold for displacement magnitude heuristic
#' @return "relative" or "absolute"
#' @export
#' @examples
#' \dontrun{
#' def_type <- detect_fnirt_def_type("warp.nii.gz")
#' }
detect_fnirt_def_type <- function(warp_path, sample_n = 200, threshold_mm = 5) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 required for FNIRT detection")
  }
  if (!file.exists(warp_path)) stop("Warp file not found: ", warp_path)

  img <- neuroim2::read_vec(warp_path)
  dim4 <- dim(img)
  if (length(dim4) < 4 || dim4[4] < 3) stop("Warp must be 4D with last dim length 3")

  aff <- neuroim2::trans(img)

  # Sample voxels uniformly (0-based indices)
  set.seed(1L)
  n_samples <- min(sample_n, prod(dim4[1:3]))
  i <- sample.int(dim4[1], n_samples, replace = TRUE) - 1L
  j <- sample.int(dim4[2], n_samples, replace = TRUE) - 1L
  k <- sample.int(dim4[3], n_samples, replace = TRUE) - 1L

  vals <- matrix(NA_real_, nrow = n_samples, ncol = 3)
  coords <- matrix(NA_real_, nrow = n_samples, ncol = 3)
  vox_to_world <- aff

  for (idx in seq_len(n_samples)) {
    vals[idx, ] <- img[i[idx] + 1, j[idx] + 1, k[idx] + 1, 1:3]
    hom <- c(i[idx], j[idx], k[idx], 1)
    coords[idx, ] <- as.numeric(vox_to_world %*% hom)[1:3]
  }

  if (all(is.na(vals))) return("relative")

  # Heuristic: absolute fields track coordinates (high correlation)
  # Relative fields are small shifts with low coord correlation
  cors <- sapply(1:3, function(k) {
    suppressWarnings(cor(vals[, k], coords[, k], use = "complete.obs"))
  })

  if (all(!is.na(cors)) && all(cors > 0.9)) {
    return("absolute")
  }

  # Fallback: magnitude heuristic
  delta <- rowMeans(abs(vals), na.rm = TRUE)
  if (median(delta, na.rm = TRUE) < threshold_mm) "relative" else "absolute"
}
