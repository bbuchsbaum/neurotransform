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
#' @param dim Optional image dimensions. Required when `affine` is
#'   right-handed because FSL's handedness swap needs the x-extent.
#' @return 4x4 scaling matrix
#' @keywords internal
fsl_vox_to_fsl <- function(affine, dim = NULL) {
  stopifnot(is.matrix(affine), all(dim(affine) == 4))
  sp <- fsl_spacing_from_affine(affine)
  swp <- diag(4)
  det_sign <- det(affine[1:3, 1:3, drop = FALSE])
  if (det_sign > 0 && is.null(dim)) {
    stop("dim must be supplied for right-handed affines to compute the FSL handedness swap")
  }
  if (det_sign > 0) {
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

# Validate one image geometry used to interpret FSL fields. A required pair
# must be present; an optional pair must be supplied completely or omitted.
.fsl_check_geometry <- function(affine, dims, label, required = TRUE,
                                why = paste0("FSL field values are scaled-voxel coordinates of the ",
                                             label, " image, which cannot be recovered from the warp file")) {
  if (is.null(affine) && is.null(dims)) {
    if (!required) return(invisible(NULL))
    stop("FSL warps require ", label, "_affine and ", label, "_dim: ", why, ".", call. = FALSE)
  }
  if (is.null(affine) || is.null(dims)) {
    stop("Supply both ", label, "_affine and ", label, "_dim, or neither.", call. = FALSE)
  }
  if (!is.matrix(affine) || !identical(dim(affine), c(4L, 4L)) ||
      any(!is.finite(affine)) || abs(det(affine[1:3, 1:3])) < 1e-12 ||
      any(abs(affine[4, ] - c(0, 0, 0, 1)) > 1e-10)) {
    stop(label, "_affine must be a finite, nonsingular 4x4 voxel-to-RAS affine", call. = FALSE)
  }
  if (!is.numeric(dims) || length(dims) != 3L || any(!is.finite(dims)) ||
      any(dims < 1 | dims != floor(dims))) {
    stop(label, "_dim must contain three positive integer dimensions ",
         "(use dim(image)[1:3] for 4D images)", call. = FALSE)
  }
  invisible(NULL)
}

# A dense FSL field maps reference FSL coordinates to source FSL coordinates.
# Normalize both representations to RAS displacements at the field lattice so
# coordinate transforms, flattened resampling plans, and Jacobians share it.
.fsl_dense_to_ras_displacement <- function(field, params) {
  validate <- .fsl_check_geometry
  sa <- params$source_affine
  sd <- params$source_dim
  # The field grid is normally the reference image grid. Explicit reference
  # geometry is needed if the dense field is sampled on a different lattice.
  ta <- params$target_affine
  td <- params$target_dim
  if (is.null(ta) && is.null(td)) {
    ta <- field$vox_to_world
    td <- field$dim
  }
  validate(sa, sd, "source")
  validate(ta, td, "target")
  world <- .grid_world_coords_matrix(field$dim, field$vox_to_world)
  native <- matrix(field$array, ncol = 3L, byrow = TRUE)
  if (!identical(params$def_type, "absolute")) {
    ref_fsl <- (cbind(world, 1) %*% t(fsl_world_to_fsl(ta, td)))[, 1:3, drop = FALSE]
    native <- native + ref_fsl
  }
  source_world <- (cbind(native, 1) %*% t(fsl_fsl_to_world(sa, sd)))[, 1:3, drop = FALSE]
  field$array <- as.numeric(t(source_world - world))
  field$def_type <- "relative"
  field
}

# FNIRT spline-coefficient files (fnirt --cout), established against native
# FSL 5.0.9 fnirt/fnirtfileutils/applywarp (inst/extdata/fsl_coef_oracle):
#   dim[1:3]      coefficient counts; dim[4] = 3 displacement components
#   pixdim[1:3]   knot spacing k in reference voxels
#   intent        2007 cubic, 2009 quadratic spline
#   intent_p1..3  reference voxel size (mm)
#   qoffset       reference dimensions
#   sform         the --aff FLIRT matrix A (source FSL -> reference FSL)
# The reference orientation and origin are not stored, so the reference
# geometry must be supplied. Displacements d are reference FSL mm, with knots
# over the internal (x-flipped for right-handed images) reference index u:
#   d(u) = sum C[a,b,c] B(ux/kx - a + o) B(uy/ky - b + o) B(uz/kz - c + o),
#   o = 1 if k > 1 else 0, with unnormalized weights, and
#   source_FSL = inv(A) ref_FSL + d.

.fsl_bspline <- function(t, order) {
  a <- abs(t)
  if (order == 3L) {
    ifelse(a < 1, 2 / 3 - a^2 + a^3 / 2, ifelse(a < 2, (2 - a)^3 / 6, 0))
  } else {
    ifelse(a < 0.5, 0.75 - a^2, ifelse(a < 1.5, (a - 1.5)^2 / 2, 0))
  }
}

# Weights mapping coefficient index (columns) to internal voxel index (rows).
.fsl_spline_weights <- function(n_vox, n_coef, knot, order) {
  offset <- if (knot > 1) 1 else 0
  outer(0:(n_vox - 1), 0:(n_coef - 1), function(u, a) .fsl_bspline(u / knot - a + offset, order))
}

# Separable evaluation: C (cx, cy, cz) -> field (Nx, Ny, Nz).
.fsl_spline_evaluate <- function(C, wx, wy, wz) {
  cd <- dim(C)
  n <- c(nrow(wx), nrow(wy), nrow(wz))
  a <- array(wx %*% matrix(C, cd[1]), c(n[1], cd[2], cd[3]))
  a <- aperm(array(wy %*% matrix(aperm(a, c(2, 1, 3)), cd[2]), c(n[2], n[1], cd[3])), c(2, 1, 3))
  array(matrix(a, ncol = cd[3]) %*% t(wz), n)
}

# Decode a coefficient file into the relative dense FSL field that
# fnirtfileutils --withaff (and fnirt --fout) would write on the reference grid.
.fsl_coef_to_dense <- function(coef, params) {
  .fsl_check_geometry(params$target_affine, params$target_dim, "target",
                      why = "the reference orientation and origin are not stored in FNIRT coefficient files")
  ref_affine <- params$target_affine
  ref_dim <- as.integer(params$target_dim)
  if (!identical(ref_dim, coef$ref_dim)) {
    stop("target_dim (", paste(ref_dim, collapse = "x"), ") does not match the reference ",
         "dimensions stored in the coefficient file (", paste(coef$ref_dim, collapse = "x"), ")",
         call. = FALSE)
  }
  voxel <- fsl_spacing_from_affine(ref_affine)[1:3]
  if (any(abs(voxel - coef$ref_voxel) > 1e-4 * coef$ref_voxel)) {
    stop("target_affine voxel sizes (", paste(signif(voxel, 6), collapse = ", "),
         ") do not match the reference voxel sizes stored in the coefficient file (",
         paste(signif(coef$ref_voxel, 6), collapse = ", "), ")", call. = FALSE)
  }
  n_coef <- dim(coef$coef)[1:3]
  if (coef$order == 3L) {
    expected <- ifelse(coef$knot > 1, ceiling((ref_dim + 1) / coef$knot) + 2, ref_dim)
    if (any(n_coef != expected)) {
      stop("Coefficient grid (", paste(n_coef, collapse = "x"), ") is inconsistent with ",
           "the reference dimensions and knot spacing (expected ",
           paste(expected, collapse = "x"), ")", call. = FALSE)
    }
  }

  w <- lapply(1:3, function(a) .fsl_spline_weights(ref_dim[a], n_coef[a], coef$knot[a], coef$order))
  d <- vapply(1:3, function(m) {
    # Keep size-1 coefficient axes (single-slice references with knot spacing 1).
    component <- array(coef$coef[, , , m, drop = FALSE], n_coef)
    as.numeric(.fsl_spline_evaluate(component, w[[1]], w[[2]], w[[3]]))
  }, numeric(prod(ref_dim)))
  if (det(ref_affine[1:3, 1:3]) > 0) {
    # Knots run over FSL's x-flipped index; reorder to stored voxel order.
    stored_x <- rep(seq.int(ref_dim[1], 1L), times = prod(ref_dim[2:3])) +
      ref_dim[1] * rep(seq_len(prod(ref_dim[2:3])) - 1L, each = ref_dim[1])
    d <- d[stored_x, , drop = FALSE]
  }

  vox <- as.matrix(expand.grid(lapply(ref_dim, function(n) 0:(n - 1))))
  ref_fsl <- (cbind(vox, 1) %*% t(fsl_vox_to_fsl(ref_affine, ref_dim)))[, 1:3, drop = FALSE]
  src_fsl <- (cbind(ref_fsl, 1) %*% t(solve(coef$affine)))[, 1:3, drop = FALSE] + d
  list(
    array = as.numeric(t(src_fsl - ref_fsl)),
    dim = ref_dim,
    vox_to_world = ref_affine,
    world_to_vox = solve(ref_affine)
  )
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
#' @param source_dim Optional source image dimensions for right-handed FSL affines
#' @param target_dim Optional target image dimensions for right-handed FSL affines
#' @param cost Path cost
#' @param method_tag Method tag
#' @return Affine3DMorphism object
#' @keywords internal
fsl_load_flirt_morphism <- function(source, target, mat_path,
                                    source_dim = NULL, target_dim = NULL,
                                    cost = 1.0, method_tag = "anatomical") {
  if (!file.exists(mat_path)) stop("FLIRT matrix not found: ", mat_path)
  mat <- as.matrix(read.table(mat_path))

  # Extract affines - caller should provide 4x4 matrices or objects with @geometry@affine
  src_aff <- if (is.matrix(source)) source else source@geometry@affine
  ref_aff <- if (is.matrix(target)) target else target@geometry@affine

  phi <- fsl_flirt_to_internal_affine(
    mat, src_aff, ref_aff,
    source_dim = source_dim,
    ref_dim = target_dim
  )

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
#' Determines whether a dense FSL warp stores relative displacements or
#' absolute source coordinates, using two kinds of evidence on a regular
#' subsample of non-zero lattice vectors:
#' \enumerate{
#'   \item Field of view (when \code{source_affine} and \code{source_dim} are
#'     given): the reading whose implied source FSL coordinates fall inside the
#'     source image (within one voxel) for a clearly larger fraction of samples
#'     (by at least 0.25) wins.
#'   \item Jacobian: field values are fitted as \eqn{v \approx M x + t} in
#'     reference FSL coordinates, so the implied reference-to-source Jacobian is
#'     \eqn{M} (absolute) or \eqn{M + I} (relative). A reading wins when the
#'     other implies a volume change at least twice as extreme
#'     (\eqn{|\log|\det J||} larger by \eqn{\log 2}).
#' }
#' The field-of-view test decides first. If neither is decisive the function
#' stops rather than guess; pass \code{def_type} explicitly in that case.
#' Scaled or axis-permuted mappings (e.g. a small or sagittally stored source)
#' can make the Jacobian test ambiguous, which the field-of-view test resolves.
#'
#' @param warp_path Path to FNIRT warp file
#' @param sample_n Approximate number of lattice points sampled
#' @param threshold_mm Deprecated and ignored; retained for compatibility.
#' @param source_affine,source_dim Optional source image voxel-to-RAS affine
#'   and dimensions, enabling the field-of-view test.
#' @return "relative" or "absolute"
#' @export
#' @examples
#' \dontrun{
#' def_type <- detect_fnirt_def_type("warp.nii.gz", source_affine = src_affine,
#'                                   source_dim = src_dim)
#' }
detect_fnirt_def_type <- function(warp_path, sample_n = 1000, threshold_mm = NULL,
                                  source_affine = NULL, source_dim = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 required for FNIRT detection")
  }
  if (!file.exists(warp_path)) stop("Warp file not found: ", warp_path)
  use_fov <- !is.null(source_affine) || !is.null(source_dim)
  if (use_fov) .fsl_check_geometry(source_affine, source_dim, "source")

  img <- neuroim2::read_vec(warp_path)
  dim4 <- dim(img)
  if (length(dim4) < 4 || dim4[4] < 3) stop("Warp must be 4D with last dim length 3")
  dims <- as.integer(dim4[1:3])

  # Regular subgrid of the lattice, 0-based.
  per_axis <- max(2L, ceiling(sample_n^(1 / 3)))
  axes <- lapply(dims, function(n) unique(round(seq(0, n - 1, length.out = min(n, per_axis)))))
  vox <- as.matrix(expand.grid(axes))
  arr <- as.array(img)
  vals <- sapply(1:3, function(k) arr[cbind(vox + 1L, k)])
  vals <- matrix(vals, ncol = 3)
  keep <- rowSums(is.finite(vals)) == 3L & rowSums(vals != 0) > 0L
  # An all-zero field is the relative identity.
  if (!any(keep)) return("relative")
  ref_fsl <- (cbind(vox, 1) %*% t(fsl_vox_to_fsl(neuroim2::trans(img), dims)))[, 1:3, drop = FALSE]
  vals <- vals[keep, , drop = FALSE]
  ref_fsl <- ref_fsl[keep, , drop = FALSE]

  if (use_fov) {
    spacing <- fsl_spacing_from_affine(source_affine)[1:3]
    upper <- (as.integer(source_dim) - 1) * spacing + spacing
    inside <- function(p) {
      mean(rowSums(sweep(p, 2, -spacing, ">=") & sweep(p, 2, upper, "<=")) == 3L)
    }
    fov <- c(absolute = inside(vals), relative = inside(vals + ref_fsl))
    if (abs(fov[["absolute"]] - fov[["relative"]]) >= 0.25) return(names(which.max(fov)))
  }

  if (all(dims >= 2L) && nrow(vals) >= 12L) {
    coef <- tryCatch(qr.solve(cbind(ref_fsl, 1), vals), error = function(e) NULL)
    if (!is.null(coef)) {
      M <- t(coef[1:3, , drop = FALSE])
      distortion <- function(J) {
        d <- abs(det(J))
        if (!is.finite(d) || d <= 0) Inf else abs(log(d))
      }
      jac <- c(absolute = distortion(M), relative = distortion(M + diag(3)))
      if (isTRUE(abs(jac[["absolute"]] - jac[["relative"]]) >= log(2))) {
        return(names(which.min(jac)))
      }
    }
  }
  stop("Cannot determine whether ", basename(warp_path), " stores relative displacements ",
       "or absolute coordinates",
       if (!use_fov) " (source_affine and source_dim enable a field-of-view test)",
       "; pass def_type explicitly.", call. = FALSE)
}
