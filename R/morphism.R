#' @title Morphism Constructors and Methods
#' @name morphism
#' @description
#' Constructors and methods for morphism classes. Morphisms represent
#' coordinate transformations with pullback semantics.
NULL

# =============================================================================
# MORPHISM HASH
# =============================================================================

#' Compute hash for a morphism
#'
#' Creates a stable hash from morphism properties for identity and caching.
#'
#' @param object A Morphism object
#' @return Character hash string
#' @export
morphism_hash <- function(object) {
  if (methods::is(object, "MorphismPath")) {
    comp_hashes <- vapply(object@morphisms, morphism_hash, character(1))
    return(compute_hash("morphism_path", comp_hashes))
  }
  if (!methods::is(object, "Morphism")) stop("object must be a Morphism")
  compute_hash(
    object@id,
    object@source,
    object@target,
    object@kind,
    object@params,
    object@method_tag
  )
}

# =============================================================================
# KIND HELPERS
# =============================================================================

#' Get the kind of a morphism
#'
#' Returns the morphism kind (identity, affine3d, warp3d, vol2surf, surf2surf).
#'
#' @param object A Morphism object
#' @return Character string indicating kind
#' @export
#' @examples
#' aff <- Affine3DMorphism("a", "b", diag(4))
#' morphism_kind(aff)  # "affine3d"
morphism_kind <- function(object) {
  if (!methods::is(object, "Morphism")) stop("object must be a Morphism")
  k <- object@kind
  if (length(k) == 1L && nzchar(k)) return(k)
  # Fallback: derive from class name
  tolower(gsub("Morphism$", "", class(object)[1]))
}

#' Check if morphism is linear (affine or identity)
#'
#' @param object A Morphism object
#' @return Logical
#' @export
is_linear_morphism <- function(object) {
  morphism_kind(object) %in% c("affine3d", "identity")
}

#' Check if morphism is a warp field
#'
#' @param object A Morphism object
#' @return Logical
#' @export
is_warp_morphism <- function(object) {
  morphism_kind(object) %in% c("warp3d")
}

#' Check if morphism is exactly invertible
#'
#' Returns TRUE only if a geometric inverse is available via \code{invert()}.
#'
#' @param object A Morphism object
#' @return Logical
#' @export
#' @examples
#' aff <- Affine3DMorphism("a", "b", diag(4))
#' is_invertible(aff)  # TRUE
is_invertible <- function(object) {
  if (!methods::is(object, "Morphism")) stop("object must be a Morphism")
  object@inverse_type %in% c("exact", "approximate")
}

#' Check if morphism has an adjoint
#'
#' Returns TRUE if \code{adjoint()} is available for this morphism.
#' For linear/invertible morphisms, the adjoint is the inverse.
#'
#' @param object A Morphism object
#' @return Logical
#' @export
has_adjoint <- function(object) {
  if (!methods::is(object, "Morphism")) stop("object must be a Morphism")
  object@inverse_type %in% c("exact", "approximate", "adjoint")
}

# =============================================================================
# CONSTRUCTORS
# =============================================================================

#' Create an identity morphism
#'
#' The identity morphism maps a domain to itself. Cost is 0.
#'
#' @param domain Domain hash string
#' @return IdentityMorphism object
#' @export
#' @examples
#' id <- IdentityMorphism("my_domain_hash")
IdentityMorphism <- function(domain) {
  if (!is.character(domain) || length(domain) != 1L || !nzchar(domain)) {
    stop("domain must be a non-empty character string")
  }
  m <- methods::new("IdentityMorphism",
           id = paste0("id_", domain),
           source = domain,
           target = domain,
           kind = "identity",
           params = list(),
           inverse_params = list(),
           coverage = 1.0,
           cost = 0.0,
           method_tag = "identity",
           hash = "",
           inverse_type = "exact",
           inverse_quality = 1.0,
           inverse_method = "identity")
  m@hash <- morphism_hash(m)
  m
}

#' Create a 3D affine morphism
#'
#' Represents an affine transformation between volume spaces.
#' The matrix maps target coords -> source coords (pullback semantics).
#'
#' @param source Source domain hash
#' @param target Target domain hash
#' @param matrix 4x4 affine matrix (target -> source)
#' @param cost Path cost (default 1.0)
#' @param method_tag Method tag (default "anatomical")
#' @return Affine3DMorphism object
#' @export
#' @examples
#' # Translation by (10, 20, 30)
#' mat <- diag(4)
#' mat[1:3, 4] <- c(10, 20, 30)
#' aff <- Affine3DMorphism("native", "mni", mat)
Affine3DMorphism <- function(source, target, matrix, cost = 1.0, method_tag = "anatomical") {
  validate_4x4_matrix(matrix, "matrix")
  if (!is.character(source) || length(source) != 1L) stop("source must be a single character")
  if (!is.character(target) || length(target) != 1L) stop("target must be a single character")

  m <- methods::new("Affine3DMorphism",
           id = paste0("affine_", digest::digest(matrix)),
           source = source,
           target = target,
           kind = "affine3d",
           params = list(),
           inverse_params = list(),
           coverage = 1.0,
           cost = cost,
           method_tag = method_tag,
           matrix = matrix,
           hash = "",
           inverse_type = "exact",
           inverse_quality = 1.0,
           inverse_method = "analytic")
  m@hash <- morphism_hash(m)
  m
}

#' Create a 3D warp field morphism
#'
#' Represents a nonlinear warp field transformation.
#'
#' @param source Source domain hash
#' @param target Target domain hash
#' @param warp_path Path to displacement field file
#' @param warp_type One of "ants", "ants_h5", "fsl", "fsl_coef", "afni", "freesurfer", "dense"
#' @param inverse_path Path to inverse warp (optional)
#' @param def_type Deformation type: "relative" (displacement) or "absolute" (coordinates)
#' @param warp_method Interpolation method for warp field lookup: "linear" or "cubic"
#' @param cost Path cost (default 1.5)
#' @param method_tag Method tag (default "anatomical")
#' @return Warp3DMorphism object
#' @export
#' @examples
#' warp <- Warp3DMorphism("native", "mni", "path/to/warp.nii.gz")
Warp3DMorphism <- function(source, target, warp_path,
                           warp_type = c("ants", "ants_h5", "fsl", "fsl_coef", "afni", "freesurfer", "dense"),
                           inverse_path = "", def_type = NULL,
                           warp_method = c("linear", "cubic"),
                           cost = 1.5, method_tag = "anatomical") {
  warp_type <- match.arg(warp_type)
  if (is.null(def_type)) {
    # ANTs H5 and other displacement fields store relative offsets (not absolute coords).
    # AFNI 3dQwarp also stores DICOM-order displacement fields.
    def_type <- "relative"
  }
  def_type <- match.arg(def_type, c("relative", "absolute"))
  warp_method <- match.arg(warp_method)
  if (!is.character(source) || length(source) != 1L) stop("source must be a single character")
  if (!is.character(target) || length(target) != 1L) stop("target must be a single character")
  if (!nzchar(warp_path)) stop("warp_path must be provided")

  # Determine inverse properties
  if (nzchar(inverse_path)) {
    inv_type <- "approximate"
    if (warp_type == "ants") {
      inv_quality <- 0.95; inv_method <- "inverse_warp"
    } else if (warp_type %in% c("fsl", "fsl_coef")) {
      inv_quality <- 0.8; inv_method <- "invwarp"
    } else {
      inv_quality <- 0.85; inv_method <- "inverse_warp"
    }
  } else {
    inv_type <- "none"; inv_quality <- NA_real_; inv_method <- ""
  }

  m <- methods::new("Warp3DMorphism",
           id = paste0("warp_", basename(warp_path)),
           source = source,
           target = target,
           kind = "warp3d",
           params = list(warp_path = warp_path, warp_type = warp_type,
                         def_type = def_type, warp_method = warp_method),
           inverse_params = list(inverse_path = inverse_path),
           coverage = 1.0,
           cost = cost,
           method_tag = method_tag,
           warp_path = warp_path,
           warp_type = warp_type,
           inverse_path = inverse_path,
           hash = "",
           inverse_type = inv_type,
           inverse_quality = inv_quality,
           inverse_method = inv_method)
  m@hash <- morphism_hash(m)
  m
}

#' Create a volume-to-surface morphism
#'
#' Maps from volume domain to surface domain using various sampling strategies.
#'
#' @param source Volume domain hash
#' @param target Surface domain hash
#' @param method Sampling method: "trilinear", "ribbon", "mid_thickness", "nearest"
#' @param ribbon_inner Path to inner (white) surface for ribbon sampling
#' @param ribbon_outer Path to outer (pial) surface for ribbon sampling
#' @param ribbon_inner_coords Pre-loaded inner surface coordinates (optional)
#' @param ribbon_outer_coords Pre-loaded outer surface coordinates (optional)
#' @param n_ribbon_samples Number of samples along ribbon normal
#' @param mid_coords Pre-computed midpoint coordinates (optional)
#' @param cost Path cost (default 2.0)
#' @param method_tag Method tag (default "anatomical")
#' @return VolToSurfMorphism object
#' @export
VolToSurfMorphism <- function(source, target,
                              method = c("trilinear", "ribbon", "mid_thickness", "nearest"),
                              ribbon_inner = "", ribbon_outer = "",
                              ribbon_inner_coords = NULL, ribbon_outer_coords = NULL,
                              n_ribbon_samples = 6L, mid_coords = NULL,
                              cost = 2.0, method_tag = "anatomical") {
  method <- match.arg(method)
  if (method == "ribbon" && (!nzchar(ribbon_inner) || !nzchar(ribbon_outer))) {
    stop("Ribbon sampling requires ribbon_inner and ribbon_outer surfaces")
  }
  if (!is.character(source) || length(source) != 1L) stop("source must be a single character")
  if (!is.character(target) || length(target) != 1L) stop("target must be a single character")

  m <- methods::new("VolToSurfMorphism",
           id = paste0("v2s_", source, "_", target),
           source = source,
           target = target,
           kind = "vol2surf",
           params = list(method = method,
                         ribbon_inner = ribbon_inner,
                         ribbon_outer = ribbon_outer,
                         ribbon_inner_coords = ribbon_inner_coords,
                         ribbon_outer_coords = ribbon_outer_coords,
                         n_ribbon_samples = n_ribbon_samples,
                         mid_coords = mid_coords),
           inverse_params = list(),
           coverage = 1.0,
           cost = cost,
           method_tag = method_tag,
           method = method,
           ribbon_inner = ribbon_inner,
           ribbon_outer = ribbon_outer,
           n_ribbon_samples = as.integer(n_ribbon_samples),
           hash = "",
           inverse_type = "adjoint",
           inverse_quality = 1.0,
           inverse_method = "backprojection")
  m@hash <- morphism_hash(m)
  m
}

#' Create a surface-to-surface morphism
#'
#' Maps between surface meshes via spherical registration.
#'
#' @param source Source surface domain hash
#' @param target Target surface domain hash
#' @param method Registration method: "sphere", "area", "sulc"
#' @param source_sphere Path to source sphere surface
#' @param target_sphere Path to target sphere surface
#' @param mapping How to map target points onto the source mesh.
#'   \code{"nearest"} maps by nearest target vertex (index correspondence);
#'   \code{"barycentric"} maps using barycentric weights on target faces.
#' @param source_vertices Optional V x 3 numeric matrix of source mesh vertices.
#'   Required for \code{mapping="nearest"} and \code{mapping="barycentric"}.
#' @param target_vertices Optional V x 3 numeric matrix of target mesh vertices.
#'   Required for \code{mapping="nearest"} and \code{mapping="barycentric"}.
#' @param faces Optional F x 3 integer matrix of triangle vertex indices (0- or 1-based).
#'   Required for \code{mapping="barycentric"}.
#' @param cost Path cost (default 1.0)
#' @param method_tag Method tag (default "anatomical")
#' @return SurfToSurfMorphism object
#' @export
SurfToSurfMorphism <- function(source, target,
                               method = c("sphere", "area", "sulc"),
                               source_sphere = "", target_sphere = "",
                               mapping = c("nearest", "barycentric"),
                               source_vertices = NULL,
                               target_vertices = NULL,
                               faces = NULL,
                               cost = 1.0, method_tag = "anatomical") {
  method <- match.arg(method)
  mapping <- match.arg(mapping)
  if (!is.character(source) || length(source) != 1L) stop("source must be a single character")
  if (!is.character(target) || length(target) != 1L) stop("target must be a single character")

  # Accept neurosurf::SurfaceGeometry and other surface-like objects.
  if (!is.null(source_vertices) && !is.matrix(source_vertices)) {
    src_ref <- coerce_surface_reference(source_vertices, require_faces = FALSE)
    source_vertices <- src_ref$vertices
  }
  if (!is.null(target_vertices) && !is.matrix(target_vertices)) {
    tgt_ref <- coerce_surface_reference(
      target_vertices,
      require_faces = identical(mapping, "barycentric")
    )
    target_vertices <- tgt_ref$vertices
    if (is.null(faces) && !is.null(tgt_ref$faces) && nrow(tgt_ref$faces) > 0L) {
      faces <- tgt_ref$faces
    }
  }

  have_src <- !is.null(source_vertices)
  have_tgt <- !is.null(target_vertices)
  have_vertices <- have_src || have_tgt
  if (have_vertices && !(have_src && have_tgt)) {
    stop("If supplying vertices, provide both source_vertices and target_vertices.")
  }

  if (have_src) {
    if (!is.matrix(source_vertices) || ncol(source_vertices) != 3) {
      stop("source_vertices must be a matrix with 3 columns")
    }
    if (!is.numeric(source_vertices) || nrow(source_vertices) < 1L) {
      stop("source_vertices must be a non-empty numeric matrix")
    }
    if (!all(is.finite(source_vertices))) {
      stop("source_vertices must contain only finite values")
    }
  }

  if (have_tgt) {
    if (!is.matrix(target_vertices) || ncol(target_vertices) != 3) {
      stop("target_vertices must be a matrix with 3 columns")
    }
    if (!is.numeric(target_vertices) || nrow(target_vertices) < 1L) {
      stop("target_vertices must be a non-empty numeric matrix")
    }
    if (!all(is.finite(target_vertices))) {
      stop("target_vertices must contain only finite values")
    }
  }

  if (have_src && have_tgt && nrow(source_vertices) != nrow(target_vertices)) {
    stop("source_vertices and target_vertices must have the same number of rows (vertex correspondence by index).")
  }

  if (identical(mapping, "barycentric") && is.null(faces)) {
    stop("SurfToSurfMorphism mapping='barycentric' requires faces")
  }

  if (!is.null(faces)) {
    if (!is.matrix(faces) || ncol(faces) != 3) {
      stop("faces must be a matrix with 3 columns (triangles)")
    }
    if (!is.numeric(faces)) {
      stop("faces must be numeric/integer vertex indices")
    }
    if (!all(is.finite(faces))) {
      stop("faces must contain only finite values")
    }
    if (any(abs(faces - round(faces)) > 1e-8)) {
      stop("faces must contain integer vertex indices")
    }
    faces_int <- matrix(as.integer(round(faces)), ncol = 3)
    if (any(is.na(faces_int))) stop("faces must not contain NA")
    if (any(faces_int < 0L)) stop("faces must not contain negative indices")
    if (any(apply(faces_int, 1L, function(x) length(unique(x)) != 3L))) {
      stop("faces rows must reference 3 distinct vertices")
    }
    if (!is.null(target_vertices)) {
      nv <- nrow(target_vertices)
      if (min(faces_int) == 0L) {
        if (max(faces_int) >= nv) stop("0-based faces indices must be in [0, n_vertices-1]")
      } else if (min(faces_int) >= 1L) {
        if (max(faces_int) > nv) stop("1-based faces indices must be in [1, n_vertices]")
      } else {
        stop("faces must be 0-based (min==0) or 1-based (min>=1)")
      }
    }
    faces <- faces_int
  }

  faces0 <- NULL
  if (!is.null(faces)) {
    faces0 <- faces
    storage.mode(faces0) <- "integer"
    if (min(faces0, na.rm = TRUE) >= 1L) {
      faces0 <- faces0 - 1L
    }
  }

  m <- methods::new("SurfToSurfMorphism",
           id = paste0("s2s_", source, "_", target),
           source = source,
           target = target,
           kind = "surf2surf",
           params = list(
             method = method,
             mapping = mapping,
             source_vertices = source_vertices,
             target_vertices = target_vertices,
             faces0 = faces0
           ),
           inverse_params = list(),
           coverage = 1.0,
           cost = cost,
           method_tag = method_tag,
           method = method,
           source_sphere = source_sphere,
           target_sphere = target_sphere,
           hash = "",
           inverse_type = "none",
           inverse_quality = NA_real_,
           inverse_method = "")
  m@hash <- morphism_hash(m)
  m
}

# =============================================================================
# SHOW METHODS
# =============================================================================

#' @rdname Morphism-class
#' @export
setMethod("show", "Morphism", function(object) {
  cls <- class(object)[1]
  cat(sprintf("<%s %s | %s -> %s | kind=%s | method=%s | inverse=%s | cost=%.3f>\n",
              cls, object@id, object@source, object@target, object@kind, object@method_tag,
              object@inverse_type, object@cost))
})

# =============================================================================
# ACCESSOR METHODS
# =============================================================================

#' @rdname source_of
#' @export
setMethod("source_of", "Morphism", function(object) object@source)

#' @rdname target_of
#' @export
setMethod("target_of", "Morphism", function(object) object@target)

# Compatibility aliases
#' @rdname source_of
#' @export
setMethod("source_domain", "Morphism", function(object) source_of(object))
#' @rdname target_of
#' @export
setMethod("target_domain", "Morphism", function(object) target_of(object))

#' @rdname source_of
#' @export
setMethod("source_of", "MorphismPath", function(object) object@source)

#' @rdname target_of
#' @export
setMethod("target_of", "MorphismPath", function(object) object@target)

# =============================================================================
# TRANSFORM METHODS (pullback)
# =============================================================================

#' @rdname transform
#' @export
setMethod("transform", "IdentityMorphism", function(morphism, coords) {
  coords
})

#' @rdname transform
#' @export
setMethod("transform", "Affine3DMorphism", function(morphism, coords) {
  apply_affine(coords, morphism@matrix)
})

#' @rdname transform
#' @export
setMethod("transform", "Warp3DMorphism", function(morphism, coords) {
  warp_transform_coords(morphism, coords)
})

#' @rdname transform
#' @export
setMethod("transform", "VolToSurfMorphism", function(morphism, coords) {
  # Surface coords are already in world; return as-is or use precomputed midpoints
  if (!is.null(morphism@params$mid_coords)) {
    return(morphism@params$mid_coords)
  }
  coords
})

#' @rdname transform
#' @export
setMethod("transform", "SurfToSurfMorphism", function(morphism, coords) {
  src_v <- morphism@params$source_vertices
  tgt_v <- morphism@params$target_vertices
  faces0 <- morphism@params$faces0
  mapping <- morphism@params$mapping %||% "nearest"

  if (is.null(src_v) || is.null(tgt_v)) {
    stop("SurfToSurfMorphism requires source_vertices and target_vertices (matrices Vx3) in the constructor.")
  }
  if (!is.matrix(coords) || ncol(coords) != 3 || !is.numeric(coords)) {
    stop("coords must be a numeric matrix with 3 columns")
  }
  if (!all(is.finite(coords))) stop("coords must contain only finite values")

  if (identical(mapping, "nearest")) {
    idx <- cpp_nearest_vertex(coords, tgt_v)
    return(src_v[idx, , drop = FALSE])
  }

  if (is.null(faces0)) {
    stop("SurfToSurfMorphism mapping='barycentric' requires faces")
  }

  w <- cpp_barycentric_weights(coords, tgt_v, faces0)
  out <- matrix(0, nrow = nrow(coords), ncol = 3)
  if (length(w$rows) == 0) return(out)

  for (k in 1:3) {
    rs <- rowsum(src_v[w$cols, k, drop = TRUE] * w$vals,
                 group = w$rows, reorder = FALSE)
    out[as.integer(rownames(rs)), k] <- rs[, 1]
  }
  out
})

#' @rdname transform
#' @export
setMethod("transform", "MorphismPath", function(morphism, coords) {
  transform_path(morphism@morphisms, coords)
})

# Compatibility alias
#' @rdname transform_coords
#' @export
setMethod("transform_coords", "Morphism", function(object, coords) {
  transform(object, coords)
})

#' @rdname transform_coords
#' @export
setMethod("transform_coords", "MorphismPath", function(object, coords) {
  transform(object, coords)
})

# =============================================================================
# COMPOSE METHODS
# =============================================================================

#' @rdname compose
#' @export
setMethod("compose", signature("Morphism", "Morphism"), function(f, g) {
  kind_f <- morphism_kind(f)
  kind_g <- morphism_kind(g)

  # Identity shortcuts
  if (kind_f == "identity") return(g)
  if (kind_g == "identity") return(f)

  # Domain compatibility check
  if (target_of(f) != source_of(g)) {
    stop("Cannot compose: target of first must equal source of second")
  }

  # Affine fusion: compose into single affine
  if (kind_f == "affine3d" && kind_g == "affine3d") {
    composed_mat <- f@matrix %*% g@matrix
    return(Affine3DMorphism(
      source = source_of(f),
      target = target_of(g),
      matrix = composed_mat,
      cost = f@cost + g@cost,
      method_tag = f@method_tag
    ))
  }

  # General case: return MorphismPath
  methods::new("MorphismPath",
      morphisms = list(f, g),
      source = source_of(f),
      target = target_of(g))
})

#' @rdname compose
#' @export
setMethod("compose", signature("MorphismPath", "Morphism"), function(f, g) {
  if (target_of(f) != source_of(g)) {
    stop("Cannot compose: target of path must equal source of morphism")
  }
  methods::new("MorphismPath",
      morphisms = c(f@morphisms, list(g)),
      source = source_of(f),
      target = target_of(g))
})

#' @rdname compose
#' @export
setMethod("compose", signature("Morphism", "MorphismPath"), function(f, g) {
  if (target_of(f) != source_of(g)) {
    stop("Cannot compose: target of morphism must equal source of path")
  }
  methods::new("MorphismPath",
      morphisms = c(list(f), g@morphisms),
      source = source_of(f),
      target = target_of(g))
})

#' @rdname compose
#' @export
setMethod("compose", signature("MorphismPath", "MorphismPath"), function(f, g) {
  if (target_of(f) != source_of(g)) {
    stop("Cannot compose: target of first path must equal source of second path")
  }
  methods::new("MorphismPath",
      morphisms = c(f@morphisms, g@morphisms),
      source = source_of(f),
      target = target_of(g))
})

# =============================================================================
# INVERT METHODS
# =============================================================================

#' @rdname invert
#' @export
setMethod("invert", "Morphism", function(object) {
  if (object@inverse_type %in% c("adjoint", "none")) {
    stop("No geometric inverse available for morphism (inverse_type=", object@inverse_type, ")")
  }
  stop("invert() not implemented for morphism type: ", class(object))
})

#' @rdname invert
#' @export
setMethod("invert", "IdentityMorphism", function(object) object)

#' @rdname invert
#' @export
setMethod("invert", "Affine3DMorphism", function(object) {
  inv <- invert_affine(object@matrix)
  Affine3DMorphism(
    source = target_of(object),
    target = source_of(object),
    matrix = inv,
    cost = object@cost,
    method_tag = object@method_tag
  )
})

#' @rdname invert
#' @export
setMethod("invert", "Warp3DMorphism", function(object) {
  if (object@inverse_type %in% c("adjoint", "none")) {
    stop("Warp3DMorphism inverse not available (inverse_type=", object@inverse_type, ")")
  }
  if (!nzchar(object@inverse_path)) {
    stop("Warp3DMorphism has no inverse_path; cannot invert")
  }
  Warp3DMorphism(
    source = target_of(object),
    target = source_of(object),
    warp_path = object@inverse_path,
    inverse_path = object@warp_path,
    warp_type = object@warp_type,
    cost = object@cost,
    method_tag = object@method_tag
  )
})

# =============================================================================
# ADJOINT METHODS
# =============================================================================

#' @rdname adjoint
#' @export
setMethod("adjoint", "Morphism", function(object, ...) {
  stop("Adjoint not implemented for morphism type: ", class(object)[1])
})

#' @rdname adjoint
#' @export
setMethod("adjoint", "IdentityMorphism", function(object, ...) {
  object
})

#' @rdname adjoint
#' @export
setMethod("adjoint", "Affine3DMorphism", function(object, ...) {
  invert(object)
})

#' @rdname adjoint
#' @export
setMethod("adjoint", "Warp3DMorphism", function(object, ...) {
  invert(object)
})

#' @rdname adjoint
#' @export
setMethod("adjoint", "VolToSurfMorphism", function(object, ...) {
  args <- list(...)
  grid <- args$grid %||% args$source_grid
  surface_coords <- args$surface_coords %||% NULL
  method <- args$method %||% "linear"
  outside <- args$outside %||% 0

  if (is.null(grid) || !inherits(grid, "Grid")) {
    stop("adjoint(VolToSurfMorphism) requires a Grid via grid= (or source_grid=).")
  }

  function(surface_values) {
    backproject_surface_to_volume(
      surface_values = surface_values,
      morphism = object,
      grid = grid,
      surface_coords = surface_coords,
      method = method,
      outside = outside
    )
  }
})

#' @rdname adjoint
#' @export
setMethod("adjoint", "SurfToSurfMorphism", function(object, ...) {
  stop("Adjoint is not implemented for SurfToSurfMorphism.")
})

# =============================================================================
# PATH VALIDATION
# =============================================================================

#' Validate a morphism path
#'
#' Checks that a path is composable (each morphism's target matches the next's source).
#' Raises an error if the path is invalid.
#'
#' @param path List of Morphism objects
#' @return Invisibly returns the path if valid
#' @export
#' @examples
#' aff1 <- Affine3DMorphism("a", "b", diag(4))
#' aff2 <- Affine3DMorphism("b", "c", diag(4))
#' validate_path(list(aff1, aff2))  # OK
validate_path <- function(path) {
  if (length(path) == 0) return(invisible(path))
  if (length(path) == 1) {
    if (!methods::is(path[[1]], "Morphism")) stop("Path elements must be Morphism objects")
    return(invisible(path))
  }

  for (i in seq_len(length(path) - 1)) {
    if (!methods::is(path[[i]], "Morphism") || !methods::is(path[[i + 1]], "Morphism")) {
      stop("Path elements must be Morphism objects")
    }
    if (target_of(path[[i]]) != source_of(path[[i + 1]])) {
      stop(sprintf("Path not composable at position %d: target '%s' != source '%s'",
                   i, target_of(path[[i]]), source_of(path[[i + 1]])))
    }
  }
  invisible(path)
}

#' Check if a morphism path is valid
#'
#' Returns TRUE if the path is composable, FALSE otherwise.
#'
#' @param path List of Morphism objects
#' @return Logical
#' @export
#' @examples
#' aff1 <- Affine3DMorphism("a", "b", diag(4))
#' aff2 <- Affine3DMorphism("b", "c", diag(4))
#' is_valid_path(list(aff1, aff2))  # TRUE
is_valid_path <- function(path) {
  tryCatch({
    validate_path(path)
    TRUE
  }, error = function(e) FALSE)
}
