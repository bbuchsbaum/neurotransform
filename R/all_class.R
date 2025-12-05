#' @title S4 Class Definitions for neurotransform
#' @name all_class
#' @description
#' S4 class definitions for the neurotransform package. This file defines the
#' morphism class hierarchy. neurotransform is a pure transform kernel; domain
#' and geometry abstractions live in higher-level packages.
#'
#' @section Class Hierarchy:
#' \preformatted{
#' Morphism (virtual)
#'   |-- IdentityMorphism
#'   |-- Affine3DMorphism
#'   |-- Warp3DMorphism
#'   |-- VolToSurfMorphism
#'   |-- SurfToSurfMorphism
#'
#' MorphismPath (lightweight list wrapper for composed paths)
#' }
#'
#' @section Direction Convention:
#' A morphism A -> B maps coordinates in B to coordinates in A (pullback semantics).
#' This is the natural convention for resampling: to get a value at target location,
#' you need to know where to sample from in the source.
NULL

# =============================================================================
# MORPHISM CLASSES
# =============================================================================

#' Virtual base class for morphisms
#'
#' A Morphism represents a coordinate transformation between domains.
#' Direction convention: a morphism A -> B maps coordinates in B to coordinates
#' in A (pullback semantics). Data flows opposite: from A to B.
#'
#' @slot id Character identifier
#' @slot source Character: source domain hash (opaque string)
#' @slot target Character: target domain hash (opaque string)
#' @slot kind Character morphism kind (identity, affine3d, warp3d, vol2surf, surf2surf)
#' @slot params List: forward params/metadata
#' @slot inverse_params List: inverse params/metadata
#' @slot method_tag Character: path selection tag (anatomical/functional/etc.)
#' @slot coverage Optional mask/metadata for valid region
#' @slot cost Numeric: cost/weight for graph traversal
#' @slot hash Stable hash of morphism
#' @slot inverse_type One of "exact", "provided", "approximate", "adjoint", "none"
#' @slot inverse_quality Numeric in range 0 to 1, optional confidence for inverse
#' @slot inverse_method Description/label for inverse strategy
#' @slot provenance List: metadata about the transform origin
#' @slot cache Environment for cached data
#' @param object A Morphism object (for show method)
#' @export
setClass(
  "Morphism",
  slots = c(
    id = "character",
    source = "character",
    target = "character",
    kind = "character",
    params = "list",
    inverse_params = "list",
    method_tag = "character",
    coverage = "ANY",
    cost = "numeric",
    hash = "character",
    inverse_type = "character",
    inverse_quality = "numeric",
    inverse_method = "character",
    provenance = "list",
    cache = "environment"
  ),
  contains = "VIRTUAL",
  prototype = list(
    id = "",
    source = "",
    target = "",
    kind = "",
    params = list(),
    inverse_params = list(),
    method_tag = "anatomical",
    coverage = NULL,
    cost = 1.0,
    hash = "",
    inverse_type = "none",
    inverse_quality = NA_real_,
    inverse_method = "",
    provenance = list(),
    cache = new.env(hash = TRUE, parent = emptyenv())
  )
)

#' Identity morphism
#'
#' The identity morphism maps a domain to itself. Required for categorical
#' structure and useful for composition.
#'
#' @export
setClass(
  "IdentityMorphism",
  contains = "Morphism"
)

#' 3D affine morphism
#'
#' Represents an affine transformation (rotation, translation, scaling, shear)
#' between volume spaces.
#'
#' @slot matrix 4x4 affine transformation matrix (maps target -> source)
#' @export
setClass(
  "Affine3DMorphism",
  slots = c(
    matrix = "matrix"
  ),
  contains = "Morphism",
  prototype = list(
    matrix = diag(4)
  )
)

#' 3D warp field morphism
#'
#' Represents a nonlinear warp field transformation (e.g., ANTs SyN, FSL FNIRT).
#' Stores metadata; the actual warp field is loaded on demand.
#'
#' @slot warp_path Path to the warp field file
#' @slot warp_type Type of warp ("ants", "ants_h5", "fsl", "afni", "freesurfer")
#' @slot inverse_path Path to inverse warp (if available)
#' @export
setClass(
  "Warp3DMorphism",
  slots = c(
    warp_path = "character",
    warp_type = "character",
    inverse_path = "character"
  ),
  contains = "Morphism",
  prototype = list(
    warp_path = "",
    warp_type = "ants",
    inverse_path = ""
  )
)

#' Volume-to-surface morphism
#'
#' Maps volume coordinates to surface vertices. Supports different sampling
#' strategies (ribbon, trilinear point sampling).
#'
#' @slot method Sampling method ("ribbon", "trilinear", "nearest", "mid_thickness")
#' @slot ribbon_inner Path to inner (white) surface for ribbon sampling
#' @slot ribbon_outer Path to outer (pial) surface for ribbon sampling
#' @slot n_ribbon_samples Number of samples along ribbon normal
#' @export
setClass(
  "VolToSurfMorphism",
  slots = c(
    method = "character",
    ribbon_inner = "character",
    ribbon_outer = "character",
    n_ribbon_samples = "integer"
  ),
  contains = "Morphism",
  prototype = list(
    method = "trilinear",
    ribbon_inner = "",
    ribbon_outer = "",
    n_ribbon_samples = 6L
  )
)

#' Surface-to-surface morphism
#'
#' Maps between surface meshes (e.g., native to fsaverage via spherical
#' registration).
#'
#' @slot method Registration method ("sphere", "area", "sulc")
#' @slot source_sphere Path to source sphere surface
#' @slot target_sphere Path to target sphere surface
#' @export
setClass(
  "SurfToSurfMorphism",
  slots = c(
    method = "character",
    source_sphere = "character",
    target_sphere = "character"
  ),
  contains = "Morphism",
  prototype = list(
    method = "sphere",
    source_sphere = "",
    target_sphere = ""
  )
)

# =============================================================================
# MORPHISM PATH CLASS
# =============================================================================

#' Morphism path (composed sequence)
#'
#' A lightweight wrapper around a list of morphisms representing a composed
#' path. When passed to transform(), applies the full sequence efficiently.
#' Order convention: path = list(f, g) means f first, then g (like g . f).
#'
#' @slot morphisms List of Morphism objects
#' @slot source Source domain hash (from first morphism)
#' @slot target Target domain hash (from last morphism)
#' @export
setClass(
  "MorphismPath",
  slots = c(
    morphisms = "list",
    source = "character",
    target = "character"
  ),
  prototype = list(
    morphisms = list(),
    source = "",
    target = ""
  )
)
