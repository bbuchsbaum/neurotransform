#' @title Generic Function Definitions
#' @name all_generic
#' @description
#' Generic function definitions for the neurotransform package.
#' These define the core API for coordinate transforms.
NULL

# =============================================================================
# CORE TRANSFORM GENERICS
# =============================================================================

#' Transform coordinates through a morphism (pullback)
#'
#' Applies the coordinate pullback: given coordinates in the target domain,
#' returns corresponding coordinates in the source domain.
#'
#' @param morphism A Morphism object or MorphismPath
#' @param coords Matrix of coordinates (N x 3) in target domain
#' @return Matrix of coordinates (N x 3) in source domain
#' @export
#' @examples
#' # Create an affine morphism
#' aff <- Affine3DMorphism("source_hash", "target_hash", diag(4))
#' coords <- matrix(c(10, 20, 30), ncol = 3)
#' transform(aff, coords)
setGeneric("transform", function(morphism, coords) {
  standardGeneric("transform")
})

#' Compose two morphisms
#'
#' Returns the composition g . f (f applied first, then g).
#' The result is a MorphismPath that can be passed to transform().
#'
#' @param f First morphism (applied first)
#' @param g Second morphism (applied second)
#' @return MorphismPath representing the composition
#' @export
#' @examples
#' aff1 <- Affine3DMorphism("a", "b", diag(4))
#' aff2 <- Affine3DMorphism("b", "c", diag(4))
#' path <- compose(aff1, aff2)  # a -> b -> c
setGeneric("compose", function(f, g) standardGeneric("compose"))

#' Invert a morphism
#'
#' Returns the inverse morphism if available (inverse_type == "exact" or "approximate").
#' For morphisms with approximate or adjoint inverses, use adjoint() instead.
#'
#' @param object A Morphism object
#' @return Inverted morphism
#' @export
#' @examples
#' aff <- Affine3DMorphism("a", "b", diag(4))
#' inv <- invert(aff)  # b -> a
setGeneric("invert", function(object) standardGeneric("invert"))

#' Get adjoint (generalized inverse)
#'
#' Returns an adjoint morphism when defined.
#'
#' For invertible morphisms, \code{adjoint()} is defined and equals \code{invert()}.
#' For other morphisms, \code{adjoint()} may not be implemented and will error.
#'
#' @param object A Morphism object
#' @param ... Additional arguments
#' @return Adjoint morphism or error if not available
#' @export
setGeneric("adjoint", function(object, ...) {
  standardGeneric("adjoint")
})

# =============================================================================
# ACCESSOR GENERICS
# =============================================================================

#' Get source domain of a morphism
#'
#' @param object A Morphism object
#' @return Source domain hash (character)
#' @export
setGeneric("source_of", function(object) standardGeneric("source_of"))

#' Get target domain of a morphism
#'
#' @param object A Morphism object
#' @return Target domain hash (character)
#' @export
setGeneric("target_of", function(object) standardGeneric("target_of"))

# Compatibility domain generics
#' @rdname source_of
#' @export
setGeneric("source_domain", function(object) standardGeneric("source_domain"))
#' @rdname target_of
#' @export
setGeneric("target_domain", function(object) standardGeneric("target_domain"))

# =============================================================================
# COMPATIBILITY ALIASES (for neurofunctor)
# =============================================================================

#' Transform coordinates (compatibility alias)
#'
#' Alias for transform() to maintain neurofunctor compatibility.
#'
#' @param object A Morphism object
#' @param coords Matrix of coordinates (N x 3)
#' @return Transformed coordinates
#' @export
setGeneric("transform_coords", function(object, coords) {
  standardGeneric("transform_coords")
})
