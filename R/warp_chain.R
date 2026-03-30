#' @title Warp Chain Utilities
#' @name warp_chain
#' @description
#' Helpers to compose morphism paths into efficient coordinate transforms.
#' Optimizations include:
#' - Affine fusion: consecutive affines are multiplied into one
#' - Affine batching: C++ applies chain of affines efficiently
#' - Warp composition: consecutive warps can be composed
#' @keywords internal
NULL

#' Build cache key for a morphism path
#'
#' @param path List of morphisms
#' @return Character hash key
#' @keywords internal
warp_chain_key <- function(path) {
  hashes <- vapply(path, morphism_hash, character(1))
  compute_hash("warp_chain", hashes)
}

#' Build a warp chain from a morphism path
#'
#' Analyzes and optimizes a morphism path for efficient application.
#'
#' @param path List of Morphism objects ordered source -> target
#' @return List describing the optimized chain
#' @keywords internal
build_warp_chain <- function(path) {
  if (length(path) == 0) {
    return(list(kind = "identity"))
  }

  # Check if all affine/identity
  all_affine <- all(vapply(path, is_linear_morphism, logical(1)))

  if (all_affine) {
    # Combine into single affine pullback: apply in reverse path order
    mat <- diag(4)
    for (m in rev(path)) {
      kind <- morphism_kind(m)
      if (kind == "affine3d") {
        mat <- m@matrix %*% mat
      }
    }
    return(list(kind = "affine", matrix = mat, path = path))
  }

  # Mixed path: preserve exact pullback order while batching consecutive affines.
  # Path morphisms are stored source -> target, so pullback evaluation applies them
  # in reverse order at target coordinates.
  segments <- list()
  current_affine <- diag(4)

  for (m in rev(path)) {
    kind <- morphism_kind(m)
    if (kind == "affine3d") {
      current_affine <- m@matrix %*% current_affine
    } else if (kind == "identity") {
      next
    } else {
      if (!all(current_affine == diag(4))) {
        segments <- c(segments, list(list(type = "affine", matrix = current_affine)))
        current_affine <- diag(4)
      }
      segments <- c(segments, list(list(type = "warp", morphism = m)))
    }
  }

  if (!all(current_affine == diag(4))) {
    segments <- c(segments, list(list(type = "affine", matrix = current_affine)))
  }

  list(kind = "mixed", segments = segments, path = path)
}

#' Apply a warp chain to coordinates
#'
#' Applies an optimized warp chain built by build_warp_chain().
#'
#' @param chain Chain object from build_warp_chain()
#' @param coords Numeric matrix (N x 3) of target world coords
#' @return Transformed coordinates in source space
#' @keywords internal
apply_warp_chain <- function(chain, coords) {
  if (chain$kind == "identity") return(coords)

  if (chain$kind == "affine") {
    return(cpp_apply_affine_chain(coords, list(chain$matrix)))
  }

  if (chain$kind == "mixed") {
    out <- coords
    pending_affines <- list()

    flush_affines <- function() {
      if (length(pending_affines) == 0) return(out)
      mats <- lapply(pending_affines, `[[`, "matrix")
      out <<- cpp_apply_affine_chain(out, mats)
      pending_affines <<- list()
    }

    for (seg in chain$segments) {
      if (seg$type == "affine") {
        pending_affines <- c(pending_affines, list(seg))
      } else {
        flush_affines()
        if (seg$type == "warp_comp") {
          cw <- compose_warps(seg$warpB, seg$warpA)
          out <- cpp_apply_warp_field(out, cw$array, cw$dim, cw$world_to_vox)
        } else if (seg$type == "warp") {
          out <- warp_transform_coords(seg$morphism, out)
        }
      }
    }
    flush_affines()
    return(out)
  }

  # Fallback: sequential application
  out <- coords
  for (m in rev(chain$path)) {
    out <- transform(m, out)
  }
  out
}
