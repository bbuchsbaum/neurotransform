#' @title Transform Path Application
#' @name transform_path
#' @description
#' Apply a morphism path to coordinates using optimized warp chain.
NULL

#' Apply a morphism path to coordinates (pullback)
#'
#' Given a path of morphisms from source to target and coordinates in the
#' target domain, returns the corresponding coordinates in the source domain.
#'
#' @section Pullback Direction:
#' Given a path `[f: A->B, g: B->C]`, transform_path applies the pullback:
#' \preformatted{
#' coords_in_C  --(g pullback)--> coords_in_B --(f pullback)--> coords_in_A
#' }
#' The path is applied from last to first (g then f).
#'
#' @param path List of Morphism objects ordered source -> target
#' @param target_coords Numeric matrix (N x 3) of target world coords
#' @return Numeric matrix (N x 3) of source world coords
#' @export
#' @examples
#' # Create affine morphisms
#' aff1 <- Affine3DMorphism("native", "mni", diag(4))
#' aff2 <- Affine3DMorphism("mni", "template", diag(4))
#' path <- list(aff1, aff2)
#'
#' # Transform template coords to native space
#' template_coords <- matrix(c(0, 0, 0, 10, 20, 30), ncol = 3, byrow = TRUE)
#' native_coords <- transform_path(path, template_coords)
transform_path <- function(path, target_coords) {
  if (length(path) == 0) return(target_coords)

  # Validate path is composable
  validate_path(path)

  # Build optimized chain and apply
  chain <- build_warp_chain(path)
  apply_warp_chain(chain, target_coords)
}
