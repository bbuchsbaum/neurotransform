#' @title Deformation fields and Jacobian maps
#' @description
#' Helpers to materialize morphisms on voxel grids for visualization or export,
#' keeping the core morphism API untouched.
NULL

setOldClass("DeformationField")

#' Create a deformation field from a morphism on a grid
#'
#' @param morphism A Morphism (affine or warp)
#' @param grid A Grid object (see grid_spec/grid_from_data)
#' @param with_jacobian Logical; attach Jacobian determinant volume
#' @param log If TRUE, store log-Jacobian
#' @return An array with dim = c(grid@dims, 3) and class "DeformationField"
#' @export
deformation_field <- function(morphism, grid, with_jacobian = TRUE, log = FALSE) {
  stopifnot(inherits(grid, "Grid"))
  coords_tgt <- grid_coords(grid)
  coords_src <- transform(morphism, coords_tgt)

  field <- array(NA_real_, dim = c(grid@dims, 3))
  field[,,,1] <- array(coords_src[, 1], dim = grid@dims)
  field[,,,2] <- array(coords_src[, 2], dim = grid@dims)
  field[,,,3] <- array(coords_src[, 3], dim = grid@dims)

  attr(field, "grid") <- grid
  attr(field, "source") <- source_of(morphism)
  attr(field, "target") <- target_of(morphism)

  if (isTRUE(with_jacobian)) {
    dets <- jacobian_det(morphism, coords_tgt, log = log)
    attr(field, "jacobian_det") <- array(dets, dim = grid@dims)
  }

  class(field) <- c("DeformationField", "array")
  field
}

#' Extract Jacobian determinant volume from a morphism on a grid
#' @export
jacobian_det_field <- function(morphism, grid, log = FALSE) {
  coords <- grid_coords(grid)
  dets <- jacobian_det(morphism, coords, log = log)
  array(dets, dim = grid@dims)
}

#' @export
setMethod("jacobian", signature("DeformationField", "missing"),
          function(morphism, coords, mode = c("pullback", "pushforward")) {
            attr(morphism, "jacobian_det")
          })
