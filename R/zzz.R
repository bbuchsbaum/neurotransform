#' @useDynLib neurotransform, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods setClass setGeneric setMethod new is slot slotNames
#' @importFrom methods validObject .hasSlot show
#' @importFrom digest digest
#' @importFrom stats cor median
#' @importFrom utils read.table
NULL

.onLoad <- function(libname, pkgname) {
  # Register default neuroim2 warp loader
  register_loader("neuroim2", load_warp_neuroim2)
  # Register ANTs H5 loader (requires hdf5r)
  register_loader("ants_h5", load_warp_ants_h5)
}
