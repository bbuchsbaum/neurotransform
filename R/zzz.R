#' @useDynLib neurotransform, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods setClass setGeneric setMethod new is slot slotNames
#' @importFrom methods validObject .hasSlot show
#' @importFrom digest digest
#' @importFrom stats cor median
#' @importFrom utils read.table
NULL

.onLoad <- function(libname, pkgname) {

  # Register default RNifti warp loader
  register_loader("rnifti", load_warp_rnifti)
}
