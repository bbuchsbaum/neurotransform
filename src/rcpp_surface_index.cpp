#include <RcppArmadillo.h>
#include "surface_index.h"

// [[Rcpp::export]]
SEXP cpp_surface_index(const Rcpp::NumericMatrix& vertices,const Rcpp::IntegerMatrix& faces) {
  Rcpp::XPtr<neurotransform::SurfaceIndex> pointer(
    new neurotransform::SurfaceIndex(Rcpp::as<arma::mat>(vertices),Rcpp::as<arma::imat>(faces)),
    true,Rf_install("neurotransform.surface_index.v1"));
  return pointer;
}

// [[Rcpp::export]]
bool cpp_surface_index_valid(SEXP pointer) {
  return TYPEOF(pointer)==EXTPTRSXP && R_ExternalPtrTag(pointer)==Rf_install("neurotransform.surface_index.v1") &&
    R_ExternalPtrAddr(pointer)!=nullptr;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_surface_index_sample(SEXP pointer,const Rcpp::NumericMatrix& coords,
                                            const Rcpp::NumericVector& data,bool closest=false) {
  if (!cpp_surface_index_valid(pointer)) Rcpp::stop("surface index is unavailable; rebuild from mesh geometry");
  Rcpp::XPtr<neurotransform::SurfaceIndex> index(pointer);
  return neurotransform::sample_surface(*index,coords,data,closest);
}
