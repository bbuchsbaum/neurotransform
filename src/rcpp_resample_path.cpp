#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Forward declarations from other compilation units
Rcpp::NumericMatrix cpp_apply_warp_field(const Rcpp::NumericMatrix& coords,
                                         const Rcpp::NumericVector& field,
                                         const Rcpp::IntegerVector& dim,
                                         const Rcpp::NumericMatrix& voxel_to_world);
Rcpp::NumericMatrix cpp_apply_warp_field_cubic(const Rcpp::NumericMatrix& coords,
                                               const Rcpp::NumericVector& field,
                                               const Rcpp::IntegerVector& dim,
                                               const Rcpp::NumericMatrix& voxel_to_world);

// Flatten a path of affines into a single matrix (3x4) applied to coords (n x 3)
inline void apply_affine_inplace(Rcpp::NumericMatrix& coords, const Rcpp::NumericMatrix& A) {
  int n = coords.nrow();
  for (int i = 0; i < n; ++i) {
    double x = coords(i,0), y = coords(i,1), z = coords(i,2);
    coords(i,0) = A(0,0)*x + A(0,1)*y + A(0,2)*z + A(0,3);
    coords(i,1) = A(1,0)*x + A(1,1)*y + A(1,2)*z + A(1,3);
    coords(i,2) = A(2,0)*x + A(2,1)*y + A(2,2)*z + A(2,3);
  }
}

//' Apply a chain of affines and then a warp to coordinates (target -> source)
//'
//' @param coords Numeric matrix (n x 3) target world coords
//' @param affines List of 4x4 matrices (applied in order)
//' @param field Numeric vector displacement field (flattened 3 x X x Y x Z)
//' @param dim Integer vector (X,Y,Z)
//' @param world_to_vox 4x4 world-to-voxel matrix for the warp
//' @param cubic Logical; if TRUE use cubic warp sampling
//' @return Numeric matrix (n x 3) source world coords
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_path_apply_affine_warp(const Rcpp::NumericMatrix& coords,
                                               const Rcpp::List& affines,
                                               const Rcpp::NumericVector& field,
                                               const Rcpp::IntegerVector& dim,
                                               const Rcpp::NumericMatrix& world_to_vox,
                                               bool cubic) {
  Rcpp::NumericMatrix out = Rcpp::clone(coords);
  // Apply affines in order
  for (int i = 0; i < affines.size(); ++i) {
    Rcpp::NumericMatrix A = affines[i];
    if (A.nrow() != 4 || A.ncol() != 4) Rcpp::stop("affine must be 4x4");
    apply_affine_inplace(out, A);
  }
  if (cubic) {
    return cpp_apply_warp_field_cubic(out, field, dim, world_to_vox);
  }
  return cpp_apply_warp_field(out, field, dim, world_to_vox);
}

//' Apply affines -> warp -> affines to coordinates (target -> source)
//'
//' @param coords n x 3 target world coords
//' @param affines_pre list of 4x4 matrices applied before warp
//' @param field displacement field
//' @param dim field dims
//' @param world_to_vox warp world_to_vox
//' @param cubic use cubic sampling
//' @param affines_post list of 4x4 matrices applied after warp
//' @return n x 3 source world coords
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_path_apply_affine_warp_affine(const Rcpp::NumericMatrix& coords,
                                                      const Rcpp::List& affines_pre,
                                                      const Rcpp::NumericVector& field,
                                                      const Rcpp::IntegerVector& dim,
                                                      const Rcpp::NumericMatrix& world_to_vox,
                                                      bool cubic,
                                                      const Rcpp::List& affines_post) {
  Rcpp::NumericMatrix out = Rcpp::clone(coords);
  for (int i = 0; i < affines_pre.size(); ++i) {
    Rcpp::NumericMatrix A = affines_pre[i];
    if (A.nrow() != 4 || A.ncol() != 4) Rcpp::stop("affine must be 4x4");
    apply_affine_inplace(out, A);
  }
  if (field.size() > 0) {
    if (cubic) {
      out = cpp_apply_warp_field_cubic(out, field, dim, world_to_vox);
    } else {
      out = cpp_apply_warp_field(out, field, dim, world_to_vox);
    }
  }
  for (int i = 0; i < affines_post.size(); ++i) {
    Rcpp::NumericMatrix A = affines_post[i];
    if (A.nrow() != 4 || A.ncol() != 4) Rcpp::stop("affine must be 4x4");
    apply_affine_inplace(out, A);
  }
  return out;
}

//' Apply arbitrary steps (affine/warp) to coordinates (target -> source)
//'
//' @param coords n x 3 target world coords
//' @param steps list of steps with kind \"affine\" or \"warp\" and associated data
//' @return n x 3 transformed coords
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_path_apply_steps(const Rcpp::NumericMatrix& coords,
                                         const Rcpp::List& steps) {
  Rcpp::NumericMatrix out = Rcpp::clone(coords);
  int nsteps = steps.size();
  for (int s = 0; s < nsteps; ++s) {
    Rcpp::List st = steps[s];
    std::string kind = Rcpp::as<std::string>(st["kind"]);
    if (kind == "affine") {
      Rcpp::NumericMatrix A = st["matrix"];
      if (A.nrow() != 4 || A.ncol() != 4) Rcpp::stop("affine must be 4x4");
      apply_affine_inplace(out, A);
    } else if (kind == "warp") {
      Rcpp::NumericVector field = st["field"];
      Rcpp::IntegerVector dim = st["dim"];
      Rcpp::NumericMatrix w2v = st["world_to_vox"];
      std::string method = Rcpp::as<std::string>(st["method"]);
      if (method == "cubic") {
        out = cpp_apply_warp_field_cubic(out, field, dim, w2v);
      } else {
        out = cpp_apply_warp_field(out, field, dim, w2v);
      }
    } else {
      Rcpp::stop("Unknown step kind");
    }
  }
  return out;
}
