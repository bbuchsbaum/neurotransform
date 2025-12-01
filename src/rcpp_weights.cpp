#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute trilinear interpolation weights
//'
//' Given voxel coordinates, compute weights for 8 neighboring voxels.
//'
//' @param vox_coords Numeric matrix (N x 3) of voxel coordinates
//' @param dims Volume dimensions
//' @return List with rows, cols, vals for sparse matrix construction
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_trilinear_weights(const Rcpp::NumericMatrix& vox_coords,
                                 const Rcpp::IntegerVector& dims) {
  int n = vox_coords.nrow();
  int nx = dims[0], ny = dims[1], nz = dims[2];
  std::vector<int> rows; rows.reserve(n*8);
  std::vector<int> cols; cols.reserve(n*8);
  std::vector<double> vals; vals.reserve(n*8);

#ifdef _OPENMP
  #pragma omp parallel
  {
    std::vector<int> local_rows; local_rows.reserve(1024);
    std::vector<int> local_cols; local_cols.reserve(1024);
    std::vector<double> local_vals; local_vals.reserve(1024);

    #pragma omp for nowait
    for (int i = 0; i < n; ++i) {
      double x = vox_coords(i,0);
      double y = vox_coords(i,1);
      double z = vox_coords(i,2);
      if (x < -0.5 || y < -0.5 || z < -0.5 || x > nx - 0.5 || y > ny - 0.5 || z > nz - 0.5) {
        continue;
      }
      int x0 = (int)std::floor(x);
      int y0 = (int)std::floor(y);
      int z0 = (int)std::floor(z);
      double fx = x - x0, fy = y - y0, fz = z - z0;
      double wx[2] = {1.0 - fx, fx};
      double wy[2] = {1.0 - fy, fy};
      double wz[2] = {1.0 - fz, fz};
      for (int iz=0; iz<2; ++iz) {
        for (int iy=0; iy<2; ++iy) {
          for (int ix=0; ix<2; ++ix) {
            double w = wx[ix]*wy[iy]*wz[iz];
            int xi = std::min(std::max(x0+ix,0), nx-1);
            int yi = std::min(std::max(y0+iy,0), ny-1);
            int zi = std::min(std::max(z0+iz,0), nz-1);
            int col = xi + yi*nx + zi*nx*ny; // 0-based linear
            local_rows.push_back(i+1); // 1-based for R
            local_cols.push_back(col+1);
            local_vals.push_back(w);
          }
        }
      }
    }
    #pragma omp critical
    {
      rows.insert(rows.end(), local_rows.begin(), local_rows.end());
      cols.insert(cols.end(), local_cols.begin(), local_cols.end());
      vals.insert(vals.end(), local_vals.begin(), local_vals.end());
    }
  }
#else
  for (int i = 0; i < n; ++i) {
    double x = vox_coords(i,0);
    double y = vox_coords(i,1);
    double z = vox_coords(i,2);
    if (x < -0.5 || y < -0.5 || z < -0.5 || x > nx - 0.5 || y > ny - 0.5 || z > nz - 0.5) {
      continue;
    }
    int x0 = (int)std::floor(x);
    int y0 = (int)std::floor(y);
    int z0 = (int)std::floor(z);
    double fx = x - x0, fy = y - y0, fz = z - z0;
    double wx[2] = {1.0 - fx, fx};
    double wy[2] = {1.0 - fy, fy};
    double wz[2] = {1.0 - fz, fz};
    for (int iz=0; iz<2; ++iz) {
      for (int iy=0; iy<2; ++iy) {
        for (int ix=0; ix<2; ++ix) {
          double w = wx[ix]*wy[iy]*wz[iz];
          int xi = std::min(std::max(x0+ix,0), nx-1);
          int yi = std::min(std::max(y0+iy,0), ny-1);
          int zi = std::min(std::max(z0+iz,0), nz-1);
          int col = xi + yi*nx + zi*nx*ny;
          rows.push_back(i+1);
          cols.push_back(col+1);
          vals.push_back(w);
        }
      }
    }
  }
#endif

  return Rcpp::List::create(
    Rcpp::Named("rows") = rows,
    Rcpp::Named("cols") = cols,
    Rcpp::Named("vals") = vals
  );
}
