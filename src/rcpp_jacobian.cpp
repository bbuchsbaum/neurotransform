#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// Cubic interpolation helpers
inline double cubic_weight(double t) {
  double at = std::abs(t);
  if (at <= 1.0) {
    return (1.5*at - 2.5)*at*at + 1.0;
  } else if (at < 2.0) {
    return ((-0.5*at + 2.5)*at - 4.0)*at + 2.0;
  }
  return 0.0;
}

inline double cubic_interp_1d(const double* vals, double t) {
  return vals[0]*cubic_weight(t+1) + vals[1]*cubic_weight(t) +
         vals[2]*cubic_weight(t-1) + vals[3]*cubic_weight(t-2);
}

// Trilinear sample displacement at voxel coords (0-based)
inline bool sample_disp_local(const double* field, int nx, int ny, int nz,
                              double x, double y, double z, double* out) {
  if (x < -0.5 || y < -0.5 || z < -0.5 || x > nx - 0.5 || y > ny - 0.5 || z > nz - 0.5) return false;
  int x0 = (int)floor(x), y0 = (int)floor(y), z0 = (int)floor(z);
  double fx = x - x0, fy = y - y0, fz = z - z0;
  double wx[2] = {1.0 - fx, fx};
  double wy[2] = {1.0 - fy, fy};
  double wz[2] = {1.0 - fz, fz};

  double acc[3] = {0,0,0};
  for (int iz=0; iz<2; ++iz) for (int iy=0; iy<2; ++iy) for (int ix=0; ix<2; ++ix) {
    double w = wx[ix]*wy[iy]*wz[iz];
    int xi = std::min(std::max(x0+ix,0), nx-1);
    int yi = std::min(std::max(y0+iy,0), ny-1);
    int zi = std::min(std::max(z0+iz,0), nz-1);
    size_t base = ((size_t(zi)*ny + yi)*nx + xi)*3;
    acc[0] += w * field[base];
    acc[1] += w * field[base + 1];
    acc[2] += w * field[base + 2];
  }
  out[0] = acc[0]; out[1] = acc[1]; out[2] = acc[2];
  return true;
}

inline bool sample_disp_local_cubic(const double* field, int nx, int ny, int nz,
                                    double x, double y, double z, double* out) {
  if (x < -0.5 || y < -0.5 || z < -0.5 || x > nx - 0.5 || y > ny - 0.5 || z > nz - 0.5) return false;
  if (x < 0.5 || y < 0.5 || z < 0.5 || x > nx - 1.5 || y > ny - 1.5 || z > nz - 1.5) {
    return sample_disp_local(field, nx, ny, nz, x, y, z, out);
  }
  int x0 = (int)floor(x) - 1;
  int y0 = (int)floor(y) - 1;
  int z0 = (int)floor(z) - 1;
  double fx = x - (x0 + 1);
  double fy = y - (y0 + 1);
  double fz = z - (z0 + 1);

  double z_slices[4][3];
  for (int iz = 0; iz < 4; ++iz) {
    double y_rows[4][3];
    for (int iy = 0; iy < 4; ++iy) {
      double x_vals[4][3];
      for (int ix = 0; ix < 4; ++ix) {
        int xi = std::min(std::max(x0 + ix, 0), nx - 1);
        int yi = std::min(std::max(y0 + iy, 0), ny - 1);
        int zi = std::min(std::max(z0 + iz, 0), nz - 1);
        size_t base = ((size_t(zi)*ny + yi)*nx + xi)*3;
        x_vals[ix][0] = field[base];
        x_vals[ix][1] = field[base + 1];
        x_vals[ix][2] = field[base + 2];
      }
      for (int c = 0; c < 3; ++c) {
        y_rows[iy][c] = cubic_interp_1d(&x_vals[0][c], fx);
      }
    }
    for (int c = 0; c < 3; ++c) {
      double tmp[4] = { y_rows[0][c], y_rows[1][c], y_rows[2][c], y_rows[3][c] };
      z_slices[iz][c] = cubic_interp_1d(tmp, fy);
    }
  }
  for (int c = 0; c < 3; ++c) {
    double tmp[4] = { z_slices[0][c], z_slices[1][c], z_slices[2][c], z_slices[3][c] };
    out[c] = cubic_interp_1d(tmp, fz);
  }
  return true;
}

//' Compute Jacobian matrices for a warp field at given coordinates
//'
//' Uses finite differences on the displacement field to compute 3x3 Jacobian
//' matrices at each query point.
//'
//' @param coords Numeric matrix (N x 3) of world coordinates
//' @param field Numeric vector of displacement values (flattened X x Y x Z x 3)
//' @param dim Integer vector of dimensions (X, Y, Z)
//' @param world_to_vox 4x4 world-to-voxel transform
//' @param vox_to_world 4x4 voxel-to-world transform (for scaling)
//' @return Numeric array (N x 3 x 3) of Jacobian matrices
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_warp_jacobian(const Rcpp::NumericMatrix& coords,
                                       const Rcpp::NumericVector& field,
                                       const Rcpp::IntegerVector& dim,
                                       const Rcpp::NumericMatrix& world_to_vox,
                                       const Rcpp::NumericMatrix& vox_to_world) {
  int n = coords.nrow();
  int nx = dim[0], ny = dim[1], nz = dim[2];
  const double* fdata = field.begin();

  // Output: n x 3 x 3 array
  Rcpp::NumericVector out(n * 9);
  out.attr("dim") = Rcpp::IntegerVector::create(n, 3, 3);

  // Finite difference step in voxel coords
  double h = 0.5;

  // Get voxel spacing from vox_to_world for proper scaling
  double vox_spacing[3];
  vox_spacing[0] = sqrt(vox_to_world(0,0)*vox_to_world(0,0) +
                        vox_to_world(1,0)*vox_to_world(1,0) +
                        vox_to_world(2,0)*vox_to_world(2,0));
  vox_spacing[1] = sqrt(vox_to_world(0,1)*vox_to_world(0,1) +
                        vox_to_world(1,1)*vox_to_world(1,1) +
                        vox_to_world(2,1)*vox_to_world(2,1));
  vox_spacing[2] = sqrt(vox_to_world(0,2)*vox_to_world(0,2) +
                        vox_to_world(1,2)*vox_to_world(1,2) +
                        vox_to_world(2,2)*vox_to_world(2,2));

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double wx = coords(i,0), wy = coords(i,1), wz = coords(i,2);

    // Convert world to voxel
    double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy + world_to_vox(0,2)*wz + world_to_vox(0,3);
    double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy + world_to_vox(1,2)*wz + world_to_vox(1,3);
    double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy + world_to_vox(2,2)*wz + world_to_vox(2,3);

    // Initialize Jacobian as identity (for the base coordinate)
    // J_ij = d(source_i) / d(target_j)
    // source = target + disp, so J = I + d(disp)/d(target)
    double J[9] = {1,0,0, 0,1,0, 0,0,1};

    // Compute partial derivatives via central differences
    double disp_plus[3], disp_minus[3];

    // d/dx
    if (sample_disp_local(fdata, nx, ny, nz, vx+h, vy, vz, disp_plus) &&
        sample_disp_local(fdata, nx, ny, nz, vx-h, vy, vz, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[0]);
      J[0] += (disp_plus[0] - disp_minus[0]) * scale;  // d(disp_x)/d(x)
      J[3] += (disp_plus[1] - disp_minus[1]) * scale;  // d(disp_y)/d(x)
      J[6] += (disp_plus[2] - disp_minus[2]) * scale;  // d(disp_z)/d(x)
    }

    // d/dy
    if (sample_disp_local(fdata, nx, ny, nz, vx, vy+h, vz, disp_plus) &&
        sample_disp_local(fdata, nx, ny, nz, vx, vy-h, vz, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[1]);
      J[1] += (disp_plus[0] - disp_minus[0]) * scale;
      J[4] += (disp_plus[1] - disp_minus[1]) * scale;
      J[7] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    // d/dz
    if (sample_disp_local(fdata, nx, ny, nz, vx, vy, vz+h, disp_plus) &&
        sample_disp_local(fdata, nx, ny, nz, vx, vy, vz-h, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[2]);
      J[2] += (disp_plus[0] - disp_minus[0]) * scale;
      J[5] += (disp_plus[1] - disp_minus[1]) * scale;
      J[8] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    // Store in column-major order for R (n x 3 x 3)
    // out[i, j, k] = out[i + n*j + n*3*k]
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        out[i + n*j + n*3*k] = J[j*3 + k];
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_warp_jacobian_cubic(const Rcpp::NumericMatrix& coords,
                                            const Rcpp::NumericVector& field,
                                            const Rcpp::IntegerVector& dim,
                                            const Rcpp::NumericMatrix& world_to_vox,
                                            const Rcpp::NumericMatrix& vox_to_world) {
  int n = coords.nrow();
  int nx = dim[0], ny = dim[1], nz = dim[2];
  const double* fdata = field.begin();

  Rcpp::NumericVector out(n * 9);
  out.attr("dim") = Rcpp::IntegerVector::create(n, 3, 3);
  double h = 0.5;

  double vox_spacing[3];
  vox_spacing[0] = sqrt(vox_to_world(0,0)*vox_to_world(0,0) +
                        vox_to_world(1,0)*vox_to_world(1,0) +
                        vox_to_world(2,0)*vox_to_world(2,0));
  vox_spacing[1] = sqrt(vox_to_world(0,1)*vox_to_world(0,1) +
                        vox_to_world(1,1)*vox_to_world(1,1) +
                        vox_to_world(2,1)*vox_to_world(2,1));
  vox_spacing[2] = sqrt(vox_to_world(0,2)*vox_to_world(0,2) +
                        vox_to_world(1,2)*vox_to_world(1,2) +
                        vox_to_world(2,2)*vox_to_world(2,2));

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double wx = coords(i,0), wy = coords(i,1), wz = coords(i,2);
    double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy + world_to_vox(0,2)*wz + world_to_vox(0,3);
    double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy + world_to_vox(1,2)*wz + world_to_vox(1,3);
    double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy + world_to_vox(2,2)*wz + world_to_vox(2,3);

    double J[9] = {1,0,0, 0,1,0, 0,0,1};
    double disp_plus[3], disp_minus[3];

    if (sample_disp_local_cubic(fdata, nx, ny, nz, vx+h, vy, vz, disp_plus) &&
        sample_disp_local_cubic(fdata, nx, ny, nz, vx-h, vy, vz, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[0]);
      J[0] += (disp_plus[0] - disp_minus[0]) * scale;
      J[3] += (disp_plus[1] - disp_minus[1]) * scale;
      J[6] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    if (sample_disp_local_cubic(fdata, nx, ny, nz, vx, vy+h, vz, disp_plus) &&
        sample_disp_local_cubic(fdata, nx, ny, nz, vx, vy-h, vz, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[1]);
      J[1] += (disp_plus[0] - disp_minus[0]) * scale;
      J[4] += (disp_plus[1] - disp_minus[1]) * scale;
      J[7] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    if (sample_disp_local_cubic(fdata, nx, ny, nz, vx, vy, vz+h, disp_plus) &&
        sample_disp_local_cubic(fdata, nx, ny, nz, vx, vy, vz-h, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[2]);
      J[2] += (disp_plus[0] - disp_minus[0]) * scale;
      J[5] += (disp_plus[1] - disp_minus[1]) * scale;
      J[8] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        out[i + n*j + n*3*k] = J[j*3 + k];
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_warp_jacobian_det_cubic(const Rcpp::NumericMatrix& coords,
                                                const Rcpp::NumericVector& field,
                                                const Rcpp::IntegerVector& dim,
                                                const Rcpp::NumericMatrix& world_to_vox,
                                                const Rcpp::NumericMatrix& vox_to_world) {
  int n = coords.nrow();
  int nx = dim[0], ny = dim[1], nz = dim[2];
  const double* fdata = field.begin();

  Rcpp::NumericVector out(n);
  double h = 0.5;

  double vox_spacing[3];
  vox_spacing[0] = sqrt(vox_to_world(0,0)*vox_to_world(0,0) +
                        vox_to_world(1,0)*vox_to_world(1,0) +
                        vox_to_world(2,0)*vox_to_world(2,0));
  vox_spacing[1] = sqrt(vox_to_world(0,1)*vox_to_world(0,1) +
                        vox_to_world(1,1)*vox_to_world(1,1) +
                        vox_to_world(2,1)*vox_to_world(2,1));
  vox_spacing[2] = sqrt(vox_to_world(0,2)*vox_to_world(0,2) +
                        vox_to_world(1,2)*vox_to_world(1,2) +
                        vox_to_world(2,2)*vox_to_world(2,2));

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double wx = coords(i,0), wy = coords(i,1), wz = coords(i,2);

    double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy + world_to_vox(0,2)*wz + world_to_vox(0,3);
    double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy + world_to_vox(1,2)*wz + world_to_vox(1,3);
    double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy + world_to_vox(2,2)*wz + world_to_vox(2,3);

    double J[9] = {1,0,0, 0,1,0, 0,0,1};
    double disp_plus[3], disp_minus[3];

    if (sample_disp_local_cubic(fdata, nx, ny, nz, vx+h, vy, vz, disp_plus) &&
        sample_disp_local_cubic(fdata, nx, ny, nz, vx-h, vy, vz, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[0]);
      J[0] += (disp_plus[0] - disp_minus[0]) * scale;
      J[3] += (disp_plus[1] - disp_minus[1]) * scale;
      J[6] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    if (sample_disp_local_cubic(fdata, nx, ny, nz, vx, vy+h, vz, disp_plus) &&
        sample_disp_local_cubic(fdata, nx, ny, nz, vx, vy-h, vz, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[1]);
      J[1] += (disp_plus[0] - disp_minus[0]) * scale;
      J[4] += (disp_plus[1] - disp_minus[1]) * scale;
      J[7] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    if (sample_disp_local_cubic(fdata, nx, ny, nz, vx, vy, vz+h, disp_plus) &&
        sample_disp_local_cubic(fdata, nx, ny, nz, vx, vy, vz-h, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[2]);
      J[2] += (disp_plus[0] - disp_minus[0]) * scale;
      J[5] += (disp_plus[1] - disp_minus[1]) * scale;
      J[8] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    out[i] = J[0]*(J[4]*J[8] - J[7]*J[5])
           - J[1]*(J[3]*J[8] - J[6]*J[5])
           + J[2]*(J[3]*J[7] - J[6]*J[4]);
  }

  return out;
}

//' Compute Jacobian determinants for a warp field at given coordinates
//'
//' Efficiently computes only the determinants without full Jacobian matrices.
//'
//' @param coords Numeric matrix (N x 3) of world coordinates
//' @param field Numeric vector of displacement values (flattened X x Y x Z x 3)
//' @param dim Integer vector of dimensions (X, Y, Z)
//' @param world_to_vox 4x4 world-to-voxel transform
//' @param vox_to_world 4x4 voxel-to-world transform (for scaling)
//' @return Numeric vector of N determinants
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_warp_jacobian_det(const Rcpp::NumericMatrix& coords,
                                           const Rcpp::NumericVector& field,
                                           const Rcpp::IntegerVector& dim,
                                           const Rcpp::NumericMatrix& world_to_vox,
                                           const Rcpp::NumericMatrix& vox_to_world) {
  int n = coords.nrow();
  int nx = dim[0], ny = dim[1], nz = dim[2];
  const double* fdata = field.begin();

  Rcpp::NumericVector out(n);

  double h = 0.5;

  double vox_spacing[3];
  vox_spacing[0] = sqrt(vox_to_world(0,0)*vox_to_world(0,0) +
                        vox_to_world(1,0)*vox_to_world(1,0) +
                        vox_to_world(2,0)*vox_to_world(2,0));
  vox_spacing[1] = sqrt(vox_to_world(0,1)*vox_to_world(0,1) +
                        vox_to_world(1,1)*vox_to_world(1,1) +
                        vox_to_world(2,1)*vox_to_world(2,1));
  vox_spacing[2] = sqrt(vox_to_world(0,2)*vox_to_world(0,2) +
                        vox_to_world(1,2)*vox_to_world(1,2) +
                        vox_to_world(2,2)*vox_to_world(2,2));

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double wx = coords(i,0), wy = coords(i,1), wz = coords(i,2);

    double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy + world_to_vox(0,2)*wz + world_to_vox(0,3);
    double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy + world_to_vox(1,2)*wz + world_to_vox(1,3);
    double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy + world_to_vox(2,2)*wz + world_to_vox(2,3);

    double J[9] = {1,0,0, 0,1,0, 0,0,1};
    double disp_plus[3], disp_minus[3];

    if (sample_disp_local(fdata, nx, ny, nz, vx+h, vy, vz, disp_plus) &&
        sample_disp_local(fdata, nx, ny, nz, vx-h, vy, vz, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[0]);
      J[0] += (disp_plus[0] - disp_minus[0]) * scale;
      J[3] += (disp_plus[1] - disp_minus[1]) * scale;
      J[6] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    if (sample_disp_local(fdata, nx, ny, nz, vx, vy+h, vz, disp_plus) &&
        sample_disp_local(fdata, nx, ny, nz, vx, vy-h, vz, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[1]);
      J[1] += (disp_plus[0] - disp_minus[0]) * scale;
      J[4] += (disp_plus[1] - disp_minus[1]) * scale;
      J[7] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    if (sample_disp_local(fdata, nx, ny, nz, vx, vy, vz+h, disp_plus) &&
        sample_disp_local(fdata, nx, ny, nz, vx, vy, vz-h, disp_minus)) {
      double scale = 1.0 / (2.0 * h * vox_spacing[2]);
      J[2] += (disp_plus[0] - disp_minus[0]) * scale;
      J[5] += (disp_plus[1] - disp_minus[1]) * scale;
      J[8] += (disp_plus[2] - disp_minus[2]) * scale;
    }

    // det(J) = J[0]*(J[4]*J[8]-J[5]*J[7]) - J[1]*(J[3]*J[8]-J[5]*J[6]) + J[2]*(J[3]*J[7]-J[4]*J[6])
    out[i] = J[0]*(J[4]*J[8] - J[7]*J[5])
           - J[1]*(J[3]*J[8] - J[6]*J[5])
           + J[2]*(J[3]*J[7] - J[6]*J[4]);
  }

  return out;
}
