#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// Cubic interpolation helpers (Catmull-Rom-style)
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

// Cubic B-spline basis (compact support [-2, 2]).
inline double bspline3_weight(double d) {
  double ad = std::abs(d);
  if (ad < 1.0) {
    return (4.0 - 6.0 * ad * ad + 3.0 * ad * ad * ad) / 6.0;
  } else if (ad < 2.0) {
    double t = 2.0 - ad;
    return (t * t * t) / 6.0;
  }
  return 0.0;
}

// Trilinear sample displacement field (x,y,z,3) given voxel coords (0-based)
inline bool sample_disp(const Rcpp::NumericVector& field, const Rcpp::IntegerVector& dim,
                        double x, double y, double z, double* out) {
  int nx = dim[0], ny = dim[1], nz = dim[2];
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
    int base = ((zi*ny + yi)*nx + xi)*3;
    acc[0] += w * field[base];
    acc[1] += w * field[base + 1];
    acc[2] += w * field[base + 2];
  }
  out[0] = acc[0]; out[1] = acc[1]; out[2] = acc[2];
  return true;
}

// Tricubic sample displacement field; falls back to trilinear near boundaries
inline bool sample_disp_cubic(const Rcpp::NumericVector& field, const Rcpp::IntegerVector& dim,
                              double x, double y, double z, double* out) {
  int nx = dim[0], ny = dim[1], nz = dim[2];
  if (x < -0.5 || y < -0.5 || z < -0.5 || x > nx - 0.5 || y > ny - 0.5 || z > nz - 0.5) return false;
  if (x < 0.5 || y < 0.5 || z < 0.5 || x > nx - 1.5 || y > ny - 1.5 || z > nz - 1.5) {
    return sample_disp(field, dim, x, y, z, out);
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
        int base = ((zi*ny + yi)*nx + xi)*3;
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

// Evaluate displacement from FNIRT B-spline coefficients at control-grid IJK.
inline bool sample_bspline_coeff_disp(const Rcpp::NumericVector& coeff,
                                      const Rcpp::IntegerVector& dim,
                                      double x, double y, double z,
                                      double* out) {
  int nx = dim[0], ny = dim[1], nz = dim[2];
  // Beyond compact support of edge knots -> effectively zero displacement.
  if (x < -2.0 || y < -2.0 || z < -2.0 ||
      x > (nx - 1) + 2.0 || y > (ny - 1) + 2.0 || z > (nz - 1) + 2.0) {
    out[0] = out[1] = out[2] = 0.0;
    return false;
  }

  int sx = static_cast<int>(std::floor(x)) - 1;
  int sy = static_cast<int>(std::floor(y)) - 1;
  int sz = static_cast<int>(std::floor(z)) - 1;

  double acc0 = 0.0, acc1 = 0.0, acc2 = 0.0;
  double wsum = 0.0;
  for (int kz = sz; kz <= sz + 3; ++kz) {
    double wz = bspline3_weight(z - kz);
    if (wz == 0.0) continue;
    if (kz < 0 || kz >= nz) continue;
    for (int ky = sy; ky <= sy + 3; ++ky) {
      double wy = bspline3_weight(y - ky);
      if (wy == 0.0) continue;
      if (ky < 0 || ky >= ny) continue;
      for (int kx = sx; kx <= sx + 3; ++kx) {
        double wx = bspline3_weight(x - kx);
        if (wx == 0.0) continue;
        if (kx < 0 || kx >= nx) continue;

        double w = wx * wy * wz;
        int base = ((kz * ny + ky) * nx + kx) * 3;
        acc0 += w * coeff[base];
        acc1 += w * coeff[base + 1];
        acc2 += w * coeff[base + 2];
        wsum += w;
      }
    }
  }

  out[0] = acc0;
  out[1] = acc1;
  out[2] = acc2;
  return wsum > 0.0;
}

//' Apply warp displacement field to coordinates
//'
//' Fast C++ implementation of trilinear interpolation for displacement fields.
//'
//' @param coords Numeric matrix (N x 3) of world coordinates
//' @param field Numeric vector of displacement values (flattened 3 x X x Y x Z)
//' @param dim Integer vector of dimensions (X, Y, Z)
//' @param voxel_to_world 4x4 matrix (actually world_to_voxel inverse)
//' @return Numeric matrix (N x 3) of transformed coordinates
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_apply_warp_field(const Rcpp::NumericMatrix& coords,
                                         const Rcpp::NumericVector& field,
                                         const Rcpp::IntegerVector& dim,
                                         const Rcpp::NumericMatrix& voxel_to_world) {
  int n = coords.nrow();
  Rcpp::NumericMatrix out(n, 3);

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double x = coords(i,0), y = coords(i,1), z = coords(i,2);
    // convert world to voxel using inverse affine
    double hom[4] = {x, y, z, 1.0};
    double vinv[4] = {0,0,0,0};
    // invert 4x4 (assumes already inverted outside); here just use provided voxel_to_world inverse
    // voxel_to_world is actually inverse (world_to_voxel)
    for (int r=0;r<4;++r) {
      for (int c=0;c<4;++c) vinv[r] += voxel_to_world(r,c) * hom[c];
    }
    double disp[3];
    if (sample_disp(field, dim, vinv[0], vinv[1], vinv[2], disp)) {
      out(i,0) = x + disp[0];
      out(i,1) = y + disp[1];
      out(i,2) = z + disp[2];
    } else {
      out(i,0) = NA_REAL;
      out(i,1) = NA_REAL;
      out(i,2) = NA_REAL;
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_apply_warp_field_cubic(const Rcpp::NumericMatrix& coords,
                                               const Rcpp::NumericVector& field,
                                               const Rcpp::IntegerVector& dim,
                                               const Rcpp::NumericMatrix& voxel_to_world) {
  int n = coords.nrow();
  Rcpp::NumericMatrix out(n, 3);

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double x = coords(i,0), y = coords(i,1), z = coords(i,2);
    double hom[4] = {x, y, z, 1.0};
    double vinv[4] = {0,0,0,0};
    for (int r=0;r<4;++r) {
      for (int c=0;c<4;++c) vinv[r] += voxel_to_world(r,c) * hom[c];
    }
    double disp[3];
    if (sample_disp_cubic(field, dim, vinv[0], vinv[1], vinv[2], disp)) {
      out(i,0) = x + disp[0];
      out(i,1) = y + disp[1];
      out(i,2) = z + disp[2];
    } else {
      out(i,0) = NA_REAL;
      out(i,1) = NA_REAL;
      out(i,2) = NA_REAL;
    }
  }
  return out;
}

//' Apply FNIRT B-spline coefficient field to coordinates
//'
//' Evaluates a cubic tensor-product B-spline coefficient field at query
//' coordinates in world space, returning transformed coordinates.
//'
//' @param coords Numeric matrix (N x 3) world coordinates
//' @param coeff Numeric vector of coefficients (flattened interleaved x,y,z)
//' @param dim Integer vector (X, Y, Z) control-grid dimensions
//' @param world_to_ctrl 4x4 world-to-control-grid affine
//' @return Numeric matrix (N x 3) transformed coordinates
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_apply_bspline_coeff_field(const Rcpp::NumericMatrix& coords,
                                                  const Rcpp::NumericVector& coeff,
                                                  const Rcpp::IntegerVector& dim,
                                                  const Rcpp::NumericMatrix& world_to_ctrl) {
  int n = coords.nrow();
  Rcpp::NumericMatrix out(n, 3);

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double x = coords(i, 0), y = coords(i, 1), z = coords(i, 2);
    double hom[4] = {x, y, z, 1.0};
    double cijk[4] = {0, 0, 0, 0};
    for (int r = 0; r < 4; ++r) {
      for (int c = 0; c < 4; ++c) cijk[r] += world_to_ctrl(r, c) * hom[c];
    }

    double disp[3];
    sample_bspline_coeff_disp(coeff, dim, cijk[0], cijk[1], cijk[2], disp);
    out(i, 0) = x + disp[0];
    out(i, 1) = y + disp[1];
    out(i, 2) = z + disp[2];
  }

  return out;
}

//' Compose two warp displacement fields
//'
//' Computes B then A composition for warp fields.
//'
//' @param fieldA Displacement field A (flattened)
//' @param dimA Dimensions of field A
//' @param world_to_voxA World-to-voxel transform for A
//' @param fieldB Displacement field B (flattened)
//' @param dimB Dimensions of field B
//' @param vox_to_worldB Voxel-to-world transform for B
//' @return List with composed field and dimensions
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_compose_warp_fields(const Rcpp::NumericVector& fieldA,
                                   const Rcpp::IntegerVector& dimA,
                                   const Rcpp::NumericMatrix& world_to_voxA,
                                   const Rcpp::NumericVector& fieldB,
                                   const Rcpp::IntegerVector& dimB,
                                   const Rcpp::NumericMatrix& vox_to_worldB) {
  // Compose B then A: resulting warp maps target coords to source coords
  int nx = dimB[0], ny = dimB[1], nz = dimB[2];
  Rcpp::NumericVector out(nx * ny * nz * 3);

  // Parallel over voxels if OpenMP available
#ifdef _OPENMP
  #pragma omp parallel for collapse(3)
#endif
  for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x) {
        // voxel center in B voxel space
        double homB[4] = {double(x), double(y), double(z), 1.0};
        double worldB[4] = {0,0,0,1};
        // Convert B voxel -> world using inverse (caller passes vox_to_worldB)
        for (int r=0;r<4;++r) {
          worldB[r] = vox_to_worldB(r,0)*homB[0] +
                      vox_to_worldB(r,1)*homB[1] +
                      vox_to_worldB(r,2)*homB[2] +
                      vox_to_worldB(r,3);
        }

        // Apply warp B displacement in world
        double dispB[3];
        if (!sample_disp(fieldB, dimB, homB[0], homB[1], homB[2], dispB)) {
          // out of bounds
          continue;
        }
        worldB[0] += dispB[0];
        worldB[1] += dispB[1];
        worldB[2] += dispB[2];

        // Map to A voxel space
        double voxA[4] = {0,0,0,0};
        for (int r=0;r<4;++r) {
          voxA[r] = world_to_voxA(r,0)*worldB[0] +
                    world_to_voxA(r,1)*worldB[1] +
                    world_to_voxA(r,2)*worldB[2] +
                    world_to_voxA(r,3)*1.0;
        }

        // Apply warp A displacement (in voxel space of A)
        double dispA[3];
        if (!sample_disp(fieldA, dimA, voxA[0], voxA[1], voxA[2], dispA)) {
          continue;
        }
        // final world after A
        double worldAfter[3] = {0,0,0};
        // convert voxA + dispA to world using inverse of world_to_voxA (i.e., vox_to_worldA)
        // vox_to_worldA = solve(world_to_voxA)
        // For simplicity, caller should pass vox_to_worldA; skip here to avoid extra inverse; instead return composite displacement in world of B
        int idx = ((z * ny + y) * nx + x) * 3;
        out[idx]     = dispB[0] + dispA[0];
        out[idx + 1] = dispB[1] + dispA[1];
        out[idx + 2] = dispB[2] + dispA[2];
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("field") = out,
    Rcpp::Named("dim") = dimB
  );
}

//' Convert absolute-coordinate warp field to displacement field
//'
//' Converts a warp field storing absolute world coordinates to displacement field.
//'
//' @param field_abs Absolute coordinate field (flattened)
//' @param dim Field dimensions
//' @param vox_to_world Voxel-to-world affine
//' @return Displacement field
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_absolute_to_displacement(const Rcpp::NumericVector& field_abs,
                                                 const Rcpp::IntegerVector& dim,
                                                 const Rcpp::NumericMatrix& vox_to_world) {
  int nx = dim[0], ny = dim[1], nz = dim[2];
  if (field_abs.size() != nx * ny * nz * 3) Rcpp::stop("absolute field has wrong length");
  Rcpp::NumericVector out(field_abs.size());

#ifdef _OPENMP
  #pragma omp parallel for collapse(3)
#endif
  for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x) {
        int base = ((z * ny + y) * nx + x) * 3;
        // voxel center in world
        double hom[4] = {double(x), double(y), double(z), 1.0};
        double world[3] = {0,0,0};
        for (int r=0;r<3;++r) {
          world[r] = vox_to_world(r,0)*hom[0] + vox_to_world(r,1)*hom[1] + vox_to_world(r,2)*hom[2] + vox_to_world(r,3);
        }
        out[base]     = field_abs[base]     - world[0];
        out[base + 1] = field_abs[base + 1] - world[1];
        out[base + 2] = field_abs[base + 2] - world[2];
      }
    }
  }
  return out;
}
