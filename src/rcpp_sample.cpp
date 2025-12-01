#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// Trilinear interpolation helper for 3D volume
inline double trilinear_interp(const double* data, int nx, int ny, int nz,
                                double vx, double vy, double vz, double outside) {
  if (vx < -0.5 || vy < -0.5 || vz < -0.5 ||
      vx > nx - 0.5 || vy > ny - 0.5 || vz > nz - 0.5) {
    return outside;
  }

  int x0 = (int)floor(vx), y0 = (int)floor(vy), z0 = (int)floor(vz);
  double fx = vx - x0, fy = vy - y0, fz = vz - z0;
  double wx[2] = {1.0 - fx, fx};
  double wy[2] = {1.0 - fy, fy};
  double wz[2] = {1.0 - fz, fz};

  double acc = 0.0;
  for (int iz = 0; iz < 2; ++iz) {
    for (int iy = 0; iy < 2; ++iy) {
      for (int ix = 0; ix < 2; ++ix) {
        int xi = std::min(std::max(x0 + ix, 0), nx - 1);
        int yi = std::min(std::max(y0 + iy, 0), ny - 1);
        int zi = std::min(std::max(z0 + iz, 0), nz - 1);
        double w = wx[ix] * wy[iy] * wz[iz];
        acc += w * data[(size_t(zi) * ny + yi) * nx + xi];
      }
    }
  }
  return acc;
}

// Nearest neighbor interpolation
inline double nearest_interp(const double* data, int nx, int ny, int nz,
                              double vx, double vy, double vz, double outside) {
  int xi = (int)round(vx);
  int yi = (int)round(vy);
  int zi = (int)round(vz);

  if (xi < 0 || yi < 0 || zi < 0 || xi >= nx || yi >= ny || zi >= nz) {
    return outside;
  }
  return data[(size_t(zi) * ny + yi) * nx + xi];
}

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

// Tricubic interpolation
inline double tricubic_interp(const double* data, int nx, int ny, int nz,
                               double vx, double vy, double vz, double outside) {
  if (vx < 0.5 || vy < 0.5 || vz < 0.5 ||
      vx > nx - 1.5 || vy > ny - 1.5 || vz > nz - 1.5) {
    // Fall back to trilinear at boundaries
    return trilinear_interp(data, nx, ny, nz, vx, vy, vz, outside);
  }

  int x0 = (int)floor(vx) - 1;
  int y0 = (int)floor(vy) - 1;
  int z0 = (int)floor(vz) - 1;
  double fx = vx - (x0 + 1);
  double fy = vy - (y0 + 1);
  double fz = vz - (z0 + 1);

  double z_slices[4];
  for (int iz = 0; iz < 4; ++iz) {
    double y_rows[4];
    for (int iy = 0; iy < 4; ++iy) {
      double x_vals[4];
      for (int ix = 0; ix < 4; ++ix) {
        int xi = std::min(std::max(x0 + ix, 0), nx - 1);
        int yi = std::min(std::max(y0 + iy, 0), ny - 1);
        int zi = std::min(std::max(z0 + iz, 0), nz - 1);
        x_vals[ix] = data[(size_t(zi) * ny + yi) * nx + xi];
      }
      y_rows[iy] = cubic_interp_1d(x_vals, fx);
    }
    z_slices[iz] = cubic_interp_1d(y_rows, fy);
  }
  return cubic_interp_1d(z_slices, fz);
}

//' Sample a 3D volume at world coordinates
//'
//' Interpolates volume data at arbitrary world coordinates using the specified
//' interpolation method.
//'
//' @param data Numeric array (3D volume)
//' @param coords Numeric matrix (N x 3) of world coordinates
//' @param world_to_vox 4x4 world-to-voxel transform
//' @param method Interpolation method: "linear", "nearest", or "cubic"
//' @param outside Value for out-of-bounds locations
//' @return Numeric vector of N interpolated values
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_sample_volume(const Rcpp::NumericVector& data,
                                       const Rcpp::NumericMatrix& coords,
                                       const Rcpp::NumericMatrix& world_to_vox,
                                       const std::string& method,
                                       double outside) {
  Rcpp::IntegerVector dims = data.attr("dim");
  if (dims.size() < 3) Rcpp::stop("data must be at least 3D");

  int nx = dims[0], ny = dims[1], nz = dims[2];
  int n = coords.nrow();
  const double* dptr = data.begin();

  // Check for 4D data
  int nvol = 1;
  if (dims.size() > 3) nvol = dims[3];

  Rcpp::NumericVector out(n * nvol);
  if (nvol > 1) {
    out.attr("dim") = Rcpp::IntegerVector::create(n, nvol);
  }

  bool use_linear = (method == "linear");
  bool use_cubic = (method == "cubic");

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double wx = coords(i,0), wy = coords(i,1), wz = coords(i,2);

    // Convert world to voxel
    double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy + world_to_vox(0,2)*wz + world_to_vox(0,3);
    double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy + world_to_vox(1,2)*wz + world_to_vox(1,3);
    double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy + world_to_vox(2,2)*wz + world_to_vox(2,3);

    for (int v = 0; v < nvol; ++v) {
      const double* vol_ptr = dptr + size_t(v) * nx * ny * nz;
      double val;

      if (use_linear) {
        val = trilinear_interp(vol_ptr, nx, ny, nz, vx, vy, vz, outside);
      } else if (use_cubic) {
        val = tricubic_interp(vol_ptr, nx, ny, nz, vx, vy, vz, outside);
      } else {
        val = nearest_interp(vol_ptr, nx, ny, nz, vx, vy, vz, outside);
      }

      if (nvol == 1) {
        out[i] = val;
      } else {
        out[i + n*v] = val;
      }
    }
  }

  return out;
}

//' Find nearest vertex for each query point
//'
//' @param coords Query coordinates (N x 3)
//' @param vertices Vertex coordinates (V x 3)
//' @return Integer vector of N vertex indices (1-based)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_nearest_vertex(const Rcpp::NumericMatrix& coords,
                                        const Rcpp::NumericMatrix& vertices) {
  int n = coords.nrow();
  int nv = vertices.nrow();
  Rcpp::IntegerVector out(n);

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double qx = coords(i,0), qy = coords(i,1), qz = coords(i,2);
    double best_dist = std::numeric_limits<double>::infinity();
    int best_idx = 0;

    for (int v = 0; v < nv; ++v) {
      double dx = vertices(v,0) - qx;
      double dy = vertices(v,1) - qy;
      double dz = vertices(v,2) - qz;
      double d2 = dx*dx + dy*dy + dz*dz;
      if (d2 < best_dist) {
        best_dist = d2;
        best_idx = v;
      }
    }
    out[i] = best_idx + 1;  // 1-based for R
  }

  return out;
}

//' Barycentric interpolation on surface mesh
//'
//' @param coords Query coordinates (N x 3)
//' @param vertices Mesh vertex coordinates (V x 3)
//' @param faces Face indices (F x 3, 1-based)
//' @param data Vertex data (V) or (V x K)
//' @return Interpolated values
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_barycentric_sample(const Rcpp::NumericMatrix& coords,
                                            const Rcpp::NumericMatrix& vertices,
                                            const Rcpp::IntegerMatrix& faces,
                                            const Rcpp::NumericVector& data) {
  int n = coords.nrow();
  int nf = faces.nrow();
  int nv = vertices.nrow();

  // Check if data is matrix (V x K) or vector (V)
  int k = 1;
  SEXP dims_sexp = data.attr("dim");
  if (dims_sexp != R_NilValue) {
    Rcpp::IntegerVector data_dims(dims_sexp);
    if (data_dims.size() == 2) {
      k = data_dims[1];
    }
  }

  Rcpp::NumericVector out(n * k);
  if (k > 1) {
    out.attr("dim") = Rcpp::IntegerVector::create(n, k);
  }

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double qx = coords(i,0), qy = coords(i,1), qz = coords(i,2);

    // Find best face (brute force - could use spatial index)
    double best_dist = std::numeric_limits<double>::infinity();
    int best_face = -1;
    double best_u = 0, best_v = 0, best_w = 0;

    for (int f = 0; f < nf; ++f) {
      // 1-based to 0-based
      int v0 = faces(f,0) - 1;
      int v1 = faces(f,1) - 1;
      int v2 = faces(f,2) - 1;

      // Triangle vertices
      double ax = vertices(v0,0), ay = vertices(v0,1), az = vertices(v0,2);
      double bx = vertices(v1,0), by = vertices(v1,1), bz = vertices(v1,2);
      double cx = vertices(v2,0), cy = vertices(v2,1), cz = vertices(v2,2);

      // Compute barycentric coordinates
      double e1x = bx-ax, e1y = by-ay, e1z = bz-az;
      double e2x = cx-ax, e2y = cy-ay, e2z = cz-az;
      double px = qx-ax, py = qy-ay, pz = qz-az;

      double d11 = e1x*e1x + e1y*e1y + e1z*e1z;
      double d12 = e1x*e2x + e1y*e2y + e1z*e2z;
      double d22 = e2x*e2x + e2y*e2y + e2z*e2z;
      double d1p = e1x*px + e1y*py + e1z*pz;
      double d2p = e2x*px + e2y*py + e2z*pz;

      double denom = d11*d22 - d12*d12 + 1e-15;
      double bv = (d22*d1p - d12*d2p) / denom;
      double bw = (d11*d2p - d12*d1p) / denom;
      double bu = 1.0 - bv - bw;

      // Check if inside or near triangle
      if (bu >= -0.01 && bv >= -0.01 && bw >= -0.01) {
        // Compute distance to centroid as tie-breaker
        double centx = (ax+bx+cx)/3, centy = (ay+by+cy)/3, centz = (az+bz+cz)/3;
        double dx = qx-centx, dy = qy-centy, dz = qz-centz;
        double dist = dx*dx + dy*dy + dz*dz;

        if (dist < best_dist) {
          best_dist = dist;
          best_face = f;
          // Clamp and renormalize
          bu = std::max(0.0, bu);
          bv = std::max(0.0, bv);
          bw = std::max(0.0, bw);
          double sum = bu + bv + bw;
          best_u = bu/sum;
          best_v = bv/sum;
          best_w = bw/sum;
        }
      }
    }

    if (best_face >= 0) {
      int v0 = faces(best_face,0) - 1;
      int v1 = faces(best_face,1) - 1;
      int v2 = faces(best_face,2) - 1;

      for (int kk = 0; kk < k; ++kk) {
        double val0 = (k == 1) ? data[v0] : data[v0 + nv*kk];
        double val1 = (k == 1) ? data[v1] : data[v1 + nv*kk];
        double val2 = (k == 1) ? data[v2] : data[v2 + nv*kk];
        double interp = best_u*val0 + best_v*val1 + best_w*val2;

        if (k == 1) {
          out[i] = interp;
        } else {
          out[i + n*kk] = interp;
        }
      }
    } else {
      // No face found - return NA
      for (int kk = 0; kk < k; ++kk) {
        if (k == 1) {
          out[i] = NA_REAL;
        } else {
          out[i + n*kk] = NA_REAL;
        }
      }
    }
  }

  return out;
}

//' Ribbon sampling from volume onto surface
//'
//' Samples along cortical ribbon (inner to outer surface) and averages.
//'
//' @param data Volume data (3D array)
//' @param inner Inner surface coordinates (N x 3)
//' @param outer Outer surface coordinates (N x 3)
//' @param world_to_vox 4x4 world-to-voxel transform
//' @param n_steps Number of samples along ribbon
//' @param method Interpolation method
//' @return Numeric vector of N averaged values
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_ribbon_sample_volume(const Rcpp::NumericVector& data,
                                              const Rcpp::NumericMatrix& inner,
                                              const Rcpp::NumericMatrix& outer,
                                              const Rcpp::NumericMatrix& world_to_vox,
                                              int n_steps,
                                              const std::string& method) {
  Rcpp::IntegerVector dims = data.attr("dim");
  int nx = dims[0], ny = dims[1], nz = dims[2];
  int n = inner.nrow();
  const double* dptr = data.begin();

  bool use_linear = (method == "linear");

  Rcpp::NumericVector out(n);

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    double x0 = inner(i,0), y0 = inner(i,1), z0 = inner(i,2);
    double x1 = outer(i,0), y1 = outer(i,1), z1 = outer(i,2);

    double sum = 0.0;
    int count = 0;

    for (int s = 0; s <= n_steps; ++s) {
      double alpha = static_cast<double>(s) / static_cast<double>(n_steps);
      double wx = x0 + (x1 - x0) * alpha;
      double wy = y0 + (y1 - y0) * alpha;
      double wz = z0 + (z1 - z0) * alpha;

      // Convert to voxel
      double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy + world_to_vox(0,2)*wz + world_to_vox(0,3);
      double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy + world_to_vox(1,2)*wz + world_to_vox(1,3);
      double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy + world_to_vox(2,2)*wz + world_to_vox(2,3);

      double val;
      if (use_linear) {
        val = trilinear_interp(dptr, nx, ny, nz, vx, vy, vz, NA_REAL);
      } else {
        val = nearest_interp(dptr, nx, ny, nz, vx, vy, vz, NA_REAL);
      }

      if (!ISNA(val)) {
        sum += val;
        count++;
      }
    }

    out[i] = (count > 0) ? sum / count : NA_REAL;
  }

  return out;
}
