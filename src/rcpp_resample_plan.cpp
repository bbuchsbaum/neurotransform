#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

struct ResamplePlan {
  int nx, ny, nz;
  int nvox;       // number of target voxels
  int nweights;   // 1 (nearest), 8 (linear), 64 (cubic)
  std::vector<int> indices;    // size nvox * nweights
  std::vector<double> weights; // size nvox * nweights
};

inline double cubic_weight(double t) {
  double at = std::abs(t);
  if (at <= 1.0) {
    return (1.5*at - 2.5)*at*at + 1.0;
  } else if (at < 2.0) {
    return ((-0.5*at + 2.5)*at - 4.0)*at + 2.0;
  }
  return 0.0;
}

static void fill_nearest(const Rcpp::NumericMatrix& coords, ResamplePlan& plan) {
  int n = coords.nrow();
  for (int i = 0; i < n; ++i) {
    double x = coords(i,0), y = coords(i,1), z = coords(i,2);
    int xi = (int)std::llround(x);
    int yi = (int)std::llround(y);
    int zi = (int)std::llround(z);
    int idx_base = i * plan.nweights;
    if (xi < 0 || yi < 0 || zi < 0 || xi >= plan.nx || yi >= plan.ny || zi >= plan.nz) {
      plan.indices[idx_base] = -1;
      plan.weights[idx_base] = 0.0;
      continue;
    }
    int linear_idx = (zi * plan.ny + yi) * plan.nx + xi;
    plan.indices[idx_base] = linear_idx;
    plan.weights[idx_base] = 1.0;
  }
}

static void fill_linear(const Rcpp::NumericMatrix& coords, ResamplePlan& plan) {
  int n = coords.nrow();
  for (int i = 0; i < n; ++i) {
    double x = coords(i,0), y = coords(i,1), z = coords(i,2);
    int xi = (int)std::floor(x);
    int yi = (int)std::floor(y);
    int zi = (int)std::floor(z);
    double fx = x - xi, fy = y - yi, fz = z - zi;

    int idx_base = i * plan.nweights;
    if (xi < 0 || yi < 0 || zi < 0 || xi >= plan.nx - 1 || yi >= plan.ny - 1 || zi >= plan.nz - 1) {
      for (int k = 0; k < plan.nweights; ++k) {
        plan.indices[idx_base + k] = -1;
        plan.weights[idx_base + k] = 0.0;
      }
      continue;
    }

    double wx[2] = {1.0 - fx, fx};
    double wy[2] = {1.0 - fy, fy};
    double wz[2] = {1.0 - fz, fz};

    int k = 0;
    for (int iz = 0; iz < 2; ++iz) {
      for (int iy = 0; iy < 2; ++iy) {
        for (int ix = 0; ix < 2; ++ix) {
          double w = wx[ix] * wy[iy] * wz[iz];
          int xi_n = xi + ix;
          int yi_n = yi + iy;
          int zi_n = zi + iz;
          int linear_idx = (zi_n * plan.ny + yi_n) * plan.nx + xi_n;
          plan.indices[idx_base + k] = linear_idx;
          plan.weights[idx_base + k] = w;
          ++k;
        }
      }
    }
  }
}

static void fill_cubic(const Rcpp::NumericMatrix& coords, ResamplePlan& plan) {
  int n = coords.nrow();
  for (int i = 0; i < n; ++i) {
    double x = coords(i,0), y = coords(i,1), z = coords(i,2);
    int xi = (int)std::floor(x);
    int yi = (int)std::floor(y);
    int zi = (int)std::floor(z);
    double fx = x - xi, fy = y - yi, fz = z - zi;

    int idx_base = i * plan.nweights;

    // If near boundary, fall back to trilinear weights into first 8 slots
    if (xi < 1 || yi < 1 || zi < 1 || xi >= plan.nx - 2 || yi >= plan.ny - 2 || zi >= plan.nz - 2) {
      // manual linear weights
      double fx = x - xi, fy = y - yi, fz = z - zi;
      if (xi < 0 || yi < 0 || zi < 0 || xi >= plan.nx - 1 || yi >= plan.ny - 1 || zi >= plan.nz - 1) {
        for (int k = 0; k < plan.nweights; ++k) {
          plan.indices[idx_base + k] = -1;
          plan.weights[idx_base + k] = 0.0;
        }
        continue;
      }
      double wx_lin[2] = {1.0 - fx, fx};
      double wy_lin[2] = {1.0 - fy, fy};
      double wz_lin[2] = {1.0 - fz, fz};
      int k = 0;
      for (int iz = 0; iz < 2; ++iz) for (int iy = 0; iy < 2; ++iy) for (int ix = 0; ix < 2; ++ix) {
        double w = wx_lin[ix] * wy_lin[iy] * wz_lin[iz];
        int xi_n = xi + ix;
        int yi_n = yi + iy;
        int zi_n = zi + iz;
        int linear_idx = (zi_n * plan.ny + yi_n) * plan.nx + xi_n;
        plan.indices[idx_base + k] = linear_idx;
        plan.weights[idx_base + k] = w;
        ++k;
      }
      for (; k < plan.nweights; ++k) {
        plan.indices[idx_base + k] = -1;
        plan.weights[idx_base + k] = 0.0;
      }
      continue;
    }

    double wx[4], wy[4], wz[4];
    for (int j = 0; j < 4; ++j) {
      wx[j] = cubic_weight(fx - (j - 1));
      wy[j] = cubic_weight(fy - (j - 1));
      wz[j] = cubic_weight(fz - (j - 1));
    }

    int k = 0;
    for (int iz = -1; iz <= 2; ++iz) {
      for (int iy = -1; iy <= 2; ++iy) {
        for (int ix = -1; ix <= 2; ++ix) {
          double w = wx[ix + 1] * wy[iy + 1] * wz[iz + 1];
          int xi_n = xi + ix;
          int yi_n = yi + iy;
          int zi_n = zi + iz;
          int linear_idx = (zi_n * plan.ny + yi_n) * plan.nx + xi_n;
          plan.indices[idx_base + k] = linear_idx;
          plan.weights[idx_base + k] = w;
          ++k;
        }
      }
    }
  }
}

// [[Rcpp::export]]
SEXP cpp_make_resample_plan(const Rcpp::IntegerVector& src_dim,
                            const Rcpp::NumericMatrix& src_coords_vox,
                            const std::string& method) {
  if (src_dim.size() < 3) Rcpp::stop("src_dim must have length >= 3");
  ResamplePlan* plan = new ResamplePlan();
  plan->nx = src_dim[0];
  plan->ny = src_dim[1];
  plan->nz = src_dim[2];
  plan->nvox = src_coords_vox.nrow();

  if (method == "nearest") {
    plan->nweights = 1;
  } else if (method == "linear") {
    plan->nweights = 8;
  } else if (method == "cubic") {
    plan->nweights = 64;
  } else {
    delete plan;
    Rcpp::stop("Unknown interpolation method");
  }

  plan->indices.resize(plan->nvox * plan->nweights, -1);
  plan->weights.resize(plan->nvox * plan->nweights, 0.0);

  if (method == "nearest") {
    fill_nearest(src_coords_vox, *plan);
  } else if (method == "linear") {
    fill_linear(src_coords_vox, *plan);
  } else {
    fill_cubic(src_coords_vox, *plan);
  }

  Rcpp::XPtr<ResamplePlan> xptr(plan, true);
  return xptr;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_apply_resample_plan(SEXP plan_xptr,
                                            const Rcpp::NumericVector& data,
                                            double outside) {
  Rcpp::XPtr<ResamplePlan> plan(plan_xptr);
  Rcpp::IntegerVector dims = data.attr("dim");
  if (dims.size() < 3) Rcpp::stop("data must be at least 3D");
  if (dims[0] != plan->nx || dims[1] != plan->ny || dims[2] != plan->nz) {
    Rcpp::stop("data dims do not match plan source dims");
  }

  const double* src = REAL(data);
  Rcpp::NumericVector out(plan->nvox);

  int nweights = plan->nweights;

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < plan->nvox; ++i) {
    int base = i * nweights;
    int idx0 = plan->indices[base];
    if (idx0 < 0) {
      out[i] = outside;
      continue;
    }
    if (nweights == 1) {
      out[i] = src[idx0];
      continue;
    }
    double acc = 0.0;
    for (int k = 0; k < nweights; ++k) {
      int idx = plan->indices[base + k];
      if (idx < 0) continue;
      acc += plan->weights[base + k] * src[idx];
    }
    out[i] = acc;
  }
  return out;
}
