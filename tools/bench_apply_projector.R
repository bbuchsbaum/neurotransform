#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this benchmark script")
}

pkgload::load_all(".", quiet = TRUE)

make_projector <- function(n_rows, n_cols, nnz_per_row) {
  Matrix::sparseMatrix(
    i = rep(seq_len(n_rows), each = nnz_per_row),
    j = sample.int(n_cols, n_rows * nnz_per_row, replace = TRUE),
    x = runif(n_rows * nnz_per_row),
    dims = c(n_rows, n_cols),
    repr = "C"
  )
}

run_shape <- function(label, n_rows, n_cols, nnz_per_row, n_time,
                      threads = 1L, check_reference = FALSE) {
  cat(sprintf("\n[%s]\n", label))
  cat(sprintf("shape=%d x %d, nnz/row=%d, T=%d, threads=%d\n",
              n_rows, n_cols, nnz_per_row, n_time, threads))

  proj <- make_projector(n_rows, n_cols, nnz_per_row)
  data <- matrix(runif(n_cols * n_time), nrow = n_cols, ncol = n_time)

  elapsed <- system.time({
    result <- neurotransform:::cpp_apply_projector(proj, data, threads = threads)
  })["elapsed"]

  cat(sprintf("cpp_apply_projector elapsed: %.3f sec\n", as.numeric(elapsed)))
  cat(sprintf("output dim: %d x %d\n", nrow(result), ncol(result)))

  if (check_reference) {
    reference <- as.matrix(proj %*% data)
    diff <- max(abs(result - reference))
    cat(sprintf("reference max abs diff: %.3e\n", diff))
    stopifnot(isTRUE(all.equal(result, reference, tolerance = 1e-10)))
  }

  invisible(result)
}

set.seed(1)

run_shape("small-balanced", 50000L, 50000L, 4L, 50L, threads = 1L, check_reference = TRUE)
run_shape("threshold-shape", 140000L, 32000L, 2L, 8L, threads = 1L, check_reference = TRUE)
run_shape("fsLR-vol-medium", 300000L, 32000L, 4L, 20L, threads = 1L, check_reference = FALSE)
run_shape("fsLR-vol-medium-4t", 300000L, 32000L, 4L, 20L, threads = 4L, check_reference = FALSE)
