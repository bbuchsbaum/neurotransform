#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(neurotransform))

# Simple micro-benchmark comparing plan vs direct for a small warp

warp_path <- system.file("extdata/ants/sample_ANTs_1Warp.nii.gz", package = "neurotransform")
morph <- Warp3DMorphism("src", "tgt", warp_path, warp_type = "ants")

src_grid <- grid_spec(dims = c(32L, 32L, 32L), affine = diag(4), domain = "src")
tgt_grid <- grid_spec(dims = c(32L, 32L, 32L), affine = diag(4), domain = "tgt")

vol <- array(runif(prod(src_grid@dims)), dim = src_grid@dims)

plan <- make_resampling_plan(morph, src_grid, tgt_grid, interpolation = "linear", reuse_count = 5L)

cat("Benchmarking resample_volume vs plan apply...\n")

t1 <- system.time({
  res1 <- resample_volume(vol, morphism = morph, target = tgt_grid, method = "linear", modulate = "none")
})
t2 <- system.time({
  res2 <- apply_resampling_plan(plan, vol, outside = NA_real_)
})

cat("resample_volume elapsed:", t1["elapsed"], "seconds\n")
cat("apply_resampling_plan elapsed:", t2["elapsed"], "seconds\n")

if (max(abs(res1 - res2), na.rm = TRUE) < 1e-6) {
  cat("Outputs match within tolerance.\n")
} else {
  cat("Outputs differ; investigate.\n")
}
