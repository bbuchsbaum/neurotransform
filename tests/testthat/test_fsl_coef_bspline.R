flatten_components <- function(arr4d) {
  nvox <- prod(dim(arr4d)[1:3])
  idx <- seq_len(nvox)
  out <- numeric(3 * nvox)
  out[3L * (idx - 1L) + 1L] <- as.numeric(arr4d[, , , 1])
  out[3L * (idx - 1L) + 2L] <- as.numeric(arr4d[, , , 2])
  out[3L * (idx - 1L) + 3L] <- as.numeric(arr4d[, , , 3])
  out
}

bs3 <- function(d) {
  d <- abs(d)
  ifelse(
    d < 1,
    (4 - 6 * d^2 + 3 * d^3) / 6,
    ifelse(d < 2, (2 - d)^3 / 6, 0)
  )
}

eval_bspline_coeff_ref <- function(coords, coeff, affine) {
  dims <- dim(coeff)[1:3]
  world_to_ctrl <- solve(affine)
  out <- matrix(NA_real_, nrow(coords), 3)

  for (i in seq_len(nrow(coords))) {
    x <- coords[i, ]
    u <- as.numeric(world_to_ctrl %*% c(x, 1))[1:3]
    start <- floor(u) - 1
    disp <- c(0, 0, 0)
    wsum <- 0
    for (kz in start[3]:(start[3] + 3)) {
      if (kz < 0 || kz >= dims[3]) next
      wz <- bs3(u[3] - kz)
      if (wz == 0) next
      for (ky in start[2]:(start[2] + 3)) {
        if (ky < 0 || ky >= dims[2]) next
        wy <- bs3(u[2] - ky)
        if (wy == 0) next
        for (kx in start[1]:(start[1] + 3)) {
          if (kx < 0 || kx >= dims[1]) next
          wx <- bs3(u[1] - kx)
          if (wx == 0) next
          w <- wx * wy * wz
          disp <- disp + w * coeff[kx + 1, ky + 1, kz + 1, ]
          wsum <- wsum + w
        }
      }
    }
    if (wsum > 0) {
      disp <- disp / wsum
    }
    out[i, ] <- x + disp
  }
  out
}

make_coef_warp_file <- function() {
  coeff <- array(0, dim = c(7, 7, 7, 3))
  # One non-zero control coefficient in X to make the expected field analytic.
  coeff[4, 4, 4, 1] <- 2
  aff <- diag(c(2, 2, 2, 1))

  space <- neuroim2::NeuroSpace(c(7, 7, 7, 3), trans = aff)
  vec <- neuroim2::DenseNeuroVec(coeff, space)
  path <- tempfile(pattern = "fnirt_coef_", fileext = ".nii.gz")
  neuroim2::write_vec(vec, path, format = "nifti")
  list(path = path, coeff = coeff, affine = aff)
}

test_that("detect_transform_type classifies coefficient warps as fsl_coef", {
  skip_if_not_installed("neuroim2")
  d <- make_coef_warp_file()
  on.exit(unlink(d$path))

  expect_equal(detect_transform_type(d$path), "fsl_coef")

  m <- read_transform(d$path, source = "src", target = "tgt")
  expect_s4_class(m, "Warp3DMorphism")
  expect_equal(m@warp_type, "fsl_coef")
})

test_that("C++ B-spline coefficient evaluator matches reference implementation", {
  skip_if_not_installed("neuroim2")
  d <- make_coef_warp_file()
  on.exit(unlink(d$path))

  coords <- rbind(
    c(6, 6, 6),
    c(6.5, 6, 6),
    c(7, 6, 6),
    c(0, 0, 0),
    c(20, 20, 20)
  )

  expected <- eval_bspline_coeff_ref(coords, d$coeff, d$affine)
  got <- cpp_apply_bspline_coeff_field(
    coords,
    flatten_components(d$coeff),
    as.integer(dim(d$coeff)[1:3]),
    solve(d$affine)
  )
  expect_equal(got, expected, tolerance = 1e-6)
})

test_that("fsl_coef morphism transform uses coefficient semantics", {
  skip_if_not_installed("neuroim2")
  d <- make_coef_warp_file()
  on.exit(unlink(d$path))

  coords <- rbind(
    c(6, 6, 6),
    c(6.5, 6, 6),
    c(7, 6, 6),
    c(0, 0, 0)
  )
  expected <- eval_bspline_coeff_ref(coords, d$coeff, d$affine)

  m <- read_transform(d$path, type = "fsl_coef", source = "src", target = "tgt")
  got <- transform(m, coords)
  expect_equal(got, expected, tolerance = 1e-6)
})

test_that("dense fsl warp path rejects likely coefficient files (no silent misread)", {
  skip_if_not_installed("neuroim2")
  d <- make_coef_warp_file()
  on.exit(unlink(d$path))

  m_bad <- Warp3DMorphism("src", "tgt", d$path, warp_type = "fsl")
  expect_error(
    transform(m_bad, matrix(c(6, 6, 6), ncol = 3)),
    "coefficient field"
  )
})

test_that("resampling plans for fsl_coef use correct pulled-back coordinates", {
  skip_if_not_installed("neuroim2")
  d <- make_coef_warp_file()
  on.exit(unlink(d$path))

  m <- read_transform(d$path, type = "fsl_coef", source = "src", target = "tgt")
  g <- grid_spec(c(7, 7, 7), d$affine)
  plan <- make_resampling_plan(m, source_grid = g, target_grid = g, reuse_count = 2L)
  expected_world <- transform(m, grid_coords(g))
  expected_vox <- (cbind(expected_world, 1) %*% t(solve(g@affine)))[, 1:3, drop = FALSE]
  expect_equal(plan$coords_src_vox, expected_vox, tolerance = 1e-6)
})
