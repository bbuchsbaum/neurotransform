# Real FLIRT/FNIRT validation against native FSL 5.0.9 outputs.
#
# The fixtures are large, so they are generated locally rather than shipped:
# run inst/extdata/fsl/register_to_mni.sh (its header shows a Docker command
# using the same pinned FSL image as the other native fixtures). Without them
# these tests skip. Every expected value below comes from FSL itself:
# applywarp and flirt resampling, convertwarp's absolute field, fnirtfileutils'
# affine-free field and Jacobian, and invwarp.

fnirt_file <- function(name) {
  system.file("extdata", "fsl", name, package = "neurotransform")
}

fnirt_fixture <- function() {
  names <- c(
    mat = "highres2standard.mat",
    warp = "highres2standard_warp.nii.gz",
    abs = "highres2standard_warp_abs.nii.gz",
    noaff = "highres2standard_warp_noaff.nii.gz",
    coef = "highres2standard_warp_coef.nii.gz",
    jac = "highres2standard_jac.nii.gz",
    inv = "standard2highres_warp.nii.gz",
    applywarp = "highres_in_mni_applywarp.nii.gz",
    flirt = "highres_in_mni_flirt.nii.gz"
  )
  paths <- vapply(names, fnirt_file, character(1))
  skip_if_not(all(nzchar(paths)),
              "Real FNIRT fixtures not generated (run inst/extdata/fsl/register_to_mni.sh)")
  src_path <- system.file("extdata/afni/ss_sub-1001_T1w.nii.gz", package = "neurotransform")
  skip_if_not(nzchar(src_path), "Source image not available")
  skip_if_not_installed("neuroim2")

  src <- neuroim2::read_vol(src_path)
  ref <- neuroim2::read_vol(paths[["applywarp"]])
  grid <- grid_spec(dim(ref)[1:3], neuroim2::trans(ref))
  list(
    path = as.list(paths), src = src, ref = ref, grid = grid, points = grid_coords(grid),
    src_affine = neuroim2::trans(src), src_dim = dim(src)[1:3],
    mni_affine = neuroim2::trans(ref), mni_dim = dim(ref)[1:3]
  )
}

fnirt_forward <- function(f, path = f$path$warp, def_type = "relative") {
  read_transform(path, type = "fsl", source = "native", target = "mni",
                 def_type = def_type,
                 source_affine = f$src_affine, source_dim = f$src_dim)
}

# Target voxels whose source sample lies at least one voxel inside the source
# grid, so boundary and padding conventions do not enter the comparison.
interior_samples <- function(f, source_points) {
  vox <- (cbind(source_points, 1) %*% t(solve(f$src_affine)))[, 1:3, drop = FALSE]
  rowSums(vox >= 1 & sweep(vox, 2, f$src_dim - 2, "<=")) == 3
}

relative_image_error <- function(actual, expected, mask) {
  actual <- as.numeric(as.array(actual))
  expected <- as.numeric(as.array(expected))
  max(abs(actual - expected)[mask]) / diff(range(expected[mask]))
}

test_that("real FNIRT outputs are detected by header and representation", {
  f <- fnirt_fixture()
  expect_equal(detect_transform_type(f$path$warp), "fsl")
  expect_equal(detect_transform_type(f$path$coef), "fsl_coef")
  expect_equal(detect_fnirt_def_type(f$path$warp), "relative")
  expect_equal(detect_fnirt_def_type(f$path$noaff), "relative")
  expect_equal(detect_fnirt_def_type(f$path$inv), "relative")
  expect_equal(detect_fnirt_def_type(f$path$abs), "absolute")

  expect_error(read_transform(f$path$warp), "Detected a dense FSL warp field")
  inferred <- read_transform(f$path$warp, source_affine = f$src_affine, source_dim = f$src_dim)
  expect_equal(inferred@warp_type, "fsl")
  expect_equal(inferred@params$def_type, "relative")
})

test_that("FNIRT resampling reproduces applywarp for relative and absolute fields", {
  f <- fnirt_fixture()
  fwd <- fnirt_forward(f)
  mask <- interior_samples(f, transform(fwd, f$points))
  expect_gt(mean(mask), 0.9)
  expected <- neuroim2::read_vol(f$path$applywarp)

  # Measured: 1.6e-5 of the intensity range (float32 field and output).
  # Convention errors produce errors of order 1e-1 or more.
  out <- resample_to(f$src, f$ref, fwd, method = "linear")
  expect_lt(relative_image_error(out, expected, mask), 1e-4)

  absolute <- fnirt_forward(f, f$path$abs, def_type = "absolute")
  out_abs <- resample_to(f$src, f$ref, absolute, method = "linear")
  expect_lt(relative_image_error(out_abs, expected, mask), 1e-4)

  # convertwarp's absolute field encodes the same mapping.
  sel <- which(mask)[seq(1, sum(mask), length.out = 20000)]
  pts <- f$points[sel, , drop = FALSE]
  expect_lt(max(abs(transform(absolute, pts) - transform(fwd, pts))), 1e-4)
})

test_that("FLIRT matrices reproduce flirt -applyxfm resampling", {
  f <- fnirt_fixture()
  aff <- read_linear_transform(f$path$mat, format = "fsl", source = "native", target = "mni",
                               source_affine = f$src_affine, source_dim = f$src_dim,
                               target_affine = f$mni_affine, target_dim = f$mni_dim)
  mask <- interior_samples(f, transform(aff, f$points))
  out <- resample_to(f$src, f$ref, aff, method = "linear")
  # Measured: 2.4e-4 of the range against flirt -noresampblur (FLIRT resamples
  # in single precision). A handedness or origin error gives about 0.5.
  expect_lt(relative_image_error(out, neuroim2::read_vol(f$path$flirt), mask), 1e-3)
})

test_that("FNIRT fields add the nonlinear part after the FLIRT affine", {
  f <- fnirt_fixture()
  fwd <- fnirt_forward(f)
  aff <- read_linear_transform(f$path$mat, format = "fsl", source = "native", target = "mni",
                               source_affine = f$src_affine, source_dim = f$src_dim,
                               target_affine = f$mni_affine, target_dim = f$mni_dim)
  mask <- interior_samples(f, transform(fwd, f$points))
  sel <- which(mask)[seq(1, sum(mask), length.out = 20000)]
  pts <- f$points[sel, , drop = FALSE]

  # fnirt --fout = FLIRT pullback plus the affine-free displacement, both in
  # source FSL coordinates: src_fsl = inv(A) ref_fsl + d(ref).
  d_fsl <- matrix(as.array(neuroim2::read_vec(f$path$noaff)), ncol = 3)[sel, , drop = FALSE]
  to_fsl <- fsl_world_to_fsl(f$src_affine, f$src_dim)
  to_world <- fsl_fsl_to_world(f$src_affine, f$src_dim)
  src_fsl <- (cbind(transform(aff, pts), 1) %*% t(to_fsl))[, 1:3] + d_fsl
  expected <- (cbind(src_fsl, 1) %*% t(to_world))[, 1:3]
  # Measured: 4.3e-6 mm. The field already contains the affine, so it must not
  # be composed with the FLIRT matrix again.
  expect_lt(max(abs(transform(fwd, pts) - expected)), 1e-4)
})

test_that("FNIRT Jacobian determinants match fnirtfileutils --jac", {
  f <- fnirt_fixture()
  fwd <- fnirt_forward(f)
  mask <- interior_samples(f, transform(fwd, f$points)) & as.numeric(as.array(f$ref)) > 0
  sel <- which(mask)[seq(1, sum(mask), length.out = 5000)]
  pts <- f$points[sel, , drop = FALSE]
  fsl_jac <- as.numeric(as.array(neuroim2::read_vol(f$path$jac)))[sel]

  ours <- jacobian_det(fwd, pts)
  rel <- abs(ours - fsl_jac) / abs(fsl_jac)
  # FSL differentiates the spline analytically; we difference the dense field
  # (measured median 0.36%, p99 2.3%).
  expect_true(all(sign(ours) == sign(fsl_jac)))
  expect_lt(median(rel), 0.01)
  expect_lt(unname(quantile(rel, 0.99)), 0.05)

  jac <- jacobian(fwd, pts[1:200, , drop = FALSE])
  expect_equal(det(jac), jacobian_det(fwd, pts[1:200, , drop = FALSE]), tolerance = 1e-8)
})

test_that("invwarp fields invert FNIRT fields, directly and through invert()", {
  f <- fnirt_fixture()
  fwd <- fnirt_forward(f)
  mask <- interior_samples(f, transform(fwd, f$points)) & as.numeric(as.array(f$ref)) > 0
  sel <- which(mask)[seq(1, sum(mask), length.out = 5000)]
  pts <- f$points[sel, , drop = FALSE]
  native_pts <- transform(fwd, pts)

  # The inverse field lies on the source grid; its values are MNI FSL
  # coordinates, so its source geometry is the MNI grid.
  inv <- read_transform(f$path$inv, type = "fsl", source = "mni", target = "native",
                        source_affine = f$mni_affine, source_dim = f$mni_dim)
  err <- sqrt(rowSums((transform(inv, native_pts) - pts)^2))
  # invwarp is an iterative approximation (measured median 0.024 mm, p99 0.14 mm).
  expect_lt(median(err), 0.1)
  expect_lt(unname(quantile(err, 0.99)), 0.5)

  paired <- Warp3DMorphism("native", "mni", f$path$warp, warp_type = "fsl",
                           inverse_path = f$path$inv,
                           source_affine = f$src_affine, source_dim = f$src_dim)
  expect_equal(transform(invert(paired), native_pts), transform(inv, native_pts),
               tolerance = 1e-10)
})

test_that("FNIRT coefficients reproduce fnirt --fout on the reference grid", {
  f <- fnirt_fixture()
  coef <- read_transform(f$path$coef, type = "fsl_coef", source = "native", target = "mni",
                         source_affine = f$src_affine, source_dim = f$src_dim,
                         target_affine = f$mni_affine, target_dim = f$mni_dim)
  # Every MNI voxel (measured 9.6e-6 mm).
  expect_lt(max(abs(transform(coef, f$points) - transform(fnirt_forward(f), f$points))), 1e-4)

  # invwarp writes a dense field, so the inverse of a coefficient warp is dense.
  paired <- Warp3DMorphism("native", "mni", f$path$coef, warp_type = "fsl_coef",
                           inverse_path = f$path$inv,
                           source_affine = f$src_affine, source_dim = f$src_dim,
                           target_affine = f$mni_affine, target_dim = f$mni_dim)
  inv <- invert(paired)
  expect_equal(inv@warp_type, "fsl")
  pts <- f$points[seq(1, nrow(f$points), length.out = 2000), , drop = FALSE]
  native_pts <- transform(coef, pts)
  direct <- read_transform(f$path$inv, type = "fsl", source_affine = f$mni_affine,
                           source_dim = f$mni_dim)
  expect_equal(transform(inv, native_pts), transform(direct, native_pts), tolerance = 1e-10)
})
