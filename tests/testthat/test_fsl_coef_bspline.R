# FNIRT spline-coefficient files against native FSL 5.0.9 (fnirt --cout,
# fnirtfileutils, applywarp --warp=coef). See tools/generate_fsl_coef_oracle.py
# and inst/extdata/fsl_coef_oracle/README.md.

fsl_coef_case <- function(id) {
  folder <- system.file("extdata", "fsl_coef_oracle", id,
                        package = "neurotransform", mustWork = TRUE)
  path <- function(name) file.path(folder, name)
  read <- function(name) neuroim2::read_vol(path(paste0(name, ".nii.gz")))
  source <- read("source")
  target <- read("target")
  ta <- neuroim2::trans(target)
  dims <- dim(target)[1:3]
  ijk <- as.matrix(expand.grid(lapply(dims, function(n) 0:(n - 1))))
  list(
    path = path, read = read, source = source, target = target,
    world = (cbind(ijk, 1) %*% t(ta))[, 1:3],
    expected = sapply(0:2, function(k) as.numeric(as.array(read(paste0("native_coord", k))))),
    mask = as.numeric(as.array(read("native_support"))) > .999,
    args = list(path = path("coef.nii.gz"), type = "fsl_coef",
                source_affine = neuroim2::trans(source), source_dim = dim(source)[1:3],
                target_affine = ta, target_dim = dims)
  )
}

fsl_coef_ids <- function() {
  grid <- expand.grid(aff = c("noaff", "aff"), ref = c("refleft", "refright"),
                      src = c("srcleft", "srcright"), stringsAsFactors = FALSE)
  paste(grid$src, grid$ref, grid$aff, sep = "_")
}

test_that("coefficient files reproduce applywarp --warp=coef for all handedness pairs", {
  skip_if_not_installed("RNifti")
  for (id in fsl_coef_ids()) {
    d <- fsl_coef_case(id)
    m <- do.call(read_transform, d$args)
    expect_gt(sum(d$mask), 200)
    # Measured: at most 1.2e-5 mm (float32 fixtures).
    expect_lt(max(abs(transform(m, d$world[d$mask, ]) - d$expected[d$mask, ])), 1e-4)

    expected <- as.numeric(as.array(d$read("native_source")))
    scale <- diff(range(expected[d$mask]))
    img <- as.numeric(as.array(resample_to(d$source, d$target, m, method = "linear")))
    expect_lt(max(abs(img - expected)[d$mask]) / scale, 1e-5)
    sg <- grid_spec(dim(d$source)[1:3], neuroim2::trans(d$source))
    tg <- grid_spec(dim(d$target)[1:3], neuroim2::trans(d$target))
    plan <- make_resampling_plan(m, sg, tg, interpolation = "linear", reuse_count = 2L)
    planned <- as.numeric(apply_resampling_plan(plan, as.array(d$source)))
    expect_lt(max(abs(planned - expected)[d$mask]) / scale, 1e-5)
  }
})

test_that("coefficient files decode to the fnirtfileutils --withaff field", {
  skip_if_not_installed("RNifti")
  for (id in fsl_coef_ids()) {
    d <- fsl_coef_case(id)
    m <- do.call(read_transform, d$args)
    dense <- read_transform(d$path("field_aff.nii.gz"), type = "fsl",
                            source_affine = d$args$source_affine, source_dim = d$args$source_dim)
    # Every reference voxel, including outside the source support (measured 2.4e-6 mm).
    expect_lt(max(abs(transform(m, d$world) - transform(dense, d$world))), 1e-5)
    expect_equal(jacobian_det(m, d$world[d$mask, ][1:20, ]),
                 jacobian_det(dense, d$world[d$mask, ][1:20, ]), tolerance = 1e-5)
  }
})

test_that("coefficient files are detected from the header", {
  skip_if_not_installed("RNifti")
  d <- fsl_coef_case("srcright_refleft_aff")
  expect_equal(detect_transform_type(d$args$path), "fsl_coef")
  renamed <- file.path(tempfile("coef_"), "field.nii.gz")
  dir.create(dirname(renamed))
  file.copy(d$args$path, renamed)
  expect_equal(detect_transform_type(renamed), "fsl_coef")
  # The dense reader refuses coefficient files even without "coef" in the name.
  as_dense <- Warp3DMorphism("s", "t", renamed, warp_type = "fsl",
                             source_affine = d$args$source_affine, source_dim = d$args$source_dim)
  expect_error(transform(as_dense, d$world[1:2, ]), "coefficient field")
})

test_that("coefficient files require and check the reference geometry", {
  skip_if_not_installed("RNifti")
  d <- fsl_coef_case("srcleft_refright_aff")
  args <- d$args
  args$target_affine <- NULL
  args$target_dim <- NULL
  expect_error(do.call(read_transform, args), "require target_affine and target_dim")

  args <- d$args
  args$target_dim <- args$target_dim + c(1, 0, 0)
  expect_error(transform(do.call(read_transform, args), d$world[1:2, ]),
               "does not match the reference dimensions")

  args <- d$args
  args$target_affine[1:3, 1:3] <- args$target_affine[1:3, 1:3] * 1.1
  expect_error(transform(do.call(read_transform, args), d$world[1:2, ]), "voxel sizes")
})

test_that("files without FSL spline headers are rejected, not guessed", {
  skip_if_not_installed("RNifti")
  d <- fsl_coef_case("srcleft_refleft_aff")
  dir <- tempfile("coef_bad_")
  dir.create(dir)

  unlabeled <- file.path(dir, "unlabeled_coef.nii.gz")
  RNifti::writeNifti(array(0, c(7L, 7L, 7L, 3L)), unlabeled)
  m <- Warp3DMorphism("s", "t", unlabeled, warp_type = "fsl_coef",
                      source_affine = diag(4), source_dim = c(7L, 7L, 7L),
                      target_affine = diag(4), target_dim = c(7L, 7L, 7L))
  expect_error(transform(m, matrix(1, 1, 3)), "Not a supported FNIRT coefficient file")

  reflected <- file.path(dir, "reflected_coef.nii.gz")
  img <- RNifti::readNifti(d$args$path)
  xform <- RNifti::xform(img, useQuaternionFirst = FALSE)
  xform[1, 1:3] <- -xform[1, 1:3]
  RNifti::sform(img) <- structure(xform, code = 1L)
  RNifti::writeNifti(img, reflected)
  args <- d$args
  args$path <- reflected
  expect_error(transform(do.call(read_transform, args), d$world[1:2, ]), "positive determinant")
})

# A coefficient file with the header fields fnirt writes (see load_warp_fsl_coef).
write_coef_file <- function(coef, knot, ref_dim, ref_voxel, path = tempfile(fileext = ".nii.gz")) {
  img <- RNifti::asNifti(coef)
  RNifti::writeNifti(img, path, datatype = "float")
  h <- RNifti::niftiHeader(path)
  h$intent_code <- 2007L
  h$pixdim[2:4] <- knot
  h$intent_p1 <- ref_voxel[1]; h$intent_p2 <- ref_voxel[2]; h$intent_p3 <- ref_voxel[3]
  h$qform_code <- 1L
  h$qoffset_x <- ref_dim[1]; h$qoffset_y <- ref_dim[2]; h$qoffset_z <- ref_dim[3]
  h$sform_code <- 1L
  h$srow_x <- c(1, 0, 0, 0); h$srow_y <- c(0, 1, 0, 0); h$srow_z <- c(0, 0, 1, 0)
  img <- RNifti::asNifti(coef, reference = h)
  RNifti::pixdim(img) <- c(knot, 1)
  RNifti::writeNifti(img, path, datatype = "float")
  path
}

test_that("single-slice coefficient files decode without dropping an axis", {
  skip_if_not_installed("RNifti")
  ref_dim <- c(10L, 12L, 1L)
  knot <- c(2L, 2L, 1L)
  counts <- ifelse(knot > 1, ceiling((ref_dim + 1) / knot) + 2, ref_dim)
  coef <- array(0, c(counts, 3L))
  coef[4, 5, 1, 1] <- 1.5
  coef[3, 4, 1, 3] <- -2
  aff <- diag(c(-2, 2, 2, 1))
  path <- write_coef_file(coef, knot, ref_dim, c(2, 2, 2))
  m <- Warp3DMorphism("s", "t", path, warp_type = "fsl_coef",
                      source_affine = aff, source_dim = ref_dim,
                      target_affine = aff, target_dim = ref_dim)

  # Independent evaluation of the FNIRT formula (left-handed, so no x flip;
  # identity --aff; same grid for source and reference).
  bs <- function(t) ifelse(abs(t) < 1, 2 / 3 - t^2 + abs(t)^3 / 2,
                           ifelse(abs(t) < 2, (2 - abs(t))^3 / 6, 0))
  ijk <- as.matrix(expand.grid(0:9, 0:11, 0L))
  d <- t(apply(ijk, 1, function(u) {
    wx <- bs(u[1] / 2 - 0:(counts[1] - 1) + 1)
    wy <- bs(u[2] / 2 - 0:(counts[2] - 1) + 1)
    wz <- bs(u[3] - 0:(counts[3] - 1))
    sapply(1:3, function(k) sum(outer(outer(wx, wy), wz) * array(coef[, , , k], counts)))
  }))
  expected_src_fsl <- ijk * 2 + d
  expected_world <- (cbind(expected_src_fsl, 1) %*% t(aff %*% solve(fsl_vox_to_fsl(aff, ref_dim))))[, 1:3]
  world <- (cbind(ijk, 1) %*% t(aff))[, 1:3]
  field <- load_warp_array(m)
  actual <- world + matrix(field$array, ncol = 3, byrow = TRUE)
  expect_equal(actual, expected_world, tolerance = 1e-6)
})

test_that("coefficient counts inconsistent with the reference are rejected", {
  skip_if_not_installed("RNifti")
  ref_dim <- c(10L, 12L, 9L)
  counts <- ceiling((ref_dim + 1) / 3) + 2
  path <- write_coef_file(array(0, c(counts - c(1L, 0L, 0L), 3L)), c(3L, 3L, 3L), ref_dim, c(2, 2, 2))
  aff <- diag(c(-2, 2, 2, 1))
  m <- Warp3DMorphism("s", "t", path, warp_type = "fsl_coef",
                      source_affine = aff, source_dim = ref_dim,
                      target_affine = aff, target_dim = ref_dim)
  expect_error(transform(m, matrix(0, 1, 3)), "inconsistent with the reference dimensions")
})

test_that("invert() takes the inverse format from the inverse file and round-trips", {
  skip_if_not_installed("RNifti")
  d <- fsl_coef_case("srcleft_refright_aff")
  geometry <- d$args[c("source_affine", "source_dim", "target_affine", "target_dim")]
  # field_aff stands in for a dense (intent 2006) inverse; only formats matter here.
  m <- do.call(Warp3DMorphism, c(list("s", "t", d$args$path, warp_type = "fsl_coef",
                                      inverse_path = d$path("field_aff.nii.gz")), geometry))
  inv <- invert(m)
  expect_equal(inv@warp_type, "fsl")
  expect_equal(inv@params$def_type, "relative")
  back <- invert(inv)
  expect_equal(back@warp_type, "fsl_coef")
  expect_equal(back@params, m@params)
  expect_equal(transform(back, d$world[1:5, ]), transform(m, d$world[1:5, ]))

  dense <- do.call(Warp3DMorphism, c(list("s", "t", d$path("field_aff.nii.gz"), warp_type = "fsl",
                                          inverse_path = d$args$path), geometry))
  expect_equal(invert(dense)@warp_type, "fsl_coef")

  expect_error(do.call(Warp3DMorphism, c(list("s", "t", d$args$path, warp_type = "fsl_coef",
                                              def_type = "absolute"), geometry)),
               "def_type must be 'relative'")
})

test_that("dense fields are identified by header intent before filename", {
  skip_if_not_installed("RNifti")
  d <- fsl_coef_case("srcleft_refright_aff")
  dir <- tempfile("coefname_")
  dir.create(dir)
  # A dense (intent 2006) field whose name mentions "coef" still loads as dense.
  named <- file.path(dir, "sub01_warpcoef_as_field.nii.gz")
  file.copy(d$path("field_aff.nii.gz"), named)
  expect_equal(detect_transform_type(named), "fsl")
  m <- read_transform(named, type = "fsl", source_affine = d$args$source_affine,
                      source_dim = d$args$source_dim)
  reference <- read_transform(d$path("field_aff.nii.gz"), type = "fsl",
                              source_affine = d$args$source_affine, source_dim = d$args$source_dim)
  expect_equal(transform(m, d$world[1:5, ]), transform(reference, d$world[1:5, ]))

  # Without an intent, the filename still guards the dense reader.
  unlabeled <- file.path(dir, "unlabeled_coef.nii.gz")
  file.copy(system.file("extdata/fsl_dense_oracle/left_left_relative/warp.nii.gz",
                        package = "neurotransform"), unlabeled)
  guarded <- Warp3DMorphism("s", "t", unlabeled, warp_type = "fsl",
                            source_affine = d$args$source_affine, source_dim = d$args$source_dim)
  expect_error(transform(guarded, d$world[1:2, ]), "coefficient field")
})
