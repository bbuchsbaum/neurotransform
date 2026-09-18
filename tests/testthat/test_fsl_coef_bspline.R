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
