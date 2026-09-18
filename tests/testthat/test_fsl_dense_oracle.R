fsl_dense_case <- function(id) {
  folder <- system.file("extdata", "fsl_dense_oracle", id,
                        package = "neurotransform", mustWork = TRUE)
  read <- function(name) neuroim2::read_vol(file.path(folder, paste0(name, ".nii.gz")))
  source <- read("source")
  target <- read("target")
  ta <- neuroim2::trans(target)
  dims <- dim(target)[1:3]
  ijk <- as.matrix(expand.grid(lapply(dims, function(n) 0:(n-1))))
  world <- (cbind(ijk, 1) %*% t(ta))[, 1:3]
  expected <- sapply(0:2, function(k) as.numeric(as.array(read(paste0("native_coord", k)))) )
  mask <- as.numeric(as.array(read("native_support"))) > .999 &
    rowSums(ijk >= 2 & sweep(ijk, 2, dims-3, "<=")) == 3
  list(folder = folder, read = read, source = source, target = target,
       world = world, expected = expected, mask = mask,
       args = list(path = file.path(folder, "warp.nii.gz"), type = "fsl",
                   def_type = if (grepl("absolute$", id)) "absolute" else "relative",
                   source_affine = neuroim2::trans(source), source_dim = dim(source)[1:3],
                   target_affine = ta, target_dim = dims))
}

test_that("dense FSL fields match native applywarp for both grid handednesses", {
  for (src in c("left", "right")) for (ref in c("left", "right")) {
    for (representation in c("relative", "absolute")) {
      d <- fsl_dense_case(paste(src, ref, representation, sep = "_"))
      m <- do.call(read_transform, d$args)
      expect_gt(sum(d$mask), 100)
      actual <- transform(m, d$world[d$mask, ])
      expect_true(all(is.finite(actual)))
      expect_lt(max(abs(actual-d$expected[d$mask, ])), 3e-5)
      img <- as.numeric(as.array(resample_to(d$source, d$target, m, method = "linear")))
      expected <- as.numeric(as.array(d$read("native_source")))
      expect_true(all(is.finite(img[d$mask])))
      expect_lt(max(abs(img[d$mask]-expected[d$mask])), 3e-6)
      # Exercise the flattened path separately; it must not interpret absolute
      # FSL coordinates as displacements, or add the target position twice.
      sg <- grid_spec(dim(d$source)[1:3], neuroim2::trans(d$source))
      tg <- grid_spec(dim(d$target)[1:3], neuroim2::trans(d$target))
      plan <- make_resampling_plan(m, sg, tg, interpolation = "linear", reuse_count = 2L)
      planned <- as.numeric(apply_resampling_plan(plan, as.array(d$source)))
      expect_lt(max(abs(planned[d$mask]-expected[d$mask])), 3e-6)
    }
  }
})

test_that("FSL geometry is required and participates in decoding and caching", {
  d <- fsl_dense_case("right_left_relative")
  # Geometry problems surface when the transform is read, not at first use.
  expect_error(read_transform(d$args$path, type = "fsl"),
               "require source_affine and source_dim")
  args <- d$args
  args$source_dim <- c(3, -1, 4)
  expect_error(do.call(read_transform, args), "positive integer")
  first <- do.call(read_transform, d$args)
  args <- d$args
  args$source_affine[1:3, 4] <- args$source_affine[1:3, 4] + c(3, -2, 1)
  second <- do.call(read_transform, args)
  second@cache <- first@cache
  points <- d$world[d$mask, ][1:5, ]
  a <- transform(first, points)
  b <- transform(second, points)
  expect_equal(b-a, matrix(rep(c(3, -2, 1), each = 5), ncol = 3), tolerance = 1e-10)
  expect_false(identical(first@hash, second@hash))
})

test_that("FSL absolute exports preserve native sampling coordinates", {
  d <- fsl_dense_case("left_right_absolute")
  m <- do.call(read_transform, d$args)
  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path))
  for (representation in c("auto", "displacements", "deformations")) {
    write_warp_field(m, path, representation)
    exported <- read_transform(path, type = "ants",
                               def_type = if (representation == "displacements") "relative" else "absolute")
    expect_lt(max(abs(transform(exported, d$world[d$mask, ])-d$expected[d$mask, ])), 3e-5)
  }
})

test_that("FSL Jacobians match finite differences of native coordinate ramps", {
  for (representation in c("relative", "absolute")) {
    d <- fsl_dense_case(paste0("right_left_", representation))
    m <- do.call(read_transform, d$args)
    dims <- dim(d$target)[1:3]
    ijk <- c(5L, 6L, 7L)
    offset <- function(v) 1L + v[1] + dims[1]*(v[2] + dims[2]*v[3])
    derivative <- sapply(1:3, function(axis) {
      step <- diag(3)[, axis]
      (d$expected[offset(ijk+step), ]-d$expected[offset(ijk-step), ])/2
    })
    expected <- derivative %*% solve(d$args$target_affine[1:3, 1:3])
    point <- d$world[offset(ijk), , drop = FALSE]
    expect_equal(unname(jacobian(m, point)@values[1, , ]), unname(expected), tolerance = 2e-5)
    expect_equal(as.numeric(jacobian_det(m, point)), det(expected), tolerance = 2e-5)
  }
})

test_that("inverting FSL fields swaps source and reference geometry", {
  d <- fsl_dense_case("right_left_relative")
  m <- Warp3DMorphism("source", "target", d$args$path, warp_type = "fsl",
                     inverse_path = d$args$path,
                     source_affine = d$args$source_affine, source_dim = d$args$source_dim,
                     target_affine = d$args$target_affine, target_dim = d$args$target_dim)
  inv <- invert(m)
  expect_equal(inv@params$source_affine, m@params$target_affine)
  expect_equal(inv@params$source_dim, m@params$target_dim)
  expect_equal(inv@params$target_affine, m@params$source_affine)
  expect_equal(inv@params$target_dim, m@params$source_dim)
  expect_equal(invert(inv)@params, m@params)
})

test_that("dense FSL fields default the reference geometry to the warp lattice", {
  for (id in c("left_left", "left_right", "right_left", "right_right")) {
    for (representation in c("relative", "absolute")) {
      d <- fsl_dense_case(paste(id, representation, sep = "_"))
      args <- d$args
      args$target_affine <- NULL
      args$target_dim <- NULL
      m <- do.call(read_transform, args)
      expect_null(m@params$target_affine)
      actual <- transform(m, d$world[d$mask, ])
      expect_lt(max(abs(actual - d$expected[d$mask, ])), 3e-5)
    }
  }
})

test_that("inverting without reference geometry reads it from the forward header", {
  d <- fsl_dense_case("left_right_relative")
  m <- Warp3DMorphism("source", "target", d$args$path, warp_type = "fsl",
                      inverse_path = d$args$path,
                      source_affine = d$args$source_affine, source_dim = d$args$source_dim)
  inv <- invert(m)
  lattice <- load_warp_neuroim2(d$args$path)
  expect_equal(inv@params$source_affine, lattice$vox_to_world, tolerance = 1e-6)
  expect_equal(inv@params$source_dim, lattice$dim)
  expect_equal(inv@params$target_affine, d$args$source_affine)
  expect_error(invert(Warp3DMorphism("s", "t", tempfile(fileext = ".nii.gz"), warp_type = "fsl",
                                     inverse_path = d$args$path,
                                     source_affine = d$args$source_affine,
                                     source_dim = d$args$source_dim)),
               "Warp file not found")
})

test_that("native FSL fields are detected as FSL and need source geometry", {
  for (id in c("left_left", "left_right", "right_left", "right_right")) {
    for (representation in c("relative", "absolute")) {
      d <- fsl_dense_case(paste(id, representation, sep = "_"))
      expect_equal(detect_transform_type(d$args$path), "fsl")
      expect_error(read_transform(d$args$path), "Detected a dense FSL warp field")
      inferred <- read_transform(d$args$path,
                                 source_affine = d$args$source_affine,
                                 source_dim = d$args$source_dim)
      expect_equal(inferred@warp_type, "fsl")
      expect_equal(inferred@params$def_type, representation)
      expect_lt(max(abs(transform(inferred, d$world[d$mask, ]) - d$expected[d$mask, ])), 3e-5)
    }
  }
})

test_that("resample_to() lends the moving image's geometry to FSL field paths", {
  for (id in c("right_left_relative", "left_right_absolute")) {
    d <- fsl_dense_case(id)
    expected <- as.numeric(as.array(d$read("native_source")))
    by_path <- as.numeric(as.array(resample_to(d$source, d$target, d$args$path, method = "linear")))
    expect_lt(max(abs(by_path[d$mask] - expected[d$mask])), 3e-6)
  }
  # Without a geometry-bearing moving image the path still needs explicit geometry.
  plain <- array(as.numeric(as.array(d$source)), dim(d$source))
  expect_error(resample_to(plain, d$target, d$args$path),
               "Detected a dense FSL warp field")
})

test_that("resample_to() does not guess the reference grid of coefficient files", {
  skip_if_not_installed("RNifti")
  folder <- system.file("extdata", "fsl_coef_oracle", "srcleft_refright_aff",
                        package = "neurotransform", mustWork = TRUE)
  source <- neuroim2::read_vol(file.path(folder, "source.nii.gz"))
  target <- neuroim2::read_vol(file.path(folder, "target.nii.gz"))
  expect_error(resample_to(source, target, file.path(folder, "coef.nii.gz")),
               "require source_affine and source_dim")
})
