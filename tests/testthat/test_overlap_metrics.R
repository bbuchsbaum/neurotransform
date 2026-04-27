test_that("dice and jaccard match analytic mask cases", {
  x <- array(c(TRUE, TRUE, FALSE, FALSE), dim = c(2, 2, 1))
  y <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))

  metrics <- overlap_metrics(x, y)

  expect_equal(metrics$n_x, 2)
  expect_equal(metrics$n_y, 2)
  expect_equal(metrics$intersection, 1)
  expect_equal(metrics$union, 3)
  expect_equal(metrics$dice, 0.5)
  expect_equal(metrics$jaccard, 1 / 3)
  expect_equal(metrics$overlap_coefficient, 0.5)
  expect_equal(metrics$volume_similarity, 1)
  expect_equal(dice_coefficient(x, y), 0.5)
  expect_equal(jaccard_index(x, y), 1 / 3)
})

test_that("overlap metrics are symmetric except directional rates", {
  x <- array(c(1, 1, 1, 0, 0, 0), dim = c(3, 2, 1))
  y <- array(c(1, 0, 0, 1, 0, 0), dim = c(3, 2, 1))

  xy <- overlap_metrics(x, y)
  yx <- overlap_metrics(y, x)

  expect_equal(xy$dice, yx$dice)
  expect_equal(xy$jaccard, yx$jaccard)
  expect_equal(xy$overlap_coefficient, yx$overlap_coefficient)
  expect_equal(xy$volume_similarity, yx$volume_similarity)
  expect_equal(xy$sensitivity, yx$precision)
  expect_equal(xy$precision, yx$sensitivity)
})

test_that("thresholds and labels define masks explicitly", {
  x <- array(c(0.2, 0.8, 1.2, 0), dim = c(2, 2, 1))
  y <- array(c(0.1, 0.9, 0, 1.5), dim = c(2, 2, 1))

  expect_equal(dice_coefficient(x, y, threshold = 0.5), 0.5)

  lx <- array(c(1, 1, 2, 0), dim = c(2, 2, 1))
  ly <- array(c(1, 2, 2, 0), dim = c(2, 2, 1))
  by_label <- overlap_metrics(lx, ly, labels = c(1, 2))

  expect_equal(by_label$label, c(1, 2))
  expect_equal(by_label$dice, c(2 / 3, 2 / 3))
  expect_equal(dice_coefficient(lx, ly, label = 1), 2 / 3)
  expect_equal(jaccard_index(lx, ly, label = 2), 1 / 2)
})

test_that("empty masks use the requested convention", {
  x <- array(0, dim = c(2, 2, 1))
  y <- array(0, dim = c(2, 2, 1))

  expect_equal(dice_coefficient(x, y), 1)
  expect_equal(dice_coefficient(x, y, empty = "zero"), 0)
  expect_true(is.na(dice_coefficient(x, y, empty = "na")))
})

test_that("missing values are either removed or rejected", {
  x <- array(c(1, NA, 0, 1), dim = c(2, 2, 1))
  y <- array(c(1, 1, 0, NA), dim = c(2, 2, 1))

  metrics <- overlap_metrics(x, y)

  expect_equal(metrics$n_x, 1)
  expect_equal(metrics$n_y, 1)
  expect_equal(metrics$dice, 1)
  expect_error(overlap_metrics(x, y, na.rm = FALSE), "missing or non-finite")
})

test_that("overlap metrics validate compatible inputs", {
  x <- array(1, dim = c(2, 2, 1))
  y <- array(1, dim = c(2, 3, 1))
  z <- array(c("a", "b", "a", "b"), dim = c(2, 2, 1))

  expect_error(overlap_metrics(x, y), "dimensions must match")
  expect_error(overlap_metrics(z, z), "numeric or logical")
  expect_error(overlap_metrics(z, z, labels = c("a", NA)), "must not contain missing")
})

test_that("image affines are checked when available", {
  skip_if_not_installed("neuroim2")

  arr <- array(1, dim = c(2, 2, 2))
  s1 <- neuroim2::NeuroSpace(dim(arr), trans = diag(4))
  s2 <- neuroim2::NeuroSpace(dim(arr), trans = diag(c(2, 1, 1, 1)))
  v1 <- neuroim2::DenseNeuroVol(arr, s1)
  v2 <- neuroim2::DenseNeuroVol(arr, s2)

  expect_error(overlap_metrics(v1, v2), "affines must match")
  expect_equal(dice_coefficient(v1, v2, check_geometry = FALSE), 1)
})
