#' Image overlap metrics
#'
#' Compute voxelwise overlap summaries for two images or masks. Inputs may be
#' plain arrays or image objects that can be coerced with `as.array()`, including
#' `neuroim2` volumes.
#'
#' By default, numeric images are converted to foreground masks with
#' `value > threshold`, and logical images use `TRUE` as foreground. Supply
#' `labels` to compute label-wise overlap instead.
#'
#' @param x,y Arrays or image objects with matching dimensions.
#' @param threshold Numeric scalar used to form foreground masks when `labels`
#'   is `NULL`.
#' @param labels Optional vector of label values to compare. When supplied,
#'   `threshold` is ignored and one row is returned per label.
#' @param na.rm Logical; if `TRUE`, paired missing or non-finite numeric voxels
#'   are removed before computing masks.
#' @param empty How to score Dice, Jaccard, overlap coefficient, and volume
#'   similarity when both masks are empty. `"one"` treats two empty masks as
#'   perfect agreement, `"zero"` returns 0, and `"na"` returns `NA_real_`.
#' @param check_geometry Logical; if `TRUE`, compare available voxel-to-world
#'   affines after confirming that dimensions match.
#' @param tolerance Numeric tolerance used for affine comparisons.
#' @return A data frame with voxel counts and overlap metrics. `sensitivity`
#'   is `intersection / n_x`; `precision` is `intersection / n_y`.
#' @export
#' @examples
#' x <- array(c(1, 1, 0, 0), dim = c(2, 2, 1))
#' y <- array(c(1, 0, 1, 0), dim = c(2, 2, 1))
#' overlap_metrics(x, y)
overlap_metrics <- function(x, y,
                            threshold = 0,
                            labels = NULL,
                            na.rm = TRUE,
                            empty = c("one", "zero", "na"),
                            check_geometry = TRUE,
                            tolerance = sqrt(.Machine$double.eps)) {
  empty <- match.arg(empty)
  empty_value <- .overlap_empty_value(empty)

  if (!isTRUE(na.rm) && !identical(na.rm, FALSE)) {
    stop("na.rm must be TRUE or FALSE")
  }
  if (!is.logical(check_geometry) || length(check_geometry) != 1L || is.na(check_geometry)) {
    stop("check_geometry must be TRUE or FALSE")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a non-negative finite numeric scalar")
  }

  xi <- .overlap_image_info(x, "x")
  yi <- .overlap_image_info(y, "y")
  .check_overlap_compatibility(xi, yi, check_geometry = check_geometry, tolerance = tolerance)

  xv <- as.vector(xi$data)
  yv <- as.vector(yi$data)
  valid <- .overlap_valid_pair(xv, yv)
  if (isTRUE(na.rm)) {
    xv <- xv[valid]
    yv <- yv[valid]
  } else if (any(!valid)) {
    stop("x and y must not contain missing or non-finite values when na.rm = FALSE")
  }

  if (!is.null(labels)) {
    labels <- .overlap_validate_labels(labels)
    rows <- lapply(labels, function(label) {
      row <- .overlap_metric_row(
        .overlap_label_mask(xv, label),
        .overlap_label_mask(yv, label),
        empty_value
      )
      cbind(data.frame(label = label, stringsAsFactors = FALSE), row)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  mx <- .overlap_threshold_mask(xv, threshold, "x")
  my <- .overlap_threshold_mask(yv, threshold, "y")
  .overlap_metric_row(mx, my, empty_value)
}

#' Dice coefficient for two images
#'
#' Convenience wrapper returning the Dice coefficient from `overlap_metrics()`.
#'
#' @inheritParams overlap_metrics
#' @param label Optional single label value. If supplied, Dice is computed for
#'   that label instead of by thresholding.
#' @return Numeric scalar Dice coefficient.
#' @export
#' @examples
#' x <- array(c(TRUE, TRUE, FALSE, FALSE), dim = c(2, 2, 1))
#' y <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))
#' dice_coefficient(x, y)
dice_coefficient <- function(x, y,
                             threshold = 0,
                             label = NULL,
                             na.rm = TRUE,
                             empty = c("one", "zero", "na"),
                             check_geometry = TRUE,
                             tolerance = sqrt(.Machine$double.eps)) {
  if (!is.null(label) && length(label) != 1L) {
    stop("label must be a single value")
  }
  overlap_metrics(
    x, y,
    threshold = threshold,
    labels = label,
    na.rm = na.rm,
    empty = empty,
    check_geometry = check_geometry,
    tolerance = tolerance
  )$dice[[1L]]
}

#' Jaccard index for two images
#'
#' Convenience wrapper returning the Jaccard index from `overlap_metrics()`.
#'
#' @inheritParams dice_coefficient
#' @return Numeric scalar Jaccard index.
#' @export
#' @examples
#' x <- array(c(1, 1, 0, 0), dim = c(2, 2, 1))
#' y <- array(c(1, 0, 1, 0), dim = c(2, 2, 1))
#' jaccard_index(x, y)
jaccard_index <- function(x, y,
                          threshold = 0,
                          label = NULL,
                          na.rm = TRUE,
                          empty = c("one", "zero", "na"),
                          check_geometry = TRUE,
                          tolerance = sqrt(.Machine$double.eps)) {
  if (!is.null(label) && length(label) != 1L) {
    stop("label must be a single value")
  }
  overlap_metrics(
    x, y,
    threshold = threshold,
    labels = label,
    na.rm = na.rm,
    empty = empty,
    check_geometry = check_geometry,
    tolerance = tolerance
  )$jaccard[[1L]]
}

.overlap_image_info <- function(x, name) {
  image_like <- inherits(x, "DenseNeuroVol") ||
    inherits(x, "DenseNeuroVec") ||
    inherits(x, "niftiImage")

  if (is.array(x) && !image_like) {
    return(list(data = x, dims = dim(x), affine = diag(4), has_geometry = FALSE))
  }

  data <- tryCatch(as.array(x), error = function(e) NULL)
  if (!is.array(data) || is.null(dim(data))) {
    stop(name, " must be an array or image object coercible with as.array()")
  }

  affine <- tryCatch(extract_affine(x), error = function(e) NULL)
  list(
    data = data,
    dims = dim(data),
    affine = affine %||% diag(4),
    has_geometry = !is.null(affine)
  )
}

.check_overlap_compatibility <- function(xi, yi, check_geometry, tolerance) {
  if (!identical(as.integer(xi$dims), as.integer(yi$dims))) {
    stop(
      "x and y dimensions must match; got ",
      paste(xi$dims, collapse = "x"),
      " and ",
      paste(yi$dims, collapse = "x")
    )
  }

  if (!isTRUE(check_geometry)) return(invisible(TRUE))
  if (!isTRUE(xi$has_geometry) && !isTRUE(yi$has_geometry)) return(invisible(TRUE))

  same_affine <- isTRUE(all.equal(
    xi$affine,
    yi$affine,
    tolerance = tolerance,
    check.attributes = FALSE
  ))
  if (!same_affine) {
    stop("x and y affines must match; set check_geometry = FALSE to compare by array index")
  }
  invisible(TRUE)
}

.overlap_valid_pair <- function(x, y) {
  .overlap_finite_or_present(x) & .overlap_finite_or_present(y)
}

.overlap_finite_or_present <- function(x) {
  if (is.numeric(x) || is.complex(x)) {
    return(is.finite(x))
  }
  !is.na(x)
}

.overlap_validate_labels <- function(labels) {
  if (!is.atomic(labels) || length(labels) < 1L) {
    stop("labels must be a non-empty atomic vector")
  }
  labels <- as.vector(labels)
  if (any(is.na(labels))) {
    stop("labels must not contain missing values")
  }
  labels
}

.overlap_threshold_mask <- function(values, threshold, name) {
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("threshold must be a finite numeric scalar")
  }
  if (is.logical(values)) {
    values <- as.integer(values)
  }
  if (!is.numeric(values)) {
    stop(name, " must be numeric or logical when labels is NULL")
  }
  values > threshold
}

.overlap_label_mask <- function(values, label) {
  values == label
}

.overlap_empty_value <- function(empty) {
  switch(empty,
    one = 1,
    zero = 0,
    na = NA_real_
  )
}

.overlap_metric_row <- function(mask_x, mask_y, empty_value) {
  n_x <- sum(mask_x)
  n_y <- sum(mask_y)
  intersection <- sum(mask_x & mask_y)
  union <- sum(mask_x | mask_y)
  both_empty <- n_x == 0L && n_y == 0L

  dice <- if (both_empty) empty_value else (2 * intersection) / (n_x + n_y)
  jaccard <- if (both_empty) empty_value else intersection / union
  overlap_coefficient <- if (both_empty) {
    empty_value
  } else if (min(n_x, n_y) == 0L) {
    0
  } else {
    intersection / min(n_x, n_y)
  }
  volume_similarity <- if (both_empty) {
    empty_value
  } else {
    1 - abs(n_x - n_y) / (n_x + n_y)
  }
  sensitivity <- if (n_x == 0L) {
    if (both_empty) empty_value else NA_real_
  } else {
    intersection / n_x
  }
  precision <- if (n_y == 0L) {
    if (both_empty) empty_value else NA_real_
  } else {
    intersection / n_y
  }

  data.frame(
    n_x = n_x,
    n_y = n_y,
    intersection = intersection,
    union = union,
    dice = dice,
    jaccard = jaccard,
    overlap_coefficient = overlap_coefficient,
    volume_similarity = volume_similarity,
    sensitivity = sensitivity,
    precision = precision
  )
}
