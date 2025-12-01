test_that("cpp_path_apply_steps matches manual affine + warp + affine", {
  coords <- matrix(c(0, 0, 0,
                     1, 1, 1), ncol = 3, byrow = TRUE)

  aff1 <- diag(4); aff1[1, 4] <- 1       # +1 in x
  aff2 <- diag(4); aff2[3, 4] <- -1      # -1 in z

  dimf <- c(2L, 2L, 2L)
  disp <- array(0, dim = c(3, dimf))
  disp[1, , , ] <- 0.5
  disp[2, , , ] <- -0.25
  disp[3, , , ] <- 0.1
  field <- as.numeric(disp)

  steps <- list(
    list(kind = "affine", matrix = aff1),
    list(kind = "warp", field = field, dim = dimf, world_to_vox = diag(4), method = "linear"),
    list(kind = "affine", matrix = aff2)
  )

  out_cpp <- cpp_path_apply_steps(coords, steps)

  expected <- rbind(
    c(1.5, -0.25, -0.9),
    c(NA_real_, NA_real_, NA_real_)
  )

  expect_equal(out_cpp, expected, tolerance = 1e-8)
})

test_that("cpp_path_apply_steps handles multiple warps", {
  coords <- matrix(c(0.2, 0.3, 0.4), ncol = 3)

  dimf <- c(2L, 2L, 2L)
  make_field <- function(dx, dy, dz) {
    arr <- array(0, dim = c(3, dimf))
    arr[1, , , ] <- dx
    arr[2, , , ] <- dy
    arr[3, , , ] <- dz
    as.numeric(arr)
  }

  steps <- list(
    list(kind = "warp", field = make_field(0.1, 0, 0), dim = dimf, world_to_vox = diag(4), method = "linear"),
    list(kind = "warp", field = make_field(0, 0.2, 0), dim = dimf, world_to_vox = diag(4), method = "linear")
  )

  out_cpp <- cpp_path_apply_steps(coords, steps)
  expect_equal(out_cpp, coords + c(0.1, 0.2, 0), tolerance = 1e-8)
})
