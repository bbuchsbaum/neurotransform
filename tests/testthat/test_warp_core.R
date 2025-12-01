test_that("apply warp field translates coords and handles bounds", {
  dim <- c(2L, 2L, 2L)
  disp <- array(0, dim = c(3, dim))
  disp[1, , , ] <- 1
  disp[2, , , ] <- -0.5
  disp[3, , , ] <- 2
  field <- as.numeric(disp)

  coords <- matrix(c(
    0, 0, 0,
    0.5, 0.5, 0.5,
    -1, 0, 0
  ), ncol = 3, byrow = TRUE)

  out <- neurotransform:::cpp_apply_warp_field(coords, field, dim, diag(4))

  expect_equal(out[1, ], c(1, -0.5, 2))
  expect_equal(out[2, ], c(1.5, 0, 2.5), tolerance = 1e-8)
  expect_true(all(is.na(out[3, ])))
})

test_that("warp composition matches sequential application for translations", {
  dim <- c(2L, 2L, 2L)
  make_field <- function(dx, dy, dz) {
    arr <- array(0, dim = c(3, dim))
    arr[1, , , ] <- dx
    arr[2, , , ] <- dy
    arr[3, , , ] <- dz
    as.numeric(arr)
  }

  fieldA <- make_field(1, 0, 0)
  fieldB <- make_field(0, 1, 0)
  vox_affine <- diag(4)

  composed <- neurotransform:::cpp_compose_warp_fields(fieldA, dim, vox_affine, fieldB, dim, vox_affine)

  coords <- matrix(c(
    0, 0, 0,
    1, 1, 1
  ), ncol = 3, byrow = TRUE)

  seq_AB <- neurotransform:::cpp_apply_warp_field(
    neurotransform:::cpp_apply_warp_field(coords, fieldB, dim, vox_affine),
    fieldA, dim, vox_affine
  )
  comp_out <- neurotransform:::cpp_apply_warp_field(coords, composed$field, composed$dim, vox_affine)

  finite <- is.finite(seq_AB)
  expect_equal(comp_out[finite], seq_AB[finite], tolerance = 1e-8)
})

test_that("warp jacobian recovers linear gradients", {
  dim <- c(3L, 3L, 3L)
  disp <- array(0, dim = c(3, dim))
  for (x in 0:(dim[1] - 1)) {
    for (y in 0:(dim[2] - 1)) {
      for (z in 0:(dim[3] - 1)) {
        disp[1, x + 1, y + 1, z + 1] <- 0.1 * x
        disp[2, x + 1, y + 1, z + 1] <- 0.2 * y
        disp[3, x + 1, y + 1, z + 1] <- -0.05 * z
      }
    }
  }
  field <- as.numeric(disp)
  coords <- matrix(c(1, 1, 1), ncol = 3)

  J <- neurotransform:::cpp_warp_jacobian(coords, field, dim, diag(4), diag(4))
  Jmat <- J[1, , ]
  expected <- diag(c(1.1, 1.2, 0.95))

  expect_equal(Jmat, expected, tolerance = 1e-6)

  dets <- neurotransform:::cpp_warp_jacobian_det(coords, field, dim, diag(4), diag(4))
  expect_equal(as.numeric(dets), prod(diag(expected)), tolerance = 1e-6)
})
