test_that("warp Jacobians differentiate physical RAS coordinates on oriented grids", {
  dims <- c(9L, 9L, 9L)
  B <- matrix(c(.1, .02, -.03, .04, -.08, .01, -.02, .03, .12), 3, byrow = TRUE)
  theta <- .37
  rotation <- matrix(c(cos(theta), -sin(theta), 0,
                       sin(theta), cos(theta), 0, 0, 0, 1), 3, byrow = TRUE)
  for (hand in c(-1, 1)) {
    affine <- diag(4)
    affine[1:3, 1:3] <- rotation %*% diag(c(hand * 1.3, -1.7, 2.1))
    affine[1:3, 3] <- affine[1:3, 3] + .1 * affine[1:3, 1]
    affine[1:3, 4] <- c(12, -8, 3)
    ijk <- as.matrix(expand.grid(0:8, 0:8, 0:8))
    world <- (cbind(ijk, 1) %*% t(affine))[, 1:3]
    displacement <- world %*% t(B)
    points <- (cbind(rbind(c(3.2, 3.7, 4.1), c(4.3, 4.2, 3.8)), 1) %*% t(affine))[, 1:3]
    for (method in c("linear", "cubic")) {
      m <- warp_from_field("source", "target", array(displacement, c(dims, 3)),
                           grid = grid_spec(dims, affine))
      m@params$warp_method <- method
      expected <- diag(3) + B
      J <- jacobian(m, points)
      for (i in seq_len(nrow(points))) {
        expect_equal(J@values[i, , ], expected, tolerance = 1e-10)
      }
      expect_equal(jacobian_det(m, points), rep(det(expected), nrow(points)), tolerance = 1e-10)
      expect_equal(transform(m, points), points %*% t(expected), ignore_attr = TRUE,
                   tolerance = 1e-10)
    }
  }
})
