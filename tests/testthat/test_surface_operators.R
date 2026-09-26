test_that('adjoint is the transpose of a frozen masked rectangular operator', {
  p <- structure(list(rows=c(1L,1L,2L,2L,3L),cols=c(1L,2L,2L,3L,1L),
    vals=c(.25,.75,.5,.5,1),n_reference=3L,n_moving=4L,method='barycentric',
    source_mask=c(TRUE,TRUE,FALSE,TRUE),target_mask=c(TRUE,FALSE,TRUE)),class='SurfaceResamplingPlan')
  x <- c(2,3,7,11); y <- c(5,7,13)
  # Independently specified W: source3 removed, target2 removed after normalization.
  operators <- list(element=rbind(c(.25,.75,0,0),c(0,0,0,0),c(1,0,0,0)),
                    none=rbind(c(.25,.75,0,0),c(0,0,0,0),c(1,0,0,0)),
                    sum=rbind(c(.2,.6,0,0),c(0,0,0,0),c(.8,0,0,0)))
  for (normal in names(operators)) {
    W <- operators[[normal]]
    forward <- apply_surface_resampling(p,x,normalize=normal)
    expect_equal(forward[c(1,3)],as.vector(W%*%x)[c(1,3)])
    adj <- apply_surface_adjoint(p,y,normalize=normal)
    expect_equal(adj,as.vector(t(W)%*%y))
    expect_equal(sum((W%*%x)*y),sum(x*adj))
  }
  expect_identical(apply_surface_adjoint(p,y,details=TRUE)$available,c(TRUE,TRUE,FALSE,FALSE))
  expect_equal(apply_surface_adjoint(p,c(5,NA,13)),apply_surface_adjoint(p,y))
  expect_error(apply_surface_adjoint(p,c(NA,7,13)),'finite')
  expect_error(apply_surface_adjoint(p,y,inverse=TRUE),'unused argument')
})

test_that('geometric reverse is rebuilt and verifies ordered identities', {
  a <- surface_mesh(rbind(c(0,0,0),c(1,0,0),c(0,1,0)),matrix(1:3,1))
  b <- surface_mesh(rbind(c(.1,.1,0),c(.8,.1,0),c(.1,.8,0)),matrix(1:3,1))
  p <- surface_resampling_plan(b,a,spherical=FALSE)
  rev <- reverse_surface_resampling_plan(p,b,a)
  direct <- surface_resampling_plan(a,b,spherical=FALSE)
  expect_equal(rev$vals,direct$vals)
  expect_equal(rev$cols,direct$cols)
  expect_equal(apply_surface_resampling(rev,1:3),c(1,2,3))
  expect_false(isTRUE(all.equal(apply_surface_resampling(rev,1:3),apply_surface_adjoint(p,1:3))))
  expect_error(reverse_surface_resampling_plan(p,a,b),'recorded ordered')
  nofaces <- surface_mesh(b@coords)
  p2 <- surface_resampling_plan(nofaces,a,spherical=FALSE)
  expect_error(reverse_surface_resampling_plan(p2,nofaces,a),'requires triangle faces')
  expect_warning(legacy <- apply_surface_resampling(p,1:3,inverse=TRUE),'deprecated')
  W <- rbind(c(.8,.1,.1),c(.1,.8,.1),c(.1,.1,.8))
  expect_equal(legacy,as.vector(t(W)%*%(1:3)))
  expect_error(apply_surface_resampling(p,1:3,inverse=TRUE,details=TRUE),'legacy inverse')
})
