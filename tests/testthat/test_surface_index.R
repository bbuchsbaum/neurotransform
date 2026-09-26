index_triplets <- function(w,n,m) {
  result <- matrix(0,n,m)
  result[cbind(w$rows,w$cols)] <- w$vals
  result
}

test_that("indexed and exhaustive triangle search agree across query regimes", {
  set.seed(70926)
  vertices <- matrix(rnorm(3*192),ncol=3)
  faces <- matrix(0:191,ncol=3,byrow=TRUE)
  query <- rbind(matrix(rnorm(900),ncol=3),vertices,
                 (vertices[faces[,1]+1,]+vertices[faces[,2]+1,])/2)
  for (closest in c(TRUE,FALSE)) for (scale in c(1e-6,1,1e6)) {
    a <- neurotransform:::cpp_barycentric_weights(query*scale,vertices*scale,faces,closest,TRUE)
    b <- neurotransform:::cpp_barycentric_weights(query*scale,vertices*scale,faces,closest,FALSE)
    expect_identical(a$face,b$face)
    expect_identical(a$distance,b$distance)
    expect_identical(a$rows,b$rows)
    expect_identical(a$cols,b$cols)
    expect_identical(a$vals,b$vals)
    expect_lt(sum(a$triangle_tests),sum(b$triangle_tests))
  }
  p <- sample.int(nrow(faces))
  a <- neurotransform:::cpp_barycentric_weights(query,vertices,faces[p,3:1],TRUE,TRUE)
  b <- neurotransform:::cpp_barycentric_weights(query,vertices,faces,TRUE,FALSE)
  expect_equal(index_triplets(a,nrow(query),nrow(vertices)),
               index_triplets(b,nrow(query),nrow(vertices)),tolerance=1e-12)
})

test_that("box pruning survives cancellation near large triangles", {
  near <- rbind(c(.5,-1,-1),c(.5,1,-1),c(.5,0,1))
  giant <- rbind(c(1e16,0,0),c(1,0,0),c(1e16,1e16,0))
  vertices <- do.call(rbind,c(rep(list(near),4),list(giant),
                              lapply(1:4,function(i) sweep(near,2,c(i*10,0,0),"+"))))
  faces <- matrix(seq_len(nrow(vertices))-1L,ncol=3,byrow=TRUE)
  query <- matrix(c(0,0,0),1)
  a <- neurotransform:::cpp_barycentric_weights(query,vertices,faces,TRUE,TRUE)
  b <- neurotransform:::cpp_barycentric_weights(query,vertices,faces,TRUE,FALSE)
  expect_equal(a$distance,.5)
  expect_identical(a$face,b$face)
  expect_identical(a$vals,b$vals)
  expect_equal(a$face,1L)
})

test_that("duplicate triangle face diagnostics have a stable tie break", {
  v <- rbind(c(0,0,0),c(1,0,0),c(0,1,0))
  f <- matrix(rep(0:2,each=20),ncol=3)
  q <- matrix(c(.2,.2,.1),1)
  for (indexed in c(FALSE,TRUE)) {
    w <- neurotransform:::cpp_barycentric_weights(q,v,f,TRUE,indexed)
    expect_equal(w$face,1L)
  }
})

test_that("surface samplers rebuild their immutable index after serialization", {
  v <- rbind(c(0,0,0),c(1,0,0),c(0,1,0))
  sampler <- surface_sampler(v,c(10,20,30),matrix(1:3,1),"barycentric")
  q <- matrix(c(.25,.25,.1),1)
  restored <- unserialize(serialize(sampler,NULL))
  expect_equal(restored@evaluate(q),17.5)
  expect_equal(restored@evaluate(q),sampler@evaluate(q))
  expect_false(neurotransform:::cpp_surface_index_valid(NULL))
})

test_that('subnormal squared distances cannot falsely prune tied planes', {
  for (s in c(1e-162,1e-200)) {
    vertices <- do.call(rbind,lapply(seq_len(10),function(i)
      cbind(rbind(c(0,0),c(1,0),c(0,1)),if(i%%2) 0 else 2*s)))
    faces <- matrix(seq_len(nrow(vertices))-1L,ncol=3,byrow=TRUE)
    query <- matrix(c(.25,.25,s),1)
    for (f in list(faces,faces[nrow(faces):1,,drop=FALSE])) {
      indexed <- cpp_barycentric_weights(query,vertices,f,indexed=TRUE)
      exhaustive <- cpp_barycentric_weights(query,vertices,f,indexed=FALSE)
      expect_identical(indexed$cols,exhaustive$cols)
      expect_identical(indexed$vals,exhaustive$vals)
      expect_identical(indexed$distance,exhaustive$distance)
      expect_identical(indexed$face,exhaustive$face)
    }
  }
})
