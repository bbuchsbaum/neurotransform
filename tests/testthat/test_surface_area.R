area_triplets <- function(mat) {
  ij <- which(mat>0,arr.ind=TRUE)
  list(rows=as.integer(ij[,1]),cols=as.integer(ij[,2]),vals=mat[ij])
}
area_dense <- function(w,n,m) {
  mat <- matrix(0,n,m); mat[cbind(w$rows,w$cols)] <- w$vals; mat
}

test_that('adaptive support chooses forward by subset, otherwise all of reverse', {
  F <- matrix(c(.5,.5,0),1)
  contained <- .surface_adaptive_weights(area_triplets(F),area_triplets(matrix(c(.75,.25,0),3)),1,3,rep(1,3),1)
  expect_false(contained$use_reverse)
  expanded <- .surface_adaptive_weights(area_triplets(F),area_triplets(matrix(c(0,.25,.75),3)),1,3,rep(1,3),1)
  expect_true(expanded$use_reverse)
  expect_identical(expanded$cols,c(2L,3L))
  absent <- .surface_adaptive_weights(area_triplets(F),area_triplets(matrix(0,3,1)),1,3,rep(1,3),1)
  expect_false(absent$use_reverse)
})

test_that('unequal areas have independently derived weights and effective measure', {
  A <- rbind(c(1,0),c(.5,.5))
  weights <- .surface_adaptive_weights(area_triplets(A),area_triplets(t(A)),2,2,c(2,3),c(5,7))
  W <- area_dense(weights,2,2)
  expect_equal(W,rbind(c(1,0),c(14/65,51/65)))
  expect_equal(weights$target_measure,c(20/17,65/17))
  x <- c(3,7)
  expect_equal(sum(weights$target_measure*as.vector(W%*%x)),sum(c(2,3)*x))
  equal <- .surface_adaptive_weights(area_triplets(A),area_triplets(t(A)),2,2,c(1,1),c(1,1))
  expect_equal(area_dense(equal,2,2),rbind(c(1,0),c(.25,.75)))
  expect_equal(sum(area_dense(equal,2,2)%*%c(1,0)),1.25) # NOT supplied-area conservation
})

test_that('ROI is applied after adaptive area correction', {
  A <- rbind(c(.5,.5,0),c(0,.5,.5))
  w <- .surface_adaptive_weights(area_triplets(A),area_triplets(t(A)),2,3,rep(1,3),rep(1,2))
  p <- structure(c(w,list(n_reference=2L,n_moving=3L,method='adaptive_bary_area',
    source_mask=c(FALSE,TRUE,TRUE),area=list(target_measure=w$target_measure))),class='SurfaceResamplingPlan')
  expect_equal(apply_surface_resampling(p,diag(3)),rbind(c(0,1,0),c(0,1/3,2/3)))
  d <- apply_surface_resampling(p,c(2,3,7),details=TRUE)
  expect_equal(sum(d$effective_target_area*d$values),10)
})

test_that('area objects require identity, positive measures and provenance', {
  m <- surface_mesh(rbind(c(0,0,0),c(1,0,0),c(0,1,0)),matrix(1:3,1))
  areas <- surface_vertex_areas(m,anatomical=m)
  expect_equal(areas$values,rep(1/6,3))
  expect_true(nzchar(areas$provenance))
  expect_error(surface_vertex_areas(m,values=1:3),'provenance')
  expect_error(surface_vertex_areas(m,values=c(0,1,1),provenance='test'),'positive')
  expect_error(surface_vertex_areas(m,values=c(NA,1,1),provenance='test'),'positive')
  expect_error(surface_vertex_areas(m,values=1:3,anatomical=m),'exactly one')
  changed <- surface_mesh(m@coords*2,m@faces+1L)
  expect_error(.surface_check_areas(areas,.surface_identity(changed),3,'test'),'matching ordered')
})

test_that('experimental adaptive plans match unequal-area Workbench basis and policy fixtures', {
  skip_if_not_installed('jsonlite')
  fixture <- jsonlite::fromJSON(system.file('extdata/adaptive_area_oracle/oracle.json',package='neurotransform'),simplifyVector=FALSE)
  mat <- function(x) matrix(unlist(x),ncol=length(x[[1]]),byrow=TRUE)
  for (case in fixture$cases) {
    moving <- surface_mesh(mat(case$vertices),mat(case$faces)+1L)
    reference <- surface_mesh(mat(case$query),mat(case$query_faces)+1L)
    source <- surface_vertex_areas(moving,values=unlist(case$source_areas),provenance='synthetic unequal-area fixture')
    target <- surface_vertex_areas(reference,values=unlist(case$target_areas),provenance='synthetic unequal-area fixture')
    p <- surface_resampling_plan(reference,moving,method='adaptive_bary_area',source_areas=source,
      target_areas=target,source_mask=unlist(case$source_mask),experimental=TRUE)
    d <- apply_surface_resampling(p,mat(case$data),details=TRUE)
    valid <- as.logical(unlist(case$valid))
    expected <- mat(case$metric)
    expect_identical(d$available[,1],valid,info=case$name)
    expect_equal(d$values[valid,,drop=FALSE],expected[valid,,drop=FALSE],tolerance=fixture$tolerance,info=case$name)
    expect_true(all(is.na(d$values[!valid,,drop=FALSE])))
    weights <- expected[,seq_len(nrow(moving@coords)),drop=FALSE]
    labels <- unlist(case$labels)
    label_sums <- vapply(sort(unique(labels)),function(key) rowSums(weights[,labels==key,drop=FALSE]),numeric(nrow(weights)))
    margin <- function(m) apply(m,1,function(row) diff(tail(sort(row),2)))
    for (policy in c('aggregate','largest')) {
      non_tie <- valid & margin(if(policy=='aggregate') label_sums else weights)>1e-5
      got <- apply_surface_resampling(p,labels,data_type='label',label_method=policy,unassigned=0L)
      expect_equal(got[non_tie],unlist(case[[policy]])[non_tie],info=paste(case$name,policy))
    }
    expect_identical(p$qualification,'experimental')
    expect_true(all(p$area$represented_source))
    expect_equal(sum(d$effective_target_area[valid]*d$values[valid,ncol(d$values)]),
                 sum(source$values[unlist(case$source_mask)>0]),tolerance=1e-12)
    expect_error(surface_resampling_plan(reference,moving,method='adaptive_bary_area'),'experimental')
    expect_error(surface_resampling_plan(reference,moving,method='adaptive_bary_area',experimental=TRUE),'SurfaceVertexAreas')
    reverse <- reverse_surface_resampling_plan(p,reference,moving)
    expect_identical(reverse$area$source$geometry,p$area$target$geometry)
  }
})

test_that('adaptive weights are invariant to area units at extreme common scales', {
  A <- rbind(c(1,0),c(.5,.5))
  for (s in c(1e-200,1e200)) {
    w <- .surface_adaptive_weights(area_triplets(A),area_triplets(t(A)),2,2,c(2,3)*s,c(5,7)*s)
    expect_equal(area_dense(w,2,2),rbind(c(1,0),c(14/65,51/65)))
    expect_equal(w$target_measure/s,c(20/17,65/17))
  }
})

test_that('adaptive effective measures follow omission and normalization semantics', {
  A <- rbind(c(.5,.5,0),c(0,.5,.5))
  w <- .surface_adaptive_weights(area_triplets(A),area_triplets(t(A)),2,3,rep(1,3),rep(1,2))
  p <- structure(c(w,list(n_reference=2L,n_moving=3L,method='adaptive_bary_area',
    area=list(target_measure=w$target_measure))),class='SurfaceResamplingPlan')
  d <- apply_surface_resampling(p,c(NA,3,7),na_policy='omit',details=TRUE)
  expect_equal(d$values,c(3,17/3))
  expect_equal(d$effective_target_area,c(.5,1.5))
  expect_equal(d$source_roi_target_area,c(1.5,1.5))
  expect_equal(sum(d$values*d$effective_target_area),10)
  dm <- apply_surface_resampling(p,cbind(c(NA,3,7),c(2,3,NA)),na_policy='omit',details=TRUE)
  expect_equal(dm$effective_target_area,cbind(c(.5,1.5),c(1.5,.5)))
  expect_equal(colSums(dm$values*dm$effective_target_area),c(10,5))
  for (normal in c('none','sum')) {
    dn <- apply_surface_resampling(p,c(2,3,7),normalize=normal,details=TRUE)
    expect_null(dn$effective_target_area)
    expect_equal(dn$source_roi_target_area,c(1.5,1.5))
  }
})
