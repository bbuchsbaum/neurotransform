policy_plan <- function(source_mask = NULL, target_mask = NULL) {
  moving <- surface_mesh(rbind(c(0,0,0), c(1,0,0), c(0,1,0)), matrix(1:3,1))
  reference <- surface_mesh(rbind(c(.3,.4,0), c(0,0,0), c(2,2,0)))
  surface_resampling_plan(reference, moving, spherical=FALSE, outside='missing',
                          source_mask=source_mask, target_mask=target_mask)
}

test_that('coverage distinguishes geometry, ROI, missingness and numerical zero', {
  p <- policy_plan(source_mask=c(1,.1,0), target_mask=c(1,0,1))
  d <- apply_surface_resampling(p, c(0,2,99), details=TRUE)
  expect_equal(d$values, c(1,NA,NA))
  expect_equal(d$source_weight_mass, c(.6,1,0))
  expect_equal(d$finite_weight_mass, c(.6,1,0))
  expect_identical(d$status, c('available','masked_target','no_geometry'))
  expect_equal(apply_surface_resampling(policy_plan(), c(0,0,0)), c(0,0,NA))
  expect_equal(apply_surface_resampling(policy_plan(source_mask=c(0,0,0)), 1:3), rep(NA_real_,3))
  expect_equal(apply_surface_resampling(policy_plan(source_mask=c(0,1,0)), 1:3), c(2,NA,NA))
  expect_error(policy_plan(source_mask=c(1,NA,0)), 'finite')
  expect_error(policy_plan(target_mask=TRUE), 'length')
})

test_that('missing policies are column specific and ignore zero contributions', {
  p <- policy_plan()
  x <- cbind(c(0,2,NA), c(0,NaN,4), c(0,2,Inf))
  expect_equal(apply_surface_resampling(p,x), matrix(c(NA,0,NA),3,3))
  d <- apply_surface_resampling(p,x,na_policy='omit',details=TRUE)
  expect_equal(d$values, rbind(c(1,16/7,1),c(0,0,0),c(NA,NA,NA)))
  expect_equal(d$finite_weight_mass, rbind(c(.6,.7,.6),c(1,1,1),c(0,0,0)))
  z <- apply_surface_resampling(p,x[,1],na_policy='omit')
  expect_equal(attr(z,'retained_weight_mass'), c(.6,1,0))
  expect_error(apply_surface_resampling(p,x,na_policy='error'), 'positive weight')
  expect_error(apply_surface_resampling(p,x,na_policy='omit',normalize='sum'), 'requires')
  expect_equal(apply_surface_resampling(p,rep(NA_real_,3)), rep(NA_real_,3))
  expect_equal(apply_surface_resampling(policy_plan(target_mask=c(0,1,0)),x,na_policy='error'),
               matrix(c(NA,0,NA),3,3))
  expect_identical(dim(apply_surface_resampling(p,matrix(1:3,3,1))), c(3L,1L))
})

test_that('categorical voting preserves keys and metadata with deterministic ties', {
  p <- policy_plan()
  table <- data.frame(key=c(0L,10L,20L), name=c('none','A','B'), red=c(0,1,0))
  aggregate <- apply_surface_resampling(p,c(10,10,20),data_type='label',label_table=table,unassigned=0L)
  expect_identical(as.vector(aggregate),c(10L,10L,0L))
  expect_identical(attr(aggregate,'label_table'),table)
  expect_identical(apply_surface_resampling(p,c(10,10,20),data_type='label',label_method='largest'),c(20L,10L,NA_integer_))
  # Exact dyadic ties avoid declaring almost-equal floating weights equal.
  p$vals <- c(.25,.25,.5,1)
  expect_equal(apply_surface_resampling(p,c(10,10,20),data_type='label')[1],10L)
  p$vals <- c(.5,.5,0,1)
  expect_equal(apply_surface_resampling(p,c(20,10,30),data_type='label',label_method='largest')[1],20L)
  expect_error(apply_surface_resampling(p,c(1.5,2,3),data_type='label'),'integer')
  expect_error(apply_surface_resampling(p,c(1,2,3),data_type='label',label_table=table),'every')
  expect_error(apply_surface_resampling(p,1:3,data_type='label',label_table=data.frame(key=c(1,1))), 'unique')
})

test_that('outside policies expose nearest fallback explicitly', {
  moving <- surface_mesh(rbind(c(0,0,0),c(1,0,0),c(0,1,0)),matrix(1:3,1))
  ref <- surface_mesh(matrix(c(2,2,0),1))
  expect_error(surface_resampling_plan(ref,moving,spherical=FALSE,outside='error'),'without valid')
  p <- surface_resampling_plan(ref,moving,spherical=FALSE,outside='nearest')
  expect_identical(p$support,'nearest_fallback')
  expect_false(apply_surface_resampling(p,1:3,details=TRUE)$geometric_support)
})

test_that('ordinary masks and non-tied labels agree with the independent Workbench fixture', {
  skip_if_not_installed('jsonlite')
  path <- system.file('extdata/surface_policy_oracle/oracle.json',package='neurotransform')
  fixture <- jsonlite::fromJSON(path,simplifyVector=TRUE)
  moving <- surface_mesh(fixture$vertices,fixture$faces+1L)
  reference <- surface_mesh(fixture$query,fixture$query_faces+1L)
  for (i in seq_len(nrow(fixture$cases))) {
    case <- fixture$cases[i,]
    p <- surface_resampling_plan(reference,moving,source_mask=case$source_mask[[1]])
    d <- apply_surface_resampling(p,fixture$data,details=TRUE)
    valid <- as.logical(case$valid[[1]])
    expect_identical(d$available[,1],valid,info=case$name)
    expect_equal(d$values[valid,,drop=FALSE],case$metric[[1]][valid,,drop=FALSE],
                 tolerance=fixture$tolerance,info=case$name)
    expect_true(all(is.na(d$values[!valid,,drop=FALSE])))
    # Exclude actual weight ties prospectively from parity; library ties have a
    # separately specified deterministic contract tested above.
    weights <- case$metric[[1]][,seq_len(nrow(moving@coords)),drop=FALSE]
    label_sums <- vapply(sort(unique(fixture$labels)),function(key)
      rowSums(weights[,fixture$labels==key,drop=FALSE]),numeric(nrow(weights)))
    margin <- function(mat) apply(mat,1,function(row) diff(tail(sort(row),2)))
    for (policy in c('aggregate','largest')) {
      non_tie <- valid & margin(if(policy=='aggregate') label_sums else weights)>1e-5
      got <- apply_surface_resampling(p,fixture$labels,data_type='label',label_method=policy,unassigned=0L)
      expect_equal(got[non_tie],case[[policy]][[1]][non_tie],info=paste(case$name,policy))
      expect_true(all(got[!valid]==0L))
    }
  }
})
