# Called by compare_adaptive_templates.py against an explicitly installed build.
library(neurotransform)
library(jsonlite)
args <- commandArgs(TRUE); folder <- args[1]; tolerance <- as.numeric(args[2])
read <- function(name) as.matrix(read.csv(file.path(folder,paste0(name,'.csv')),header=FALSE))
config <- fromJSON(file.path(folder,'config.json'))
v <- read('vertices'); f <- read('faces')+1L; q <- read('query'); qf <- read('query_faces')+1L
x <- read('data'); expected <- read('expected'); labels <- read('labels')[,1]
source_mask <- read('roi')[,1]; expected_valid <- read('valid')[,1]>0
expected_labels <- read('expected_labels')
reference <- surface_mesh(q,qf)
make <- function(faces) {
  moving <- surface_mesh(v,faces)
  source <- surface_vertex_areas(moving,values=read('source_area')[,1],units='mm^2',provenance=config$source_area_sha256)
  target <- surface_vertex_areas(reference,values=read('target_area')[,1],units='mm^2',provenance=config$target_area_sha256)
  surface_resampling_plan(reference,moving,method='adaptive_bary_area',source_areas=source,
    target_areas=target,source_mask=source_mask,experimental=TRUE)
}
elapsed <- system.time(p <- make(f))['elapsed']
application <- system.time(d <- apply_surface_resampling(p,x,details=TRUE))['elapsed']
valid <- d$available[,1]
coverage_equal <- identical(valid,expected_valid)
comparable <- valid & expected_valid
errors <- abs(d$values[comparable,,drop=FALSE]-expected[comparable,,drop=FALSE])
worst <- arrayInd(order(errors,decreasing=TRUE)[seq_len(min(20,length(errors)))],dim(errors))
worst[,1] <- which(comparable)[worst[,1]]
reordered <- make(f[nrow(f):1,,drop=FALSE])
y2 <- apply_surface_resampling(reordered,x)
permutation_error <- max(abs(y2[valid,,drop=FALSE]-d$values[valid,,drop=FALSE]))
label_results <- lapply(c('aggregate','largest'),function(method) apply_surface_resampling(p,labels,
  data_type='label',label_method=method,unassigned=0L))
label_mismatch <- vapply(seq_len(2),function(i) sum((label_results[[i]]!=expected_labels[,i]) & valid),integer(1))
label_worst <- lapply(seq_len(2),function(i) head(which((label_results[[i]]!=expected_labels[,i]) & valid),30))
# This identity uses the effective target measure h, not supplied target area b.
left <- colSums(d$values[valid,,drop=FALSE]*d$effective_target_area[valid])
right <- colSums(x*(p$area$source$values*(source_mask>0)))
scale <- pmax(1,colSums(abs(x)*(p$area$source$values*(source_mask>0))))
area_error <- max(abs(left-right)/scale)
row_sums <- numeric(nrow(q)); grouped <- rowsum(p$vals,p$rows,reorder=FALSE)
row_sums[as.integer(rownames(grouped))] <- grouped[,1]
geometry_valid <- all(is.finite(p$vals)) && all(p$vals>0) && max(abs(row_sums-1))<1e-12 && all(p$area$represented_source)
result <- list(max_absolute_error=max(errors),rms_error=sqrt(mean(errors^2)),
  error_quantiles=quantile(errors,c(.5,.95,.99,1)),field_max_errors=apply(errors,2,max),
  worst=worst, coverage_equal=coverage_equal, unsupported_outputs=sum(!valid),
  label_mismatches=setNames(as.list(label_mismatch),c('aggregate','largest')),
  unsupported_workbench_label_keys=lapply(seq_len(2),function(i) sort(unique(expected_labels[!valid,i]))),
  native_unassigned_label_key=0L,
  label_mismatch_examples=setNames(label_worst,c('aggregate','largest')),
  effective_area_identity_scaled_error=area_error,geometry_valid=geometry_valid,
  max_row_sum_error=max(abs(row_sums-1)),permutation_error=permutation_error,
  elapsed_seconds=unname(elapsed),application_seconds=unname(application),native_timing=p$timing,
  targets=nrow(q),nnz=length(p$vals),reverse_selected_rows=sum(p$area$use_reverse),
  passed=coverage_equal && geometry_valid && max(errors)<tolerance && all(label_mismatch==0) &&
    permutation_error<1e-12 && area_error<1e-12,
  dll_sha256=digest::digest(file=getLoadedDLLs()[['neurotransform']][['path']],algo='sha256'),
  session=capture.output(sessionInfo()))
write_json(result,if(length(args)>2L) args[3] else file.path(folder,'result.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
quit(status=if(result$passed) 0L else 1L)
