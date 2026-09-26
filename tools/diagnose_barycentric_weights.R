library(neurotransform)
library(jsonlite)
base <- commandArgs(TRUE)[1]
for (folder in list.dirs(base, recursive=FALSE)) {
 read <- function(name) as.matrix(read.csv(file.path(folder,paste0(name,'.csv')),header=FALSE))
 v <- read('vertices'); f <- read('faces'); q <- read('query'); x <- read('data'); expected <- read('expected')
 p <- surface_resampling_plan(surface_mesh(q), surface_mesh(v,f+1L))
 got <- apply_surface_resampling(p,x,normalize='none')
 bad <- which(apply(abs(got-expected),1,max)>5e-5)
 if(!length(bad)) bad <- which.max(apply(abs(got-expected),1,max))
 write_json(list(indices=bad, cols=p$cols[p$rows %in% bad], rows=p$rows[p$rows %in% bad], weights=p$vals[p$rows %in% bad], got=got[bad,,drop=FALSE]),file.path(folder,'diagnose.json'),auto_unbox=FALSE,digits=NA)
}
