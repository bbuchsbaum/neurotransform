base <- '/tmp/neurotransform-full-barycentric-measured-20260926'
for(folder in list.dirs(base,recursive=FALSE)) {
  read <- function(name) as.matrix(read.csv(file.path(folder,paste0(name,'.csv')),header=FALSE))
  saveRDS(list(vertices=read('vertices'),faces=read('faces'),query=read('query')),
          file.path('/tmp',paste0('neurotransform-benchmark-',basename(folder),'.rds')))
  gc()
}
