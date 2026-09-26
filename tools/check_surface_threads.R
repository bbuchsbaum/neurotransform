# Run in separate R processes with OMP_NUM_THREADS=1 and 4, then compare RDS values.
library(neurotransform)
args <- commandArgs(TRUE)
fixture <- jsonlite::fromJSON(system.file('extdata/barycentric_oracle/oracle.json',package='neurotransform'))
v <- fixture$cases$vertices[[3]]; f <- fixture$cases$faces[[3]]
set.seed(20926)
q <- rbind(matrix(rnorm(900),ncol=3)*100,v)
weights <- neurotransform:::cpp_barycentric_weights(q,v,f)
exhaustive <- neurotransform:::cpp_barycentric_weights(q,v,f,indexed=FALSE)
fields <- c('rows','cols','vals','distance','face')
stopifnot(identical(weights[fields],exhaustive[fields]))
sampler <- surface_sampler(v,v,faces=f+1L,method='barycentric',projection='closest')
saveRDS(list(weights=weights[fields],samples=sampler@evaluate(q)),args[1])
