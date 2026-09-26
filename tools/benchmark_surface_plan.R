# Isolate engine construction/application from CSV parsing and oracle comparisons.
library(neurotransform)
args <- commandArgs(TRUE); input <- readRDS(args[1]); output <- args[2]
source <- surface_mesh(input$vertices,input$faces+1L)
reference <- surface_mesh(input$query)
set.seed(20926)
p <- input$vertices/sqrt(rowSums(input$vertices^2))
data <- cbind(p,p[,1]*p[,2],sin(seq_len(nrow(p))*.017),cos(seq_len(nrow(p))*.031),
              runif(nrow(p),-1,1),matrix(0,nrow(p),4))
data[cbind(sample.int(nrow(p),4),8:11)] <- 1
gc()
baseline <- as.numeric(system2('ps',c('-o','rss=','-p',Sys.getpid()),stdout=TRUE))*1024
construction <- system.time(plan <- surface_resampling_plan(reference,source))['elapsed']
application <- system.time(values <- apply_surface_resampling(plan,data,normalize='none'))['elapsed']
scaled <- mesh_set_radius(source,100)
queries <- mesh_set_radius(reference,100)@coords
sampler_construction <- system.time(sampler <- surface_sampler(scaled,data,method='barycentric',projection='closest'))['elapsed']
sampling <- system.time(sampled <- sampler@evaluate(queries))['elapsed']
error <- max(abs(sampled-values))
receipt <- list(vertices=nrow(p),targets=nrow(input$query),columns=ncol(data),
  construction_seconds=unname(construction),application_seconds=unname(application),
  native_timing=plan$timing,sampler_construction_seconds=unname(sampler_construction),
  cached_sampling_seconds=unname(sampling),sampler_plan_max_error=error,
  baseline_resident_bytes=baseline,passed=is.finite(error) && error<1e-12 && construction<60,
  dll_sha256=digest::digest(file=getLoadedDLLs()[['neurotransform']][['path']],algo='sha256'),
  session=capture.output(sessionInfo()))
jsonlite::write_json(receipt,output,auto_unbox=TRUE,pretty=TRUE,digits=NA)
quit(status=if(receipt$passed) 0L else 1L)
