args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(args[[1]], mustWork = TRUE)
pkgload::load_all(quiet = TRUE)
manifest <- jsonlite::fromJSON(file.path(out, "manifest.json"), simplifyVector = FALSE)

write_array <- function(x, path, affine) {
  im <- RNifti::asNifti(x, reference = list(xyzt_units = 2L))
  im <- RNifti::`pixdim<-`(im, c(sqrt(colSums(affine[1:3, 1:3]^2)), rep(1, max(0, length(dim(x))-3))))
  im <- RNifti::`qform<-`(im, structure(affine, code = 1L))
  im <- RNifti::`sform<-`(im, structure(affine, code = 1L))
  RNifti::writeNifti(im, path, datatype = "float")
}

results <- list()
for (case in manifest$cases) {
  if (!is.null(case$missing)) next
  cat("Candidate:", case$id, "\n")
  result <- tryCatch({
    folder <- file.path(out, case$id)
    src <- neuroim2::read_vol(file.path(folder, "source.nii.gz"))
    target <- neuroim2::read_vol(file.path(folder, "target.nii.gz"))
    sa <- neuroim2::trans(src)
    ta <- neuroim2::trans(target)
    dims <- as.integer(dim(target)[1:3])
    path <- file.path(folder, case$transform)
    if (case$kind == "affine") {
      m <- read_linear_transform(path, format = switch(case$family, ANTs = "itk", AFNI = "afni", FSL = "fsl"),
                                  source_affine = sa, target_affine = ta,
                                  source_dim = dim(src)[1:3], target_dim = dims)
    } else if (case$kind == "h5") {
      m <- read_transform(path, type = "ants_h5")
    } else if (case$kind == "composite") {
      # AFNI -nwarp evaluates its list left-to-right; MorphismPath stores
      # source-to-target components and evaluates right-to-left. Thus the
      # native string "affine warp" becomes the path list(warp, affine).
      # https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dNwarpApply.html
      a <- read_linear_transform(file.path(folder, "affine.aff12.1D"), format = "afni",
                                 source = "middle", target = "target")
      w <- read_transform(file.path(folder, "warp.nii.gz"), type = "afni",
                          source = "source", target = "middle")
      m <- methods::new("MorphismPath", morphisms = list(w, a), source = "source", target = "target")
    } else {
      m <- read_transform(path, type = tolower(case$family),
                          source_affine = sa, source_dim = dim(src)[1:3],
                          target_affine = ta, target_dim = dims,
                          def_type = if (identical(case$representation, "absolute")) "absolute" else "relative")
    }
    ijk <- as.matrix(expand.grid(lapply(dims, function(n) 0:(n-1))))
    points <- (cbind(ijk, 1) %*% t(ta))[, 1:3, drop = FALSE]
    mapped <- transform(m, points)
    write_array(array(mapped, c(dims, 3)), file.path(folder, "ours_coords.nii.gz"), ta)
    # Exercise the public resampling path separately from coordinate evaluation.
    img <- resample_to(src, target, m, method = "linear")
    write_array(as.array(img), file.path(folder, "ours_source.nii.gz"), ta)
    sampler <- volume_sampler(src, method = "linear", outside = 0)
    bad_sign <- points + sweep(mapped - points, 2, c(-1, -1, 1), "*")
    variants <- list(list(id = "wrong_sign", label = "Negate X/Y displacement", coords = bad_sign))
    if (is(m, "Affine3DMorphism")) {
      variants[[length(variants)+1L]] <- list(id = "wrong_direction", label = "Invert affine direction",
                                             coords = transform(invert(m), points))
    } else if (is(m, "MorphismPath") && length(m@morphisms) == 2L) {
      parts <- rev(m@morphisms)
      ids <- c("source", "middle", "target")
      for (i in 1:2) {
        parts[[i]]@source <- ids[[i]]
        parts[[i]]@target <- ids[[i+1]]
      }
      wrong <- methods::new("MorphismPath", morphisms = parts, source = ids[[1]], target = ids[[3]])
      variants[[length(variants)+1L]] <- list(id = "wrong_order", label = "Swap affine/warp order",
                                             coords = transform(wrong, points))
    }
    for (v in variants) {
      write_array(array(sampler@evaluate(v$coords), dims), file.path(folder, paste0(v$id, ".nii.gz")), ta)
      write_array(array(v$coords, c(dims, 3)), file.path(folder, paste0(v$id, "_coords.nii.gz")), ta)
    }
    # Jacobians are diagnostics; invalid/outside samples remain visible as NA.
    jac_error <- NULL
    jac <- tryCatch(jacobian_det(m, points), error = function(e) {
      jac_error <<- conditionMessage(e)
      rep(NA_real_, nrow(points))
    })
    write_array(array(jac, dims), file.path(folder, "ours_jacobian.nii.gz"), ta)
    list(id = case$id, status = "evaluated", jacobian_error = jac_error,
         variants = lapply(variants, function(v) v[c("id", "label")]))
  }, error = function(e) list(id = case$id, status = "error", message = conditionMessage(e)))
  results[[length(results)+1L]] <- result
}
jsonlite::write_json(list(cases = results, session = capture.output(sessionInfo())),
                      file.path(out, "candidate-results.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
