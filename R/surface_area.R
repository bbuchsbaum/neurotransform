#' Bind vertex areas to a registered surface geometry
#'
#' Supply either positive area values with explicit provenance, or an anatomical
#' mesh with the same ordered faces and vertex identities as the registration
#' geometry. Mesh areas allocate one third of each Euclidean triangle's area to
#' each corner. Registration spheres are never used implicitly as anatomy.
#'
#' @param geometry Registration `SurfaceMesh` to which the areas correspond.
#' @param anatomical Optional anatomical `SurfaceMesh` with matching topology.
#' @param values Optional finite positive numeric area vector.
#' @param units Area units, for example `"mm^2"`; both sides must use the same units.
#' @param provenance Nonempty source description, required for supplied values.
#' @return A `SurfaceVertexAreas` object with values, geometry identity, units and provenance.
#' @export
surface_vertex_areas <- function(geometry, anatomical = NULL, values = NULL,
                                 units = "mm^2", provenance = NULL) {
  geometry <- if (inherits(geometry, "SurfaceMesh")) geometry else surface_mesh(geometry)
  if (is.null(anatomical) == is.null(values)) stop("supply exactly one of anatomical or values")
  if (!is.character(units) || length(units) != 1L || is.na(units) || !nzchar(trimws(units)))
    stop("units must be one nonempty string")
  anatomy_hash <- NULL
  if (!is.null(anatomical)) {
    anatomy <- if (inherits(anatomical, "SurfaceMesh")) anatomical else surface_mesh(anatomical)
    if (nrow(anatomy@coords) != nrow(geometry@coords) ||
        !identical(anatomy@faces, geometry@faces) || !nrow(anatomy@faces))
      stop("anatomical mesh must have the same ordered faces and vertex count as geometry")
    validate_surface_mesh(anatomy, spherical = FALSE, error = TRUE)
    f <- anatomy@faces + 1L
    a <- anatomy@coords[f[,2],,drop=FALSE] - anatomy@coords[f[,1],,drop=FALSE]
    b <- anatomy@coords[f[,3],,drop=FALSE] - anatomy@coords[f[,1],,drop=FALSE]
    scale <- pmax(abs(a[,1]),abs(a[,2]),abs(a[,3]),abs(b[,1]),abs(b[,2]),abs(b[,3]))
    a <- a/scale; b <- b/scale
    cross <- cbind(a[,2]*b[,3]-a[,3]*b[,2],a[,3]*b[,1]-a[,1]*b[,3],a[,1]*b[,2]-a[,2]*b[,1])
    area <- .5 * scale^2 * sqrt(rowSums(cross^2))
    values <- .surface_rowsum(rep(area/3,3),as.vector(f),nrow(anatomy@coords))
    anatomy_hash <- .surface_identity(anatomy)
    if (is.null(provenance)) provenance <- paste("Euclidean triangle areas from anatomical mesh",anatomy_hash)
  }
  if (!is.numeric(values) || !is.null(dim(values)) || length(values) != nrow(geometry@coords) ||
      any(!is.finite(values) | values <= 0)) stop("areas must be finite, positive and have one value per vertex")
  if (!is.character(provenance) || length(provenance) != 1L || is.na(provenance) || !nzchar(trimws(provenance)))
    stop("provenance must be one nonempty source description")
  structure(list(values=as.numeric(values), geometry=.surface_identity(geometry),
    units=units, provenance=provenance, anatomical_identity=anatomy_hash),class="SurfaceVertexAreas")
}

.surface_check_areas <- function(areas, geometry, n, name) {
  if (!inherits(areas,"SurfaceVertexAreas") || !identical(areas$geometry,geometry))
    stop(name," must be SurfaceVertexAreas bound to the matching ordered geometry")
  if (!is.numeric(areas$values) || length(areas$values)!=n || any(!is.finite(areas$values) | areas$values<=0))
    stop(name," must contain finite positive areas for every vertex")
  if (!is.character(areas$provenance) || length(areas$provenance)!=1L ||
      is.na(areas$provenance) || !nzchar(trimws(areas$provenance))) stop(name," lacks provenance")
  if (!is.character(areas$units) || length(areas$units)!=1L ||
      is.na(areas$units) || !nzchar(trimws(areas$units))) stop(name," lacks units")
  areas
}

# Positive sparse support is significant: an extra reverse contributor changes
# the selected row wholesale. Do not threshold small weights for oracle parity.
.surface_adaptive_weights <- function(forward, reverse, n_reference, n_moving,
                                      source_area, target_area) {
  f <- forward$vals>0; r <- reverse$vals>0
  fr <- forward$rows[f]; fc <- forward$cols[f]; fv <- forward$vals[f]
  gr <- reverse$cols[r]; gc <- reverse$rows[r]; gv <- reverse$vals[r]
  fkey <- as.double(fr) + as.double(n_reference)*(fc-1)
  gkey <- as.double(gr) + as.double(n_reference)*(gc-1)
  use_reverse <- rep(FALSE,n_reference)
  use_reverse[gr[!gkey %in% fkey]] <- TRUE
  keep_f <- !use_reverse[fr]; keep_g <- use_reverse[gr]
  rows <- c(fr[keep_f],gr[keep_g]); cols <- c(fc[keep_f],gc[keep_g])
  selected <- c(fv[keep_f],gv[keep_g])
  # A common target-area scale cancels in column correction. Normalizing it
  # first avoids overflow in column sums for otherwise valid area units.
  weighted <- selected*(target_area/max(target_area))[rows]
  correction <- .surface_rowsum(weighted,cols,n_moving)
  corrected <- (weighted/correction[cols])*source_area[cols]
  measure <- .surface_rowsum(corrected,rows,n_reference)
  vals <- corrected/measure[rows]
  if (any(!is.finite(vals) | vals<=0) || any(!is.finite(measure)))
    stop("adaptive area arithmetic produced nonfinite or nonpositive weights")
  ord <- order(rows,cols)
  list(rows=as.integer(rows[ord]),cols=as.integer(cols[ord]),vals=vals[ord],
       target_measure=measure,represented_source=correction>0,use_reverse=use_reverse)
}
