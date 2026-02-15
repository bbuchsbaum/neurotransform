#' @title X5 Transform IO
#' @name x5_io
#' @description
#' Read and write X5 transform files for linear and nonlinear transforms.
NULL

#' Create an X5 domain descriptor
#'
#' @param grid Logical
#' @param size Integer vector
#' @param mapping Optional matrix
#' @param coordinates Optional coordinate kind string
#' @return An `X5Domain` object
#' @export
x5_domain <- function(grid, size, mapping = NULL, coordinates = NULL) {
  structure(
    list(
      grid = isTRUE(grid),
      size = as.integer(size),
      mapping = mapping,
      coordinates = coordinates
    ),
    class = "X5Domain"
  )
}

#' Create an X5 transform descriptor
#'
#' @param type `"linear"`, `"nonlinear"`, or `"composite"`
#' @param transform Numeric array of transform parameters
#' @param subtype Optional subtype
#' @param representation Optional representation label
#' @param metadata Optional named list or JSON string
#' @param dimension_kinds Optional character vector
#' @param domain Optional `X5Domain`
#' @param inverse Optional inverse array
#' @param jacobian Optional Jacobian array
#' @param additional_parameters Optional array
#' @param array_length Integer array length
#' @return An `X5Transform` object
#' @export
x5_transform <- function(type,
                         transform,
                         subtype = NULL,
                         representation = NULL,
                         metadata = NULL,
                         dimension_kinds = NULL,
                         domain = NULL,
                         inverse = NULL,
                         jacobian = NULL,
                         additional_parameters = NULL,
                         array_length = 1L) {
  structure(
    list(
      type = as.character(type),
      transform = transform,
      subtype = subtype,
      representation = representation,
      metadata = metadata,
      dimension_kinds = dimension_kinds,
      domain = domain,
      inverse = inverse,
      jacobian = jacobian,
      additional_parameters = additional_parameters,
      array_length = as.integer(array_length)[1]
    ),
    class = "X5Transform"
  )
}

.x5_encode_metadata <- function(metadata) {
  if (is.null(metadata)) return(NULL)
  if (is.character(metadata) && length(metadata) == 1L) return(metadata)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::toJSON(metadata, auto_unbox = TRUE))
  }
  .stop_transform_io("jsonlite is required to encode non-string X5 metadata.")
}

.x5_decode_metadata <- function(value) {
  if (is.null(value) || !length(value)) return(NULL)
  text <- as.character(value[[1]])
  if (!nzchar(text)) return(NULL)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    out <- tryCatch(jsonlite::fromJSON(text), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  text
}

.x5_sorted_keys <- function(g) {
  keys <- names(g)
  ord <- order(suppressWarnings(as.integer(keys)), na.last = TRUE)
  keys[ord]
}

.x5_write_group <- function(g, node) {
  g$create_attr("Type", robj = as.character(node$type))
  g$create_attr("ArrayLength", robj = as.integer(node$array_length %||% 1L))
  if (!is.null(node$subtype)) g$create_attr("SubType", robj = as.character(node$subtype))
  if (!is.null(node$representation)) g$create_attr("Representation", robj = as.character(node$representation))
  if (!is.null(node$metadata)) g$create_attr("Metadata", robj = .x5_encode_metadata(node$metadata))

  g$create_dataset("Transform", robj = node$transform)
  if (!is.null(node$dimension_kinds)) {
    g$create_dataset("DimensionKinds", robj = as.character(node$dimension_kinds))
  }
  if (!is.null(node$domain)) {
    d <- g$create_group("Domain")
    d$create_dataset("Grid", robj = as.integer(if (isTRUE(node$domain$grid)) 1L else 0L))
    d$create_dataset("Size", robj = as.integer(node$domain$size))
    if (!is.null(node$domain$mapping)) d$create_dataset("Mapping", robj = node$domain$mapping)
    if (!is.null(node$domain$coordinates)) d$create_attr("Coordinates", robj = as.character(node$domain$coordinates))
  }
  if (!is.null(node$inverse)) g$create_dataset("Inverse", robj = node$inverse)
  if (!is.null(node$jacobian)) g$create_dataset("Jacobian", robj = node$jacobian)
  if (!is.null(node$additional_parameters)) g$create_dataset("AdditionalParameters", robj = node$additional_parameters)
}

.x5_read_group <- function(g) {
  domain <- NULL
  if ("Domain" %in% names(g)) {
    d <- g[["Domain"]]
    dattrs <- hdf5r::h5attributes(d)
    domain <- x5_domain(
      grid = as.logical(as.integer(d[["Grid"]]$read())[1]),
      size = as.integer(d[["Size"]]$read()),
      mapping = if ("Mapping" %in% names(d)) d[["Mapping"]]$read() else NULL,
      coordinates = if ("Coordinates" %in% names(dattrs)) as.character(dattrs$Coordinates) else NULL
    )
  }

  gattrs <- hdf5r::h5attributes(g)
  x5_transform(
    type = as.character(gattrs$Type),
    transform = g[["Transform"]]$read(),
    subtype = if ("SubType" %in% names(gattrs)) as.character(gattrs$SubType) else NULL,
    representation = if ("Representation" %in% names(gattrs)) as.character(gattrs$Representation) else NULL,
    metadata = if ("Metadata" %in% names(gattrs)) .x5_decode_metadata(gattrs$Metadata) else NULL,
    dimension_kinds = if ("DimensionKinds" %in% names(g)) as.character(g[["DimensionKinds"]]$read()) else NULL,
    domain = domain,
    inverse = if ("Inverse" %in% names(g)) g[["Inverse"]]$read() else NULL,
    jacobian = if ("Jacobian" %in% names(g)) g[["Jacobian"]]$read() else NULL,
    additional_parameters = if ("AdditionalParameters" %in% names(g)) g[["AdditionalParameters"]]$read() else NULL,
    array_length = if ("ArrayLength" %in% names(gattrs)) as.integer(gattrs$ArrayLength)[1] else 1L
  )
}

#' Write X5 transform file
#'
#' @param path Output file path
#' @param x5_list List of `X5Transform` objects
#' @return Invisibly returns `path`
#' @export
write_x5 <- function(path, x5_list) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    .stop_transform_io("hdf5r is required to write X5 files.")
  }
  if (!is.list(x5_list) || !length(x5_list)) {
    .stop_transform_io("x5_list must be a non-empty list of X5Transform objects.")
  }
  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all())
  h5$create_attr("Format", robj = "X5")
  h5$create_attr("Version", robj = as.integer(1L))
  tg <- h5$create_group("TransformGroup")
  for (i in seq_along(x5_list)) {
    g <- tg$create_group(as.character(i - 1L))
    .x5_write_group(g, x5_list[[i]])
  }
  invisible(path)
}

#' Read X5 transform file
#'
#' @param path Input file path
#' @return List of `X5Transform` objects
#' @export
read_x5 <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    .stop_transform_io("hdf5r is required to read X5 files.")
  }
  if (!file.exists(path)) .stop_transform_file(sprintf("X5 file not found: %s", path))
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all())

  attrs <- hdf5r::h5attributes(h5)
  if (!("Format" %in% names(attrs)) || as.character(attrs$Format) != "X5") {
    .stop_transform_file("Input file is not in X5 format.")
  }
  tg <- h5[["TransformGroup"]]
  if (is.null(tg)) .stop_transform_file("X5 file missing TransformGroup.")
  lapply(.x5_sorted_keys(tg), function(k) .x5_read_group(tg[[k]]))
}

.x5_nodes_from_transform <- function(x) {
  if (inherits(x, "MorphismPath")) {
    return(unlist(lapply(x@morphisms, .x5_nodes_from_transform), recursive = FALSE))
  }
  if (inherits(x, "Affine3DMorphism")) {
    return(list(x5_transform(
      type = "linear",
      transform = x@matrix,
      dimension_kinds = c("space", "space")
    )))
  }
  if (inherits(x, "Warp3DMorphism")) {
    w <- load_warp_array(x)
    field <- .unflatten_warp_components(w$array, w$dim)
    repr <- if (identical(x@params$def_type %||% "relative", "absolute")) "deformations" else "displacements"
    dom <- x5_domain(
      grid = TRUE,
      size = as.integer(w$dim),
      mapping = w$vox_to_world,
      coordinates = "cartesian"
    )
    return(list(x5_transform(
      type = "nonlinear",
      subtype = "densefield",
      representation = repr,
      transform = field,
      dimension_kinds = c("space", "space", "space", "vector"),
      domain = dom
    )))
  }
  if (is_affine_matrix(x)) {
    return(list(x5_transform(
      type = "linear",
      transform = x,
      dimension_kinds = c("space", "space")
    )))
  }
  .stop_transform_io("Unsupported transform type for X5 conversion.")
}

.transform_from_x5_node <- function(node, source, target) {
  if (identical(node$type, "linear")) {
    mat <- as.matrix(node$transform)
    validate_4x4_matrix(mat, "x5_linear_transform")
    return(Affine3DMorphism(source = source, target = target, matrix = mat))
  }
  if (identical(node$type, "nonlinear")) {
    if (is.null(node$domain) || !isTRUE(node$domain$grid)) {
      .stop_transform_file("X5 nonlinear transform requires a grid Domain.")
    }
    dims <- as.integer(node$domain$size)
    field <- node$transform
    if (length(dim(field)) == 2L && ncol(field) == 3L && nrow(field) == prod(dims)) {
      field <- array(field, dim = c(dims, 3L))
    }
    if (length(dim(field)) != 4L || dim(field)[4] != 3L) {
      .stop_transform_file("Unsupported X5 nonlinear transform shape; expected c(X,Y,Z,3).")
    }
    repn <- tolower(as.character(node$representation %||% "displacements"))
    representation <- if (identical(repn, "deformations")) "deformations" else "displacements"
    grid <- grid_spec(dims = dims, affine = as.matrix(node$domain$mapping))
    return(warp_from_field(source, target, field = field, grid = grid, representation = representation))
  }
  .stop_transform_file(sprintf("Unsupported X5 transform type '%s'.", node$type))
}

.transform_from_x5 <- function(nodes, source = "source", target = "target") {
  if (!length(nodes)) .stop_transform_file("X5 file contains no transforms.")
  if (length(nodes) == 1L) return(.transform_from_x5_node(nodes[[1]], source, target))

  morphs <- vector("list", length(nodes))
  curr_source <- source
  for (i in seq_along(nodes)) {
    curr_target <- if (i == length(nodes)) target else sprintf(".x5_stage_%d", i)
    morphs[[i]] <- .transform_from_x5_node(nodes[[i]], curr_source, curr_target)
    curr_source <- curr_target
  }
  methods::new("MorphismPath", morphisms = morphs, source = source, target = target)
}
