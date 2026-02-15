#' @title FreeSurfer LTA IO
#' @name lta_io
#' @description
#' Read/write FreeSurfer LTA linear transform files (single or array forms).
NULL

.lta_drop_comments <- function(lines) {
  out <- trimws(sub("#.*$", "", lines))
  out[nzchar(out)]
}

.lta_parse_kv <- function(line) {
  m <- regexec("^([[:alnum:]_]+)\\s*=\\s*(.*)$", line)
  r <- regmatches(line, m)[[1]]
  if (length(r) < 3L) return(NULL)
  list(key = r[2], value = r[3])
}

.lta_default_volgeom <- function() {
  list(
    valid = 0L,
    filename = "unknown.nii.gz",
    volume = c(1, 1, 1),
    voxelsize = c(1, 1, 1),
    xras = c(1, 0, 0),
    yras = c(0, 1, 0),
    zras = c(0, 0, 1),
    cras = c(0, 0, 0)
  )
}

.lta_parse_volgeom <- function(lines, i) {
  vg <- .lta_default_volgeom()
  keys <- c("valid", "filename", "volume", "voxelsize", "xras", "yras", "zras", "cras")
  for (k in keys) {
    if (i > length(lines)) .stop_transform_file("Unexpected end of file while reading volume geometry.")
    kv <- .lta_parse_kv(lines[[i]])
    if (is.null(kv) || kv$key != k) .stop_transform_file("Malformed LTA volume geometry block.")
    if (identical(k, "filename")) {
      vg$filename <- kv$value
    } else {
      vals <- scan(text = kv$value, quiet = TRUE)
      if (!length(vals)) .stop_transform_file(sprintf("Malformed volume geometry key '%s'.", k))
      vg[[k]] <- vals
    }
    i <- i + 1L
  }
  vg$valid <- as.integer(vg$valid[1])
  vg$volume <- as.integer(vg$volume[1:3])
  vg$voxelsize <- as.numeric(vg$voxelsize[1:3])
  vg$xras <- as.numeric(vg$xras[1:3])
  vg$yras <- as.numeric(vg$yras[1:3])
  vg$zras <- as.numeric(vg$zras[1:3])
  vg$cras <- as.numeric(vg$cras[1:3])
  list(vg = vg, i = i)
}

.lta_volgeom_affine <- function(vg) {
  A <- cbind(vg$xras, vg$yras, vg$zras)
  A <- A %*% diag(vg$voxelsize)
  b <- vg$cras - as.numeric(A %*% (vg$volume / 2))
  M <- diag(4)
  M[1:3, 1:3] <- A
  M[1:3, 4] <- b
  M
}

.lta_matrix_to_internal <- function(m_L, type, src_vg = NULL, dst_vg = NULL) {
  M <- m_L
  if (identical(type, 0L)) {
    if (!is.null(src_vg) && !is.null(dst_vg)) {
      M <- .lta_volgeom_affine(dst_vg) %*% m_L %*% solve(.lta_volgeom_affine(src_vg))
    }
  } else if (!identical(type, 1L)) {
    warning(sprintf("LTA type %s is not fully supported; interpreting matrix as RAS-to-RAS.", type))
  }
  solve(M)
}

.lta_internal_to_matrix <- function(mat_internal) {
  solve(mat_internal)
}

.read_lta <- function(path) {
  if (!file.exists(path)) .stop_transform_file(sprintf("LTA file not found: %s", path))
  lines <- .lta_drop_comments(readLines(path, warn = FALSE))
  if (!length(lines)) .stop_transform_file("Empty LTA file.")

  i <- 1L
  type <- NA_integer_
  nxforms <- NA_integer_

  while (i <= length(lines)) {
    kv <- .lta_parse_kv(lines[[i]])
    if (!is.null(kv) && kv$key == "type") {
      type <- as.integer(scan(text = kv$value, quiet = TRUE)[1])
      i <- i + 1L
      break
    }
    i <- i + 1L
  }
  if (is.na(type)) .stop_transform_file("LTA missing type header.")

  while (i <= length(lines)) {
    kv <- .lta_parse_kv(lines[[i]])
    if (!is.null(kv) && kv$key == "nxforms") {
      nxforms <- as.integer(scan(text = kv$value, quiet = TRUE)[1])
      i <- i + 1L
      break
    }
    i <- i + 1L
  }
  if (is.na(nxforms) || nxforms < 1L) .stop_transform_file("LTA missing nxforms header.")

  xforms <- vector("list", nxforms)
  for (k in seq_len(nxforms)) {
    while (i <= length(lines) && !startsWith(lines[[i]], "mean")) i <- i + 1L
    if (i > length(lines)) .stop_transform_file("Malformed LTA transform block (missing mean).")
    i <- i + 1L

    if (i > length(lines) || !startsWith(lines[[i]], "sigma")) {
      .stop_transform_file("Malformed LTA transform block (missing sigma).")
    }
    i <- i + 1L

    if (i > length(lines) || !grepl("^1\\s+4\\s+4$", lines[[i]])) {
      .stop_transform_file("Malformed LTA transform block (missing matrix header '1 4 4').")
    }
    i <- i + 1L
    if (i + 3L > length(lines)) .stop_transform_file("Malformed LTA transform block (truncated matrix).")
    mvals <- scan(text = paste(lines[i:(i + 3L)], collapse = "\n"), quiet = TRUE)
    if (length(mvals) < 16L) .stop_transform_file("Malformed LTA matrix block.")
    m_L <- matrix(mvals[1:16], nrow = 4, byrow = TRUE)
    i <- i + 4L

    src_vg <- dst_vg <- NULL
    if (i <= length(lines) && startsWith(lines[[i]], "src volume info")) {
      i <- i + 1L
      parsed <- .lta_parse_volgeom(lines, i)
      src_vg <- parsed$vg
      i <- parsed$i
    }
    if (i <= length(lines) && startsWith(lines[[i]], "dst volume info")) {
      i <- i + 1L
      parsed <- .lta_parse_volgeom(lines, i)
      dst_vg <- parsed$vg
      i <- parsed$i
    }

    xforms[[k]] <- list(m_L = m_L, src = src_vg, dst = dst_vg)
  }

  list(type = type, nxforms = nxforms, xforms = xforms)
}

.format_volgeom_line <- function(key, value) {
  if (identical(key, "filename")) {
    return(sprintf("filename = %s", value))
  }
  if (identical(key, "valid")) {
    return(sprintf("valid = %d  # volume info %svalid", as.integer(value[1]), if (as.integer(value[1]) == 1L) "" else "in"))
  }
  if (identical(key, "volume")) {
    return(sprintf("volume = %d %d %d", as.integer(value[1]), as.integer(value[2]), as.integer(value[3])))
  }
  if (identical(key, "voxelsize")) {
    return(sprintf("voxelsize = %.15e %.15e %.15e", value[1], value[2], value[3]))
  }
  if (identical(key, "xras")) {
    return(sprintf("xras   = %.15e %.15e %.15e", value[1], value[2], value[3]))
  }
  if (identical(key, "yras")) {
    return(sprintf("yras   = %.15e %.15e %.15e", value[1], value[2], value[3]))
  }
  if (identical(key, "zras")) {
    return(sprintf("zras   = %.15e %.15e %.15e", value[1], value[2], value[3]))
  }
  sprintf("cras   = %.15e %.15e %.15e", value[1], value[2], value[3])
}

.write_lta <- function(path, mats_internal, src_vg = NULL, dst_vg = NULL, type = 1L) {
  if (!is.list(mats_internal)) mats_internal <- list(mats_internal)
  if (!length(mats_internal)) .stop_transform_io("No matrices provided for LTA writing.")
  for (m in mats_internal) validate_4x4_matrix(m, "lta_matrix")

  src_vg <- src_vg %||% .lta_default_volgeom()
  dst_vg <- dst_vg %||% .lta_default_volgeom()

  lines <- c(
    "# LTA-array file created by neurotransform",
    sprintf("type      = %d", as.integer(type)),
    sprintf("nxforms   = %d", length(mats_internal))
  )

  for (mat in mats_internal) {
    m_L <- .lta_internal_to_matrix(mat)
    lines <- c(
      lines,
      "mean      = 0.0000 0.0000 0.0000",
      "sigma     = 1.0000",
      "1 4 4",
      apply(m_L, 1, function(r) paste(sprintf("%18.15e", r), collapse = " ")),
      "src volume info",
      .format_volgeom_line("valid", src_vg$valid),
      .format_volgeom_line("filename", src_vg$filename),
      .format_volgeom_line("volume", src_vg$volume),
      .format_volgeom_line("voxelsize", src_vg$voxelsize),
      .format_volgeom_line("xras", src_vg$xras),
      .format_volgeom_line("yras", src_vg$yras),
      .format_volgeom_line("zras", src_vg$zras),
      .format_volgeom_line("cras", src_vg$cras),
      "dst volume info",
      .format_volgeom_line("valid", dst_vg$valid),
      .format_volgeom_line("filename", dst_vg$filename),
      .format_volgeom_line("volume", dst_vg$volume),
      .format_volgeom_line("voxelsize", dst_vg$voxelsize),
      .format_volgeom_line("xras", dst_vg$xras),
      .format_volgeom_line("yras", dst_vg$yras),
      .format_volgeom_line("zras", dst_vg$zras),
      .format_volgeom_line("cras", dst_vg$cras)
    )
  }

  lines <- c(lines, "subject unknown", "fscale 0.000000", "")
  writeLines(lines, con = path)
  invisible(path)
}
