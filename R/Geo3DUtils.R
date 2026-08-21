####################################################
### Internal coordinate-dimension validation
####################################################

.GeoCoordinateDimension <- function(coordx = NULL, coordy = NULL,
                                    coordz = NULL, coordx_dyn = NULL) {
  if (!is.null(coordx_dyn)) {
    if (!is.list(coordx_dyn) || !length(coordx_dyn)) {
      stop("coordx_dyn must be a non-empty list of coordinate matrices.",
           call. = FALSE)
    }
    dims <- vapply(seq_along(coordx_dyn), function(i) {
      z <- as.matrix(coordx_dyn[[i]])
      if (!(ncol(z) %in% c(2L, 3L))) {
        stop("Each element of coordx_dyn must have two or three columns.",
             call. = FALSE)
      }
      ncol(z)
    }, integer(1L))
    if (length(unique(dims)) != 1L) {
      stop("All elements of coordx_dyn must have the same number of columns.",
           call. = FALSE)
    }
    return(unname(dims[[1L]]))
  }

  if ((is.matrix(coordx) || is.data.frame(coordx)) &&
      is.null(coordy) && is.null(coordz)) {
    return(ncol(as.matrix(coordx)))
  }
  if (!is.null(coordz)) return(3L)
  if (!is.null(coordy)) return(2L)
  if (is.matrix(coordx) || is.data.frame(coordx)) {
    return(ncol(as.matrix(coordx)))
  }
  NA_integer_
}

.GeoNormalizeDynamicCoordinates <- function(
    coordx_dyn, coordt = NULL, context = "This function") {
  if (!is.list(coordx_dyn) || !length(coordx_dyn)) {
    stop(context, ": coordx_dyn must be a non-empty list of coordinate matrices.",
         call. = FALSE)
  }

  blocks <- lapply(seq_along(coordx_dyn), function(i) {
    z <- as.matrix(coordx_dyn[[i]])
    if (!is.numeric(z)) {
      stop(context, ": coordx_dyn[[", i, "]] must be numeric.",
           call. = FALSE)
    }
    if (length(dim(z)) != 2L || !(ncol(z) %in% c(2L, 3L))) {
      stop(context, ": each element of coordx_dyn must have two or three columns.",
           call. = FALSE)
    }
    if (nrow(z) < 1L) {
      stop(context, ": each element of coordx_dyn must contain at least one location.",
           call. = FALSE)
    }
    if (any(!is.finite(z))) {
      stop(context, ": coordx_dyn contains non-finite coordinates.",
           call. = FALSE)
    }
    storage.mode(z) <- "double"
    unname(z)
  })

  dims <- vapply(blocks, ncol, integer(1L))
  if (length(unique(dims)) != 1L) {
    stop(context, ": all elements of coordx_dyn must have the same number of columns.",
         call. = FALSE)
  }

  if (!is.null(coordt)) {
    if (!is.numeric(coordt) || any(!is.finite(coordt))) {
      stop(context, ": coordt must be a finite numeric vector.", call. = FALSE)
    }
    if (length(coordt) != length(blocks)) {
      stop(context, ": coordx_dyn must contain one coordinate matrix for each element of coordt.",
           call. = FALSE)
    }
  }

  blocks
}

.GeoValidateDynamicData <- function(data, coordx_dyn,
                                    context = "This function",
                                    require_list = FALSE) {
  sizes <- vapply(coordx_dyn, nrow, integer(1L))
  total <- sum(sizes)

  if (isTRUE(require_list) && !is.list(data)) {
    stop(context, ": data must be a list with one component per temporal instant when coordx_dyn is used.",
         call. = FALSE)
  }

  if (is.list(data)) {
    if (length(data) != length(coordx_dyn)) {
      stop(context, ": data and coordx_dyn must have the same number of temporal components.",
           call. = FALSE)
    }
    for (i in seq_along(data)) {
      if (!is.numeric(data[[i]]) || any(!is.finite(data[[i]]))) {
        stop(context, ": data[[", i, "]] must be finite and numeric.",
             call. = FALSE)
      }
      if (length(data[[i]]) != sizes[i]) {
        stop(context, ": data[[", i, "]] must contain one value per row of coordx_dyn[[", i, "]].",
             call. = FALSE)
      }
    }
  } else {
    if (!is.numeric(data) || any(!is.finite(data))) {
      stop(context, ": data must be finite and numeric.", call. = FALSE)
    }
    if (length(data) != total) {
      stop(context, ": data must contain one value per dynamic observation (",
           total, " values expected).", call. = FALSE)
    }
  }

  invisible(sizes)
}

.GeoValidateDynamicRows <- function(x, coordx_dyn, name = "X",
                                    context = "This function") {
  if (is.null(x)) return(invisible(NULL))

  sizes <- vapply(coordx_dyn, nrow, integer(1L))
  total <- sum(sizes)

  if (is.list(x)) {
    if (length(x) != length(coordx_dyn)) {
      stop(context, ": ", name,
           " must have one component per temporal instant.", call. = FALSE)
    }
    for (i in seq_along(x)) {
      xi <- as.matrix(x[[i]])
      if (!is.numeric(xi) || any(!is.finite(xi))) {
        stop(context, ": ", name, "[[", i, "]] must be finite and numeric.",
             call. = FALSE)
      }
      if (nrow(xi) != sizes[i]) {
        stop(context, ": ", name, "[[", i,
             "]] must have one row per row of coordx_dyn[[", i, "]].",
             call. = FALSE)
      }
    }
  } else {
    xi <- as.matrix(x)
    if (!is.numeric(xi) || any(!is.finite(xi))) {
      stop(context, ": ", name, " must be finite and numeric.",
           call. = FALSE)
    }
    if (nrow(xi) != total) {
      stop(context, ": ", name,
           " must have one row per dynamic observation (", total,
           " rows expected).", call. = FALSE)
    }
  }

  invisible(sizes)
}

.GeoValidateCoordinateCompatibility <- function(
    coordx = NULL, coordy = NULL, coordz = NULL, coordx_dyn = NULL,
    loc, context = "Kriging") {
  p_obs <- .GeoCoordinateDimension(coordx, coordy, coordz, coordx_dyn)
  loc <- as.matrix(loc)
  p_loc <- ncol(loc)

  if (!(p_loc %in% c(2L, 3L))) {
    stop(context, ": prediction locations must have two or three columns.",
         call. = FALSE)
  }
  if (is.na(p_obs)) return(invisible(c(observed = p_obs, prediction = p_loc)))
  if (!(p_obs %in% c(2L, 3L))) {
    stop(context, ": observed coordinates must have two or three columns.",
         call. = FALSE)
  }
  if (p_obs != p_loc) {
    stop(context, ": observed coordinates are ", p_obs,
         "-dimensional, but prediction locations are ", p_loc,
         "-dimensional. Supply coordinates with the same spatial dimension.",
         call. = FALSE)
  }

  invisible(c(observed = p_obs, prediction = p_loc))
}

.GeoValidate3D <- function(coordx = NULL, coordy = NULL, coordz = NULL,
                           coordx_dyn = NULL, distance = "Eucl",
                           grid = FALSE, approximate = FALSE,
                           method = NULL, context = "This function") {
  p <- .GeoCoordinateDimension(coordx, coordy, coordz, coordx_dyn)
  if (is.na(p)) return(invisible(p))
  if (!(p %in% c(2L, 3L))) {
    stop(context, " requires coordinates with two or three columns.",
         call. = FALSE)
  }

  distance_key <- tolower(gsub("[[:blank:]]", "", as.character(distance)[1L]))
  if (p == 3L && !identical(distance_key, "eucl")) {
    stop(context, ": three-dimensional coordinates require distance = 'Eucl'.",
         call. = FALSE)
  }
  if (p == 3L && isTRUE(grid)) {
    stop(context, ": regular three-dimensional grids are not implemented; ",
         "supply their Cartesian product as an explicit N x 3 matrix with ",
         "grid = FALSE.", call. = FALSE)
  }
  if (p == 3L && isTRUE(approximate)) {
    method_key <- if (is.null(method)) "" else
      toupper(gsub("[[:blank:]]", "", as.character(method)[1L]))
    if (!identical(method_key, "TB")) {
      if (identical(method_key, "CE")) {
        stop(context, ": CE approximate simulation is implemented only for ",
             "two-dimensional regular grids. Use method = 'TB' for a purely ",
             "spatial univariate irregular three-dimensional field, or use ",
             "method = 'cholesky'.", call. = FALSE)
      }
      stop(context, ": approximate three-dimensional simulation is supported ",
           "only by method = 'TB' for purely spatial univariate fields. Use ",
           "method = 'cholesky' otherwise.", call. = FALSE)
    }
  }
  invisible(p)
}
