## Fast exact duplicate-location validation.
## The native checker uses an open-addressing hash table and stops at the
## first duplicate, so the expected cost is O(n) and no distance matrix is
## created. Internal resampling/refitting paths can suppress a repeated check
## after the original design has already been validated.

.GeoDuplicateCheckOption <- "GeoModels.skip_duplicate_coordinate_check"

.GeoWithDuplicateCheckSuppressed <- function(expr) {
  old <- getOption(.GeoDuplicateCheckOption, FALSE)
  options(structure(list(TRUE), names = .GeoDuplicateCheckOption))
  on.exit(options(structure(list(old), names = .GeoDuplicateCheckOption)), add = TRUE)
  force(expr)
}

.GeoDuplicateCheckEnabled <- function() {
  !isTRUE(getOption(.GeoDuplicateCheckOption, FALSE))
}

.GeoSuppressDuplicateChecksForNestedCalls <- function() {
  old <- getOption(.GeoDuplicateCheckOption, FALSE)
  options(structure(list(TRUE), names = .GeoDuplicateCheckOption))
  old
}

.GeoRestoreDuplicateCheckOption <- function(old) {
  options(structure(list(old), names = .GeoDuplicateCheckOption))
  invisible(NULL)
}

.GeoCoordMatrixForDuplicateCheck <- function(coordx, coordy = NULL, coordz = NULL,
                                              context = "GeoModels") {
  if (!is.null(coordy)) {
    z <- if (is.null(coordz)) cbind(coordx, coordy) else cbind(coordx, coordy, coordz)
  } else {
    z <- as.matrix(coordx)
  }
  z <- as.matrix(z)
  if (!is.numeric(z) || !nrow(z) || !(ncol(z) %in% 1:3) || any(!is.finite(z))) {
    stop(context, ": coordinates must be a finite numeric matrix.", call. = FALSE)
  }
  storage.mode(z) <- "double"
  unname(z)
}

.GeoCheckDuplicateMatrix <- function(z, context, label = "observation locations") {
  hit <- .Call("GeoDuplicateCoordinates", z, PACKAGE = "GeoModels")
  if (hit[2L] > 0) {
    stop(
      context, ": duplicated ", label, " are not allowed. First duplicate: observations ",
      format(hit[1L], scientific = FALSE, trim = TRUE), " and ",
      format(hit[2L], scientific = FALSE, trim = TRUE), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.GeoCheckDuplicateLocations <- function(coordx = NULL, coordy = NULL, coordz = NULL,
                                         coordt = NULL, coordx_dyn = NULL,
                                         spacetime = FALSE, bivariate = FALSE,
                                         grid = FALSE, context = "GeoModels",
                                         check.duplicates = FALSE) {
  if (!is.logical(check.duplicates) || length(check.duplicates) != 1L || is.na(check.duplicates)) {
    stop(context, ": check.duplicates must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(check.duplicates) || !.GeoDuplicateCheckEnabled()) return(invisible(TRUE))

  check_time <- function(tt) {
    if (is.null(tt)) return(invisible(TRUE))
    if (!is.numeric(tt) || any(!is.finite(tt)))
      stop(context, ": coordt must be finite and numeric.", call. = FALSE)
    jj <- anyDuplicated(as.numeric(tt))
    if (jj) {
      ii <- match(tt[jj], tt[seq_len(jj - 1L)])
      stop(context, ": duplicated time coordinates are not allowed. First duplicate: times ",
           ii, " and ", jj, ".", call. = FALSE)
    }
    invisible(TRUE)
  }

  if (isTRUE(grid) && is.null(coordx_dyn)) {
    axes <- Filter(Negate(is.null), list(coordx, coordy, coordz))
    labs <- c("x", "y", "z")[seq_along(axes)]
    for (j in seq_along(axes)) {
      a <- as.numeric(axes[[j]])
      if (any(!is.finite(a))) stop(context, ": grid coordinates must be finite.", call. = FALSE)
      k <- anyDuplicated(a)
      if (k) stop(context, ": duplicated values in the ", labs[j],
                  " grid axis would create duplicated locations.", call. = FALSE)
    }
    if (isTRUE(spacetime)) check_time(coordt)
    return(invisible(TRUE))
  }

  if (isTRUE(bivariate)) {
    if (!is.null(coordx_dyn)) {
      if (!is.list(coordx_dyn) || length(coordx_dyn) != 2L)
        stop(context, ": bivariate coordx_dyn must contain two coordinate blocks.", call. = FALSE)
      for (j in 1:2) {
        z <- .GeoCoordMatrixForDuplicateCheck(coordx_dyn[[j]], context = context)
        .GeoCheckDuplicateMatrix(z, context,
                                 paste0("observation locations within variable ", j))
      }
    } else {
      z <- .GeoCoordMatrixForDuplicateCheck(coordx, coordy, coordz, context)
      .GeoCheckDuplicateMatrix(z, context, "observation locations")
    }
    return(invisible(TRUE))
  }

  if (isTRUE(spacetime) && !is.null(coordx_dyn)) {
    check_time(coordt)
    if (!is.list(coordx_dyn))
      stop(context, ": coordx_dyn must be a list of coordinate matrices.", call. = FALSE)
    for (j in seq_along(coordx_dyn)) {
      z <- .GeoCoordMatrixForDuplicateCheck(coordx_dyn[[j]], context = context)
      .GeoCheckDuplicateMatrix(z, context,
                               paste0("space-time observation locations within time slice ", j))
    }
    return(invisible(TRUE))
  }

  z <- .GeoCoordMatrixForDuplicateCheck(coordx, coordy, coordz, context)
  .GeoCheckDuplicateMatrix(z, context,
                           if (isTRUE(spacetime)) "spatial locations" else "observation locations")
  if (isTRUE(spacetime)) check_time(coordt)
  invisible(TRUE)
}
