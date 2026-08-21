####################################################
### File name: GeoDistances.R
####################################################

GeoDistances <- function(coordx = NULL, coordy = NULL, coordz = NULL,
                         distance = c("Eucl", "Chor", "Geod"),
                         radius = 1) {

  distance <- match.arg(distance)

  ##########################################################################
  ## Build and validate the coordinate matrix
  ##########################################################################

  matrix_input <-
    !is.null(coordx) &&
    is.null(coordy) &&
    is.null(coordz) &&
    (is.matrix(coordx) || is.data.frame(coordx))

  if (matrix_input) {

    coords <- as.matrix(coordx)

  } else {

    if (is.matrix(coordx) || is.data.frame(coordx)) {
      stop(
        "When 'coordx' is a matrix or data frame, 'coordy' and 'coordz' ",
        "must be NULL.",
        call. = FALSE
      )
    }

    if (is.null(coordx) || is.null(coordy)) {
      stop(
        "Provide either a coordinate matrix in 'coordx' or vectors ",
        "'coordx' and 'coordy'.",
        call. = FALSE
      )
    }

    if (length(coordx) != length(coordy)) {
      stop(
        "'coordx' and 'coordy' must have the same length.",
        call. = FALSE
      )
    }

    if (!is.null(coordz) && length(coordz) != length(coordx)) {
      stop(
        "'coordz' must have the same length as 'coordx' and 'coordy'.",
        call. = FALSE
      )
    }

    coords <- if (is.null(coordz)) {
      cbind(coordx, coordy)
    } else {
      cbind(coordx, coordy, coordz)
    }
  }

  if (!is.numeric(coords)) {
    stop("Coordinates must be numeric.", call. = FALSE)
  }

  if (length(dim(coords)) != 2L) {
    stop("Coordinates must be a two-dimensional matrix.", call. = FALSE)
  }

  n <- nrow(coords)
  p <- ncol(coords)

  if (n < 1L) {
    stop("At least one coordinate is required.", call. = FALSE)
  }

  if (!(p %in% c(2L, 3L))) {
    stop("Coordinates must have two or three columns.", call. = FALSE)
  }

  if (any(!is.finite(coords))) {
    stop("Coordinates must contain only finite values.", call. = FALSE)
  }

  ##########################################################################
  ## Metric-specific checks
  ##########################################################################

  if (!is.numeric(radius) ||
      length(radius) != 1L ||
      !is.finite(radius) ||
      radius <= 0) {
    stop("'radius' must be a finite positive scalar.", call. = FALSE)
  }

  if (distance %in% c("Geod", "Chor")) {

    if (p != 2L) {
      stop(
        "'Geod' and 'Chor' require exactly two columns: ",
        "longitude and latitude in decimal degrees.",
        call. = FALSE
      )
    }

    latitude <- coords[, 2L]

    if (any(latitude < -90 | latitude > 90)) {
      stop(
        "Latitude values must lie between -90 and 90 degrees.",
        call. = FALSE
      )
    }
  }

  ##########################################################################
  ## Memory and compiled-code index checks
  ##########################################################################

  ## geo_distances() currently uses signed int indices i + j*n in C.
  max_n_c <- floor(sqrt(.Machine$integer.max))

  if (n > max_n_c) {
    stop(
      sprintf(
        paste0(
          "The compiled routine uses 32-bit integer matrix indices and ",
          "cannot safely construct a full distance matrix for n > %d."
        ),
        max_n_c
      ),
      call. = FALSE
    )
  }

  n2 <- as.double(n) * as.double(n)
  output_gb <- 8 * n2 / 1024^3

  ## During the call and matrix construction, more than one copy may exist.
  if (output_gb > 4) {
    warning(
      sprintf(
        paste0(
          "The n x n distance matrix alone requires approximately %.2f GB; ",
          "peak memory usage can be substantially larger."
        ),
        output_gb
      ),
      call. = FALSE
    )
  }

  ##########################################################################
  ## Call the compiled routine
  ##########################################################################

  type_dist <- switch(
    distance,
    Eucl = 0L,
    Chor = 1L,
    Geod = 2L
  )

  out <- numeric(n2)

  res <- dotCall64::.C64(
    "geo_distances",
    SIGNATURE = c(
      "double", "integer", "integer",
      "integer", "double", "double"
    ),
    INTENT = c("r", "r", "r", "r", "r", "w"),
    coords = as.double(coords),
    ncoords = as.integer(n),
    p = as.integer(p),
    type_dist = as.integer(type_dist),
    radius = as.double(radius),
    out = out,
    NAOK = FALSE,
    PACKAGE = "GeoModels",
    VERBOSE = 0
  )$out

  if (length(res) != n2 || any(!is.finite(res))) {
    stop(
      "The compiled distance routine returned an invalid result.",
      call. = FALSE
    )
  }

  ans <- matrix(res, nrow = n, ncol = n, byrow = FALSE)

  ## Preserve observation labels when available.
  rn <- rownames(coords)
  if (!is.null(rn)) {
    dimnames(ans) <- list(rn, rn)
  }

  ans
}