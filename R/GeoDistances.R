GeoDistances <- function(coordx = NULL, coordy = NULL, coordz = NULL,
                         distance = c("Eucl", "Chor", "Geod"), radius = 1) {
  distance <- match.arg(distance)

  # Case 1: coordx is already a matrix of coordinates
  # (and coordy/coordz are not provided)
  if (!is.null(coordx) && is.matrix(coordx) && is.null(coordy) && is.null(coordz)) {
    coords <- coordx
  } else {
    # Case 2: at least coordx and coordy must be provided as vectors
    if (is.null(coordx) || is.null(coordy)) {
      stop("You must provide either a matrix in 'coordx' or the vectors 'coordx' and 'coordy'")
    }
    # Build matrix with optional coordz
    coords <- if (!is.null(coordz)) cbind(coordx, coordy, coordz) else cbind(coordx, coordy)
  }
  coords <- as.matrix(coords)
  n <- nrow(coords); p <- ncol(coords)
  if (!(p %in% c(2,3))) stop("Coordinates must have 2 or 3 columns")
  # Select distance type
  type_dist <- switch(distance, Eucl = 0L, Chor = 1L, Geod = 2L)
  out <- numeric(n * n)

  # Call compiled C routine
  res <- dotCall64::.C64(
    "geo_distances",
    SIGNATURE = c("double","integer","integer","integer","double","double"),
    INTENT    = c("r","r","r","r","r","w"),
    coords    = as.double(coords),
    ncoords   = as.integer(n),
    p         = as.integer(p),
    type_dist = as.integer(type_dist),
    radius    = as.double(radius),
    out       = out,
    NAOK=TRUE, PACKAGE="GeoModels"
  )$out

  # Reshape output into n x n distance matrix
  matrix(res, nrow = n, ncol = n, byrow = FALSE)
}
