GeoAniso <- function(coords,
                     anisopars = c(0, 1),
                     inverse = FALSE) {

  tol = sqrt(.Machine$double.eps)
  if (!is.matrix(coords) || !is.numeric(coords)) {
    stop("coords must be a numeric matrix.")
  }

  ncc <- ncol(coords)

  if (!(ncc %in% c(2L, 3L))) {
    stop("coords must have two or three columns.")
  }

  if (any(!is.finite(coords))) {
    stop("coords contains non-finite values.")
  }

  if (!is.numeric(anisopars) ||
      length(anisopars) != 2L ||
      any(!is.finite(anisopars))) {
    stop(
      "anisopars must contain two finite numeric values: ",
      "angle and anisotropy ratio."
    )
  }

  if (!is.logical(inverse) ||
      length(inverse) != 1L ||
      is.na(inverse)) {
    stop("inverse must be TRUE or FALSE.")
  }

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      !is.finite(tol) ||
      tol < 0) {
    stop("tol must be a finite non-negative scalar.")
  }

  angle <- anisopars[1L]
  stretch <- anisopars[2L]

  if (angle < 0 || angle > pi) {
    stop("The anisotropy angle must be between 0 and pi radians.")
  }

  if (stretch < 1 - tol) {
    stop("The anisotropy ratio must be greater than or equal to 1.")
  }

  # Accept harmless floating-point deviations slightly below one.
  stretch <- max(stretch, 1)

  rotation <- matrix(
    c(
      cos(angle), -sin(angle),
      sin(angle),  cos(angle)
    ),
    nrow = 2L,
    ncol = 2L
  )

  transform_2d <- function(x) {
    if (inverse) {
      # Inverse of rotation %*% diag(c(1, 1 / stretch)).
      x %*% diag(c(1, stretch)) %*% t(rotation)
    } else {
      x %*% rotation %*% diag(c(1, 1 / stretch))
    }
  }

  if (ncc == 2L) {
    coordstransf <- transform_2d(coords)
  } else {
    if (!all(abs(coords[, 3L]) <= tol)) {
      stop(
        "Anisotropy for genuine three-dimensional coordinates ",
        "is not implemented."
      )
    }

    coordstransf <- coords
    coordstransf[, 1:2] <- transform_2d(
      coords[, 1:2, drop = FALSE]
    )
  }

  dimnames(coordstransf) <- dimnames(coords)
  coordstransf
}