GeoAniso <- function(coords, anisopars = c(0, 1), inverse = FALSE)
{
  if (!is.matrix(coords)) stop(" coords is not a matrix ")

  ncc <- ncol(coords)
  if (ncc != 2 && ncc != 3) stop(" coords dimension is wrong ")

  if (length(anisopars) != 2)
    stop(" anisopars parameters must be a vector with 2 elements (anisotropy angle and anisotropy ratio) ")

  angle <- anisopars[1]
  stretch <- anisopars[2]

  if (stretch < 1) {
    stretch <- round(stretch, digits = 6)
    if (stretch < 1) stop("anisotropy ratio must be greater or equal  than 1 ")
  }

  if (!(angle >= 0 && angle <= pi))
    stop("the anisotropy angle  must be between 0 and pi ")

  ############################
  ### 2D anisotropy
  ############################
  if (ncc == 2) {
    rm <- matrix(c(cos(angle), -sin(angle),
                   sin(angle),  cos(angle)), ncol = 2)
    tm <- diag(c(1, 1 / stretch))

    if (inverse) coordstransf <- coords %*% solve(rm %*% tm)
    else         coordstransf <- coords %*% rm %*% tm
  }

  ############################
  ### 3D coordinates
  ### If z is identically zero, apply the usual 2D anisotropy to x-y and
  ### leave z unchanged. True 3D anisotropy is not implemented.
  ############################
  if (ncc == 3) {
    if (all(coords[, 3] == 0)) {
      rm2 <- matrix(c(cos(angle), -sin(angle),
                      sin(angle),  cos(angle)), ncol = 2)
      tm2 <- diag(c(1, 1 / stretch))

      rm <- rbind(cbind(rm2, c(0, 0)), c(0, 0, 1))
      tm <- rbind(cbind(tm2, c(0, 0)), c(0, 0, 1))

      if (inverse) coordstransf <- coords %*% solve(rm %*% tm)
      else         coordstransf <- coords %*% rm %*% tm
    } else {
      warning("3D anisotropy is not implemented; returning original coordinates.")
      coordstransf <- coords
    }
  }

  return(coordstransf)
}