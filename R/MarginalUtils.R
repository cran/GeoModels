# Internal marginal-distribution utilities shared by diagnostics and simulation.

.GeoMarginalInvLambertH <- function(x, tail) {
  tail <- as.numeric(tail)[1L]
  if (!is.finite(tail)) stop("The tail parameter must be finite.", call. = FALSE)
  x <- as.numeric(x)
  if (abs(tail) < .Machine$double.eps) return(x)
  value <- sqrt(VGAM::lambertW(tail * x * x) / tail)
  sign(x) * value
}

.GeoMarginalTukeyTail <- function(tail, name = "tail") {
  tail <- as.numeric(tail)
  if (length(tail) != 1L || !is.finite(tail) || tail < 0 || tail >= 0.5) {
    stop(name, " must lie in [0, 0.5).", call. = FALSE)
  }
  tail
}

.GeoMarginalTukeyH <- function(z, tail) {
  z <- as.numeric(z)
  z * exp(0.5 * as.numeric(tail)[1L] * z^2)
}

.GeoMarginalTukeyH2 <- function(z, tail1, tail2) {
  z <- as.numeric(z)
  out <- numeric(length(z))
  pos <- z >= 0
  if (any(pos)) out[pos] <- .GeoMarginalTukeyH(z[pos], tail1)
  if (any(!pos)) out[!pos] <- .GeoMarginalTukeyH(z[!pos], tail2)
  out
}

.GeoMarginalQTukeyH <- function(p, tail) {
  .GeoMarginalTukeyH(stats::qnorm(p), tail)
}

.GeoMarginalQTukeyH2 <- function(p, tail1, tail2) {
  .GeoMarginalTukeyH2(stats::qnorm(p), tail1, tail2)
}

.GeoMarginalDTukeyH <- function(z, tail) {
  tail <- as.numeric(tail)[1L]
  z <- as.numeric(z)
  if (!is.finite(tail) || abs(tail) < .Machine$double.eps) {
    return(stats::dnorm(z))
  }
  b <- .GeoMarginalInvLambertH(z, tail)
  stats::dnorm(b) / (exp(0.5 * tail * b^2) * (1 + tail * b^2))
}

.GeoMarginalQSkewLaplace <- function(p, skew) {
  p <- as.numeric(p)
  skew <- as.numeric(skew)[1L]
  out <- rep(NA_real_, length(p))
  right <- p >= skew
  if (any(right)) {
    out[right] <- -(log1p(-p[right]) - log1p(-skew)) / skew
  }
  if (any(!right)) {
    out[!right] <- (log(p[!right]) - log(skew)) / (1 - skew)
  }
  out
}

.GeoMarginalPSkewLaplace <- function(x, skew) {
  x <- as.numeric(x)
  skew <- as.numeric(skew)[1L]
  out <- numeric(length(x))
  left <- x < 0
  if (any(left)) out[left] <- skew * exp((1 - skew) * x[left])
  if (any(!left)) out[!left] <- 1 - (1 - skew) * exp(-skew * x[!left])
  out
}

.GeoMarginalQTwoPieceGaussian <- function(p, skew) {
  p <- as.numeric(p); skew <- as.numeric(skew)[1L]
  out <- numeric(length(p))
  split <- 0.5 * (1 + skew)
  left <- p > 0 & p < split
  right <- p >= split & p <= 1
  if (any(left)) out[left] <- (1 + skew) * stats::qnorm(p[left] / (1 + skew))
  if (any(right)) out[right] <- (1 - skew) * stats::qnorm((p[right] - skew) / (1 - skew))
  out
}

.GeoMarginalQTwoPieceStudentT <- function(p, skew, df) {
  p <- as.numeric(p); skew <- as.numeric(skew)[1L]; df <- as.numeric(df)[1L]
  out <- numeric(length(p))
  split <- 0.5 * (1 + skew)
  left <- p > 0 & p < split
  right <- p >= split & p <= 1
  if (any(left)) out[left] <- (1 + skew) * stats::qt(p[left] / (1 + skew), df = df)
  if (any(right)) out[right] <- (1 - skew) * stats::qt((p[right] - skew) / (1 - skew), df = df)
  out
}

.GeoMarginalQTwoPieceTukeyH <- function(p, skew, tail) {
  p <- as.numeric(p); skew <- as.numeric(skew)[1L]
  out <- numeric(length(p))
  split <- 0.5 * (1 + skew)
  left <- p > 0 & p < split
  right <- p >= split & p <= 1
  if (any(left)) {
    out[left] <- (1 + skew) * .GeoMarginalQTukeyH(p[left] / (1 + skew), tail)
  }
  if (any(right)) {
    out[right] <- (1 - skew) * .GeoMarginalQTukeyH((p[right] - skew) / (1 - skew), tail)
  }
  out
}

.GeoMarginalInvTukeyH2 <- function(x, tail1, tail2) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  pos <- x >= 0
  if (any(pos)) out[pos] <- .GeoMarginalInvLambertH(x[pos], tail1)
  if (any(!pos)) out[!pos] <- .GeoMarginalInvLambertH(x[!pos], tail2)
  out
}

.GeoMarginalQLogLogistic <- function(p, shape, rate = 1, scale = 1 / rate,
                                     lower.tail = TRUE, log.p = FALSE) {
  if (missing(shape)) stop("argument 'shape' is missing, with no default")
  if (shape <= 0) stop("'shape' must be positive")
  if (scale <= 0) stop("'scale' must be positive")
  if (log.p) p <- exp(p)
  if (any(p < 0 | p > 1)) stop("probabilities must be between 0 and 1")
  if (!lower.tail) p <- 1 - p
  out <- numeric(length(p))
  out[p == 0] <- 0
  out[p == 1] <- Inf
  valid <- p > 0 & p < 1
  if (any(valid)) {
    pv <- p[valid]
    out[valid] <- scale * (pv / (1 - pv))^(1 / shape)
  }
  out
}

.GeoMarginalDLogLogistic <- function(x, shape, rate = 1, scale = 1 / rate) {
  if (missing(shape)) stop("argument 'shape' is missing, with no default")
  if (shape <= 0) stop("'shape' must be positive")
  if (scale <= 0) stop("'scale' must be positive")
  out <- numeric(length(x))
  pos <- x > 0 & is.finite(x)
  z <- x[pos] / scale
  out[pos] <- (shape / scale) * z^(shape - 1) / (1 + z^shape)^2
  out
}

# Canonical aliases shared by marginal diagnostics.
.GeoMarginalModelAliases <- c(
  Gauss = "Gaussian",
  SkewGauss = "SkewGaussian",
  LogGauss = "LogGaussian",
  Loglogistic = "LogLogistic",
  TwoPieceGauss = "TwoPieceGaussian",
  tukeyh = "Tukeyh",
  tukeyh2 = "Tukeyh2",
  gamma = "Gamma",
  weibull = "Weibull"
)

.GeoMarginalCanonicalModel <- function(model) {
  model <- as.character(model)[1L]
  if (model %in% names(.GeoMarginalModelAliases)) {
    model <- unname(.GeoMarginalModelAliases[[model]])
  }
  model
}

.GeoMarginalSplitCDF <- function(x, fun_left, fun_right) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  right <- x >= 0
  if (any(right)) out[right] <- fun_right(x[right])
  if (any(!right)) out[!right] <- fun_left(x[!right])
  out
}

.GeoMarginalPTukeyH <- function(x, tail) {
  stats::pnorm(.GeoMarginalInvLambertH(x, tail))
}

.GeoMarginalPTukeyH2 <- function(x, tail1, tail2) {
  stats::pnorm(.GeoMarginalInvTukeyH2(x, tail1, tail2))
}

.GeoMarginalDSkewLaplace <- function(x, skew) {
  x <- as.numeric(x)
  skew <- as.numeric(skew)[1L]
  skew * (1 - skew) * ifelse(
    x < 0,
    exp((1 - skew) * x),
    exp(-skew * x)
  )
}

.GeoMarginalPTwoPieceGaussian <- function(x, skew) {
  skew <- as.numeric(skew)[1L]
  .GeoMarginalSplitCDF(
    x,
    fun_left = function(z) (1 + skew) * stats::pnorm(z / (1 + skew)),
    fun_right = function(z) skew + (1 - skew) * stats::pnorm(z / (1 - skew))
  )
}

.GeoMarginalPTwoPieceStudentT <- function(x, skew, df) {
  skew <- as.numeric(skew)[1L]
  df <- as.numeric(df)[1L]
  .GeoMarginalSplitCDF(
    x,
    fun_left = function(z) (1 + skew) * stats::pt(z / (1 + skew), df = df),
    fun_right = function(z) skew + (1 - skew) * stats::pt(z / (1 - skew), df = df)
  )
}

.GeoMarginalPTwoPieceTukeyH <- function(x, skew, tail) {
  skew <- as.numeric(skew)[1L]
  .GeoMarginalSplitCDF(
    x,
    fun_left = function(z) (1 + skew) * .GeoMarginalPTukeyH(z / (1 + skew), tail),
    fun_right = function(z) skew + (1 - skew) * .GeoMarginalPTukeyH(z / (1 - skew), tail)
  )
}

.GeoMarginalPLogLogistic <- function(x, shape, rate = 1, scale = 1 / rate) {
  if (missing(shape)) stop("argument 'shape' is missing, with no default")
  if (length(shape) != 1L || !is.finite(shape) || shape <= 0) {
    stop("'shape' must be a positive finite scalar")
  }
  scale <- as.numeric(scale)
  if (any(scale <= 0, na.rm = TRUE)) stop("'scale' must be positive")
  z <- (x / scale)^shape
  z / (1 + z)
}
