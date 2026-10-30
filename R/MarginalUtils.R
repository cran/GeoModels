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

.GeoSkewStudentTAlpha <- function(skew) {
  delta <- as.numeric(skew)
  if (length(delta) != 1L || !is.finite(delta) || abs(delta) >= 1) {
    stop("For SkewStudentT, skew must lie strictly between -1 and 1.",
         call. = FALSE)
  }
  ## GeoModels uses the latent skewness delta; sn uses the shape alpha.
  delta / sqrt(1 - delta^2)
}

.GeoMarginalPLogLogistic <- function(x, shape, rate = 1, scale = 1 / rate) {
  if (missing(shape)) stop("argument 'shape' is missing, with no default")
  if (length(shape) != 1L || !is.finite(shape) || shape <= 0) {
    stop("'shape' must be a positive finite scalar")
  }
  scale <- as.numeric(scale)
  if (any(scale <= 0, na.rm = TRUE)) stop("'scale' must be positive")
  ## Work on the log scale to avoid overflow in (x / scale)^shape.
  ## pmax maps the nonpositive support to log(0), while preserving NA/NaN.
  stats::plogis(shape * (log(pmax(x, 0)) - log(scale)))
}


## E|T_nu| under GeoModels' convention nu = 1/df.  The ratio expansion
## prevents cancellation/overflow for large nu; nu=Inf is its exact limit.
.GeoStudentAbsMean <- function(nu) {
  nu <- as.numeric(nu)
  if (length(nu) != 1L || is.na(nu) || nu <= 1)
    stop("The Student marginal mean requires nu > 1.", call. = FALSE)
  if (is.infinite(nu)) return(sqrt(2 / pi))
  x <- 0.5 * (nu - 2)
  if (x >= 1000) {
    r <- 1 / x
    log_a <- r * (-0.25 + r^2 * (1 / 96 - r^2 / 320))
    return(exp(0.5 * (log(2 / pi) - log1p(-2 / nu) + log_a)))
  }
  exp(0.5 * (log(nu) - log(pi)) +
        lgamma(0.5 * (nu - 1)) - lgamma(0.5 * nu))
}

## Observed-scale marginal mean for linear-prediction diagnostics.
## Keep direct and copula parametrizations separate (notably Kumaraswamy).
## Unsupported models must not silently be treated as having zero mean.
.GeoMarginalMean <- function(model_type, params, eta, n = 1, copula = NULL) {
  model_type <- as.integer(model_type)
  eta <- as.numeric(eta)
  if (length(model_type) != 1L || is.na(model_type) ||
      !length(eta) || any(!is.finite(eta)))
    stop("Invalid model or marginal location vector.", call. = FALSE)
  scalar <- function(name) {
    value <- suppressWarnings(as.numeric(params[[name]]))
    if (length(value) != 1L || !is.finite(value))
      stop("A finite scalar parameter '", name, "' is required.", call. = FALSE)
    value
  }
  sizes <- function() {
    nn <- as.numeric(n)
    if (length(nn) == 1L) nn <- rep(nn, length(eta))
    if (length(nn) != length(eta) || any(!is.finite(nn)) || any(nn < 0))
      stop("n must contain one valid size or one size per observation.", call. = FALSE)
    nn
  }
  bounded <- function(mean01) {
    lo <- scalar("min"); hi <- scalar("max")
    if (hi <= lo) stop("Bounded margins require min < max.", call. = FALSE)
    lo + (hi - lo) * mean01
  }
  if (model_type == 13L)
    stop("GeoDoScores does not define circular scores for Wrapped; use circular prediction diagnostics, not Euclidean scores on angles.",
         call. = FALSE)

  # Count means; the links do not depend on the copula, and n may vary by site.
  if (model_type %in% c(2L, 11L, 51L)) return(sizes() * pnorm(eta))
  if (model_type == 49L) return(sizes() * plogis(eta))
  if (model_type %in% c(14L, 16L, 45L, 52L, 56L)) {
    nn <- if (model_type == 14L) 1 else sizes()
    mm <- nn * pnorm(eta, lower.tail = FALSE) / pnorm(eta)
    if (model_type %in% c(45L, 56L)) mm <- mm * pnorm(scalar("pmu"), lower.tail = FALSE)
    return(as.numeric(mm))
  }
  if (model_type %in% c(30L, 36L, 46L, 47L)) return(exp(eta))
  if (model_type %in% c(43L, 44L, 53L, 57L, 58L))
    return(pnorm(scalar("pmu"), lower.tail = FALSE) * exp(eta))

  if (!is.null(copula))
    return(as.numeric(.GeoCopulaMarginalMoments(model_type, params, eta)$mean))

  if (model_type %in% c(1L, 12L, 35L, 34L, 25L)) return(eta)
  if (model_type %in% c(21L, 22L, 24L, 26L)) return(exp(eta))
  if (model_type == 10L) return(eta + scalar("skew") * sqrt(2 / pi))
  if (model_type %in% c(18L, 27L, 37L)) {
    df_param <- scalar("df")
    if (df_param < 0) stop("df cannot be negative.", call. = FALSE)
    nu <- if (df_param == 0) Inf else 1 / df_param
    shift <- sqrt(scalar("sill")) * scalar("skew") * .GeoStudentAbsMean(nu)
    return(eta + if (model_type == 27L) -2 * shift else shift)
  }
  if (model_type == 29L)
    return(eta - 2 * scalar("skew") * sqrt(2 * scalar("sill") / pi))
  if (model_type == 38L)
    return(eta - 2 * scalar("skew") * sqrt(2 * scalar("sill") / pi) /
             (1 - scalar("tail")))
  if (model_type %in% c(20L, 40L))
    return(as.numeric(.GeoPlainTransformedMarginalMoments(model_type, params, eta)$mean))
  if (model_type %in% c(9L, 41L)) {
    g <- scalar("skew"); h <- scalar("tail")
    if (h < 0 || h >= 1) stop("The Tukeygh mean requires 0 <= tail < 1.", call. = FALSE)
    shift <- if (g == 0) 0 else expm1(g^2 / (2 * (1 - h))) / (g * sqrt(1 - h))
    return(eta + sqrt(scalar("sill")) * shift)
  }
  if (model_type == 39L) {
    nu <- scalar("df"); shape <- scalar("shape")
    if (nu <= 0 || shape < 0)
      stop("TwoPieceBimodal requires df > 0 and shape >= 0.", call. = FALSE)
    alpha <- 2 * (shape + 1) / nu
    abs_mean <- sqrt(2) * exp(lgamma(nu / 2 + 1 / alpha) - lgamma(nu / 2))
    return(eta - 2 * scalar("skew") * sqrt(scalar("sill")) * abs_mean)
  }
  if (model_type == 59L) {
    sk <- scalar("skew")
    if (sk <= 0 || sk >= 1) stop("SkewLaplace skew must lie in (0,1).", call. = FALSE)
    return(eta + sqrt(scalar("sill")) * (1 / sk - 1 / (1 - sk)))
  }
  if (model_type == 28L) {
    a <- scalar("shape1"); b <- scalar("shape2")
    return(rep(bounded(a / (a + b)), length(eta)))
  }
  if (model_type %in% c(33L, 42L)) {
    if (model_type == 33L) {
      # Direct model: F(x)=1-(1-x^shape2)^shape1, unlike copula margins.
      a <- scalar("shape2"); b <- scalar("shape1")
    } else {
      a <- scalar("shape")
      b <- log(0.5) / log1p(-plogis(eta)^a)
    }
    if (a <= 0 || any(!is.finite(b)) || any(b <= 0))
      stop("Invalid Kumaraswamy shape parameters.", call. = FALSE)
    mm <- bounded(exp(log(b) + lbeta(1 + 1 / a, b)))
    return(if (length(mm) == 1L) rep(mm, length(eta)) else as.numeric(mm))
  }
  stop("GeoDoScores marginal centering is not implemented for model id ",
       model_type, ".", call. = FALSE)
}
