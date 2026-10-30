GeoPit <- function(object, type = c("Uniform", "Gaussian"), data_to_pred = NULL)
{
  type <- match.arg(type)

  valid_class <- inherits(object, "GeoFit") ||
    inherits(object, "GeoKrig") ||
    inherits(object, "GeoKrigloc")

  if (!valid_class) {
    stop("A GeoFit, GeoKrig, or GeoKrigloc object is required.")
  }

  if (isTRUE(object$bivariate)) {
    stop("The spatial bivariate case is not implemented in GeoPit().")
  }

  model <- .GeoMarginalCanonicalModel(object$model)
  EST <- inherits(object, "GeoFit")
  PRED <- inherits(object, "GeoKrig") || inherits(object, "GeoKrigloc")

  if (EST && identical(object$residual_type, "Pearson")) {
    stop(
      "GeoPit() is not defined for Pearson-residual objects. Use the original discrete GeoFit object for PIT diagnostics.",
      call. = FALSE
    )
  }

  if (PRED && !(model %in% c("Gaussian", "Gauss"))) {
    stop(
      "Predictive PIT from GeoKrig/GeoKrigloc is currently supported only for Gaussian models, because non-Gaussian linear kriging does not provide a full conditional predictive distribution.",
      call. = FALSE
    )
  }

  pp <- c(unlist(object$param), unlist(object$fixed))

  getpar <- function(name, default = NA_real_) {
    val <- suppressWarnings(as.numeric(pp[name][1L]))
    if (length(val) == 0L || !is.finite(val)) default else val
  }

  require_par <- function(name) {
    val <- getpar(name)
    if (!is.finite(val)) {
      stop("Parameter '", name, "' is missing in the fitted object.")
    }
    val
  }

  scaled_to_unit <- function(x, lower, upper) {
    if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
      stop("The support parameters must satisfy min < max.")
    }
    (x - lower) / (upper - lower)
  }

  if (EST) {
    if (is.list(object$data)) {
      dd <- unlist(object$data, use.names = FALSE)
    } else {
      dd <- c(t(object$data))
    }
    dd <- as.numeric(dd)

    ## A site-specific external mean vector takes precedence over X beta.
    external_mean <- .GeoMean_external(object$fixed, n = length(dd))

    if (!is.null(external_mean)) {
      allmean <- external_mean
    } else if (!is.null(object$X)) {
      X <- object$X
      if (is.list(X)) {
        X <- as.matrix(do.call(rbind, X))
      } else {
        X <- as.matrix(X)
      }

      if (nrow(X) != length(dd)) {
        stop("The number of rows of X does not match the number of observations.")
      }

      beta <- .GeoMean_beta(pp, ncol(X))
      allmean <- as.numeric(X %*% beta)
    } else {
      allmean <- rep(getpar("mean", 0), length(dd))
    }

    ## Models without a sill parameter use the structural value one.
    allsill <- rep(getpar("sill", 1), length(dd))
  } else {
    target <- data_to_pred
    if (is.null(target) && !is.null(object$data_to_pred)) {
      target <- object$data_to_pred
    }
    if (is.null(target)) {
      stop(
        "For a GeoKrig/GeoKrigloc object, supply the observed prediction targets through data_to_pred.",
        call. = FALSE
      )
    }
    dd <- as.numeric(c(target))
    allmean <- as.numeric(c(object$pred))
    allsill <- as.numeric(c(object$mse))
    if (!length(allsill)) {
      stop("Predictive PIT requires mse=TRUE in GeoKrig/GeoKrigloc.", call. = FALSE)
    }
  }

  dd <- as.numeric(dd)
  allmean <- as.numeric(allmean)
  allsill <- as.numeric(allsill)

  if (length(allmean) == 1L) allmean <- rep(allmean, length(dd))
  if (length(allsill) == 1L) allsill <- rep(allsill, length(dd))

  if (length(allmean) != length(dd)) {
    stop("The length of the mean vector does not match the data length.")
  }
  if (length(allsill) != length(dd)) {
    stop("The length of the variance vector does not match the data length.")
  }
  if (any(!is.finite(allsill)) || any(allsill <= 0)) {
    stop("All variance values used by GeoPit() must be positive and finite.")
  }

  data <- NULL

  #######################################
  ## Gaussian marginal models
  if (model %in% c("Gaussian", "Gaussian_misp_Binomial",
                   "Gaussian_misp_Poisson", "Gaussian_misp_BinomialNeg")) {
    data <- pnorm(dd, mean = allmean, sd = sqrt(allsill))
  }

  #######################################
  ## Positive-support marginal models
  if (model == "Weibull") {
    shape <- require_par("shape")
    if (shape <= 0) stop("The Weibull shape parameter must be positive.")
    data <- pweibull(
      dd,
      shape = shape,
      scale = exp(allmean) / gamma(1 + 1 / shape)
    )
  }

  if (model == "Gamma") {
    shape <- require_par("shape")
    if (shape <= 0) stop("The Gamma shape parameter must be positive.")
    data <- pgamma(
      dd,
      shape = shape / 2,
      rate = shape / (2 * exp(allmean))
    )
  }

  ## GeoModels LogGaussian parametrization:
  ## Y = exp(MM) * exp(sqrt(VV) * Z) / exp(VV / 2).
  if (model == "LogGaussian") {
    data <- suppressWarnings(
      pnorm((log(dd) - allmean + allsill / 2) / sqrt(allsill))
    )
  }

  if (model == "LogLogistic") {
    shape <- require_par("shape")
    if (shape <= 1) {
      stop("The LogLogistic shape parameter must be greater than one.")
    }
    cc <- gamma(1 + 1 / shape) * gamma(1 - 1 / shape)
    scale <- exp(allmean) / cc
    data <- .GeoMarginalPLogLogistic(dd, shape = shape, scale = scale)
  }

  #######################################
  ## Bounded-support marginal models
  if (model == "Beta") {
    shape1 <- require_par("shape1")
    shape2 <- require_par("shape2")
    pmin <- require_par("min")
    pmax <- require_par("max")
    if (shape1 <= 0 || shape2 <= 0) {
      stop("The Beta shape parameters must be positive.")
    }
    zz <- scaled_to_unit(dd, pmin, pmax)
    data <- pbeta(zz, shape1 = shape1, shape2 = shape2)
  }

  if (model == "Beta2") {
    shape <- require_par("shape")
    pmin <- require_par("min")
    pmax <- require_par("max")
    if (shape <= 0) stop("The Beta2 shape parameter must be positive.")
    mu <- plogis(allmean)
    zz <- scaled_to_unit(dd, pmin, pmax)
    data <- pbeta(
      zz,
      shape1 = mu * shape,
      shape2 = (1 - mu) * shape
    )
  }

  if (model == "Kumaraswamy") {
    shape1 <- require_par("shape1")
    shape2 <- require_par("shape2")
    pmin <- require_par("min")
    pmax <- require_par("max")
    if (shape1 <= 0 || shape2 <= 0) {
      stop("The Kumaraswamy shape parameters must be positive.")
    }
    zz <- scaled_to_unit(dd, pmin, pmax)
    data <- numeric(length(zz))
    data[zz >= 1] <- 1
    inside <- which(zz > 0 & zz < 1)
    if (length(inside)) {
      data[inside] <- 1 - (1 - zz[inside]^shape1)^shape2
    }
  }

  if (model == "Kumaraswamy2") {
    shape <- require_par("shape")
    pmin <- require_par("min")
    pmax <- require_par("max")
    if (shape <= 0) {
      stop("The Kumaraswamy2 shape parameter must be positive.")
    }
    mu <- plogis(allmean)
    zz <- scaled_to_unit(dd, pmin, pmax)
    second_shape <- log(0.5) / log1p(-mu^shape)
    if (any(!is.finite(second_shape)) || any(second_shape <= 0)) {
      stop("Kumaraswamy2 produced invalid location-specific shape parameters.")
    }
    data <- numeric(length(zz))
    data[zz >= 1] <- 1
    inside <- which(zz > 0 & zz < 1)
    if (length(inside)) {
      data[inside] <- 1 - (1 - zz[inside]^shape)^second_shape[inside]
    }
  }

  #######################################
  ## Real-valued non-Gaussian marginal models
  if (model == "SkewGaussian") {
    skew <- require_par("skew")
    sill <- getpar("sill", 1)
    omega <- sqrt((skew^2 + allsill) / allsill)
    alpha <- skew / sqrt(sill)
    data <- sn::psn(
      (dd - allmean) / sqrt(allsill),
      xi = 0,
      omega = omega,
      alpha = alpha
    )
  }

  if (model %in% c("StudentT", "Gaussian_misp_StudentT")) {
    inv_df <- require_par("df")
    if (inv_df <= 0) stop("The inverse degrees-of-freedom parameter must be positive.")
    df <- 1 / inv_df
    data <- pt((dd - allmean) / sqrt(allsill), df = df)
  }

  if (model %in% c("SkewStudentT", "Gaussian_misp_SkewStudentT")) {
    inv_df <- require_par("df")
    if (inv_df <= 0) stop("The inverse degrees-of-freedom parameter must be positive.")
    alpha <- .GeoSkewStudentTAlpha(require_par("skew"))
    nu <- 1 / inv_df
    data <- sn::pst(
      (dd - allmean) / sqrt(allsill),
      xi = 0,
      omega = 1,
      alpha = alpha,
      nu = nu
    )
  }

  if (model == "Logistic") {
    data <- plogis((dd - allmean) / sqrt(allsill))
  }

  if (model == "SkewLaplace") {
    skew <- require_par("skew")
    if (skew <= 0 || skew >= 1) {
      stop("The SkewLaplace skew parameter must lie in (0,1).")
    }
    z <- (dd - allmean) / sqrt(allsill)
    data <- .GeoMarginalPSkewLaplace(z, skew)
  }

  if (model == "SinhAsinh") {
    tail <- require_par("tail")
    skew <- require_par("skew")
    if (tail <= 0) stop("The SinhAsinh tail parameter must be positive.")
    x <- (dd - allmean) / sqrt(allsill)
    data <- pnorm(sinh(tail * asinh(x) - skew))
  }

  if (model == "Tukeyh") {
    tail <- require_par("tail")
    if (tail < 0 || tail >= 0.5) stop("The Tukeyh tail parameter must lie in [0,0.5).")
    x <- (dd - allmean) / sqrt(allsill)
    data <- .GeoMarginalPTukeyH(x, tail)
  }

  if (model == "Tukeyh2") {
    tail1 <- require_par("tail1")
    tail2 <- require_par("tail2")
    if (tail1 < 0 || tail1 >= 0.5 || tail2 < 0 || tail2 >= 0.5)
      stop("The Tukeyh2 tail1 and tail2 parameters must lie in [0,0.5).")
    x <- (dd - allmean) / sqrt(allsill)
    data <- .GeoMarginalPTukeyH2(x, tail1, tail2)
  }

  if (model == "TwoPieceGaussian") {
    skew <- require_par("skew")
    x <- (dd - allmean) / sqrt(allsill)
    data <- .GeoMarginalPTwoPieceGaussian(x, skew)
  }

  if (model == "TwoPieceStudentT") {
    skew <- require_par("skew")
    inv_df <- require_par("df")
    if (inv_df <= 0) stop("The inverse degrees-of-freedom parameter must be positive.")
    df <- 1 / inv_df
    x <- (dd - allmean) / sqrt(allsill)
    data <- .GeoMarginalPTwoPieceStudentT(x, skew, df)
  }

  if (model == "TwoPieceTukeyh") {
    skew <- require_par("skew")
    tail <- require_par("tail")
    x <- (dd - allmean) / sqrt(allsill)
    data <- .GeoMarginalPTwoPieceTukeyH(x, skew, tail)
  }

  if (model == "TwoPieceBimodal") {
    skew <- require_par("skew")
    df <- require_par("df")
    delta <- require_par("shape")
    alpha <- 2 * (delta + 1) / df
    nn <- 2^(1 - alpha / 2)
    x <- (dd - allmean) / sqrt(allsill)
    data <- .GeoMarginalSplitCDF(
      x,
      fun_left = function(z) {
        0.5 * (1 + skew) *
          (1 - pgamma((-z / (1 + skew))^alpha,
                      shape = df, scale = 1 / nn))
      },
      fun_right = function(z) {
        0.5 * (1 + skew) +
          0.5 * (1 - skew) *
          pgamma((z / (1 - skew))^alpha,
                 shape = df, scale = 1 / nn)
      }
    )
  }

  #######################################
  ## Discrete marginal models
  if (model == "Binomial") {
    data <- pbinom(dd, size = object$n, prob = pnorm(allmean))
  }

  if (model == "BinomialNeg") {
    data <- pnbinom(dd, size = object$n, prob = pnorm(allmean))
  }

  if (model == "Poisson") {
    data <- ppois(dd, lambda = exp(allmean))
  }

  if (model == "PoissonGamma") {
    shape <- require_par("shape")
    if (shape <= 0) {
      stop("The PoissonGamma shape parameter must be positive.")
    }
    data <- pnbinom(dd, size = shape, mu = exp(allmean))
  }

  if (is.null(data)) {
    stop("PIT is not implemented for model '", model, "' in GeoPit().")
  }

  ## Numerical protection. The raw Uniform PIT is left on [0, 1], while
  ## the Gaussian transform is clipped only to avoid infinite qnorm values.
  data <- pmin(1, pmax(0, data))

  if (type == "Gaussian") {
    eps <- .Machine$double.eps
    data <- qnorm(pmin(1 - eps, pmax(eps, data)))
  }

  object$data <- data
  object
}
