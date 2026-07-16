GeoPit <- function(object, type = "Uniform")
{
  if (!(type == "Uniform" || type == "Gaussian"))
    stop("The type parameter can be Uniform or Gaussian")

  if (!(inherits(object, "GeoFit") || inherits(object, "GeoKrig") || inherits(object, "GeoKrigloc")))
    stop("A GeoFit, GeoKrig or GeoKrigloc object is needed as input\n")

  model <- object$model
  EST <- inherits(object, "GeoFit")

  pp <- c(unlist(object$param), unlist(object$fixed))

  ## GeoPit is intended to be applied after GeoResiduals() when
  ## a fitted GeoFit object contains covariates.  Prediction objects
  ## are allowed because pred/mse are already site-specific.
  if (EST) {
    mean_pars <- names(pp)[substr(names(pp), 1, 4) == "mean"]
    if (!is.null(object$X) && length(mean_pars) > 1) {
      stop("GeoPit() is intended for residual GeoFit objects when covariates are present. ",
           "Please use GeoResiduals(object) before calling GeoPit().")
    }
  }

  getpar <- function(name, default = NA_real_) {
    val <- suppressWarnings(as.numeric(pp[name][1]))
    if (length(val) == 0 || is.na(val)) default else val
  }

  inv_lamb <- function(x, tail) {
    if (!is.finite(tail) || abs(tail) < .Machine$double.eps)
      return(x)
    value <- sqrt(VGAM::lambertW(tail * x * x) / tail)
    sign(x) * value
  }

  split_cdf <- function(x, fun_left, fun_right) {
    res <- numeric(length(x))
    id_right <- which(x >= 0)
    id_left  <- which(x < 0)
    if (length(id_right) > 0) res[id_right] <- fun_right(x[id_right])
    if (length(id_left)  > 0) res[id_left]  <- fun_left(x[id_left])
    res
  }

  if (EST) {
    if (is.list(object$data)) dd <- unlist(object$data)
    else dd <- c(t(object$data))

    ## If an external fixed mean vector was used, keep it as the PIT mean.
    if (!is.null(object$fixed$mean) && length(object$fixed$mean) > 1) {
      allmean <- as.numeric(object$fixed$mean)
    } else if (!is.null(object$X)) {
      X <- object$X
      if (is.list(X)) X <- as.matrix(do.call(rbind, X))
      else X <- as.matrix(X)

      namesnuis <- NuisParam(object$model, object$bivariate, object$numbetas)
      nuisance <- pp[namesnuis]
      sel <- substr(names(nuisance), 1, 4) == "mean"
      beta <- as.numeric(nuisance[sel])
      beta <- beta[is.finite(beta)]

      if (length(beta) == ncol(X)) {
        allmean <- as.numeric(X %*% beta)
      } else {
        allmean <- rep(getpar("mean", 0), length(dd))
      }
    } else {
      allmean <- rep(getpar("mean", 0), length(dd))
    }

    allsill <- rep(getpar("sill", 1), length(dd))
  } else {
    dd <- as.numeric(c(object$data_to_pred))
    allmean <- as.numeric(c(object$pred))
    allsill <- as.numeric(c(object$mse))
  }

  dd <- as.numeric(c(dd))
  allmean <- as.numeric(c(allmean))
  allsill <- as.numeric(c(allsill))

  if (length(allmean) == 1) allmean <- rep(allmean, length(dd))
  if (length(allsill) == 1) allsill <- rep(allsill, length(dd))

  if (length(allmean) != length(dd))
    stop("length of mean vector does not match data length in GeoPit()")
  if (length(allsill) != length(dd))
    stop("length of variance vector does not match data length in GeoPit()")

  if (object$bivariate)
    stop("The spatial bivariate case is not implemented yet")

  data <- NULL

  #######################################
  if (model %in% c("Gaussian", "Gaussian_misp_Binomial",
                   "Gaussian_misp_Poisson", "Gaussian_misp_BinomialNeg")) {
    MM <- allmean
    VV <- allsill
    data <- pnorm(dd, mean = MM, sd = sqrt(VV))
  }

  #######################################
  if (model %in% c("Weibull")) {
    MM <- allmean
    sh <- getpar("shape")
    data <- pweibull(dd, shape = sh, scale = exp(MM) / gamma(1 + 1 / sh))
  }

  #######################################
  if (model %in% c("Beta2")) {
    MM <- allmean
    mm <- 1 / (1 + exp(-MM))
    sh <- getpar("shape")
    pmin <- getpar("min")
    pmax <- getpar("max")
    data <- pbeta((dd - pmin) / (pmax - pmin), shape1 = mm * sh, shape2 = (1 - mm) * sh)
  }

  #######################################
  if (model %in% c("Kumaraswamy2")) {
    MM <- allmean
    mm <- 1 / (1 + exp(-MM))
    sh <- getpar("shape")
    pmin <- getpar("min")
    pmax <- getpar("max")
    shape1 <- log(0.5) / log1p(-mm^sh)
    data <- 1 - (1 - ((dd - pmin) / (pmax - pmin))^sh)^shape1
  }

  #######################################
  if (model %in% c("SkewGaussian")) {
    MM <- allmean
    VV <- allsill
    sk <- getpar("skew")
    omega <- sqrt((sk^2 + VV) / VV)
    alpha <- sk / sqrt(getpar("sill", 1))
    data <- sn::psn((dd - MM) / sqrt(VV), xi = 0, omega = omega, alpha = alpha)
  }

  #######################################
  if (model %in% c("StudentT", "Gaussian_misp_StudentT")) {
    MM <- allmean
    VV <- allsill
    data <- pt((dd - MM) / sqrt(VV), df = round(1 / getpar("df")))
  }

  #######################################
  if (model %in% c("SkewStudentT", "Gaussian_misp_SkewStudentT")) {
    MM <- allmean
    VV <- allsill
    alpha <- getpar("skew")
    nu <- round(1 / getpar("df"))
    data <- sn::pst((dd - MM) / sqrt(VV), xi = 0, omega = 1, alpha = alpha, nu = nu)
  }

  #######################################
  if (model %in% c("Gamma")) {
    MM <- allmean
    shape <- getpar("shape")
    data <- pgamma(dd, shape = shape / 2, rate = shape / (2 * exp(MM)))
  }

  #######################################
  ## GeoModels LogGaussian parametrization from GeoSim:
  ## Y = exp(MM) * exp(sqrt(VV) * Z) / exp(VV / 2), Z ~ N(0,1).
  if (model %in% c("LogGaussian", "LogGauss")) {
    MM <- allmean
    VV <- allsill
    data <- suppressWarnings(pnorm((log(dd) - MM + VV / 2) / sqrt(VV)))
  }

  #######################################
  if (model %in% c("LogLogistic")) {
    pllogis1 <- function(q, shape, rate = 1, scale = 1 / rate,
                         lower.tail = TRUE, log.p = FALSE) {
      x <- (q / scale)^shape
      p <- x / (1 + x)
      if (!lower.tail) p <- 1 - p
      if (log.p) p <- log(p)
      p
    }
    MM <- allmean
    shape <- getpar("shape")
    cc <- gamma(1 + 1 / shape) * gamma(1 - 1 / shape)
    data <- pllogis1(dd, shape = shape, scale = exp(MM) / cc)
  }

  #######################################
  if (model %in% c("Logistic")) {
    MM <- allmean
    VV <- allsill
    data <- (1 + exp(-(dd - MM) / sqrt(VV)))^(-1)
  }

  #######################################
  if (model %in% c("SinhAsinh")) {
    MM <- allmean
    VV <- allsill
    tail <- getpar("tail")
    skew <- getpar("skew")
    data <- pnorm(sinh(tail * asinh((dd - MM) / sqrt(VV)) - skew))
  }

  #######################################
  if (model %in% c("Tukeyh")) {
    MM <- allmean
    VV <- allsill
    tail <- getpar("tail")
    x <- (dd - MM) / sqrt(VV)
    data <- pnorm(inv_lamb(x, tail))
  }

  #######################################
  if (model %in% c("Tukeyh2")) {
    MM <- allmean
    VV <- allsill
    tail1 <- getpar("tail1")
    tail2 <- getpar("tail2")
    x <- (dd - MM) / sqrt(VV)
    data <- numeric(length(x))
    id_pos <- which(x > 0)
    id_neg <- which(x <= 0)
    if (length(id_pos) > 0) data[id_pos] <- pnorm(inv_lamb(x[id_pos], tail1))
    if (length(id_neg) > 0) data[id_neg] <- pnorm(inv_lamb(x[id_neg], tail2))
  }

  #######################################
  if (model %in% c("TwoPieceGaussian")) {
    MM <- allmean
    VV <- allsill
    skew <- getpar("skew")
    x <- (dd - MM) / sqrt(VV)
    data <- split_cdf(
      x,
      fun_left  = function(z) (1 + skew) * pnorm(z / (1 + skew)),
      fun_right = function(z) skew + (1 - skew) * pnorm(z / (1 - skew))
    )
  }

  #######################################
  if (model %in% c("TwoPieceStudentT")) {
    MM <- allmean
    VV <- allsill
    skew <- getpar("skew")
    df <- 1 / getpar("df")
    x <- (dd - MM) / sqrt(VV)
    data <- split_cdf(
      x,
      fun_left  = function(z) (1 + skew) * pt(z / (1 + skew), df = df),
      fun_right = function(z) skew + (1 - skew) * pt(z / (1 - skew), df = df)
    )
  }

  #######################################
  if (model %in% c("TwoPieceTukeyh")) {
    MM <- allmean
    VV <- allsill
    skew <- getpar("skew")
    tail <- getpar("tail")
    pTukeyh <- function(z) pnorm(inv_lamb(z, tail))
    x <- (dd - MM) / sqrt(VV)
    data <- split_cdf(
      x,
      fun_left  = function(z) (1 + skew) * pTukeyh(z / (1 + skew)),
      fun_right = function(z) skew + (1 - skew) * pTukeyh(z / (1 - skew))
    )
  }

  #######################################
  if (model %in% c("TwoPieceBimodal")) {
    MM <- allmean
    VV <- allsill
    skew <- getpar("skew")
    df <- getpar("df")
    delta <- getpar("shape")
    alpha <- 2 * (delta + 1) / df
    nn <- 2^(1 - alpha / 2)
    x <- (dd - MM) / sqrt(VV)
    data <- split_cdf(
      x,
      fun_left = function(z) {
        0.5 * (1 + skew) *
          (1 - pgamma((-z / (1 + skew))^alpha, shape = df, scale = 1 / nn))
      },
      fun_right = function(z) {
        0.5 * (1 + skew) +
          0.5 * (1 - skew) * pgamma((z / (1 - skew))^alpha, shape = df, scale = 1 / nn)
      }
    )
  }

  #######################################
  ## discrete
  if (model %in% c("Binomial")) {
    MM <- allmean
    data <- pbinom(dd, size = object$n, prob = pnorm(MM))
  }

  if (model %in% c("BinomialNeg")) {
    MM <- allmean
    data <- pnbinom(dd, size = object$n, prob = pnorm(MM))
  }

  if (model %in% c("Poisson")) {
    MM <- allmean
    data <- ppois(dd, lambda = exp(MM))
  }

  if (is.null(data))
    stop("PIT is not implemented for this model in GeoPit()")

  if (type == "Gaussian") {
    eps <- .Machine$double.eps
    data <- pmin(1 - eps, pmax(eps, data))
    data <- qnorm(data)
  }

  object$data <- data
  return(object)
}