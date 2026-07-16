####################################################
### File name: GeoCorrFct_Cop.R
####################################################

GeoCorrFct_Cop <- function(x, t = NULL, corrmodel, model = "Gaussian", copula = "Gaussian",
                           distance = "Eucl", param, radius = 6371, n = 1,
                           covariance = FALSE, variogram = FALSE)
{
  ############################################################################
  ## C wrappers
  ############################################################################
  biv_unif_CopulaClayton <- function(a, b, c, d) {
    sol <- .C("biv_unif_CopulaClayton_call",
              as.double(a), as.double(b), as.double(c), as.double(d),
              ress = double(1),
              PACKAGE = "GeoModels", DUP = TRUE, NAOK = TRUE)
    exp(sol$ress)
  }

  biv_unif_CopulaGauss <- function(a, b, c) {
    sol <- .C("biv_unif_CopulaGauss_call",
              as.double(a), as.double(b), as.double(c),
              ress = double(1),
              PACKAGE = "GeoModels", DUP = TRUE, NAOK = TRUE)
    sol$ress
  }

  biv_unif_CopulaClayton <- Vectorize(biv_unif_CopulaClayton,
                                      vectorize.args = c("a", "b"))
  biv_unif_CopulaGauss <- Vectorize(biv_unif_CopulaGauss,
                                    vectorize.args = c("a", "b"))

  arg_clayton <- function(x, y, rho, nu, q1, q2) {
    q1(x) * q2(y) * biv_unif_CopulaClayton(x, y, rho, nu)
  }

  arg_gaussian <- function(x, y, rho, q1, q2) {
    q1(x) * q2(y) * biv_unif_CopulaGauss(x, y, rho)
  }

  ############################################################################
  ## Copula-induced correlation
  ############################################################################
  corr_copula <- function(rho, copula, q1, q2, e1, e2, v1, v2, nu) {

    sqrt_v1v2 <- sqrt(v1 * v2)
    e1_e2 <- e1 * e2

    if (!is.finite(sqrt_v1v2) || sqrt_v1v2 <= 0) {
      stop("Invalid marginal variances in corr_copula().")
    }

    if (copula == "Clayton") {
      res <- numeric(length(rho))

      for (i in seq_along(rho)) {
        if (!is.finite(rho[i])) {
          res[i] <- NA_real_
        } else if (rho[i] > 0.9999 && rho[i] <= 1) {
          res[i] <- 1
        } else {
          integral_result <- pracma::integral2(arg_clayton, 0, 1, 0, 1,
                                               rho = rho[i], nu = nu,
                                               q1 = q1, q2 = q2)$Q
          res[i] <- (integral_result - e1_e2) / sqrt_v1v2
        }
      }
      return(res)
    }

    if (copula == "Gaussian") {
      res <- numeric(length(rho))

      for (i in seq_along(rho)) {
        if (!is.finite(rho[i])) {
          res[i] <- NA_real_
        } else if (rho[i] > 0.9999 && rho[i] <= 1) {
          res[i] <- 1
        } else {
          integral_result <- pracma::integral2(arg_gaussian, 0, 1, 0, 1,
                                               rho = rho[i], q1 = q1, q2 = q2)$Q
          res[i] <- (integral_result - e1_e2) / sqrt_v1v2
        }
      }
      return(res)
    }

    stop("Unsupported copula in corr_copula().")
  }

  ############################################################################
  ## Latent correlation from C code
  ############################################################################
  CorrelationFct <- function(bivariate, corrmodel, lags, lagt, numlags, numlagt,
                             mu, model, nuisance, param, N) {
    if (!bivariate) {
      p <- dotCall64::.C64(
        "VectCorrelation",
        SIGNATURE = c("double", "integer", "double", "integer", "integer",
                      "double", "integer", "double", "double", "double", "integer"),
        corr = dotCall64::vector_dc("double", numlags * numlagt),
        corrmodel, lags, numlags, numlagt, mu, model, nuisance, param, lagt, N,
        INTENT = c("rw", rep("r", 10)),
        PACKAGE = "GeoModels", VERBOSE = 0, NAOK = TRUE
      )
      return(p$corr)
    }

    nn <- numlags * 4L
    p <- dotCall64::.C64(
      "VectCorrelation_biv",
      SIGNATURE = c("double", "double", "integer", "double", "integer", "integer",
                    "double", "integer", "double", "double", "double", "integer"),
      corr = dotCall64::vector_dc("double", nn),
      vario = dotCall64::vector_dc("double", nn),
      corrmodel, lags, numlags, numlagt, mu, model, nuisance, param, lagt, N,
      INTENT = c("rw", "rw", rep("r", 10)),
      PACKAGE = "GeoModels", VERBOSE = 0, NAOK = TRUE
    )
    p$corr
  }

  ############################################################################
  ## Checks
  ############################################################################
  call <- match.call()

  if (is.null(CkCorrModel(corrmodel)))
    stop("The name of the correlation model is not correct\n")
  if (is.null(CkModel(model)))
    stop("The name of the model is not correct\n")
  if (!is.numeric(x))
    stop("Distances must be numeric\n")
  if (any(x < 0))
    stop("Distances must be positive\n")
  if (!is.list(param))
    stop("param must be a list\n")

  supported_models <- c("Gaussian", "Beta2")
  if (!(model %in% supported_models)) {
    stop("GeoCorrFct_Cop currently supports only Gaussian and Beta2 margins.")
  }

  supported_copulas <- c("Gaussian", "Clayton")
  if (!(copula %in% supported_copulas)) {
    stop("GeoCorrFct_Cop currently supports only Gaussian and Clayton copulas.")
  }

  spacetime <- CheckST(CkCorrModel(corrmodel))
  bivariate <- CheckBiv(CkCorrModel(corrmodel))

  if (bivariate) {
    stop("GeoCorrFct_Cop is currently implemented only for univariate margins.")
  }

  ############################################################################
  ## Setup
  ############################################################################
  nx <- length(x)

  if (is.null(t)) t <- 0

  if (spacetime) {
    nt <- length(t)
    A <- as.matrix(expand.grid(x, t))
  } else {
    t <- 0
    nt <- 1L
  }

  num_betas <- c(1, 1)
  if (!any(names(param) == "mean")) param$mean <- 0

  mu <- as.numeric(param$mean)

  parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
  nuisance <- c(param)[NuisParam(model, FALSE, num_betas)]

  sel <- substr(names(nuisance), 1, 4) == "mean"
  mm <- as.numeric(nuisance[sel])
  nuisance <- nuisance[!sel]

  if (!any(names(nuisance) == "nugget") || !is.finite(as.numeric(nuisance["nugget"]))) {
    nuisance <- c(nuisance, nugget = 0)
  }

  ############################################################################
  ## Latent correlation
  ############################################################################
  correlation <- CorrelationFct(FALSE, CkCorrModel(corrmodel), x, t, nx, nt, mu,
                                CkModel(model), nuisance, parcorr, n)

  nugget_val <- as.numeric(nuisance["nugget"])

  ## This follows the original implementation: cc is the latent correlation
  ## used inside the copula integral, without adding the point-mass nugget at h=0.
  cc <- correlation * (1 - nugget_val)

  if (length(t) > 1L) {
    correlation <- correlation * (1 - nugget_val) +
      nugget_val * I(A[, 1] == 0 & A[, 2] == 0)
  } else {
    correlation <- correlation * (1 - nugget_val) +
      nugget_val * I(x == 0)
  }

  ############################################################################
  ## Marginal quantiles and moments
  ############################################################################
  if (model == "Beta2") {

    delta <- as.numeric(nuisance["shape"])
    if (!is.finite(delta) || delta <= 0)
      stop("Invalid shape parameter for Beta2 margin.")

    mu1 <- mu2 <- 1 / (1 + exp(-mm))

    a1 <- a2 <- mu1 * delta
    b1 <- b2 <- (1 - mu1) * delta

    q1 <- function(x) qbeta(x, a1, b1)
    q2 <- function(x) qbeta(x, a2, b2)

    e1 <- a1 / (a1 + b1)
    e2 <- a2 / (a2 + b2)
    v1 <- a1 * b1 / ((a1 + b1)^2 * (a1 + b1 + 1))
    v2 <- a2 * b2 / ((a2 + b2)^2 * (a2 + b2 + 1))

    min_val <- suppressWarnings(as.numeric(param["min"]))
    max_val <- suppressWarnings(as.numeric(param["max"]))
    if (!is.finite(min_val)) min_val <- 0
    if (!is.finite(max_val)) max_val <- 1

    dd <- max_val - min_val
    if (!is.finite(dd) || dd <= 0)
      stop("Invalid min/max parameters for Beta2 margin.")

    vs <- sqrt(v1 * v2) * dd^2
  }

  if (model == "Gaussian") {

    sill <- as.numeric(nuisance["sill"])
    if (!is.finite(sill) || sill <= 0)
      stop("Invalid sill parameter for Gaussian margin.")

    v1 <- v2 <- sill
    e1 <- e2 <- mm

    q1 <- function(x) qnorm(x, mean = e1, sd = sqrt(sill))
    q2 <- function(x) qnorm(x, mean = e2, sd = sqrt(sill))

    vs <- sill
  }

  ############################################################################
  ## Copula covariance/correlation
  ############################################################################
  if (copula == "Clayton") {
    nu <- suppressWarnings(as.numeric(param["nu"]))
    if (!is.finite(nu) || nu <= 0)
      stop("Clayton copula requires a positive parameter 'nu'.")

    cova <- corr_copula(cc, "Clayton", q1, q2, e1, e2, v1, v2, nu)

  } else if (copula == "Gaussian") {

    cova <- corr_copula(cc, "Gaussian", q1, q2, e1, e2, v1, v2, 0)

  } else {
    stop("Unsupported copula in GeoCorrFct_Cop().")
  }

  if (!covariance) vs <- 1

  res <- cova * vs
  if (variogram) res <- vs * (1 - cova)

  GeoCorrFct <- list(
    corr = res,
    distances = x,
    times = t,
    model = model,
    distance = distance,
    param = param,
    radius = radius,
    n = n,
    covariance = covariance,
    variogram = variogram,
    spacetime = spacetime,
    bivariate = bivariate
  )

  structure(c(GeoCorrFct, call = call), class = c("GeoCorrFct"))
}