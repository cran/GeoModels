GeoSimcond <- function(estobj = NULL, data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
                       coordx_dyn = NULL, corrmodel, distance = "Eucl",
                       grid = FALSE, loc, maxdist = NULL, maxtime = NULL, method = "Cholesky",
                       model = "Gaussian", n = 1, nrep = 1, local = FALSE, L = 1000,
                       neighb = NULL, param, anisopars = NULL, radius = 1, sparse = FALSE,
                       time = NULL, copula = NULL, X = NULL, Xloc = NULL, Mloc = NULL,
                       parallel = FALSE, ncores = 6, progress = FALSE, n_iter = 25L,
                       check.duplicates = FALSE, nloc = NULL, mcmc_thin = 1L) {

################################################################################# 
### Helper function for conditional skew gaussian simulation 
################################################################################# 
  # Model-specific conditional engines are implemented in SimcondUtils.R.


###############################################################
############### end helper functions ##########################
###############################################################

###############################################################
############################ MAIN FUNCTION ####################
###############################################################
  call <- match.call()

  ###### Check GeoFit object ######
  if (!is.null(estobj)) {
    if (!inherits(estobj, "GeoFit"))
      stop("need a 'GeoFit' object as input\n")

    data <- estobj$data
    if (!estobj$grid) {
      if (!estobj$bivariate) {
        if (is.null(estobj$coordx_dyn)) coordx <- cbind(estobj$coordx, estobj$coordy, estobj$coordz)
        else coordx <- NULL
      } else {
        if (is.null(estobj$coordx_dyn)) {
          coordx <- estobj$coordx[1:estobj$ns[1]]
          coordy <- estobj$coordy[1:estobj$ns[2]]
        } else {
          coordx <- NULL
        }
      }
    } else {
      coordx <- estobj$coordx; coordy <- estobj$coordy; coordz <- estobj$coordz
    }
    if (length(estobj$coordt) == 1) coordt <- NULL else coordt <- estobj$coordt
    coordx_dyn <- estobj$coordx_dyn
    corrmodel  <- estobj$corrmodel
    model      <- estobj$model
    distance   <- estobj$distance
    grid       <- estobj$grid
    n          <- estobj$n
    param      <- append(estobj$param, estobj$fixed)
    radius     <- estobj$radius
    copula     <- estobj$copula
    anisopars  <- estobj$anisopars
    X          <- estobj$X
  }
  #### end check

  ### Additional checks
  if (is.null(CkModel(model))) stop("The name of the model is not correct\n")
  if (!is.character(corrmodel) || is.null(CkCorrModel(corrmodel)))
    stop("The name of the correlation model is wrong\n")
  if (!(method %in% c("Cholesky", "TB", "CE"))) {
    stop("The method of unconditional simulation is not correct\n")
  }
  corrmodel <- gsub("[[:blank:]]", "", corrmodel)
  model     <- gsub("[[:blank:]]", "", model)
  distance  <- gsub("[[:blank:]]", "", distance)

  # Canonical aliases used throughout GeoModels.
  if (model == "Gauss") model <- "Gaussian"
  if (model == "LogGauss") model <- "LogGaussian"
  if (model == "SkewGauss") model <- "SkewGaussian"
  if (model == "poisson") model <- "Poisson"
  if (model %in% c("Binary", "Bernoulli")) {
    if (!is.null(nloc) && any(as.numeric(nloc) != 1))
      stop("Binary/Bernoulli conditional simulation requires nloc = 1 when supplied.", call. = FALSE)
    model <- "Binomial"
    n <- 1L
    if (is.null(nloc)) nloc <- 1L
  }
  if (model %in% c("Geom", "Geometric")) {
    if (!is.null(nloc) && any(as.numeric(nloc) != 1))
      stop("Geometric conditional simulation requires nloc = 1 when supplied.", call. = FALSE)
    model <- "BinomialNeg"
    n <- 1L
    if (is.null(nloc)) nloc <- 1L
  }

  # The count-constrained sampler mixes rapidly under the validated direct
  # Binomial/Negative-Binomial paths, so its user-facing default is 25 sweeps.
  # Other latent Gibbs engines retain their historical 1000-sweep default when
  # n_iter is omitted; an explicitly supplied n_iter is always respected.
  direct_count_model <- is.null(copula) && model %in% c("Binomial", "BinomialNeg")
  if (missing(n_iter) && !direct_count_model) n_iter <- 1000L

  #### checking cores ######
  ncores <- .GeoResolveWorkers(parallel, ncores)
  parallel <- isTRUE(parallel) && ncores > 1L

  # Determine model type (space / spacetime / bivariate)
  model_ck  <- CkCorrModel(corrmodel)
  bivariate <- CheckBiv(model_ck)
  spacetime <- CheckST(model_ck)
  space     <- !spacetime && !bivariate

  # Phase-1 support is intentionally conservative.  Unsupported domains are
  # rejected before covariance matrices or unconditional simulations are built.
  if (bivariate) {
    stop("GeoSimcond currently supports univariate correlation models only.", call. = FALSE)
  }
  if (spacetime || !is.null(coordt) || !is.null(time)) {
    stop("GeoSimcond space-time conditional simulation is not implemented yet.", call. = FALSE)
  }
  if (!is.null(coordx_dyn)) {
    stop("GeoSimcond dynamic coordinates are not implemented yet.", call. = FALSE)
  }
  if (isTRUE(grid)) {
    stop("GeoSimcond grid = TRUE is not implemented yet; supply explicit coordinates.", call. = FALSE)
  }
  if (!is.logical(local) || length(local) != 1L || is.na(local)) {
    stop("local must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(sparse)) {
    stop("GeoSimcond sparse = TRUE is not implemented yet.", call. = FALSE)
  }
  if (identical(method, "CE")) {
    stop(paste0(
      "GeoSimcond method = 'CE' is not available for the current explicit-coordinate ",
      "conditional simulation workflow; use method = 'Cholesky' or method = 'TB'."
    ), call. = FALSE)
  }

  .GeoCheckDuplicateLocations(
    coordx = coordx, coordy = coordy, coordz = coordz,
    spacetime = FALSE, bivariate = FALSE, grid = FALSE,
    context = "GeoSimcond", check.duplicates = check.duplicates
  )
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

  if (!is.numeric(nrep) || length(nrep) != 1L || !is.finite(nrep) ||
      nrep < 1 || nrep != as.integer(nrep)) {
    stop("nrep must be a positive integer.", call. = FALSE)
  }
  nrep <- as.integer(nrep)
  if (!is.numeric(n_iter) || length(n_iter) != 1L || !is.finite(n_iter) ||
      n_iter < 1 || n_iter != as.integer(n_iter))
    stop("n_iter must be a positive integer.", call. = FALSE)
  n_iter <- as.integer(n_iter)
  if (!is.numeric(mcmc_thin) || length(mcmc_thin) != 1L || !is.finite(mcmc_thin) ||
      mcmc_thin < 1 || mcmc_thin != as.integer(mcmc_thin))
    stop("mcmc_thin must be a positive integer.", call. = FALSE)
  mcmc_thin <- as.integer(mcmc_thin)

  if (is.vector(loc)) loc <- matrix(loc, nrow = 1L)
  loc <- as.matrix(loc)
  if (!is.numeric(loc) || nrow(loc) < 1L || any(!is.finite(loc))) {
    stop("loc must be a non-empty finite numeric matrix.", call. = FALSE)
  }

  # Explicit model/copula support matrix.  This prevents an unsupported model
  # from falling through to an undefined 'res' object.
  supported_plain <- c(
    "Gaussian", "LogGaussian", "Tukeyh", "Tukeyh2", "SinhAsinh",
    "SkewGaussian", "Gamma", "Weibull", "Binomial", "BinomialNeg"
  )
  continuous_copula_margins <- c(
    "Gaussian", "StudentT", "LogGaussian", "Gamma", "Weibull",
    "Beta", "Beta2", "Kumaraswamy", "Kumaraswamy2",
    "Logistic", "SkewLaplace"
  )
  gaussian_copula_margins <- unique(c(
    continuous_copula_margins, "Tukeyh", "Tukeyh2", "SinhAsinh"
  ))
  supported_copula <- list(
    Gaussian = gaussian_copula_margins,
    Clayton = continuous_copula_margins,
    SkewGaussian = continuous_copula_margins
  )
  if (is.null(copula)) {
    if (!(model %in% supported_plain)) {
      stop(sprintf("GeoSimcond does not currently support model = '%s' without a copula.", model),
           call. = FALSE)
    }
  } else {
    copula <- gsub("[[:blank:]]", "", as.character(copula)[1L])
    if (!(copula %in% names(supported_copula))) {
      stop(sprintf("GeoSimcond does not currently support copula = '%s'.", copula),
           call. = FALSE)
    }
    if (!(model %in% supported_copula[[copula]])) {
      stop(sprintf(
        paste0(
          "GeoSimcond implements validated copula conditional simulation only for continuous margins. ",
          "Model '%s' is not supported with copula='%s'. Supported continuous margins are: %s."
        ),
        model, copula, paste(supported_copula[[copula]], collapse = ", ")
      ), call. = FALSE)
    }
  }

  ## Local conditional simulation is implemented through local Gaussian
  ## substitution.  It is therefore available for Gaussian fields, monotone
  ## one-Gaussian transformations, and all continuous margins under a Gaussian
  ## copula.  Gibbs-based latent constructions still require global conditioning.
  if (isTRUE(local)) {
    if (is.null(neighb) && is.null(maxdist)) {
      stop("For local conditional simulation, specify neighb or maxdist.",
           call. = FALSE)
    }
    if (!is.null(neighb)) {
      if (!is.numeric(neighb) || length(neighb) != 1L || !is.finite(neighb) ||
          neighb < 1 || abs(neighb - round(neighb)) > sqrt(.Machine$double.eps)) {
        stop("neighb must be a positive integer for local conditional simulation.",
             call. = FALSE)
      }
      neighb <- as.integer(round(neighb))
    }
    if (!is.null(maxdist) &&
        (!is.numeric(maxdist) || length(maxdist) != 1L || !is.finite(maxdist) || maxdist < 0)) {
      stop("maxdist must be a finite non-negative scalar.", call. = FALSE)
    }

    local_plain <- c("Gaussian", "LogGaussian", "Tukeyh", "Tukeyh2", "SinhAsinh")
    local_supported <- if (is.null(copula)) {
      model %in% local_plain
    } else {
      identical(copula, "Gaussian")
    }
    if (!local_supported) {
      stop(
        paste0(
          "GeoSimcond local=TRUE is currently supported for Gaussian fields, ",
          "LogGaussian/Tukeyh/Tukeyh2/SinhAsinh transformations, and continuous ",
          "margins with copula='Gaussian'. Gibbs-based SkewGaussian, Gamma/Weibull, ",
          "direct Binomial/BinomialNeg models, Clayton and SkewGaussian copulas still require local=FALSE."
        ),
        call. = FALSE
      )
    }
  }

  # Observed coordinates matrix
  if (is.matrix(coordx) || is.data.frame(coordx)) {
    coord_obs <- as.matrix(coordx)
  } else {
    coord_obs <- data.frame(x = coordx)
    if (!is.null(coordy)) coord_obs$y <- coordy
    if (!is.null(coordz)) coord_obs$z <- coordz
    coord_obs <- as.matrix(coord_obs)
  }
  if (!is.numeric(coord_obs) || nrow(coord_obs) < 1L || any(!is.finite(coord_obs))) {
    stop("Observed coordinates must be a non-empty finite numeric matrix.", call. = FALSE)
  }
  if (ncol(loc) != ncol(coord_obs)) {
    stop("Observed and prediction coordinates must have the same dimension.", call. = FALSE)
  }
  data <- as.numeric(data)
  if (length(data) != nrow(coord_obs) || any(!is.finite(data))) {
    stop("data must contain one finite value per observed location.", call. = FALSE)
  }
  if (model == "LogGaussian" && any(data <= 0)) {
    stop("LogGaussian conditional simulation requires strictly positive data.", call. = FALSE)
  }
  if (model %in% c("Gamma", "Weibull") && any(data < 0)) {
    stop(sprintf("%s conditional simulation requires non-negative data.", model),
         call. = FALSE)
  }
  if (!is.null(copula) && model %in% c("Beta", "Beta2")) {
    lower_support <- as.numeric(param$min)
    upper_support <- as.numeric(param$max)
    if (length(lower_support) != 1L || length(upper_support) != 1L ||
        !is.finite(lower_support) || !is.finite(upper_support) ||
        lower_support >= upper_support) {
      stop(sprintf("%s requires finite param$min < param$max.", model), call. = FALSE)
    }
    if (any(data < lower_support | data > upper_support)) {
      stop(sprintf("%s observations must lie inside [param$min, param$max].", model), call. = FALSE)
    }
  }
  if (!is.null(copula) && model %in% c("Poisson", "BinomialNeg", "Binomial")) {
    if (any(data < 0) || any(abs(data - round(data)) > sqrt(.Machine$double.eps))) {
      stop(sprintf("%s observations must be non-negative integers.", model),
           call. = FALSE)
    }
  }
  if (!is.null(copula) && model == "Binomial") {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n < 1 || n != as.integer(n)) {
      stop("n must be a positive integer for the Binomial model.", call. = FALSE)
    }
    n <- as.integer(n)
    if (any(data > n)) stop("Binomial observations cannot exceed n.", call. = FALSE)
  }
  if (!is.null(copula) && model == "BinomialNeg") {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0) {
      stop("n must be strictly positive for the BinomialNeg model.", call. = FALSE)
    }
  }

  # Mean and covariate specification.  The observation and prediction
  # means are selected independently:
  #   observations: param$mean vector, or X %*% beta, or intercept mean;
  #   predictions:  Mloc, or Xloc %*% beta, or intercept mean.
  if (spacetime && is.null(time))
    stop("At least one temporal instant is needed for space-time conditional simulation\n")
  n_obs_mean <- length(data)
  n_loc_mean <- nrow(as.matrix(loc)) * if(spacetime) length(time) else 1L
  external_mean <- .GeoMean_external(param, n = n_obs_mean)
  if (!is.null(external_mean) && !is.null(X))
    stop("A site-specific observation mean and X cannot be used together.\n")
  if (!is.null(Mloc) && !is.null(Xloc))
    stop("Mloc and Xloc are mutually exclusive.\n")
  if (!is.null(external_mean) && is.null(Mloc))
    stop("Mloc is required when the observation mean is site-specific.\n")

  if (!is.null(X)) {
    X <- .GeoMean_design(X, n = n_obs_mean, name = "X")
    p_mean <- ncol(X)
    invisible(.GeoMean_beta(param, p_mean))
    if (is.null(Xloc) && is.null(Mloc))
      stop("Supply either Xloc or Mloc for prediction locations.\n")
  } else {
    if (!is.null(Xloc)) stop("Xloc requires an observation design matrix X.\n")
    p_mean <- 1L
    if (is.null(external_mean)) {
      if (length(setdiff(.GeoMean_supplied_names(param), "mean")))
        stop("Multiple mean coefficients require an explicit design matrix X.\n")
      invisible(.GeoMean_beta(param, 1L))
    }
  }
  if (!is.null(Xloc))
    Xloc <- .GeoMean_design(Xloc, n = n_loc_mean, p = p_mean, name = "Xloc")
  if (!is.null(Mloc)) {
    Mloc <- as.numeric(Mloc)
    if (length(Mloc) == 1L) Mloc <- rep(Mloc, n_loc_mean)
    if (length(Mloc) != n_loc_mean || any(!is.finite(Mloc)))
      stop("Mloc must contain one finite value per prediction point.\n")
  }
  # Do not anisotropize loc here.  Anisotropy is propagated to
  # GeoSim/GeoKrig/GeoKrigWeights/GeoKriglocWeights, which handle it consistently.
  loc_orig <- loc
  # Temporal coordinates to use
  coordt_use <- if (spacetime || bivariate) coordt else NULL

  # Helper for multiplicative models: eta(s) = log(mu(s)) = X(s)' beta.
  # If X/Xloc are supplied, they are used; otherwise a constant mean is assumed.
  # Mloc, when supplied, is interpreted on the eta = log(mu) scale.
  clip_probability <- function(p, context = "probability transform") {
    p <- as.numeric(p)
    if (any(!is.finite(p))) {
      stop(sprintf("%s produced non-finite probabilities.", context), call. = FALSE)
    }
    eps <- .Machine$double.eps
    pmin(pmax(p, eps), 1 - eps)
  }

  get_eta_obs_loc <- function(param, X, Xloc, Mloc, n_obs, n_loc) {
    external <- .GeoMean_external(param, n = n_obs)
    eta_obs <- .GeoMean_vector(param, X = X, n = n_obs, external = external)

    if (!is.null(Mloc)) {
      eta_loc <- as.numeric(Mloc)
      if (length(eta_loc) == 1L) eta_loc <- rep(eta_loc, n_loc)
      if (length(eta_loc) != n_loc)
        stop("Mloc must contain one value per prediction point.\n")
    } else if (!is.null(Xloc)) {
      Xloc_mat <- .GeoMean_design(Xloc, n = n_loc, name = "Xloc")
      beta <- .GeoMean_beta(param, ncol(Xloc_mat))
      eta_loc <- as.numeric(Xloc_mat %*% beta)
    } else {
      if (!is.null(X)) {
        stop(
          "Xloc (or Mloc) must be supplied for prediction when X is used at the observed locations.",
          call. = FALSE
        )
      }
      eta_loc <- rep(.GeoMean_beta(param, 1L)[1L], n_loc)
    }

    list(obs = eta_obs, loc = eta_loc)
  }

############################################################################
######################### not copula models ################################
############################################################################
res <- NULL
if (is.null(copula)) {
  if (model == "Gaussian") {
    res <- .GeoSimcondGaussian(data, corrmodel, nrep, method, L,
                    param,
                    coord_obs, loc, coordt_use, time, X, Xloc, Mloc, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
  }

  ################################################################################
  ############### monotone transformations of one Gaussian RF ####################
  ################################################################################
  if (model == "LogGaussian") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    vv <- param$sill

    # Multiplicative parametrization:
    # Y(s) = exp(eta(s)) * exp(sqrt(vv) Z(s) - vv/2), hence E[Y(s)] = exp(eta(s)).
    datanorm <- (log(data) - eta$obs + vv/2)/sqrt(vv)   # gaussian scale

    param_gauss <- .GeoMean_intercept_param(param, 0)
    param_gauss$sill <- 1

    res <- .GeoSimcondGaussian(datanorm, corrmodel, nrep, method, L, param_gauss,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)

    res <- lapply(res, function(r) exp(eta$loc) * exp(sqrt(vv) * r - vv/2)) # back-transformation
  }

  if (model == "Tukeyh") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    vv <- param$sill
    tail <- .GeoMarginalTukeyTail(param$tail, "tail")
    datanorm <- .GeoMarginalInvLambertH((data - eta$obs)/sqrt(vv), tail)
    param_latent <- .GeoMean_intercept_param(param, 0)
    param_latent$sill <- 1; param_latent$tail <- NULL
    res <- .GeoSimcondGaussian(datanorm, corrmodel, nrep, method, L, param_latent,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) eta$loc + sqrt(vv) * .GeoMarginalTukeyH(r, tail))
  }

  if (model == "Tukeyh2") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    vv <- param$sill
    tail1 <- .GeoMarginalTukeyTail(param$tail1, "tail1")
    tail2 <- .GeoMarginalTukeyTail(param$tail2, "tail2")
    datanorm <- .GeoMarginalInvTukeyH2((data - eta$obs)/sqrt(vv), tail1, tail2)
    param_latent <- .GeoMean_intercept_param(param, 0)
    param_latent$sill <- 1; param_latent$tail1 <- NULL; param_latent$tail2 <- NULL
    res <- .GeoSimcondGaussian(datanorm, corrmodel, nrep, method, L, param_latent,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) {
      eta$loc + sqrt(vv) * .GeoMarginalTukeyH2(r, tail1, tail2)
    })
  }

  if (model == "SinhAsinh") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    vv <- as.numeric(param$sill)
    skew <- as.numeric(param$skew)
    tail <- as.numeric(param$tail)
    if (length(vv) != 1L || !is.finite(vv) || vv <= 0)
      stop("SinhAsinh sill must be a finite positive scalar.", call. = FALSE)
    if (length(skew) != 1L || !is.finite(skew))
      stop("SinhAsinh skew must be a finite scalar.", call. = FALSE)
    if (length(tail) != 1L || !is.finite(tail) || tail <= 0)
      stop("SinhAsinh tail must be a finite positive scalar.", call. = FALSE)

    ## GeoSim and the SAS construction use
    ## Y = eta + sqrt(vv) * sinh((asinh(Z) + skew) / tail),
    ## whose exact inverse is
    ## Z = sinh(tail * asinh((Y-eta)/sqrt(vv)) - skew).
    datanorm <- sinh(tail * asinh((data - eta$obs)/sqrt(vv)) - skew)
    param_latent <- .GeoMean_intercept_param(param, 0)
    param_latent$sill <- 1; param_latent$skew <- NULL; param_latent$tail <- NULL
    res <- .GeoSimcondGaussian(datanorm, corrmodel, nrep, method, L, param_latent,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) {
      eta$loc + sqrt(vv) * sinh((asinh(r) + skew) / tail)
    })
  }
  
  ################################################################################
  ############### non monotone transformations of more than one Gaussian RF ######
  ################################################################################
  if (model == "SkewGaussian") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    param_skew <- .GeoMean_intercept_param(param, 0)
    res <- .GeoSimcondSkewGaussian(coord_obs, loc, data, param_skew, corrmodel, local = local,
                               neighb = neighb, L = L, nrep = nrep, method = method,
                               distance = distance, radius = radius, anisopars = anisopars,
                               parallel = parallel, ncores = ncores, progress = progress,
                               n_iter = n_iter,
                               mean_obs = eta$obs, mean_loc = eta$loc)
  }

  if (model == "Gamma") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    res <- .GeoSimcondGamma(coord_obs, loc, data, param, corrmodel, local = local,
                        neighb = neighb, L = L, nrep = nrep, method = method,
                        distance = distance, radius = radius, anisopars = anisopars,
                        parallel = parallel, ncores = ncores, progress = progress,
                        mean_obs = eta$obs, mean_loc = eta$loc)
  }
  
  if (model == "Weibull") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    param.weibull <- .GeoMean_intercept_param(param, 0)
    param.weibull$shape <- 2
    ck <- gamma((param$shape + 1)/param$shape)
    data.t <- (exp(-eta$obs) * data * ck)^param$shape
    res <- .GeoSimcondGamma(coord_obs, loc, data.t, param.weibull, corrmodel, local = local,
                        neighb = neighb, L = L, nrep = nrep, method = method,
                        distance = distance, radius = radius, anisopars = anisopars,
                        parallel = parallel, ncores = ncores, progress = progress,
                        mean_obs = rep(0, n_obs_mean), mean_loc = rep(0, n_loc_mean))
    res <- lapply(res, function(x) exp(eta$loc) * (x^(1/param$shape))/ck)
  }
  
  if (model == "Poisson")  { stop("Poisson not implemented...") }

  if (model %in% c("Binomial", "BinomialNeg")) {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    res <- .GeoSimcondCountRF(
      model = model, coord_obs = coord_obs, loc = loc, data = data,
      n = n, nloc = nloc, param = param, corrmodel = corrmodel,
      mean_obs = eta$obs, mean_loc = eta$loc,
      distance = distance, radius = radius, anisopars = anisopars,
      nrep = nrep, n_iter = n_iter, mcmc_thin = mcmc_thin,
      method = method, local = local, parallel = parallel, progress = progress
    )
  }

} else { 
###########################################################################
################## copula models ###########################################
###########################################################################

  ## Continuous margins are first mapped one-to-one to their copula uniforms.
  ## Gaussian and skew-Gaussian copulas then condition their latent fields;
  ## the Clayton-like copula conditions the Gaussian-square representation
  ## subject to the observed latent ratio.  Discrete margins are intentionally
  ## not routed through the old mid-PIT approximation.

  .GeoCopula_cond_margin <- function(model, param, data) {
    scalar <- function(x, name, lower = -Inf, upper = Inf,
                       lower_open = FALSE, upper_open = FALSE) {
      x <- suppressWarnings(as.numeric(x))
      bad <- length(x) != 1L || !is.finite(x)
      if (!bad) {
        bad <- if (lower_open) x <= lower else x < lower
        bad <- bad || if (upper_open) x >= upper else x > upper
      }
      if (bad) stop(sprintf("Invalid %s for model '%s'.", name, model), call. = FALSE)
      x
    }

    ## Only margins whose marginal scale is parameterized by sill require it.
    ## Positive/bounded margins such as Gamma, Weibull, Beta/Beta2 and
    ## Kumaraswamy/Kumaraswamy2 determine their marginal variance from their
    ## shape parameters, so GeoFit objects for those models legitimately do not
    ## contain param$sill.  Do not reject them before their PIT is constructed.
    sill_models <- c(
      "Gaussian", "StudentT", "LogGaussian", "Logistic", "SkewLaplace",
      "SinhAsinh", "Tukeyh", "Tukeyh2"
    )
    sill <- NULL
    if (model %in% sill_models) {
      sill <- scalar(param$sill, "sill", lower = 0, lower_open = TRUE)
    }
    nugget <- if (is.null(param$nugget)) 0 else
      scalar(param$nugget, "nugget", lower = 0, upper = 1, upper_open = TRUE)

    eta <- NULL
    eta_needed <- model %in% c(
      "Gaussian", "StudentT", "LogGaussian", "Gamma", "Weibull",
      "Beta2", "Kumaraswamy2", "Logistic", "SkewLaplace",
      "SinhAsinh", "Tukeyh", "Tukeyh2"
    )
    if (eta_needed) {
      eta <- get_eta_obs_loc(param, X, Xloc, Mloc, n_obs_mean, n_loc_mean)
    }

    if (model == "Gaussian") {
      u <- stats::pnorm(data, mean = eta$obs, sd = sqrt(sill))
      qloc <- function(p) stats::qnorm(p, mean = eta$loc, sd = sqrt(sill))
    } else if (model == "StudentT") {
      eps_df <- scalar(param$df, "df", lower = 0, lower_open = TRUE)
      df_t <- 1 / eps_df
      u <- stats::pt((data - eta$obs) / sqrt(sill), df = df_t)
      qloc <- function(p) eta$loc + sqrt(sill) * stats::qt(p, df = df_t)
    } else if (model == "LogGaussian") {
      if (any(data <= 0))
        stop("LogGaussian conditional simulation requires strictly positive data.", call. = FALSE)
      u <- stats::plnorm(data, meanlog = eta$obs - sill / 2, sdlog = sqrt(sill))
      qloc <- function(p) stats::qlnorm(p, meanlog = eta$loc - sill / 2, sdlog = sqrt(sill))
    } else if (model == "Gamma") {
      shape <- scalar(param$shape, "shape", lower = 0, lower_open = TRUE)
      if (any(data <= 0))
        stop("Gamma conditional simulation requires strictly positive data.", call. = FALSE)
      mu_obs <- exp(eta$obs); mu_loc <- exp(eta$loc)
      u <- stats::pgamma(data, shape = shape / 2, rate = shape / (2 * mu_obs))
      qloc <- function(p) stats::qgamma(p, shape = shape / 2, rate = shape / (2 * mu_loc))
    } else if (model == "Weibull") {
      shape <- scalar(param$shape, "shape", lower = 0, lower_open = TRUE)
      if (any(data < 0))
        stop("Weibull conditional simulation requires non-negative data.", call. = FALSE)
      mu_obs <- exp(eta$obs); mu_loc <- exp(eta$loc)
      denom <- gamma(1 + 1 / shape)
      u <- stats::pweibull(data, shape = shape, scale = mu_obs / denom)
      qloc <- function(p) stats::qweibull(p, shape = shape, scale = mu_loc / denom)
    } else if (model == "Logistic") {
      u <- stats::plogis(data, location = eta$obs, scale = sqrt(sill))
      qloc <- function(p) stats::qlogis(p, location = eta$loc, scale = sqrt(sill))
    } else if (model == "SinhAsinh") {
      tail <- scalar(param$tail, "tail", lower = 0, lower_open = TRUE)
      skew <- scalar(param$skew, "skew")
      z <- sinh(tail * asinh((data - eta$obs) / sqrt(sill)) - skew)
      u <- stats::pnorm(z)
      qloc <- function(p) {
        zq <- stats::qnorm(p)
        eta$loc + sqrt(sill) * sinh((asinh(zq) + skew) / tail)
      }
    } else if (model == "Tukeyh") {
      tail <- scalar(param$tail, "tail", lower = 0, upper = 0.5, upper_open = TRUE)
      x <- (data - eta$obs) / sqrt(sill)
      u <- stats::pnorm(.GeoMarginalInvLambertH(x, tail))
      qloc <- function(p) {
        eta$loc + sqrt(sill) * .GeoMarginalQTukeyH(p, tail)
      }
    } else if (model == "Tukeyh2") {
      tail1 <- scalar(param$tail1, "tail1", lower = 0, upper = 0.5, upper_open = TRUE)
      tail2 <- scalar(param$tail2, "tail2", lower = 0, upper = 0.5, upper_open = TRUE)
      x <- (data - eta$obs) / sqrt(sill)
      z <- numeric(length(x))
      pos <- x >= 0
      if (any(pos)) z[pos] <- .GeoMarginalInvLambertH(x[pos], tail1)
      if (any(!pos)) z[!pos] <- .GeoMarginalInvLambertH(x[!pos], tail2)
      u <- stats::pnorm(z)
      qloc <- function(p) {
        eta$loc + sqrt(sill) * .GeoMarginalQTukeyH2(p, tail1, tail2)
      }
    } else if (model == "SkewLaplace") {
      skew <- scalar(param$skew, "skew", lower = 0, upper = 1,
                     lower_open = TRUE, upper_open = TRUE)
      u <- .GeoMarginalPSkewLaplace((data - eta$obs) / sqrt(sill), skew)
      qloc <- function(p) eta$loc + sqrt(sill) * .GeoMarginalQSkewLaplace(p, skew)
    } else if (model %in% c("Beta", "Beta2", "Kumaraswamy", "Kumaraswamy2")) {
      pmin0 <- scalar(param$min, "min")
      pmax0 <- scalar(param$max, "max")
      if (pmin0 >= pmax0)
        stop(sprintf("%s requires min < max.", model), call. = FALSE)
      if (any(data < pmin0 | data > pmax0))
        stop(sprintf("%s observations must lie inside [min, max].", model), call. = FALSE)
      y <- (data - pmin0) / (pmax0 - pmin0)

      if (model == "Beta") {
        a <- scalar(param$shape1, "shape1", lower = 0, lower_open = TRUE)
        b <- scalar(param$shape2, "shape2", lower = 0, lower_open = TRUE)
        u <- stats::pbeta(y, shape1 = a, shape2 = b)
        qloc <- function(p) pmin0 + (pmax0 - pmin0) * stats::qbeta(p, shape1 = a, shape2 = b)
      } else if (model == "Beta2") {
        shape <- scalar(param$shape, "shape", lower = 0, lower_open = TRUE)
        mu_obs <- stats::plogis(eta$obs); mu_loc <- stats::plogis(eta$loc)
        u <- stats::pbeta(y, shape1 = mu_obs * shape, shape2 = (1 - mu_obs) * shape)
        qloc <- function(p) pmin0 + (pmax0 - pmin0) *
          stats::qbeta(p, shape1 = mu_loc * shape, shape2 = (1 - mu_loc) * shape)
      } else if (model == "Kumaraswamy") {
        a <- scalar(param$shape1, "shape1", lower = 0, lower_open = TRUE)
        b <- scalar(param$shape2, "shape2", lower = 0, lower_open = TRUE)
        pkuma <- function(x, aa, bb) 1 - (1 - x^aa)^bb
        qkuma <- function(p, aa, bb) (1 - (1 - p)^(1 / bb))^(1 / aa)
        u <- pkuma(y, a, b)
        qloc <- function(p) pmin0 + (pmax0 - pmin0) * qkuma(p, a, b)
      } else {
        a <- scalar(param$shape, "shape", lower = 0, lower_open = TRUE)
        mu_obs <- stats::plogis(eta$obs); mu_loc <- stats::plogis(eta$loc)
        b_obs <- log(0.5) / log1p(-mu_obs^a)
        b_loc <- log(0.5) / log1p(-mu_loc^a)
        if (any(!is.finite(c(b_obs, b_loc))) || any(c(b_obs, b_loc) <= 0))
          stop("Kumaraswamy2 produced invalid location-specific shape parameters.", call. = FALSE)
        pkuma <- function(x, aa, bb) 1 - (1 - x^aa)^bb
        qkuma <- function(p, aa, bb) (1 - (1 - p)^(1 / bb))^(1 / aa)
        u <- pkuma(y, a, b_obs)
        qloc <- function(p) pmin0 + (pmax0 - pmin0) * qkuma(p, a, b_loc)
      }
    } else {
      stop(sprintf("Internal error: unsupported continuous copula margin '%s'.", model),
           call. = FALSE)
    }

    list(
      u_obs = clip_probability(u, "marginal PIT"),
      q_loc = function(p) as.numeric(qloc(clip_probability(p, "predictive PIT"))),
      nugget = nugget
    )
  }

  .GeoCopula_latent_param <- function(param, corrmodel, nugget, skew = NULL) {
    corr_names <- CorrParam(corrmodel)
    out <- as.list(param[corr_names])
    out$mean <- 0
    out$nugget <- nugget
    out$sill <- 1
    if (!is.null(skew)) out$skew <- as.numeric(skew)
    out
  }

  margin <- .GeoCopula_cond_margin(model, param, data)

  if (copula == "Gaussian") {
    datanorm <- stats::qnorm(margin$u_obs)
    param_latent <- .GeoCopula_latent_param(param, corrmodel, margin$nugget)

    res <- .GeoSimcondGaussian(
      datanorm, corrmodel, nrep, method, L, param_latent,
      coord_obs, loc, coordt_use, time, NULL, NULL, NULL,
      distance, radius, anisopars, local, neighb, maxdist, maxtime,
      space, spacetime, bivariate, parallel, ncores, progress
    )
    res <- lapply(res, function(z) margin$q_loc(stats::pnorm(z)))
  }

  if (copula == "Clayton") {
    if (!space) {
      stop("Clayton copula conditional simulation is currently implemented only for purely spatial data.",
           call. = FALSE)
    }
    if (is.null(names(param)) || !("nu" %in% names(param))) {
      stop("For copula='Clayton', param$nu is required.", call. = FALSE)
    }
    nu_clayton <- suppressWarnings(as.numeric(param[["nu", exact = TRUE]]))
    if (length(nu_clayton) != 1L || !is.finite(nu_clayton) || nu_clayton <= 0 ||
        abs(nu_clayton - round(nu_clayton)) > sqrt(.Machine$double.eps)) {
      stop("For copula='Clayton', param$nu must be a positive integer.", call. = FALSE)
    }
    nu_clayton <- as.integer(round(nu_clayton))
    param_latent <- .GeoCopula_latent_param(param, corrmodel, margin$nugget)

    u_res <- .GeoSimcondClaytonCopula(
      coord_obs = coord_obs, loc = loc, u_obs = margin$u_obs,
      nu = nu_clayton, latent_param = param_latent, corrmodel = corrmodel,
      nrep = nrep, method = method, L = L,
      distance = distance, radius = radius, anisopars = anisopars,
      parallel = parallel, ncores = ncores, progress = progress,
      n_iter = n_iter
    )
    res <- lapply(u_res, margin$q_loc)
  }

  if (copula == "SkewGaussian") {
    if (!space) {
      stop("Skew-Gaussian copula conditional simulation is currently implemented only for purely spatial data.",
           call. = FALSE)
    }
    if (is.null(names(param)) || !("nu" %in% names(param))) {
      stop("For copula='SkewGaussian', param$nu (eta) is required.", call. = FALSE)
    }
    eta_cop <- suppressWarnings(as.numeric(param[["nu", exact = TRUE]]))
    if (length(eta_cop) != 1L || !is.finite(eta_cop) || abs(eta_cop) >= 1) {
      stop("For copula='SkewGaussian', param$nu is eta and must lie strictly between -1 and 1.",
           call. = FALSE)
    }
    gamma_eta <- eta_cop / sqrt(1 - eta_cop^2)
    omega_eta <- sqrt(1 + gamma_eta^2)

    datanorm <- sn::qsn(
      margin$u_obs, xi = 0, omega = omega_eta, alpha = gamma_eta
    )
    param_latent <- .GeoCopula_latent_param(
      param, corrmodel, margin$nugget, skew = gamma_eta
    )

    res <- .GeoSimcondSkewGaussian(
      coord_obs, loc, datanorm, param_latent, corrmodel,
      local = local, neighb = neighb, L = L, nrep = nrep,
      method = method, distance = distance, radius = radius,
      anisopars = anisopars, parallel = parallel, ncores = ncores,
      progress = progress, n_iter = n_iter
    )
    res <- lapply(res, function(z) {
      u <- sn::psn(z, xi = 0, omega = omega_eta, alpha = gamma_eta)
      margin$q_loc(u)
    })
  }

###########################################################################
}


if (is.null(res) || !is.list(res) || length(res) < 1L) {
  stop("GeoSimcond did not produce any conditional simulations for the requested model.",
       call. = FALSE)
}
expected_nloc <- nrow(loc)
valid_lengths <- vapply(res, length, integer(1L))
if (any(valid_lengths != expected_nloc)) {
  stop("A conditional simulation has an unexpected number of prediction values.",
       call. = FALSE)
}
if (any(!vapply(res, function(z) all(is.finite(z)), logical(1L)))) {
  stop("At least one conditional simulation contains non-finite values.",
       call. = FALSE)
}

## Summarize without materializing a second nrep x nloc dense matrix.  The
## conditional simulations are already retained in res, so use a streaming
## Welford update for the mean and sample variance.
if (length(res) == 1L) {
  cond_mean <- as.numeric(res[[1L]])
  cond_var <- rep(NA_real_, length(cond_mean))
} else {
  cond_mean <- numeric(expected_nloc)
  cond_M2 <- numeric(expected_nloc)
  k <- 0L
  for (z in res) {
    k <- k + 1L
    z <- as.numeric(z)
    delta <- z - cond_mean
    cond_mean <- cond_mean + delta / k
    cond_M2 <- cond_M2 + delta * (z - cond_mean)
  }
  cond_var <- cond_M2 / (k - 1L)
}
############################################################################
############################################################################
############################ OUTPUT ########################################
GeoSimcond_out <- list(
    bivariate = bivariate,
    coordx = coordx,
    coordy = coordy,
    coordz = coordz,
    coordt = coordt,
    corrmodel = corrmodel,
    data = data,
    distance = distance,
    grid = grid,
    loc = loc,
    copula = copula,
    numcoord = nrow(coord_obs),
    numloc = nrow(loc),
    numtime = length(coordt),
    numt = length(time),
    maxdist = maxdist,
    maxtime = maxtime,
    model = model,
    n = n,
    nloc = nloc,
    param = param,
    method = method,
    local = local,
    neighb = neighb,
    L = L,
    mcmc_thin = mcmc_thin,
    condsim = res,            # list of nrep vectors, each of length numloc
    cond_mean = cond_mean,    # conditional mean across replications
    cond_var = cond_var,      # conditional variance across replications
    radius = radius,
    spacetime = spacetime,
    time = time
  )
  structure(c(GeoSimcond_out, call = call), class = c("GeoSimcond"))
}
