####################################################
### Internal simulation validation helpers
####################################################

.GeoSimulationCanonicalModel <- function(model) {
  if(length(model) != 1L || is.na(model))
    stop("model must be a single non-missing character string.", call. = FALSE)
  model <- gsub("[[:blank:]]", "", as.character(model))
  aliases <- c(
    Gauss = "Gaussian",
    SkewGauss = "SkewGaussian",
    LogGauss = "LogGaussian",
    TwoPieceGauss = "TwoPieceGaussian",
    Binary = "Binomial",
    Bernoulli = "Binomial",
    Geom = "BinomialNeg",
    Geometric = "BinomialNeg",
    tukeyh = "Tukeyh",
    tukeyh2 = "Tukeyh2",
    gamma = "Gamma",
    weibull = "Weibull",
    Loglogistic = "LogLogistic",
    Binomiallogistic = "BinomialLogistic",
    poisson = "Poisson",
    poissongamma = "PoissonGamma"
  )
  if(model %in% names(aliases)) unname(aliases[[model]]) else model
}

.GeoSimulationEffectiveN <- function(model_requested, model, n) {
  requested <- gsub("[[:blank:]]", "", as.character(model_requested)[1L])
  if(requested %in% c("Binary", "Bernoulli")) return(1L)
  if(requested %in% c("Geom", "Geometric")) {
    nn <- .GeoValidateIntegerVector(n, "n", nobs = NULL, minimum = 1L)
    if(length(nn) != 1L || nn[1L] != 1L) {
      stop("For Geometric simulation, n is fixed to 1 (the r = 1 Negative-Binomial case).",
           call. = FALSE)
    }
    return(1L)
  }
  n
}

.GeoValidateNrep <- function(nrep) {
  if(!is.numeric(nrep) || length(nrep) != 1L || !is.finite(nrep) ||
     nrep < 1 || abs(nrep - round(nrep)) > sqrt(.Machine$double.eps)) {
    stop("nrep must be a positive integer.", call. = FALSE)
  }
  as.integer(round(nrep))
}

.GeoValidateIntegerVector <- function(x, name, nobs = NULL, minimum = 1L) {
  x <- as.numeric(x)
  tol <- sqrt(.Machine$double.eps)
  if(!length(x) || any(!is.finite(x)) || any(x < minimum) ||
     any(abs(x - round(x)) > tol)) {
    stop(sprintf("%s must contain integer values greater than or equal to %d.",
                 name, as.integer(minimum)), call. = FALSE)
  }
  if(!is.null(nobs) && !(length(x) %in% c(1L, as.integer(nobs)))) {
    stop(sprintf("%s must have length 1 or the number of simulated observations (%d).",
                 name, as.integer(nobs)), call. = FALSE)
  }
  as.integer(round(x))
}

.GeoValidateLatentIntegerParameters <- function(model, param) {
  tol <- sqrt(.Machine$double.eps)
  pos_integer <- function(x, name, minimum = 1L) {
    x <- as.numeric(x)
    if(length(x) != 1L || !is.finite(x) || x < minimum ||
       abs(x - round(x)) > tol) {
      stop(sprintf("For direct %s simulation, %s must be an integer >= %d because it determines the number of latent Gaussian fields.",
                   model, name, as.integer(minimum)), call. = FALSE)
    }
    as.integer(round(x))
  }

  if(model == "Gamma") {
    param$shape <- pos_integer(param$shape, "param$shape")
  }
  if(model == "Beta") {
    param$shape1 <- pos_integer(param$shape1, "param$shape1")
    param$shape2 <- pos_integer(param$shape2, "param$shape2")
  }
  if(model %in% c("StudentT", "SkewStudentT", "TwoPieceStudentT")) {
    eps_df <- as.numeric(param$df)
    if(length(eps_df) != 1L || !is.finite(eps_df) || eps_df <= 0)
      stop(sprintf("For direct %s simulation, param$df must be the positive reciprocal degrees-of-freedom parameter.", model), call. = FALSE)
    nu <- 1 / eps_df
    if(!is.finite(nu) || nu < 3 || abs(nu - round(nu)) > tol) {
      stop(sprintf("For direct %s simulation, 1/param$df must be an integer >= 3 because the model is constructed from that many latent Gaussian fields.", model), call. = FALSE)
    }
  }
  if(model == "TwoPieceBimodal") {
    param$df <- pos_integer(param$df, "param$df")
  }
  if(model %in% c("PoissonGamma", "PoissonGammaZIP")) {
    shape <- as.numeric(param$shape)
    m <- 2 * shape
    if(length(shape) != 1L || !is.finite(shape) || shape <= 0 ||
       abs(m - round(m)) > tol) {
      stop(sprintf("For direct %s simulation, 2*param$shape must be a positive integer because the Gamma mixing field is generated from squared latent Gaussian fields.", model), call. = FALSE)
    }
  }
  param
}

.GeoValidateSimulationMarginal <- function(model, param) {
  scalar_positive <- function(x, name, strict = TRUE) {
    x <- as.numeric(x)
    if(length(x) != 1L || !is.finite(x) || (strict && x <= 0) || (!strict && x < 0))
      stop(sprintf("%s must be a finite %s scalar.", name,
                   if(strict) "positive" else "non-negative"), call. = FALSE)
    x
  }
  scalar_interval <- function(x, name, lo, hi) {
    x <- as.numeric(x)
    if(length(x) != 1L || !is.finite(x) || x <= lo || x >= hi)
      stop(sprintf("%s must lie strictly between %g and %g.", name, lo, hi), call. = FALSE)
    x
  }

  if(!is.null(param$sill)) scalar_positive(param$sill, "param$sill")
  if(!is.null(param$nugget)) {
    nug <- as.numeric(param$nugget)
    if(length(nug) != 1L || !is.finite(nug) || nug < 0 || nug >= 1)
      stop("param$nugget must be a finite scalar in [0,1).", call. = FALSE)
  }
  if(model %in% c("Weibull", "Kumaraswamy2"))
    scalar_positive(param$shape, "param$shape")
  if(model == "LogLogistic") {
    sh <- as.numeric(param$shape)
    if(length(sh) != 1L || !is.finite(sh) || sh <= 1)
      stop("param$shape must be finite and greater than 1 for LogLogistic simulation under the mean parametrization.", call. = FALSE)
  }
  if(model == "TwoPieceBimodal") {
    sh <- as.numeric(param$shape)
    if(length(sh) != 1L || !is.finite(sh) || sh < 0)
      stop("param$shape must be a finite non-negative scalar for TwoPieceBimodal simulation.", call. = FALSE)
  }
  if(model == "Kumaraswamy") {
    scalar_positive(param$shape1, "param$shape1")
    scalar_positive(param$shape2, "param$shape2")
  }
  if(model %in% c("PoissonZIP", "PoissonGammaZIP", "BinomialNegZINB")) {
    pmu <- as.numeric(param$pmu)
    if(length(pmu) != 1L || !is.finite(pmu))
      stop("param$pmu must be a finite scalar for zero-inflated simulation.", call. = FALSE)
    for(nm in c("nugget1", "nugget2")) {
      ng <- as.numeric(param[[nm]])
      if(length(ng) != 1L || !is.finite(ng) || ng < 0 || ng >= 1)
        stop(sprintf("param$%s must be a finite scalar in [0,1) for zero-inflated simulation.", nm), call. = FALSE)
    }
  }
  if(model %in% c("Beta", "Kumaraswamy", "Kumaraswamy2")) {
    lo <- as.numeric(param$min); hi <- as.numeric(param$max)
    if(length(lo) != 1L || length(hi) != 1L || !is.finite(lo) ||
       !is.finite(hi) || lo >= hi)
      stop(sprintf("Direct %s simulation requires finite param$min < param$max.", model), call. = FALSE)
  }
  if(model == "SinhAsinh") scalar_positive(param$tail, "param$tail")
  if(model == "SkewLaplace") scalar_interval(param$skew, "param$skew", 0, 1)
  if(model %in% c("SkewStudentT", "TwoPieceGaussian", "TwoPieceStudentT",
                 "TwoPieceTukeyh", "TwoPieceBimodal"))
    scalar_interval(param$skew, "param$skew", -1, 1)
  if(model == "Tukeyh") {
    tail <- as.numeric(param$tail)
    if(length(tail) != 1L || !is.finite(tail) || tail < 0 || tail >= 0.5)
      stop("param$tail must lie in [0, 0.5) for Tukeyh simulation.", call. = FALSE)
  }
  if(model == "TwoPieceTukeyh") {
    tail <- as.numeric(param$tail)
    if(length(tail) != 1L || !is.finite(tail) || tail < 0)
      stop("param$tail must be a finite non-negative scalar.", call. = FALSE)
  }
  if(model == "Tukeyh2") {
    for(nm in c("tail1", "tail2")) {
      x <- as.numeric(param[[nm]])
      if(length(x) != 1L || !is.finite(x) || x < 0 || x >= 0.5)
        stop(sprintf("param$%s must lie in [0, 0.5) for Tukeyh2 simulation.", nm), call. = FALSE)
    }
  }
  if(model == "Tukeygh") {
    tail <- as.numeric(param$tail)
    skew <- as.numeric(param$skew)
    if(length(tail) != 1L || !is.finite(tail) || tail < 0)
      stop("param$tail must be a finite non-negative scalar.", call. = FALSE)
    if(length(skew) != 1L || !is.finite(skew))
      stop("param$skew must be a finite scalar.", call. = FALSE)
  }
  param
}

.GeoDirectSimulationModels <- function() {
  c(
    "Gaussian", "Tukeygh", "SkewGaussian", "Binomial", "StudentT",
    "Wrapped", "BinomialNeg", "SkewStudentT", "SinhAsinh", "Gamma",
    "LogGaussian", "LogLogistic", "Logistic", "SkewLaplace", "Weibull",
    "TwoPieceStudentT", "Beta", "TwoPieceGaussian", "Poisson",
    "Kumaraswamy", "Tukeyh", "TwoPieceTukeyh", "TwoPieceBimodal",
    "Tukeyh2", "Kumaraswamy2", "PoissonZIP", "BinomialNegZINB",
    "PoissonGamma", "BinomialLogistic", "PoissonGammaZIP"
  )
}

.GeoValidateDirectSimulationModel <- function(model, simulator = "GeoSim") {
  if(model == "Beta2") {
    stop(sprintf("Direct Beta2 simulation is not implemented in %s(); use GeoSimCopula(..., model = 'Beta2') with an explicit copula.", simulator),
         call. = FALSE)
  }
  if(grepl("_misp_", model, fixed = TRUE) || startsWith(model, "Gaussian_misp_") ||
     startsWith(model, "Binary_misp_")) {
    stop(sprintf("%s is an inferential misspecification model and is not a data-generating model in %s(). Simulate from the corresponding non-misspecified model instead.",
                 model, simulator), call. = FALSE)
  }
  if(!(model %in% .GeoDirectSimulationModels())) {
    stop(sprintf("%s does not implement direct simulation for model = '%s'. Use a supported direct model or GeoSimCopula() when a copula construction is intended.",
                 simulator, model), call. = FALSE)
  }
  invisible(TRUE)
}

.GeoTBCorrelationModels <- function() {
  c("Matern", "GenWend", "Genwend", "Hypergeometric", "HyperGeometric",
    "hypergeometric", "GenWend_Matern", "Genwend_Matern",
    "Hypergeometric_Matern", "HyperGeometric_Matern", "hypergeometric_Matern",
    "Kummer", "Kummer_Matern", "Kummer_matern")
}

.GeoValidateTBCorrelationModel <- function(corrmodel, bivariate = FALSE) {
  if(isTRUE(bivariate)) {
    if(!(corrmodel %in% c("Bi_matern", "Bi_Matern")))
      stop(sprintf("GeoSimapprox(method='TB') does not implement bivariate correlation model '%s'.", corrmodel), call. = FALSE)
  } else if(!(corrmodel %in% .GeoTBCorrelationModels())) {
    stop(sprintf("GeoSimapprox(method='TB') does not implement correlation model '%s'. Use GeoSim() or a supported TB correlation family.", corrmodel), call. = FALSE)
  }
  invisible(TRUE)
}

.GeoValidateCopulaMarginal <- function(model, param) {
  ## Copula margins use quantile transforms; integer restrictions from the
  ## direct latent-field constructions do not apply here.
  if(!is.null(param$sill)) {
    ss <- as.numeric(param$sill)
    if(length(ss) != 1L || !is.finite(ss) || ss <= 0)
      stop("param$sill must be a finite positive scalar.", call. = FALSE)
  }
  if(model == "SinhAsinh") {
    tail <- as.numeric(param$tail); skew <- as.numeric(param$skew)
    if(length(tail) != 1L || !is.finite(tail) || tail <= 0)
      stop("SinhAsinh tail must be a finite positive scalar.", call. = FALSE)
    if(length(skew) != 1L || !is.finite(skew))
      stop("SinhAsinh skew must be a finite scalar.", call. = FALSE)
  }
  if(model == "Tukeyh") {
    tail <- as.numeric(param$tail)
    if(length(tail) != 1L || !is.finite(tail) || tail < 0 || tail >= 0.5)
      stop("Tukeyh tail must lie in [0, 0.5).", call. = FALSE)
  }
  if(model == "TwoPieceTukeyh") {
    tail <- as.numeric(param$tail)
    if(length(tail) != 1L || !is.finite(tail) || tail < 0)
      stop("TwoPieceTukeyh tail must be a finite non-negative scalar.", call. = FALSE)
  }
  if(model == "Tukeyh2") {
    vals <- c(as.numeric(param$tail1), as.numeric(param$tail2))
    if(length(vals) != 2L || any(!is.finite(vals)) || any(vals < 0) || any(vals >= 0.5))
      stop("Tukeyh2 tail1 and tail2 must lie in [0, 0.5).", call. = FALSE)
  }
  if(model == "Tukeygh") {
    vals <- c(as.numeric(param$tail), as.numeric(param$skew))
    if(length(vals) != 2L || any(!is.finite(vals)) || vals[1L] < 0)
      stop("Tukeygh requires finite skew and a finite non-negative tail.", call. = FALSE)
  }
  if(model %in% c("TwoPieceGaussian", "TwoPieceStudentT", "TwoPieceTukeyh")) {
    sk <- as.numeric(param$skew)
    if(length(sk) != 1L || !is.finite(sk) || abs(sk) >= 1)
      stop(sprintf("%s skew must lie strictly between -1 and 1.", model), call. = FALSE)
  }
  if(model == "SkewGaussian") {
    sk <- as.numeric(param$skew)
    if(length(sk) != 1L || !is.finite(sk))
      stop("SkewGaussian skew must be a finite scalar.", call. = FALSE)
  }
  invisible(TRUE)
}

.GeoCopulaSimulationModels <- function() {
  c(
    "Gaussian", "Binomial", "BinomialNeg", "BinomialNegZINB", "Poisson",
    "PoissonZIP", "Gamma", "Weibull", "LogLogistic", "LogGaussian",
    "Logistic", "StudentT", "SkewGaussian", "SkewStudentT", "SinhAsinh",
    "Tukeyh", "Tukeyh2", "SkewLaplace", "Tukeygh", "TwoPieceGaussian",
    "TwoPieceStudentT", "TwoPieceTukeyh", "Kumaraswamy", "Kumaraswamy2",
    "Beta", "Beta2"
  )
}

.GeoValidateCopulaSimulationModel <- function(model) {
  if(!(model %in% .GeoCopulaSimulationModels())) {
    stop(sprintf("GeoSimCopula does not implement marginal model '%s'.", model),
         call. = FALSE)
  }
  invisible(TRUE)
}

.GeoSimulationCanonicalCorrModel <- function(corrmodel) {
  corrmodel <- gsub("[[:blank:]]", "", as.character(corrmodel)[1L])
  aliases <- c(
    matern = "Matern",
    Genwend = "GenWend",
    HyperGeometric = "Hypergeometric",
    hypergeometric = "Hypergeometric",
    Genwend_Matern = "GenWend_Matern",
    HyperGeometric_Matern = "Hypergeometric_Matern",
    hypergeometric_Matern = "Hypergeometric_Matern",
    Kummer_matern = "Kummer_Matern",
    Bi_Matern = "Bi_matern",
    matern_matern = "Matern_Matern",
    Genwend_Genwend = "GenWend_GenWend"
  )
  if(corrmodel %in% names(aliases)) unname(aliases[[corrmodel]]) else corrmodel
}

.GeoValidateExactSimulationMethod <- function(method) {
  method <- gsub("[[:blank:]]", "", as.character(method)[1L])
  if(!(method %in% c("cholesky", "svd")))
    stop("GeoSim method must be either 'cholesky' or 'svd'.", call. = FALSE)
  method
}

.GeoValidateCopulaSimulationMethod <- function(method) {
  method <- gsub("[[:blank:]]", "", as.character(method)[1L])
  key <- tolower(method)
  if (key == "tb") return("TB")
  if (key %in% c("cholesky", "svd")) return(key)
  stop("GeoSimCopula method must be one of 'cholesky', 'svd', or 'TB'.",
       call. = FALSE)
}

.GeoValidateTBDistance <- function(distance, context = "GeoSimapprox") {
  key <- tolower(gsub("[[:blank:]]", "", as.character(distance)[1L]))
  if (!identical(key, "eucl")) {
    stop(sprintf("%s(method='TB') requires distance = 'Eucl' because turning-bands simulation operates on Cartesian coordinates.",
                 context), call. = FALSE)
  }
  invisible(TRUE)
}

.GeoValidateTBControls <- function(L) {
  if(!is.numeric(L) || length(L) != 1L || !is.finite(L) || L < 1 ||
     abs(L - round(L)) > sqrt(.Machine$double.eps))
    stop("L must be a positive integer for turning-bands simulation.", call. = FALSE)
  as.integer(round(L))
}

.GeoTBParallelMode <- function(dime, L, nrep, parallel = FALSE,
                               ncores = NULL, force_inner = FALSE) {
  ## Empirical TB scheduling rule calibrated on the two large-data regimes
  ## exercised by GeoSimapprox/GeoSimcond.  The decision deliberately uses
  ## both the work in one field (N * L) and the total work across replicates
  ## (N * L * nrep): the former determines whether spatial chunking can pay
  ## for worker overhead, while nrep determines when distributing complete
  ## independent fields is preferable.
  if (!isTRUE(parallel)) return("serial")

  if (!is.null(ncores)) {
    nc <- suppressWarnings(as.numeric(ncores)[1L])
    if (is.finite(nc) && nc <= 1) return("serial")
  }

  if (isTRUE(force_inner)) return("inner")

  N <- suppressWarnings(as.numeric(dime)[1L])
  LL <- suppressWarnings(as.numeric(L)[1L])
  R <- suppressWarnings(as.numeric(nrep)[1L])

  if (!is.finite(N) || !is.finite(LL) || !is.finite(R) ||
      N <= 0 || LL <= 0 || R <= 0)
    return("serial")

  get_positive_option <- function(name, default) {
    value <- suppressWarnings(as.numeric(getOption(name, default))[1L])
    if (!is.finite(value) || value <= 0) default else value
  }

  min_one_work <- get_positive_option(
    "GeoModels.tb_parallel_min_one_work", 1e7
  )
  min_total_work <- get_positive_option(
    "GeoModels.tb_parallel_min_total_work", 5e8
  )
  large_N <- get_positive_option(
    "GeoModels.tb_parallel_large_N", 1e5
  )
  heavy_one_work <- get_positive_option(
    "GeoModels.tb_parallel_heavy_one_work", 2e8
  )
  outer_nrep_large <- get_positive_option(
    "GeoModels.tb_parallel_outer_nrep_large", 100
  )
  outer_nrep_small <- get_positive_option(
    "GeoModels.tb_parallel_outer_nrep_small", 200
  )

  work_one <- N * LL
  work_all <- work_one * R

  ## For small individual fields, or too little aggregate work, starting a
  ## multisession plan costs more than it saves.
  if (work_one < min_one_work || work_all < min_total_work)
    return("serial")

  ## Large spatial payloads and very expensive individual fields reach the
  ## outer-replicate crossover sooner.  Otherwise keep using the inner TB
  ## path until there are enough independent replicates to amortize export
  ## and worker-start costs.
  outer_nrep <- if (N >= large_N || work_one >= heavy_one_work)
    outer_nrep_large else outer_nrep_small

  if (R > 1 && R >= outer_nrep) return("outer")

  "inner"
}
.GeoSimulationLatentSpec <- function(model, param, n) {
  k <- 1L
  npoi <- 1
  if(model %in% c("SkewGaussian", "LogGaussian", "TwoPieceGaussian", "TwoPieceTukeyh")) k <- 1L
  if(model %in% c("Weibull", "Wrapped")) k <- 2L
  if(model %in% c("LogLogistic", "Logistic", "SkewLaplace")) k <- 4L
  if(model == "Binomial") k <- max(n)
  if(model == "BinomialLogistic") k <- 2L * max(n)
  if(model %in% c("BinomialNeg", "BinomialNegZINB")) k <- 99999L
  if(model %in% c("Poisson", "PoissonZIP")) {
    k <- 2L
    npoi <- 999999999
  }
  if(model %in% c("PoissonGamma", "PoissonGammaZIP")) {
    k <- 2L + as.integer(round(2 * param$shape))
    npoi <- 999999999
  }
  if(model == "PoissonWeibull") {
    k <- 4L
    npoi <- 999999999
  }
  if(model %in% c("PoissonZIP", "BinomialNegZINB", "PoissonGammaZIP")) {
    param$nugget <- param$nugget1
  }
  if(model == "Gamma") k <- as.integer(param$shape)
  if(model == "Beta") k <- as.integer(param$shape1) + as.integer(param$shape2)
  if(model %in% c("Kumaraswamy", "Kumaraswamy2")) k <- 4L
  if(model == "StudentT") k <- as.integer(round(1 / param$df)) + 1L
  if(model == "TwoPieceBimodal") k <- as.integer(param$df) + 1L
  if(model %in% c("SkewStudentT", "TwoPieceStudentT")) {
    k <- as.integer(round(1 / param$df)) + 2L
  }
  list(k = as.integer(k), npoi = npoi, param = param)
}

.GeoSimulationSquareSums <- function(dd, index, component = 1L) {
  index <- as.integer(index)
  if(!length(index)) return(numeric(dim(dd)[1L]))
  z <- matrix(dd[, component, index, drop = FALSE], nrow = dim(dd)[1L])
  rowSums(z * z)
}

.GeoSimulationBinomialCounts <- function(dd, n) {
  dime <- dim(dd)[1L]
  k <- dim(dd)[3L]
  nn <- if(length(n) == 1L) rep.int(as.integer(n), dime) else as.integer(n)
  out <- numeric(dime)
  for(i in seq_len(k)) out <- out + dd[, 1L, i] * (nn >= i)
  out
}

.GeoSimulationBinomialLogisticCounts <- function(dd, n, mean) {
  dime <- dim(dd)[1L]
  k <- dim(dd)[3L]
  half_k <- k %/% 2L
  nn <- if(length(n) == 1L) rep.int(as.integer(n), dime) else as.integer(n)
  out <- numeric(dime)
  for(j in seq_len(half_k)) {
    i <- 2L * j - 1L
    ee <- 0.5 * (dd[, 1L, i]^2 + dd[, 1L, i + 1L]^2)
    out <- out + as.numeric(mean + log(exp(ee) - 1) > 0) * (nn >= j)
  }
  out
}

.GeoSimulationFormat <- function(sim, grid, space, numxgrid = NULL, numygrid = NULL,
                                 numtime = 1L, numcoord = length(sim), byrow = TRUE) {
  if(!grid) {
    if(space) return(c(sim))
    return(matrix(sim, nrow = numtime, ncol = numcoord, byrow = byrow))
  }
  if(space) return(array(sim, c(numxgrid, numygrid)))
  array(sim, c(numxgrid, numygrid, numtime))
}

.GeoSimulationSplitDynamic <- function(sim, ns) {
  ns <- as.integer(ns)
  ends <- cumsum(ns)
  starts <- c(1L, head(ends, -1L) + 1L)
  values <- c(sim)
  Map(function(first, last) values[first:last], starts, ends)
}

.GeoSimulationTransformUnivariate <- function(model, sim, dd = NULL, simDD = NULL,
                                               mm = 0, vv = 1, sk = 0, tl = 0,
                                               t1l = 0, t2l = 0, bimo = NULL,
                                               param = NULL, k = 1L,
                                               spacetime = FALSE) {
  byrow <- TRUE
  needs_format <- TRUE
  handled <- TRUE

  if(model %in% c("SkewGaussian", "SkewGauss")) {
    sim <- mm + sk * c(abs(dd[, , 1L])) + sqrt(vv) * c(t(simDD))
  } else if(model == "SkewStudentT") {
    denom <- sqrt(.GeoSimulationSquareSums(dd, seq_len(k - 2L)) / (k - 2L))
    bb <- sk * abs(dd[, , k - 1L]) + sqrt(1 - sk^2) * t(simDD)
    sim <- mm + sqrt(vv) * (bb / denom)
  } else if(model == "StudentT") {
    denom <- sqrt(.GeoSimulationSquareSums(dd, seq_len(k - 1L)) / (k - 1L))
    sim <- mm + sqrt(vv) * (c(t(simDD)) / denom)
  } else if(model %in% c("TwoPieceGaussian", "TwoPieceGauss")) {
    base <- t(simDD)
    discrete <- dd[, , 1L]
    pp <- qnorm((1 - sk) / 2)
    sel <- discrete <= pp
    discrete[sel] <- 1 - sk
    discrete[!sel] <- -1 - sk
    sim <- mm + c(sqrt(vv) * (abs(base) * discrete))
  } else if(model == "TwoPieceTukeyh") {
    base <- t(simDD)
    base <- base * exp(tl * base^2 / 2)
    discrete <- dd[, , 1L]
    pp <- qnorm((1 - sk) / 2)
    sel <- discrete <= pp
    discrete[sel] <- 1 - sk
    discrete[!sel] <- -1 - sk
    sim <- mm + c(sqrt(vv) * (abs(base) * discrete))
  } else if(model == "TwoPieceBimodal") {
    sq <- .GeoSimulationSquareSums(dd, seq_len(k - 1L))
    alpha <- 2 * (bimo + 1) / (k - 1L)
    base <- sq / 2^(1 - alpha / 2)
    discrete <- dd[, , k]
    pp <- qnorm((1 - sk) / 2)
    sel <- discrete <= pp
    discrete[sel] <- 1 - sk
    discrete[!sel] <- -1 - sk
    sim <- mm + c(sqrt(vv) * base^(1 / alpha) * discrete)
  } else if(model == "TwoPieceStudentT") {
    denom <- sqrt(.GeoSimulationSquareSums(dd, seq_len(k - 2L)) / (k - 2L))
    base <- c(t(simDD)) / denom
    discrete <- dd[, , k]
    pp <- qnorm((1 - sk) / 2)
    sel <- discrete <= pp
    discrete[sel] <- 1 - sk
    discrete[!sel] <- -1 - sk
    sim <- mm + c(sqrt(vv) * (abs(base) * discrete))
  } else if(model %in% c("LogLogistic", "Logistic", "SkewLaplace")) {
    sim1 <- .GeoSimulationSquareSums(dd, 1:2) / 2
    sim2 <- .GeoSimulationSquareSums(dd, 3:4) / 2
    if(model == "LogLogistic") {
      sim <- exp(mm) * (sim1 / sim2)^(1 / param$shape) /
        (gamma(1 + 1 / param$shape) * gamma(1 - 1 / param$shape))
    } else if(model == "Logistic") {
      sim <- mm + log(sim1 / sim2) * (vv)^(0.5)
    } else {
      sim <- mm + (sim1 / param$skew - sim2 / (1 - param$skew)) * (vv)^(0.5)
    }
  } else if(model %in% c("Gamma", "Weibull")) {
    sq <- .GeoSimulationSquareSums(dd, seq_len(k))
    if(model == "Weibull") {
      sim <- exp(mm) * (sq / 2)^(1 / param$shape) / gamma(1 + 1 / param$shape)
    } else {
      sim <- exp(mm) * sq / k
    }
  } else if(model %in% c("Beta", "Kumaraswamy", "Kumaraswamy2")) {
    if(model == "Beta") {
      sim1 <- .GeoSimulationSquareSums(dd, seq_len(as.integer(param$shape1)))
      sim2 <- .GeoSimulationSquareSums(
        dd,
        seq.int(as.integer(param$shape1) + 1L,
                as.integer(param$shape1) + as.integer(param$shape2))
      )
      sim <- param$min + (param$max - param$min) * sim1 / (sim1 + sim2)
    } else {
      sim1 <- .GeoSimulationSquareSums(dd, 1:2)
      sim2 <- .GeoSimulationSquareSums(dd, 3:4)
      ratio <- sim1 / (sim1 + sim2)
      if(model == "Kumaraswamy") {
        sim <- param$min + (param$max - param$min) *
          (1 - ratio^(1 / param$shape1))^(1 / param$shape2)
      } else {
        med <- plogis(mm)
        shape1 <- as.numeric(param$shape)
        shape2 <- log(0.5) / log1p(-(med^shape1))
        if(any(!is.finite(shape2)) || any(shape2 <= 0)) {
          stop("Kumaraswamy2 produced invalid location-specific shape parameters.", call. = FALSE)
        }
        sim <- param$min + (param$max - param$min) *
          (1 - ratio^(1 / shape2))^(1 / shape1)
      }
    }
  } else if(model == "Wrapped") {
    if(spacetime) mm <- matrix(mm, nrow = nrow(sim), ncol = ncol(sim), byrow = TRUE)
    sim <- (sim + mm) %% (2 * pi)
    needs_format <- FALSE
  } else if(model == "Gaussian") {
    sim <- c(sim)
    byrow <- FALSE
  } else if(model %in% c("LogGaussian", "LogGauss")) {
    base <- c(t(sim))
    sim <- exp(mm) * (exp(sqrt(vv) * base) / exp(vv / 2))
  } else if(model == "Tukeygh") {
    base <- c(t(sim))
    if(!sk && !tl) sim <- mm + sqrt(vv) * base
    if(!sk && tl) sim <- mm + sqrt(vv) * base * exp(tl * base^2 / 2)
    if(!tl && sk) sim <- mm + sqrt(vv) * (exp(sk * base) - 1) / sk
    if(tl && sk) sim <- mm + sqrt(vv) * (exp(sk * base) - 1) * exp(0.5 * tl * base^2) / sk
  } else if(model == "Tukeyh") {
    base <- c(t(sim))
    if(!tl) sim <- mm + sqrt(vv) * base
    if(tl) sim <- mm + sqrt(vv) * base * exp(tl * base^2 / 2)
  } else if(model == "Tukeyh2") {
    base <- c(t(sim))
    sel <- base >= 0
    transformed <- numeric(length(base))
    transformed[sel] <- base[sel] * exp(t1l * base[sel]^2 / 2)
    transformed[!sel] <- base[!sel] * exp(t2l * base[!sel]^2 / 2)
    sim <- mm + sqrt(vv) * transformed
  } else if(model == "SinhAsinh") {
    base <- c(t(sim))
    sim <- mm + sqrt(vv) * sinh((1 / tl) * (asinh(base) + sk))
  } else {
    handled <- FALSE
  }

  list(sim = sim, byrow = byrow, needs_format = needs_format, handled = handled)
}

