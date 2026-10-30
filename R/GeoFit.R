####################################################
### File name: GeoFit.r
####################################################

.GeoFit_validate_likelihood_type <- function(likelihood, type) {
  if (!is.character(likelihood) || length(likelihood) != 1L || is.na(likelihood))
    stop("likelihood must be a single character string.", call. = FALSE)
  if (!is.character(type) || length(type) != 1L || is.na(type))
    stop("type must be a single character string.", call. = FALSE)

  if (identical(likelihood, "Difference") || identical(type, "Difference")) {
    stop("The Difference composite likelihood is no longer supported.", call. = FALSE)
  }

  allowed <- switch(likelihood,
    Full = c("Standard", "Tapering", "Tapering1", "Tapering2", "CV"),
    Marginal = c("Pairwise", "Independence"),
    Conditional = "Pairwise",
    Marginal_2 = "Pairwise",
    NULL
  )
  if (is.null(allowed)) {
    stop("The setting name of the likelihood objects is not correct.", call. = FALSE)
  }
  if (!(type %in% allowed)) {
    stop(sprintf(
      "likelihood='%s' is not compatible with type='%s'. Allowed type value(s): %s.",
      likelihood, type, paste(allowed, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}


.GeoFit_parameter_suggestions <- function(bad, valid) {
  valid <- unique(valid[nzchar(valid)])
  if (!length(valid)) return(character(0))
  out <- character(0)
  for (nm in bad) {
    dd <- as.numeric(utils::adist(nm, valid, ignore.case = TRUE))
    if (!length(dd) || all(!is.finite(dd))) next
    jj <- which.min(dd)
    lim <- max(1, min(3, ceiling(nchar(nm) / 2)))
    if (is.finite(dd[jj]) && dd[jj] <= lim)
      out <- c(out, sprintf("'%s' -> '%s'", nm, valid[jj]))
  }
  out
}

.GeoFit_unit_sill_models <- function() {
  c(
    "Weibull", "Poisson", "Binomial", "Gamma", "LogLogistic",
    "BinomialNeg", "Binary", "Bernoulli", "Geom", "Geometric",
    "Gaussian_misp_Poisson", "Binary_misp_BinomialNeg",
    "PoissonZIP", "Gaussian_misp_PoissonZIP", "BinomialNegZINB", "BinomNeg",
    "PoissonZIP1", "Gaussian_misp_PoissonZIP1", "BinomialNegZINB1",
    "Gaussian_misp_PoissonGamma", "PoissonGamma", "PoissonGammaZIP",
    "PoissonGammaZIP1", "Beta2", "Kumaraswamy2", "Beta", "Kumaraswamy"
  )
}

.GeoFit_normalize_unit_sill <- function(start, fixed, model, bivariate) {
  if (bivariate || !(model %in% .GeoFit_unit_sill_models())) return(fixed)

  if (!is.null(start$sill))
    stop("'sill' is fixed at 1 for this model and must not be supplied in 'start'.",
         call. = FALSE)

  if (is.null(fixed)) fixed <- list()
  if (!is.null(fixed$sill)) {
    fsill <- fixed$sill
    ok <- is.numeric(fsill) && length(fsill) == 1L && !is.na(fsill) &&
      is.finite(fsill) && isTRUE(all.equal(as.numeric(fsill), 1,
                                           tolerance = sqrt(.Machine$double.eps)))
    if (!ok)
      stop("'sill' is fixed at 1 for this model; if supplied in 'fixed', it must equal 1.",
           call. = FALSE)
  }
  fixed$sill <- 1
  fixed
}

.GeoFit_validate_parameter_names <- function(start, fixed, model, corrmodel,
                                             copula, bivariate, X) {
  num_betas <- if (!bivariate) {
    if (is.null(X)) 1L
    else if (is.list(X)) ncol(X[[1L]])
    else ncol(as.matrix(X))
  } else {
    if (is.null(X)) c(1L, 1L)
    else if (is.list(X)) c(ncol(X[[1L]]), ncol(X[[2L]]))
    else rep(ncol(as.matrix(X)), 2L)
  }

  nuisance_start <- unique(NuisParam(model, bivariate, num_betas, copula = copula))
  nuisance_fixed <- unique(NuisParam2(model, bivariate, num_betas, copula = copula))
  corr <- unique(CorrParam(corrmodel))

  check_one <- function(x, arg, nuisance) {
    valid <- unique(c(nuisance, corr))
    if (is.null(x) || !length(x)) return(invisible(TRUE))
    if (!is.list(x)) return(invisible(TRUE))  # CkInput reports the type error.

    nm <- names(x)
    if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) {
      stop(sprintf("GeoFit: every entry of '%s' must have a parameter name.", arg),
           call. = FALSE)
    }
    dup <- unique(nm[duplicated(nm)])
    if (length(dup)) {
      stop(sprintf("GeoFit: duplicated parameter name(s) in '%s': %s.",
                   arg, paste(dup, collapse = ", ")), call. = FALSE)
    }

    bad <- setdiff(nm, valid)
    if (!length(bad)) return(invisible(TRUE))

    cp <- if (is.null(copula)) "" else sprintf(" with copula='%s'", copula)
    msg <- sprintf(
      "GeoFit: invalid parameter name(s) in '%s': %s.\nValid nuisance/marginal parameter names for model='%s'%s: %s.\nValid correlation parameter names for corrmodel='%s': %s.",
      arg, paste(bad, collapse = ", "), model, cp,
      if (length(nuisance)) paste(nuisance, collapse = ", ") else "none",
      corrmodel, if (length(corr)) paste(corr, collapse = ", ") else "none"
    )
    sug <- .GeoFit_parameter_suggestions(bad, valid)
    if (length(sug))
      msg <- paste0(msg, "\nPossible correction(s): ", paste(sug, collapse = "; "), ".")
    stop(msg, call. = FALSE)
  }

  check_one(start, "start", nuisance_start)
  check_one(fixed, "fixed", nuisance_fixed)

  if (!is.null(start) && !is.null(fixed) && is.list(start) && is.list(fixed)) {
    overlap <- intersect(names(start), names(fixed))
    if (length(overlap)) {
      stop(sprintf(
        "GeoFit: parameter(s) cannot be supplied in both 'start' and 'fixed': %s.",
        paste(overlap, collapse = ", ")
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

.GeoFit_check_st_pairwise_support <- function(model, likelihood) {
  code <- CkModel(model)
  ## These sets mirror the space-time native routines registered in
  ## GeoModels_init.c.  Spatial-only kernels are deliberately not considered
  ## here: this guard exists only to prevent an unavailable _st2mem dispatch.
  marginal_codes <- c(
    1L, 2L, 10L, 11L, 12L, 13L, 14L, 16L, 20L, 21L, 22L, 24L, 25L,
    26L, 27L, 28L, 29L, 30L, 33L, 34L, 35L, 36L, 38L, 39L, 40L, 42L,
    43L, 44L, 45L, 46L, 47L, 49L, 59L
  )
  conditional_codes <- marginal_codes
  allowed <- if (identical(likelihood, "Marginal")) marginal_codes else conditional_codes

  if (!(code %in% allowed)) {
    stop(sprintf(
      "Model '%s' is not implemented for spatio-temporal %s pairwise composite likelihood.",
      model, tolower(likelihood)
    ), call. = FALSE)
  }
  invisible(TRUE)
}

GeoFit <- function(data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
                   coordx_dyn = NULL, copula = NULL, corrmodel = NULL,
                   distance = "Eucl", fixed = NULL, anisopars = NULL,
                   est.aniso = c(FALSE, FALSE), grid = FALSE,
                   likelihood = "Marginal", lower = NULL, maxdist = Inf,
                   neighb = NULL, p_neighb = 1, maxtime = Inf, memdist = TRUE,
                   method = "cholesky", model = "Gaussian", n = 1,
                   onlyvar = FALSE, optimizer = "Nelder-Mead", radius = 1,
                   score = FALSE, sensitivity = FALSE, sparse = FALSE,
                   start = NULL, thin_method = "bernoulli", type = "Pairwise",
                   upper = NULL, varest = FALSE, weighted = FALSE, X = NULL,
                    spobj = NULL, spdata = NULL, check.duplicates = FALSE)
{
  ###########  first preliminary check  ###############
  call <- match.call()

  if (!isTRUE(memdist)) {
    warning(
      "'memdist = FALSE' is deprecated and ignored; GeoFit now always precomputes the selected pair structure.",
      call. = FALSE
    )
  }
  memdist <- TRUE

  {

    ## start is genuinely optional.  WlsStart/StartParam provide conservative
    ## internal initial values for parameters omitted by the user.
    if (is.null(start)) start <- list()
    grid_requested <- isTRUE(grid)

    ## normalize strings early (before checks)
    if (!is.null(corrmodel)) corrmodel <- gsub("[[:blank:]]", "", corrmodel)
    if (!is.null(model))     model     <- gsub("[[:blank:]]", "", model)
    if (!is.null(distance))  distance  <- gsub("[[:blank:]]", "", distance)
    if (!is.null(optimizer)) optimizer <- gsub("[[:blank:]]", "", optimizer)
    if (!is.null(likelihood))likelihood<- gsub("[[:blank:]]", "", likelihood)
    if (!is.null(type))      type      <- gsub("[[:blank:]]", "", type)

    .GeoFit_validate_likelihood_type(likelihood, type)

    ## default corrmodel for Independence+Marginal
    if (is.null(corrmodel) && likelihood == "Marginal" && type == "Independence") {
      if (is.null(coordt)) {
        corrmodel <- "Exponential"
        tlist <- list(nugget = 0, scale = 1)
      } else {
        corrmodel <- "Exp_Exp"
        tlist <- list(nugget = 0, scale_s = 1, scale_t = 1)
      }
      fixed <- append(fixed, tlist)
    }
    nosym = FALSE
    ## checks on corrmodel/model/copula
    if (!is.character(corrmodel) || is.null(CkCorrModel(corrmodel)))
      stop("the name of the correlation model is wrong\n")
    model_code <- CkModel(model)
    if (is.null(model_code))
      stop("The name of the model is not correct\n")

    ## GeoLik has explicit full-likelihood objectives only for this set of
    ## model codes.  Older calls for other models (notably Gamma/Weibull)
    ## reached GeoLik with no objective name selected and failed with internal
    ## errors such as `object 'lname' not found`.  Stop here with a stable,
    ## user-facing capability message instead.
    if (likelihood == "Full") {
      full_model_codes <- c(1L, 20L, 22L, 34L, 35L, 36L, 37L, 40L, 47L)
      if (!(model_code %in% full_model_codes)) {
        stop(
          sprintf(
            "GeoFit: Full likelihood is not implemented for model='%s'. Use pairwise likelihood where available.",
            model
          ),
          call. = FALSE
        )
      }
    }

    ## Preserve the inferential model separately from the data-generating
    ## model.  For misspecified likelihoods these are intentionally different
    ## and bootstrap routines must be able to recover the original criterion.
    estimation_model <- model

    if (!is.null(copula)) {
      if (!copula %in% c("Clayton", "Gaussian", "SkewGaussian"))
        stop("the type of copula is wrong\n")
    }

    ## Independence requires Marginal likelihood
    if (type == "Independence" && likelihood != "Marginal")
      stop("Independence likelihood must be coupled with Marginal likelihood\n")

    ## Pairwise requirements (early)
    if (type == "Pairwise") {
      if (is.null(coordt)) {
        if (is.null(neighb) && is.null(maxdist))
          stop("neighb or maxdist must be fixed\n")
      } else {
        ## space-time: require (neighb OR maxdist) AND maxtime
        if ((is.null(neighb) && is.null(maxdist)) || is.null(maxtime))
          stop("neighb or maxdist and maxtime must be fixed\n")
      }
    }

    taper <- NULL
    tapsep <- NULL

    if (likelihood == "Marginal" && type == "Independence") {
      anisopars <- NULL
      est.aniso <- c(FALSE, FALSE)
    }

    if (!is.null(X) && is.null(coordx_dyn)) X <- as.matrix(X)

    if (!is.null(neighb)) {
      if (!is.numeric(neighb) || length(neighb) != 1L || !is.finite(neighb) ||
          neighb < 1 || abs(neighb - round(neighb)) > sqrt(.Machine$double.eps)) {
        stop("neighb must be a finite positive integer", call. = FALSE)
      }
      neighb <- as.integer(round(neighb))
    }

    if (!is.null(anisopars) && !is.list(anisopars))
      stop("anisopars must be a list with two elements\n")

    if (!is.character(optimizer)) stop("invalid optimizer\n")
    if (!is.character(distance))  stop("invalid distance\n")

  
    bivariate <- CheckBiv(CkCorrModel(corrmodel))
    spacetime <- CheckST(CkCorrModel(corrmodel))
    space <- !spacetime && !bivariate

    if (bivariate && type == "Pairwise") {
      if (CkModel(model) != 1L || likelihood != "Marginal") {
        stop(
          "Bivariate pairwise fitting is currently implemented only for model='Gaussian' with likelihood='Marginal'.",
          call. = FALSE
        )
      }
    }

    if (spacetime && type == "Pairwise" &&
        likelihood %in% c("Marginal", "Conditional", "Marginal_2")) {
      .GeoFit_check_st_pairwise_support(
        model, if (likelihood == "Marginal_2") "Conditional" else likelihood
      )
    }

    ## Pairwise copula likelihood is currently a univariate spatial implementation.
    ## Keep the accepted marginal set synchronized with biv_cop() and the
    ## native Comp_Pair_*Cop2mem wrappers.
    if (!is.null(copula) && type == "Pairwise") {
      if (!space)
        stop("Pairwise copula fitting is currently implemented only for univariate purely spatial data.",
             call. = FALSE)
      if (!(likelihood %in% c("Marginal", "Conditional")))
        stop("Copula fitting with type='Pairwise' requires likelihood='Marginal' or 'Conditional'.",
             call. = FALSE)

      model_cop <- switch(model,
        Gauss = "Gaussian",
        LogGauss = "LogGaussian",
        weibull = "Weibull",
        gamma = "Gamma",
        model
      )
      continuous_copula_margins <- c(
        "Gaussian", "StudentT", "LogGaussian", "Gamma", "Weibull",
        "Beta", "Beta2", "Kumaraswamy", "Kumaraswamy2",
        "Logistic", "SkewLaplace"
      )
      discrete_copula_margins <- c("Poisson", "Binomial", "BinomialNeg")
      supported_copula_margins <- if (copula %in% c("Gaussian", "SkewGaussian"))
        c(continuous_copula_margins, discrete_copula_margins)
      else
        continuous_copula_margins

      if (!(model_cop %in% supported_copula_margins)) {
        stop(sprintf(
          "Marginal model '%s' is not implemented for pairwise %s-copula fitting. Supported margins are: %s.",
          model, copula, paste(supported_copula_margins, collapse = ", ")
        ), call. = FALSE)
      }

      ## The Clayton option in GeoModels is the constructive Clayton-like RF
      ## of Bevilacqua et al.; its nu parameter is a positive integer. Because
      ## nu is discrete, it must be fixed rather than optimized continuously.
      if (copula == "Clayton") {
        ## Use exact name matching: '$nu' can partially match 'nugget' in R.
        has_nu_fixed <- !is.null(fixed) && "nu" %in% names(fixed)
        has_nu_start <- !is.null(start) && "nu" %in% names(start)
        if (!has_nu_fixed || has_nu_start)
          stop("For copula='Clayton', nu must be supplied only in fixed and must be a positive integer.",
               call. = FALSE)
        nu_clay <- as.numeric(fixed[["nu", exact = TRUE]])
        if (length(nu_clay) != 1L || !is.finite(nu_clay) || nu_clay <= 0 ||
            abs(nu_clay - round(nu_clay)) > sqrt(.Machine$double.eps))
          stop("For copula='Clayton', fixed$nu must be a finite positive integer.",
               call. = FALSE)
      }

      ## For the skew-Gaussian copula, nu is the bounded eta parameter in the
      ## paper. It may be estimated (start$nu) or fixed (fixed$nu), but it must
      ## be supplied in exactly one of them. Exact name matching is essential:
      ## in R, fixed$nu can otherwise partially match fixed$nugget.
      ## The objective also guards every optimizer proposal so Nelder-Mead
      ## cannot leave (-1,1).
      if (copula == "SkewGaussian") {
        has_nu_start <- !is.null(start) && "nu" %in% names(start)
        has_nu_fixed <- !is.null(fixed) && "nu" %in% names(fixed)
        if (identical(has_nu_start, has_nu_fixed))
          stop("For copula='SkewGaussian', nu (eta) must be supplied exactly once: either in start (estimated) or in fixed.",
               call. = FALSE)
        nu0 <- if (has_nu_start) start[["nu", exact = TRUE]] else fixed[["nu", exact = TRUE]]
        nu0 <- suppressWarnings(as.numeric(nu0))
        if (length(nu0) != 1L || !is.finite(nu0) || abs(nu0) >= 1)
          stop("For copula='SkewGaussian', nu (eta) must lie strictly between -1 and 1.",
               call. = FALSE)
      }
    }

    ###### checking if neighb or maxdist or maxtime has been specified when using cl
    if (space) {
      if (type == "Pairwise" && (likelihood == "Marginal" || likelihood == "Conditional")) {
        if (is.null(neighb) && isTRUE(maxdist == Inf))
          stop("neighb or maxdist must be specificed when using marginal or conditional pairwise composite likelihood\n")
      }
    }

    if (bivariate) {
      if (type == "Pairwise" && (likelihood == "Marginal" || likelihood == "Conditional")) {
        ## keep your original intent, but avoid sum(is.null(...)) which is odd
        if (is.null(neighb) && isTRUE(maxdist == Inf))
          stop("neighb or maxdist must be specificed when using marginal or conditional pairwise composite likelihood\n")
      }
    }

    if (spacetime) {
      if (type == "Pairwise" && (likelihood == "Marginal" || likelihood == "Conditional")) {
        if ((is.null(neighb) && isTRUE(maxdist == Inf)) && isTRUE(maxtime == Inf))
          stop("neighb or maxdist and maxtime must be specificed when using marginal or conditional pairwise composite likelihood\n")
        if ((is.null(neighb) && isTRUE(maxdist == Inf)) && isTRUE(maxtime < Inf))
          stop("neighb or maxdist must be specificed when using marginal or conditional pairwise composite likelihood\n")
        if ((!is.null(neighb) || isTRUE(maxdist < Inf)) && isTRUE(maxtime == Inf))
          stop("maxtime must be specificed when using marginal or conditional pairwise composite likelihood\n")
      }
    }
    ########

    ##############################################################################
    ###### extracting sp object informations if necessary              ###########
    ##############################################################################
    if (!is.null(spobj)) {
      if (space || bivariate) {
        a <- sp2Geo(spobj, spdata)
        coordx <- a$coords
        if (!a$pj) { if (distance != "Chor") distance <- "Geod" }
      }
      if (spacetime) {
        a <- sp2Geo(spobj, spdata)
        coordx <- a$coords
        coordt <- a$coordt
        if (!a$pj) { if (distance != "Chor") distance <- "Geod" }
      }
      if (!is.null(a$Y)) data <- a$Y
      if (!is.null(a$X)) X <- a$X
    }

    if (!is.null(anisopars) && tolower(distance) != "eucl") {
      stop("anisopars requires distance='Eucl'", call. = FALSE)
    }

    ## Canonicalize regular-grid inputs before generic input checking and
    ## likelihood setup.  This avoids dimension-dependent grid branches and
    ## makes spatial, space-time and bivariate grids use the same explicit
    ## coordinate/data layouts as irregular designs.
    if (grid_requested) {
      if (!is.null(coordx_dyn))
        stop("grid=TRUE cannot be combined with coordx_dyn", call. = FALSE)
      gg <- .GeoVariogram_normalize_grid(
        data = data, coordx = coordx, coordy = coordy, coordz = coordz,
        coordt = coordt, bivariate = bivariate
      )
      data <- gg$data
      coordx <- gg$coordx
      coordy <- gg$coordy
      coordz <- gg$coordz
      grid <- FALSE
    }

    if (spacetime && !is.null(coordx_dyn)) {
      coordx_dyn <- .GeoNormalizeDynamicCoordinates(
        coordx_dyn, coordt = coordt, context = "GeoFit"
      )
      .GeoValidateDynamicData(
        data, coordx_dyn, context = "GeoFit", require_list = TRUE
      )
      .GeoValidateDynamicRows(X, coordx_dyn, name = "X", context = "GeoFit")
    }
    if (bivariate) {
      biv_info <- .GeoBivariateCoordinateInfo(
        coordx = coordx, coordy = coordy, coordz = coordz,
        coordx_dyn = coordx_dyn, grid = grid, context = "GeoFit"
      )
      if (biv_info$dynamic) coordx_dyn <- biv_info$blocks
      biv_data <- .GeoBivariateDataVector(
        data, biv_info$ns, context = "GeoFit", allow_na = TRUE
      )
      ## Downstream fitting code historically expects a 2 x n matrix for
      ## common locations and a two-element list for heterotopic locations.
      data <- if (biv_info$dynamic) {
        list(
          biv_data[seq_len(biv_info$ns[1L])],
          biv_data[biv_info$ns[1L] + seq_len(biv_info$ns[2L])]
        )
      } else {
        matrix(biv_data, nrow = 2L, byrow = TRUE)
      }
      if (!is.null(X)) {
        X <- .GeoBivariateDesign(X, biv_info$ns, context = "GeoFit")
      }
      .GeoBivariateValidateFitMeans(
        start, fixed, X, biv_info$ns, context = "GeoFit"
      )
    }

    .GeoCheckDuplicateLocations(
      coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
      coordx_dyn = coordx_dyn, spacetime = spacetime, bivariate = bivariate,
      grid = grid, context = "GeoFit", check.duplicates = check.duplicates
    )
    .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
    on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

    ##############################################################################
    ###### Normalize and validate the univariate mean specification
    ##
    ## A vector fixed$mean is an externally supplied mean at the observations.
    ## It must be removed before CkInput/WlsStart, otherwise unlist() turns it
    ## into apparent coefficients mean1, mean2, ... .  A scalar fixed$mean is
    ## instead the intercept coefficient beta[1].
    MM <- NULL
    if (!bivariate) {
      nobs_mean <- length(unlist(data, use.names = FALSE))
      mean_check <- .GeoMean_validate_fit(start, fixed, X, nobs_mean)
      X <- mean_check$X
      MM <- mean_check$external
      if (!is.null(MM)) fixed$mean <- 0
    }

    ## Validate user-supplied parameter names before internal defaults (such as
    ## nugget=0 or the fixed unit sill used by constructive count/positive
    ## models) are added.  Exact matching is intentional: partial names such as
    ## `shap` must not silently match `shape`.
    .GeoFit_validate_parameter_names(
      start = start, fixed = fixed, model = model, corrmodel = corrmodel,
      copula = copula, bivariate = bivariate, X = X
    )

    ##############################################################################
    ###### setting nugget if missing
    if (!bivariate) {
      if (!sum(substr(names(unlist(append(start, fixed))), 1, 6) == "nugget"))
        fixed$nugget <- 0
    }

    ###############################################################
    ###############################################################
    fixed <- .GeoFit_normalize_unit_sill(start, fixed, model, bivariate)

    ##### all parameters are estimated
    allest <- FALSE
    if (!bivariate) {
      if ((length(c(CorrParam(corrmodel), NuisParam2(model, bivariate, 2, copula = copula))) == length(start)) &&
          is.null(fixed)) {
        fixed <- list(nugget = 0)
        tempstart <- start
        start$nugget <- NULL
        allest <- TRUE
      }
    }

    #############################################################################
  .st_old <- .GeoSTEnter(corrmodel,
    coordx = coordx, coordy = coordy,
    coordz = coordz, coordx_dyn = coordx_dyn,
    distance = distance, context = "GeoFit")
  on.exit(.GeoSTRestore(.st_old), add = TRUE)

    checkinput <- CkInput(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data,
                          distance, "Fitting", fixed, grid, likelihood, maxdist, maxtime,
                          model, n, optimizer, NULL, radius, start, taper, tapsep, type,
                          varest, weighted, copula, X)

    if (!is.null(checkinput$error))
      stop(checkinput$error)

    ### Initialization global variables:
    GeoFit <- NULL
    sensmat <- varcov <- varimat <- parscale <- NULL

    ### Initialization parameters:
    coordt <- unname(coordt)
    if (is.null(coordx_dyn)) {
      coordx <- unname(coordx)
      coordy <- unname(coordy)
      coordz <- unname(coordz)
    }

    .fit_context <- getOption("GeoModels.fit_context_cache", NULL)
    .fit_success <- FALSE
    .fit_uses_globals <- !(likelihood == "Marginal" && type == "Independence")
    on.exit({
      if (.fit_uses_globals && (!.fit_success || !is.environment(.fit_context))) {
        .C("DeleteGlobalVar2", PACKAGE = "GeoModels", NAOK = TRUE)
        if (is.environment(.fit_context) &&
            exists("wls_initparam", envir = .fit_context, inherits = FALSE))
          rm("wls_initparam", envir = .fit_context)
      }
    }, add = TRUE)

    initparam <- WlsStart(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data,
                          distance, "Fitting", fixed, grid, likelihood, maxdist, neighb,
                          maxtime, model, n, NULL, parscale, optimizer == "L-BFGS-B",
                          radius, start, taper, tapsep, type, varest, weighted, copula,
                          X, memdist, nosym, p_neighb, thin_method)

    if (identical(type, "Pairwise") &&
        (!is.numeric(initparam$numpairs) || length(initparam$numpairs) != 1L ||
         !is.finite(initparam$numpairs) || initparam$numpairs < 1L)) {
      stop("No pairwise contributions were retained. Increase 'p_neighb', 'neighb', 'maxdist', or 'maxtime'.",
           call. = FALSE)
    }

    ### fixing initparam if all parameters are estimated
    if (allest) {
      if (!bivariate) {
        bb <- initparam$num_betas
        ccpp <- CorrParam(corrmodel)
        aabb <- NuisParam2(model, bivariate, bb, copula = copula)
        initparam$namescorr <- ccpp
        initparam$namesnuis <- aabb
        initparam$param <- unlist(tempstart)
        initparam$namesparam <- names(initparam$param)
        initparam$fixed <- NULL
        oo <- double(length(ccpp)); names(oo) <- ccpp; oo[names(oo) != "kl"] <- 1
        initparam$flagcorr <- oo
        oo <- double(length(aabb)); names(oo) <- aabb; oo[names(oo) != "kl"] <- 1
        initparam$flagnuis <- oo
      }
    }

    if (type == "Independence") {
      marginal_names <- NuisParam(
        model, bivariate = initparam$bivariate,
        num_betas = initparam$num_betas - 1
      )
      if (!any(marginal_names %in% names(initparam$param)))
        stop("No marginal parameters to estimate")
    }

    ## moving sill from starting to fixed parameters if necessary (in some model sill must be 1)
    if (sum(initparam$namesparam == "sill") == 1) {
      if (initparam$model %in% c(2,14,16,21,28,33,42,50,26,24,30,46,43,11,54)) {
        initparam$param <- initparam$param[initparam$namesparam != "sill"]
        initparam$namesparam <- names(initparam$param)
        a <- 1; names(a) <- "sill"
        initparam$fixed <- c(initparam$fixed, a)
      }
    }

    if (!is.null(initparam$error)) stop(initparam$error)

    ## After WlsStart() has initialized any missing regression coefficients,
    ## enforce the canonical mapping beta = (mean, mean1, ...) to the columns
    ## of X.  This check is deliberately post-initialization so callers may
    ## omit starting values and let WlsStart() provide them.
    if (!bivariate && !is.null(X) && is.null(MM)) {
      p_mean <- ncol(X)
      if (length(initparam$num_betas) != 1L ||
          is.na(initparam$num_betas) || initparam$num_betas != p_mean) {
        stop(
          "The number of mean coefficients initialized by WlsStart() does not match ncol(X).",
          call. = FALSE
        )
      }
      expected_mean <- .GeoMean_names(p_mean)
      available_mean <- unique(c(names(initparam$param), names(initparam$fixed)))
      missing_mean <- setdiff(expected_mean, available_mean)
      unexpected_mean <- setdiff(
        .GeoMean_supplied_names(c(initparam$param, initparam$fixed)),
        expected_mean
      )
      if (length(missing_mean) || length(unexpected_mean)) {
        msg <- paste0(
          "The columns of X must correspond, in order, to coefficients ",
          paste(expected_mean, collapse = ", "), "."
        )
        if (length(missing_mean))
          msg <- paste0(msg, " Missing: ", paste(missing_mean, collapse = ", "), ".")
        if (length(unexpected_mean))
          msg <- paste0(msg, " Unexpected: ", paste(unexpected_mean, collapse = ", "), ".")
        stop(msg, call. = FALSE)
      }
    }

    if (bivariate && !is.null(X)) {
      p_mean <- ncol(X)
      if (length(initparam$num_betas) != 2L ||
          any(is.na(initparam$num_betas)) ||
          any(initparam$num_betas != c(p_mean, p_mean))) {
        stop(
          "The number of bivariate mean coefficients initialized by WlsStart() does not match ncol(X).",
          call. = FALSE
        )
      }
      expected_mean <- c(
        .GeoMean_names(p_mean, "mean_1"),
        .GeoMean_names(p_mean, "mean_2")
      )
      available_mean <- unique(c(names(initparam$param), names(initparam$fixed)))
      supplied_mean <- unique(c(
        .GeoMean_supplied_names(c(initparam$param, initparam$fixed), "mean_1"),
        .GeoMean_supplied_names(c(initparam$param, initparam$fixed), "mean_2")
      ))
      missing_mean <- setdiff(expected_mean, available_mean)
      unexpected_mean <- setdiff(supplied_mean, expected_mean)
      if (length(missing_mean) || length(unexpected_mean)) {
        msg <- paste0(
          "The columns of X must correspond, in order, to bivariate coefficients ",
          paste(expected_mean, collapse = ", "), "."
        )
        if (length(missing_mean))
          msg <- paste0(msg, " Missing: ", paste(missing_mean, collapse = ", "), ".")
        if (length(unexpected_mean))
          msg <- paste0(msg, " Unexpected: ", paste(unexpected_mean, collapse = ", "), ".")
        stop(msg, call. = FALSE)
      }
    }

    ## Defensive fallback: some StartParam/WlsStart paths may not carry this
    ## flag into initparam.  Downstream methods such as GeoCovariogram()
    ## expect every GeoFit object to contain a non-NULL logical bivariate field.
    if (is.null(initparam$bivariate)) initparam$bivariate <- bivariate

    if (!(optimizer %in% c("L-BFGS-B","nlminb","nlm","nmkb","nmk","multiNelder-Mead",
                           "multinlminb","BFGS","Nelder-Mead","optimize","SANN",
                           "bobyqa","sbplx")))
      stop("optimizer is not correct\n")

    if ((optimizer %in% c("L-BFGS-B","nlminb","nmkb","multinlminb","bobyqa","sbplx")) &&
        is.null(lower) && is.null(upper))
      stop("lower and upper bound are missing\n")

    ######################## handling lower and upper bound parameters####################################
    if (optimizer %in% c("L-BFGS-B","nlminb","nmkb","multinlminb","multiNelder-Mead","bobyqa","sbplx") ||
        length(initparam$param) == 1) {

      if (!is.null(lower) || !is.null(upper)) {
        if (!is.list(lower) || !is.list(upper)) stop("lower and upper bound must be a list\n")
        if (sum(unlist(lower) > unlist(upper)) > 0)
          stop("some values of the lower bound is greater of the upper bound\n")

        if (sum(names(lower) == "sill") == 1) {
          if (initparam$model %in% c(2,14,16,21,42,50,26,24,30,46,43,11,54)) {
            lower <- lower[names(lower) != "sill"]
            upper <- upper[names(upper) != "sill"]
          }
        }

        ## setting alphabetic order
        lower <- lower[order(names(lower))]
        upper <- upper[order(names(upper))]

        npar <- length(initparam$param)
        ll <- as.numeric(unlist(lower, use.names = TRUE))
        uu <- as.numeric(unlist(upper, use.names = TRUE))
        names(ll) <- names(lower)
        names(uu) <- names(upper)

        if (length(ll) != npar || length(uu) != npar)
          stop("lower and upper bound must be of the same length of starting values\n")

        if (sum(sort(names(initparam$param)) == sort(names(upper))) < npar ||
            sum(sort(names(initparam$param)) == sort(names(lower))) < npar)
          stop("the names of parameters in the lower and/or upper bounds do not match with starting parameters names.\n")

        replace_zero <- ll == 0
        if (!is.na(.GeoSTCode(corrmodel))) replace_zero[names(ll) == "sep"] <- FALSE
        ll[replace_zero] <- .Machine$double.eps
        uu[uu == Inf] <- 1e+12
        initparam$upper <- uu
        initparam$lower <- ll
      }
    }

    ## When only one parameter is estimated, the composite-likelihood code
    ## switches internally to stats::optimize().  optimize() requires finite
    ## scalar bounds; automatic bounds from WlsStart() can be 0 and Inf
    ## for parameters such as scale.  Sanitize only the internal bounds here.
    if (!onlyvar && length(initparam$param) == 1) {
      pname <- names(initparam$param)

      if (is.null(initparam$lower) || is.null(initparam$upper))
        stop("lower and upper bounds are required when only one parameter is estimated\n")

      lb <- initparam$lower
      ub <- initparam$upper

      if (!is.null(names(lb)) && pname %in% names(lb)) lb <- lb[pname] else lb <- lb[1]
      if (!is.null(names(ub)) && pname %in% names(ub)) ub <- ub[pname] else ub <- ub[1]

      lb <- as.numeric(lb)
      ub <- as.numeric(ub)

      if (length(lb) != 1 || length(ub) != 1 || is.na(lb) || is.na(ub))
        stop("invalid lower/upper bounds for optimize()\n")

      p0 <- as.numeric(initparam$param[1L])
      if (!is.finite(p0))
        stop("invalid starting value for one-parameter optimization\n")

      if (lb == 0 && !(pname == "sep" && !is.na(.GeoSTCode(corrmodel))))
        lb <- .Machine$double.eps

      ## Positive parameters keep the broad automatic upper bound; the
      ## one-dimensional composite-likelihood optimizer now works on the log
      ## scale, so [eps, 1e12] no longer destroys resolution near the start.
      ## For genuinely unbounded parameters (for example a single mean), use
      ## a finite start-centred interval because stats::optimize() requires it.
      if (!is.finite(lb) && !is.finite(ub)) {
        span <- 20 * (1 + abs(p0))
        lb <- p0 - span
        ub <- p0 + span
      } else {
        if (!is.finite(lb)) lb <- p0 - 20 * (1 + abs(p0))
        if (!is.finite(ub)) ub <- 1e+12
      }

      if (!is.finite(lb) || !is.finite(ub) || lb >= ub)
        stop("invalid lower/upper bounds for optimize()\n")

      names(lb) <- pname
      names(ub) <- pname
      initparam$lower <- lb
      initparam$upper <- ub
    }
    ############################################################

    ## MM was extracted before CkInput/WlsStart.  Keep the internal scalar
    ## placeholder fixed at zero; all likelihood paths receive MM explicitly.
    if (!is.null(MM)) initparam$mean <- 0

    ###################################################################################
    ####################################################################
    #############  handling anisotropy parameters ######################
    ####################################################################
    update.aniso <- function(param, namesparam, fixed, namesfixed, lower, upper,
                             anisopars, estimate_aniso)
    {
      un_anisopars <- unlist(anisopars)
      namesaniso <- names(un_anisopars)

      anisostart <- unlist(anisopars)[estimate_aniso]
      anisofixed <- unlist(anisopars)[!estimate_aniso]

      if (length(anisostart) == 0) anisostart <- NULL
      if (length(anisofixed) == 0) anisofixed <- NULL

      ll <- c(0, 1)
      uu <- c(pi, 1e+25)
      lwr <- c(lower, ll[estimate_aniso])
      upr <- c(upper, uu[estimate_aniso])

      param <- c(param, anisostart)
      fixed <- c(fixed, anisofixed)
      namesparam <- names(param)
      namesfixed <- names(fixed)

      if (sum(!is.na(fixed[namesaniso]))) { # updating fixed values
        if (!estimate_aniso[2] && estimate_aniso[1]) fixed["ratio"] <- un_anisopars["ratio"]
        if (!estimate_aniso[1] && estimate_aniso[2]) fixed["angle"] <- un_anisopars["angle"]
        if (!estimate_aniso[1] && !estimate_aniso[2]) {
          fixed["angle"] <- un_anisopars["angle"]
          fixed["ratio"] <- un_anisopars["ratio"]
        }
      }

      list(param = param, namesparam = namesparam,
           fixed = fixed, namesfixed = namesfixed,
           lower = lwr, upper = upr)
    }

    aniso <- FALSE
    if (!is.null(anisopars)) {
      aniso <- TRUE
      namesaniso <- c("angle","ratio")
      qq <- update.aniso(initparam$param, initparam$namesparam,
                         initparam$fixed, initparam$namesfixed,
                         initparam$lower, initparam$upper,
                         anisopars, est.aniso)
      initparam$param <- qq$param
      initparam$fixed <- qq$fixed
      initparam$namesparam <- qq$namesparam
      initparam$namesfixed <- qq$namesfixed
      initparam$lower <- qq$lower
      initparam$upper <- qq$upper
    }

    ###################################################################################
    ###################################################################################

    ## Full likelihood:
    if (likelihood == "Full")
      fitted <- Lik(copula, initparam$bivariate, initparam$coordx, initparam$coordy,
                    initparam$coordz, initparam$coordt, coordx_dyn, initparam$corrmodel,
                    unname(initparam$data), initparam$fixed, initparam$flagcorr,
                    initparam$flagnuis, grid, initparam$lower, method, initparam$model,
                    initparam$namescorr, initparam$namesnuis, initparam$namesparam,
                    initparam$numcoord, initparam$numpairs, initparam$numparamcorr,
                    initparam$numtime, optimizer, onlyvar, initparam$param, initparam$radius,
                    initparam$setup, initparam$spacetime, sparse, varest, taper, initparam$type,
                    initparam$upper, initparam$ns, unname(initparam$X), initparam$neighb,
                    MM, aniso, score)

    ## Composite likelihood.  The selected pair structure is always
    ## precomputed and reused during optimization, Hessian and score
    ## evaluation.  The legacy CompLik() path is intentionally retired.
    if ((likelihood %in% c("Marginal","Conditional","Marginal_2")) && type == "Pairwise") {
      fitted <- CompLik2(copula, initparam$bivariate, initparam$coordx, initparam$coordy,
                         initparam$coordz, initparam$coordt, coordx_dyn, initparam$corrmodel,
                         unname(initparam$data), initparam$distance, initparam$flagcorr,
                         initparam$flagnuis, initparam$fixed, grid, initparam$likelihood,
                         initparam$lower, initparam$model, initparam$n, initparam$namescorr,
                         initparam$namesnuis, initparam$namesparam, initparam$numparam,
                         initparam$numparamcorr, optimizer, onlyvar, initparam$param,
                         initparam$spacetime, initparam$type, initparam$upper, varest,
                         initparam$weighted, initparam$ns, unname(initparam$X), sensitivity,
                         initparam$colidx, initparam$rowidx, initparam$neighb, MM, aniso, score)
    }

    if (likelihood == "Marginal" && type == "Independence")
      fitted <- CompIndLik2(initparam$bivariate, initparam$coordx, initparam$coordy, initparam$coordz,
                            initparam$coordt, coordx_dyn, unname(initparam$data),
                            initparam$flagcorr, initparam$flagnuis, initparam$fixed, grid,
                            initparam$lower, initparam$model, initparam$n, initparam$namescorr,
                            initparam$namesnuis, initparam$namesparam, initparam$numparam,
                            optimizer, onlyvar, initparam$param, initparam$spacetime, initparam$type,
                            initparam$upper, names(upper), varest, initparam$ns, unname(initparam$X),
                            sensitivity, copula, MM, score)

    ## misspecified models
    missp <- FALSE
    if (model == "Gaussian_misp_Tukeygh")      { model <- "Tukeygh";      missp <- TRUE }
    if (model == "Gaussian_misp_Poisson")      { model <- "Poisson";      missp <- TRUE }
    if (model == "Gaussian_misp_Binomial")     { model <- "Binomial";     missp <- TRUE }
    if (model == "Gaussian_misp_BinomialNeg")  { model <- "BinomialNeg";  missp <- TRUE }
    if (model == "Gaussian_misp_PoissonGamma") { model <- "PoissonGamma"; missp <- TRUE }
    if (model == "Gaussian_misp_PoissonZIP")   { model <- "PoissonZIP";   missp <- TRUE }
    if (model == "Gaussian_misp_StudentT")     { model <- "StudentT";     missp <- TRUE }
    if (model == "Gaussian_misp_SkewStudentT") { model <- "SkewStudentT"; missp <- TRUE }
    if (model == "Binary_misp_BinomialNeg")    { model <- "BinomialNeg";  missp <- TRUE }

    numtime <- 1
    if (initparam$spacetime) numtime <- length(coordt)
    if (initparam$bivariate) numtime <- 2
    dimat <- initparam$numcoord

    ## Native context cleanup is registered before initialization and runs
    ## on normal return, errors and user interrupts.

    ff <- as.list(initparam$fixed)
    if (!is.null(MM)) ff$mean <- MM

    if (length(initparam$param) == 1) optimizer <- "optimize"

    if (aniso) anisopars <- as.list(c(fitted$par, ff)[namesaniso])

    if (is.null(unlist(ff))) ff <- NULL

    ## special case: maxdist and neighb = NULL
   # if (likelihood != "Full") {
   #   if (is.null(neighb) && is.numeric(maxdist) && likelihood == "Marginal") {
   #     fitted$value <- 2 * fitted$value
   #     initparam$numpairs <- 2 * initparam$numpairs
   #   }
   # }

    if (!is.null(coordt) && is.null(coordx_dyn)) {
      if (is.null(coordz)) {
        initparam$coordx <- initparam$coordx[1:(length(initparam$coordx) / length(initparam$coordt))]
        initparam$coordy <- initparam$coordy[1:(length(initparam$coordy) / length(initparam$coordt))]
      } else {
        initparam$coordx <- initparam$coordx[1:(length(initparam$coordx) / length(initparam$coordt))]
        initparam$coordy <- initparam$coordy[1:(length(initparam$coordy) / length(initparam$coordt))]
        initparam$coordz <- initparam$coordz[1:(length(initparam$coordz) / length(initparam$coordt))]
      }
    }

    if (model %in% c("Weibull","Poisson","Binomial","Gamma","PoissonGamma","Gaussian_misp_PoissonGamma",
                     "LogLogistic","BinomialNeg","BinomNeg","Binary","Bernoulli","Geom","Geometric","Binary_misp_BinomialNeg",
                     "Gaussian_misp_Poisson","PoissonZIP","PoissonGammaZIP","PoissonGammaZIP1",
                     "Gaussian_misp_PoissonZIP","BinomialNegZINB","PoissonZIP1","Gaussian_misp_PoissonZIP1",
                     "BinomialNegZINB1","Beta2","Kumaraswamy2","Beta","Kumaraswamy")) {
      if (!is.null(ff$sill)) ff$sill <- NULL
    }

    conf.int <- NULL
    pvalues <- NULL

    if (!is.null(initparam$coordz) && all(initparam$coordz == 0) &&
        (is.na(.GeoSTCode(corrmodel)) || .GeoSTCurrentDimension() != 3L))
      initparam$coordz <- NULL

    if (likelihood == "Full" && type == "Standard") {
      if (varest) {
        alpha <- 0.95
        conf.int <- pvalues <- NULL
        if (is.numeric(fitted$stderr)) {
          aa <- qnorm(1 - (1 - alpha) / 2) * fitted$stderr
          pp <- as.numeric(fitted$par)
          low <- pp - aa
          upp <- pp + aa
          conf.int <- rbind(low, upp)
          pvalues <- 2 * pnorm(-abs(pp / fitted$stderr))
        }
      }
    }

    if (bivariate && is.null(coordx_dyn)) {
      initparam$coordx <- initparam$coordx[1:dimat]
      initparam$coordy <- initparam$coordy[1:dimat]
    }

    ### Set the output object:
    GeoFit <- list(
      anisopars = anisopars,
      bivariate = if (is.null(initparam$bivariate)) bivariate else initparam$bivariate,
      claic = fitted$claic,
      clbic = fitted$clbic,
      colidx = as.integer(initparam$colidx),
      coordx = initparam$coordx,
      coordy = initparam$coordy,
      coordz = initparam$coordz,
      coordt = initparam$coordt,
      coordx_dyn = coordx_dyn,
      duplicate_locations_checked = isTRUE(check.duplicates),
      conf.int = conf.int,
      convergence = fitted$convergence,
      copula = copula,
      corrmodel = corrmodel,
      data = initparam$data,
      distance = distance,
      est.aniso = est.aniso,
      fixed = ff,
      grid = grid_requested,
      iterations = fitted$counts,
      likelihood = likelihood,
      logCompLik = fitted$value,
      lower = lower,
      message = fitted$message,
      model = model,
      estimation_model = estimation_model,
      n = initparam$n,
      ns = initparam$ns,
      numbetas = initparam$num_betas,
      numcoord = initparam$numcoord,
      numtime = initparam$numtime,
      optimizer = optimizer,
      param = as.list(fitted$par),
      p_neighb = p_neighb,
      thin_method = thin_method,
      memdist = TRUE,
      nozero = initparam$setup$nozero,
      score = fitted$score,
      maxdist = maxdist,
      maxtime = maxtime,
      neighb = initparam$neighb,
      numpairs = initparam$numpairs,
      missp = missp,
      pvalues = pvalues,
      radius = radius,
      rowidx = as.integer(initparam$rowidx),
      spacetime = initparam$spacetime,
      stderr = fitted$stderr,
      sensmat = fitted$sensmat,
      upper = upper,
      varcov = fitted$varcov,
      varimat = fitted$varimat,
      type = type,
      weighted = initparam$weighted,
      X = X
    )

    .fit_result <- structure(c(GeoFit, call = call), class = c("GeoFit"))
    .fit_success <- TRUE
    .fit_result
  }
}

################################################################################

print.GeoFit <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if (x$likelihood == "Full") {
    method <- "Likelihood"
    if (x$type == "Tapering") { claic <- "CLAIC"; clbic <- "CLBIC" }
    else { claic <- "AIC"; clbic <- "BIC" }
  } else {
    method <- "Composite-Likelihood"
    claic <- "CLAIC"
    clbic <- "CLBIC"
  }

  missp <- ""
  if (isTRUE(x$missp)) missp <- "misspecified"

  if (x$model %in% c("Gaussian","Gauss")) { process <- "Gaussian"; model <- "Gaussian" }
  if (x$model == "Gamma") { process <- "Gamma"; model <- "Gamma" }
  if (x$model == "TwoPieceBimodal") { process <- "TwoPieceBimodal"; model <- "TwoPieceBimodal" }
  if (x$model == "LogLogistic") { process <- "LogLogistic"; model <- "LogLogistic" }
  if (x$model == "Gaussian_misp_Poisson") { process <- "Poisson"; model <- "Misspecified Gaussian Poisson" }
  if (x$model == "Gaussian_misp_Binomial") { process <- "Binomial"; model <- "Misspecified Gaussian Binomial" }
  if (x$model == "Gaussian_misp_PoissonZIP") { process <- "PoissonZIP"; model <- "Misspecified Gaussian Poisson Inflated" }
  if (x$model == "Poisson") { process <- "Poisson"; model <- "Poisson" }
  if (x$model == "PoissonGamma") { process <- "PoissonGamma"; model <- "PoissonGamma" }
  if (x$model == "Gaussian_misp_PoissonGamma") { process <- "PoissonGamma"; model <- "Misspecified Gaussian PoissonGamma" }
  if (x$model == "PoissonZIP") { process <- "PoissonZIP"; model <- "PoissonZIP" }
  if (x$model == "PoissonGammaZIP") { process <- "PoissonGammaZIP"; model <- "PoissonGammaZIP" }
  if (x$model == "Beta2") { process <- "Beta2"; model <- "Beta2" }
  if (x$model == "Gaussian_misp_StudentT") { process <- "StudentT"; model <- "Misspecified Gaussian StudentT" }
  if (x$model == "StudentT") { process <- "StudentT"; model <- "StudentT" }
  if (x$model == "Gaussian_misp_Tukeygh") { process <- "Tukeygh"; model <- "Misspecified Gaussian Tukeygh" }
  if (x$model == "Tukeygh") { process <- "Tukeygh"; model <- "Tukeygh" }
  if (x$model == "Gaussian_misp_SkewStudentT") { process <- "SkewStudentT"; model <- "Misspecified Gaussian SkewStudentT" }
  if (x$model == "SkewStudentT") { process <- "SkewStudentT"; model <- "SkewStudentT" }
  if (x$model == "Logistic") { process <- "Logistic"; model <- "Logistic" }
  if (x$model == "Tukeyh") { process <- "Tukeyh"; model <- "Tukeyh" }
  if (x$model == "Tukeyh2") { process <- "Tukeyh2"; model <- "Tukeyh2" }
  if (x$model == "Gamma2") { process <- "Gamma2"; model <- "Gamma2" }
  if (x$model %in% c("LogGauss","LogGaussian")) { process <- "Log Gaussian"; model <- "LogGaussian" }
  if (x$model %in% c("SkewGauss","SkewGaussian")) { process <- "Skew Gaussian"; model <- "SkewGaussian" }
  if (x$model == "SkewLaplace") { process <- "Skew Laplace"; model <- "SkewLaplace" }
  if (x$model == "TwoPieceStudentT") { process <- "TwoPiece StudentT"; model <- "TwoPieceStudentT" }
  if (x$model == "TwoPieceTukeyh") { process <- "TwoPiece Tukeyh"; model <- "TwoPieceTukeyh" }
  if (x$model %in% c("TwoPieceGaussian","TwoPieceGauss")) { process <- "TwoPiece Gaussian"; model <- "TwoPieceGaussian" }
  if (x$model == "SinhAsinh") { process <- "SinhAsinh"; model <- "SinhAsinh" }
  if (x$model == "Wrapped") { process <- "Wrapped"; model <- "Wrapped" }
  if (x$model == "Weibull") { process <- "Weibull"; model <- "Weibull" }
  if (x$model == "Binomial") { process <- "Binomial"; model <- "Binomial" }
  if (x$model == "BinomialLogistic") { process <- "BinomialLogistic"; model <- "BinomialLogistic" }
  if (x$model == "Kumaraswamy") { process <- "Kumaraswamy"; model <- "Kumaraswamy" }
  if (x$model == "Kumaraswamy2") { process <- "Kumaraswamy2"; model <- "Kumaraswamy2" }
  if (x$model == "Beta") { process <- "Beta"; model <- "Beta" }
  if (x$model %in% c("Binomial_TwoPieceGaussian","Binomial_TwoPieceGauss")) { process <- "Binomial TwoPiece Gaussian"; model <- "Binomial_TwoPieceGauss" }
  if (x$model %in% c("BinomialNeg_TwoPieceGaussian","BinomialNeg_TwoPieceGauss")) { process <- "Negative Binomial TwoPiece Gaussian"; model <- "BinomialNeg_TwoPieceGauss" }
  if (x$model == "Binomial2") { process <- "Binomial"; model <- "Binomial2" }
  if (x$model == "BinomialNeg") { process <- "BinomialNeg"; model <- "BinomialNeg" }
  if (x$model == "Binary_misp_BinomialNeg") { process <- "BinomialNeg"; model <- "Misspecified Binary BinomialNeg" }
  if (x$model == "BinomialNegLogistic") { process <- "BinomialNegLogistic"; model <- "BinomialNegLogistic" }
  if (x$model == "BinomialNegZINB") { process <- "BinomialNegZINB"; model <- "BinomialNegZINB" }
  if (x$model %in% c("Geom","Geometric")) { process <- "Geometric"; model <- "Geometric" }
  if (x$model == "PoisBin") { process <- "Poisson Binomial"; model <- "PoisBin" }
  if (x$model == "PoisBinNeg") { process <- "Poisson NegBinomial"; model <- "PoisBinNeg" }

  if (x$bivariate) { biv <- "bivariate"; x$numtime <- 1 }
  else { biv <- "univariate" }

  cat("\n##################################################################")
  cat("\nMaximum", missp, method, "Fitting of", process, "Random Fields\n")
  if (!is.null(x$copula)) cat("\nCopula:", x$copula, "\n")
  cat("\nSetting:", x$likelihood, method, "\n")
  cat("\nModel:", model, "\n")
  cat("\nDistance:", x$distance, "\n")
  cat("\nType of the likelihood objects:", x$type, x$method, "\n")

  ## Pair-thinning information, when available.  For backward
  ## compatibility with older GeoFit objects, retrieve thin_method from
  ## the object if present, otherwise from the stored call; if it is not
  ## available but thinning was used, report the default.
  p_retained <- NULL
  if (!is.null(x$p_neighb)) {
    p_retained <- suppressWarnings(as.numeric(x$p_neighb)[1])
  } else if (!is.null(x$call) && !is.null(as.list(x$call)$p_neighb)) {
    p_retained <- suppressWarnings(as.numeric(as.character(as.list(x$call)$p_neighb))[1])
  }

  thin_method <- NULL
  if (!is.null(x$thin_method)) {
    thin_method <- as.character(x$thin_method)[1]
  } else if (!is.null(x$call) && !is.null(as.list(x$call)$thin_method)) {
    thin_method <- paste(deparse(as.list(x$call)$thin_method), collapse = "")
    thin_method <- gsub('^"|"$', "", thin_method)
    thin_method <- gsub("^'|'$", "", thin_method)
  }

  if (!is.null(p_retained) && is.finite(p_retained) && p_retained < 1) {
    if (is.null(thin_method) || !nzchar(thin_method)) thin_method <- "bernoulli"
    cat("Pair thinning method:", thin_method, "\n")
    cat("(Expected) Retained pairs (%):",
        format(100 * p_retained, digits = digits, nsmall = 2), "\n")
  }

  cat("\nCovariance model:", x$corrmodel, "\n")
  cat("\nOptimizer:", x$optimizer, "\n")
  cat("\nNumber of spatial coordinates:", x$numcoord, "\n")
  if (x$spacetime) cat("Number of dependent temporal realisations:", x$numtime, "\n")
  cat("Type of the random field:", biv, "\n")
  cat("Number of estimated parameters:", length(x$param), "\n")
  cat("\nType of convergence:", x$convergence, "")
  cat("\nMaximum log-", method, " value: ",
      format(x$logCompLik, digits = digits, nsmall = 2), "\n", sep = "")

  if (!is.null(x$claic))
    cat(claic, ":", format(x$claic, digits = digits), "\n")

  if (!is.null(x$clbic))
    cat(clbic, ":", format(x$clbic, digits = digits), "\n")

  cat("\nEstimated parameters:\n")
  print.default(unlist(x$param), digits = digits, print.gap = 2, quote = FALSE)

  if (!is.null(x$stderr)) {
    cat("\nStandard errors:\n")
    print.default(x$stderr, digits = digits, print.gap = 2, quote = FALSE)
  }

  cat("\n##################################################################\n")
  invisible(x)
}


summary.GeoFit <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  print.GeoFit(object, digits = digits, ...)
}