####################################################
### File name: GeoFit.r
####################################################
GeoFit2 <- function(data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
                    coordx_dyn = NULL, copula = NULL, corrmodel,
                    distance = "Eucl", fixed = NULL, anisopars = NULL,
                    est.aniso = c(FALSE, FALSE), grid = FALSE,
                    likelihood = "Marginal", lower = NULL, maxdist = Inf,
                    neighb = NULL, p_neighb = 1, maxtime = Inf, memdist = TRUE,
                    method = "cholesky", model = "Gaussian", n = 1, onlyvar = FALSE,
                    optimizer = "Nelder-Mead", radius = 1, score = FALSE,
                    sensitivity = FALSE, sparse = FALSE, start = NULL,
                    thin_method = "bernoulli", type = "Pairwise", upper = NULL,
                    varest = FALSE, weighted = FALSE, X = NULL, 
                    spobj = NULL, spdata = NULL)
{
  ###########  first preliminary check  ###############
  call <- match.call()

  suppressWarnings({

    ## normalize strings early
    corrmodel  <- gsub("[[:blank:]]", "", corrmodel)
    model      <- gsub("[[:blank:]]", "", model)
    distance   <- gsub("[[:blank:]]", "", distance)
    optimizer  <- gsub("[[:blank:]]", "", optimizer)
    likelihood <- gsub("[[:blank:]]", "", likelihood)
    type       <- gsub("[[:blank:]]", "", type)
    nosym = FALSE
    if (is.null(CkModel(model)))
      stop("The name of the model is not correct\n")

    if (!is.null(copula)) {
      if (!((copula == "Clayton") || (copula == "Gaussian") || (copula == "SkewGaussian")))
        stop("the type of copula is wrong\n")
    }

    if (type == "Independence")
      stop("use GeoFit for Independence composite likelihood\n")

    ## corrmodel check after trimming
    if (!is.character(corrmodel) || is.null(CkCorrModel(corrmodel)))
      stop("the name of the correlation model is wrong\n")

    if (!is.logical(memdist)) memdist <- FALSE
    if (!is.null(X)) X <- as.matrix(X)

    if (is.numeric(neighb)) {
      neighb <- round(neighb)
      if (any(neighb < 1)) stop("neighb must be an integer >= 1\n")
    }

    if (type == "Pairwise") {
      ## minimal requirement: at least one of neighb/maxdist present
      if (is.null(neighb) && is.null(maxdist))
        stop("neighb or maxdist must be fixed for Pairwise composite likelihood\n")
    }

    if (!is.null(anisopars) && !is.list(anisopars))
      stop("anisopars must be a list with two elements\n")

    if (!is.character(optimizer)) stop("invalid optimizer\n")
    if (!is.character(distance))  stop("invalid distance\n")

    bivariate <- CheckBiv(CkCorrModel(corrmodel))
    spacetime <- CheckST(CkCorrModel(corrmodel))
    space <- !spacetime && !bivariate

    taper <- NULL
    tapsep <- NULL

    ###### checking if neighb or maxdist or maxtime has been specified when using cl
    if (space || bivariate) {
      if (type == "Pairwise" && (likelihood == "Marginal" || likelihood == "Conditional")) {
        if (is.null(neighb) && isTRUE(maxdist == Inf))
          stop("neighb or maxdist must be specified when using marginal or conditional pairwise composite likelihood\n")
      }
    }

    if (spacetime) {
      if (type == "Pairwise" && (likelihood == "Marginal" || likelihood == "Conditional")) {
        if ((is.null(neighb) && isTRUE(maxdist == Inf)) && isTRUE(maxtime == Inf))
          stop("neighb or maxdist and maxtime must be specified when using marginal or conditional pairwise composite likelihood\n")
        if ((is.null(neighb) && isTRUE(maxdist == Inf)) && isTRUE(maxtime < Inf))
          stop("neighb or maxdist must be specified when using marginal or conditional pairwise composite likelihood\n")
        if ((!is.null(neighb) || isTRUE(maxdist < Inf)) && isTRUE(maxtime == Inf))
          stop("maxtime must be specified when using marginal or conditional pairwise composite likelihood\n")
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
      if (!is.null(a$Y) && !is.null(a$X)) { data <- a$Y; X <- a$X }
    }

    ###############################################################
    ## setting nugget if missing (univariate only)
    if (!bivariate) {
      if (!sum(names(unlist(append(start, fixed))) == "nugget"))
        fixed$nugget <- 0
    }

    ###############################################################
    ## models where sill must be 1
    if (!bivariate) {
      if (model %in% c("Weibull","Poisson","Binomial","Gamma","LogLogistic",
                       "PoissonGamma","PoissonGammaZIP","PoissonGammaZIP1",
                       "BinomialNeg","Bernoulli","Geometric","Gaussian_misp_Poisson",
                       "PoissonZIP","Gaussian_misp_PoissonZIP","BinomialNegZINB",
                       "PoissonZIP1","Gaussian_misp_PoissonZIP1","BinomialNegZINB1",
                       "Beta2","Kumaraswamy2","Beta","Kumaraswamy")) {
        if (!is.null(start$sill)) stop("sill parameter must not be considered for this model\n")
        fixed$sill <- 1
      }
    }

    #############################################################################
    checkinput <- CkInput(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data,
                          distance, "Fitting", fixed, grid, likelihood, maxdist, maxtime,
                          model, n, optimizer, NULL, radius, start, taper, tapsep,
                          type, varest, weighted, copula, X)

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

    initparam <- WlsStart(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data,
                          distance, "Fitting", fixed, grid, likelihood, maxdist, neighb,
                          maxtime, model, n, NULL, parscale, optimizer == "L-BFGS-B",
                          radius, start, taper, tapsep, type, varest, weighted, copula,
                          X, memdist, nosym, p_neighb, thin_method)

    ## in the case on external fixed mean
    MM <- NULL
    if (!is.null(fixed)) {
      if (length(fixed$mean) > 1) { MM <- as.numeric(fixed$mean); initparam$mean <- 1e-07 }
    }

    if (!is.null(initparam$error)) stop(initparam$error)

    ## checking for upper and lower bound for methods with bounds (+ optimize case)
    if ((optimizer %in% c("L-BFGS-B","nlminb","nmkb","multinlminb","bobyqa","sbplx")) &&
        is.null(lower) && is.null(upper))
      stop("lower and upper bound are missing\n")

    if (!(optimizer %in% c("L-BFGS-B","nlminb","nlm","nmkb","nmk","multiNelder-Mead",
                           "multinlminb","BFGS","Nelder-Mead","optimize","SANN",
                           "bobyqa","sbplx")))
      stop("optimizer is not correct\n")

    #### bounds handling
    if (optimizer %in% c("L-BFGS-B","nlminb","nmkb","multinlminb","multiNelder-Mead","bobyqa","sbplx") ||
        length(initparam$param) == 1) {

      if (!is.null(lower) || !is.null(upper)) {

        if (!is.list(lower) || !is.list(upper))
          stop("lower and upper bound must be a list\n")

        if (sum(unlist(lower) > unlist(upper)) > 0)
          stop("some values of the lower bound is greater than the upper bound\n")

        if (sum(names(lower) == "sill") == 1) {
          if (initparam$model %in% c(2,14,16,21,42,50,26,24,30,46,43,11)) {
            lower <- lower[names(lower) != "sill"]
            upper <- upper[names(upper) != "sill"]
          }
        }

        lower <- lower[order(names(lower))]
        upper <- upper[order(names(upper))]

        npar <- length(initparam$param)
        ll <- as.numeric(lower)
        uu <- as.numeric(upper)

        if (length(ll) != npar || length(uu) != npar)
          stop("lower and upper bound must be of the same length of starting values\n")

        ## safer name matching (order independent)
        if (sum(sort(names(initparam$param)) == sort(names(upper))) < npar ||
            sum(sort(names(initparam$param)) == sort(names(lower))) < npar)
          stop("the names in lower/upper bounds do not match starting parameter names\n")

        ll[ll == 0] <- .Machine$double.eps
        uu[uu == Inf] <- 1e+12

        initparam$upper <- uu
        initparam$lower <- ll
      }
    }

    ###############################################################################################
    ## first pass (your original code)
    fitted_ini <- CompIndLik2(initparam$bivariate, initparam$coordx, initparam$coordy, initparam$coordz, initparam$coordt,
                              coordx_dyn, unname(initparam$data),
                              initparam$flagcorr, initparam$flagnuis, initparam$fixed, grid,
                              initparam$lower, initparam$model, initparam$n,
                              initparam$namescorr, initparam$namesnuis,
                              initparam$namesparam, initparam$numparam, optimizer, onlyvar,
                              initparam$param, initparam$spacetime, initparam$type,
                              initparam$upper, names(upper), varest, initparam$ns, unname(initparam$X),
                              sensitivity, copula, MM, FALSE)

    ######################################################
    ###### updating starting and names parameters
    ######################################################
    namespp <- names(fitted_ini$par)
    aa <- append(initparam$param, initparam$fixed)
    sel <- match(namespp, names(aa))
    sel <- sel[!is.na(sel)]
    aa[sel] <- fitted_ini$par

    sel <- match(names(aa), initparam$namesparam)
    sel <- sel[!is.na(sel)]
    initparam$param <- aa[sel]

    ######################################################
    ## updating with aniso parameters
    update.aniso <- function(param, namesparam, fixed, namesfixed, lower, upper, anisopars, estimate_aniso)
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

      if (sum(!is.na(fixed[namesaniso]))) {
        if (!estimate_aniso[2] && estimate_aniso[1]) fixed["ratio"] <- un_anisopars["ratio"]
        if (!estimate_aniso[1] && estimate_aniso[2]) fixed["angle"] <- un_anisopars["angle"]
        if (!estimate_aniso[1] && !estimate_aniso[2]) {
          fixed["angle"] <- un_anisopars["angle"]
          fixed["ratio"] <- un_anisopars["ratio"]
        }
      }

      list(param = param, fixed = fixed,
           namesparam = namesparam, namesfixed = namesfixed,
           lower = lwr, upper = upr)
    }

    aniso <- FALSE
    if (!is.null(anisopars)) {
      aniso <- TRUE
      namesaniso <- c("angle","ratio")
      qq <- update.aniso(initparam$param, initparam$namesparam, initparam$fixed, initparam$namesfixed,
                         initparam$lower, initparam$upper, anisopars, est.aniso)
      initparam$param <- qq$param
      initparam$fixed <- qq$fixed
      initparam$namesparam <- qq$namesparam
      initparam$namesfixed <- qq$namesfixed
      initparam$lower <- qq$lower
      initparam$upper <- qq$upper
    }

    ###################################################################################
    ## Full likelihood:
    if (likelihood == "Full")
      fitted <- Lik(copula, initparam$bivariate, initparam$coordx, initparam$coordy, initparam$coordz,
                    initparam$coordt, coordx_dyn, initparam$corrmodel, unname(initparam$data),
                    initparam$fixed, initparam$flagcorr, initparam$flagnuis, grid, initparam$lower,
                    method, initparam$model, initparam$namescorr, initparam$namesnuis, initparam$namesparam,
                    initparam$numcoord, initparam$numpairs, initparam$numparamcorr, initparam$numtime,
                    optimizer, onlyvar, initparam$param, initparam$radius, initparam$setup, initparam$spacetime,
                    sparse, varest, taper, initparam$type, initparam$upper, initparam$ns, unname(initparam$X),
                    initparam$neighb, MM, aniso, score)

    ## Composite likelihood:
    if ((likelihood %in% c("Marginal","Conditional","Marginal_2")) && type == "Pairwise") {

      if (!memdist) {
        fitted <- CompLik(copula, initparam$bivariate, initparam$coordx, initparam$coordy, initparam$coordz,
                          initparam$coordt, coordx_dyn, initparam$corrmodel, unname(initparam$data),
                          initparam$distance, initparam$flagcorr, initparam$flagnuis, initparam$fixed, grid,
                          initparam$likelihood, initparam$lower, initparam$model, initparam$n,
                          initparam$namescorr, initparam$namesnuis, initparam$namesparam, initparam$numparam,
                          initparam$numparamcorr, optimizer, onlyvar, initparam$param, initparam$spacetime,
                          initparam$type, initparam$upper, varest, initparam$weighted, initparam$ns,
                          unname(initparam$X), sensitivity, MM, aniso, score)
      }

      if (memdist) {
        fitted <- CompLik2(copula, initparam$bivariate, initparam$coordx, initparam$coordy, initparam$coordz,
                           initparam$coordt, coordx_dyn, initparam$corrmodel, unname(initparam$data),
                           initparam$distance, initparam$flagcorr, initparam$flagnuis, initparam$fixed, grid,
                           initparam$likelihood, initparam$lower, initparam$model, initparam$n,
                           initparam$namescorr, initparam$namesnuis, initparam$namesparam, initparam$numparam,
                           initparam$numparamcorr, optimizer, onlyvar, initparam$param, initparam$spacetime,
                           initparam$type, initparam$upper, varest, initparam$weighted, initparam$ns,
                           unname(initparam$X), sensitivity, initparam$colidx, initparam$rowidx, initparam$neighb,
                           MM, aniso, score)
      }
    }

    ## misspecified models
    missp <- FALSE
    if (model == "Gaussian_misp_Tukeygh")      { model <- "Tukeygh";      missp <- TRUE }
    if (model == "Gaussian_misp_Poisson")      { model <- "Poisson";      missp <- TRUE }
    if (model == "Gaussian_misp_Binomial")     { model <- "Binomial";     missp <- TRUE }
    if (model == "Gaussian_misp_PoissonGamma") { model <- "PoissonGamma"; missp <- TRUE }
    if (model == "Gaussian_misp_PoissonZIP")   { model <- "PoissonZIP";   missp <- TRUE }
    if (model == "Gaussian_misp_StudentT")     { model <- "StudentT";     missp <- TRUE }
    if (model == "Gaussian_misp_SkewStudentT") { model <- "SkewStudentT"; missp <- TRUE }

    ## housekeeping global vars
    if (!(likelihood == "Marginal" && type == "Independence")) {
      if (memdist) .C("DeleteGlobalVar2", PACKAGE = "GeoModels", DUP = TRUE, NAOK = TRUE)
      else         .C("DeleteGlobalVar" , PACKAGE = "GeoModels", DUP = TRUE, NAOK = TRUE)
    }

    ## special case: maxdist and neighb = NULL and likelihood="Marginal"
    if (likelihood != "Full") {
      if (is.null(neighb) && is.numeric(maxdist) && likelihood == "Marginal") {
        fitted$value <- 2 * fitted$value
        initparam$numpairs <- 2 * initparam$numpairs
      }
    }

    ff <- as.list(initparam$fixed)
    if (!is.null(MM)) ff$mean <- MM
    if (is.null(unlist(ff))) ff <- NULL

    if (length(initparam$param) == 1) optimizer <- "optimize"
    if (aniso) anisopars <- as.list(c(fitted$par, ff)[namesaniso])

    if (!is.null(coordt) && is.null(coordx_dyn)) {
      initparam$coordx <- initparam$coordx[1:(length(initparam$coordx) / length(initparam$coordt))]
      initparam$coordy <- initparam$coordy[1:(length(initparam$coordy) / length(initparam$coordt))]
    }

    if (all(initparam$coordz == 0)) initparam$coordz <- NULL

    conf.int <- NULL
    pvalues <- NULL

    if (likelihood == "Full" && type == "Standard" && varest) {
      alpha <- 0.95
      aa <- qnorm(1 - (1 - alpha) / 2) * fitted$stderr
      pp <- as.numeric(fitted$par)
      low <- pp - aa
      upp <- pp + aa
      conf.int <- rbind(low, upp)
      pvalues <- 2 * pnorm(-abs(pp / fitted$stderr))
    }

    if (model %in% c("Weibull","Poisson","Binomial","Gamma",
                     "LogLogistic","BinomialNeg","Bernoulli","Geometric",
                     "Gaussian_misp_Poisson","PoissonZIP","Gaussian_misp_PoissonZIP",
                     "BinomialNegZINB","PoissonZIP1","Gaussian_misp_PoissonZIP1",
                     "BinomialNegZINB1","Beta2","Kumaraswamy2","Beta","Kumaraswamy")) {
      if (!is.null(ff$sill)) ff$sill <- NULL
    }

    ### Set the output object:
    GeoFit <- list(
      anisopars = anisopars,
      bivariate = initparam$bivariate,
      claic = fitted$claic,
      clbic = fitted$clbic,
      coordx = initparam$coordx,
      coordy = initparam$coordy,
      coordz = initparam$coordz,
      coordt = initparam$coordt,
      coordx_dyn = coordx_dyn,
      conf.int = conf.int,
      convergence = fitted$convergence,
      copula = copula,
      corrmodel = corrmodel,
      data = initparam$data,
      distance = distance,
      est.aniso = est.aniso,
      fixed = ff,
      grid = grid,
      iterations = fitted$counts,
      likelihood = likelihood,
      logCompLik = fitted$value,
      lower = lower,
      message = fitted$message,
      model = model,
      n = initparam$n,
      ns = initparam$ns,
      numbetas = initparam$num_betas,
      numcoord = initparam$numcoord,
      numtime = initparam$numtime,
      optimizer = optimizer,
      param = as.list(fitted$par),
      p_neighb = p_neighb,
      nozero = initparam$setup$nozero,
      score = fitted$score,
      maxdist = maxdist,
      maxtime = maxtime,
      neighb = initparam$neighb,
      numpairs = initparam$numpairs,
      missp = missp,
      pvalues = pvalues,
      radius = radius,
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

    structure(c(GeoFit, call = call), class = c("GeoFit"))
  })
}
