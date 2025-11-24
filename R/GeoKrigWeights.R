GeoKrigWeights <- function(
  coordx, coordy = NULL, coordz = NULL, coordt = NULL,
  coordx_dyn = NULL, corrmodel, distance = "Eucl", grid = FALSE, loc,
  method = "cholesky", model = "Gaussian", n = 1, nloc = NULL,
  param, anisopars = NULL, radius = 1, sparse = FALSE, time = NULL,
  which = 1, copula = NULL, X = NULL, Xloc = NULL, Mloc=NULL
) {
  ############################################################################
  ## START - Preprocessing
  ############################################################################
  call <- match.call()

  ## Validazioni base
  if (is.null(CkModel(model))) stop("The name of the model is not correct\n")
  if (!is.character(corrmodel) || is.null(CkCorrModel(corrmodel))) stop("the name of the correlation model is wrong")
  corrmodel <- gsub("[[:blank:]]", "", corrmodel)
  model     <- gsub("[[:blank:]]", "", model)
  distance  <- gsub("[[:blank:]]", "", distance)
  method    <- gsub("[[:blank:]]", "", method)

  bivariate <- CheckBiv(CkCorrModel(corrmodel))
  spacetime <- CheckST(CkCorrModel(corrmodel))
  space <- !spacetime && !bivariate

  ## Checks su loc
  if (is.vector(loc)) loc <- matrix(loc, nrow = 1) else if (!is.matrix(loc)) loc <- as.matrix(loc)
  if (!(ncol(loc) == 2 || ncol(loc) == 3)) stop("loc parameter must be a matrix N X 2 or N X 3")

  if (!is.null(Xloc)) { if (is.vector(Xloc)) Xloc <- matrix(Xloc, nrow = 1) else if (!is.matrix(Xloc)) Xloc <- as.matrix(Xloc) }
  if (!is.null(X) && !is.matrix(X)) X <- as.matrix(X)

  if (CheckST(CkCorrModel(corrmodel))) { if (is.null(time)) stop("At least one temporal instants is needed for space-time kriging\n") }

  loc_orig <- loc
  if (!is.null(anisopars)) loc <- GeoAniso(loc, c(anisopars$angle, anisopars$ratio))

  ## Setup dimensioni
  if (is.null(time)) time <- 0
  numloc <- nrow(loc)
  tloc   <- length(time); if (!tloc) tloc <- 1

  if (ncol(loc) == 2) { locx <- loc[, 1]; locy <- loc[, 2]; locz <- double(length(loc[, 1])) }
  if (ncol(loc) == 3) { locx <- loc[, 1]; locy <- loc[, 2]; locz <- loc[, 3] }

  ############################################################################
  ## CALCOLO MATRICE DI COVARIANZA E PESI
  ############################################################################

  ## Gestione mean parametrizzate speciali
  logGausstemp <- SinhAsinhtemp <- Tukeyh2temp <- Tukeyhtemp <- FALSE
  Mtemp <- NULL
  if (model %in% c("Weibull", "Gamma", "LogLogistic", "LogGaussian")) {
    paramtemp <- param
    sel <- substr(names(param), 1, 4) == "mean"
    meantemp <- names(param[sel])
    param <- param[!sel]
    if (length(paramtemp$mean) > 1) Mtemp <- paramtemp$mean
    param$mean <- 0
    Xtemp <- X; X <- NULL
  }

  ## Costruzione matrice di covarianza
  covmatrix <- GeoCovmatrix(
    coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
    coordx_dyn = coordx_dyn, corrmodel = corrmodel, distance = distance,
    grid = grid, model = model, n = n, param = param, anisopars = anisopars,
    radius = radius, sparse = sparse, copula = copula, X = X
  )
  covmatrix$param <- unlist(covmatrix$param)
  
  if (bivariate) tloc <- 1
  spacetime_dyn <- !is.null(covmatrix$coordx_dyn)
  if (!spacetime_dyn) dimat <- covmatrix$numcoord * covmatrix$numtime else dimat <- sum(covmatrix$ns)
  dimat2 <- numloc * tloc

  if (is.null(X)) XX <- rep(1, dimat) else XX <- X
  MM <- NULL
  if (!is.null(Mtemp)) { MM <- Mtemp; param$mean <- 0 }
  else { if (length(param$mean) > 1) { MM <- param$mean; param$mean <- 0 } }
  
  if (model %in% c("Weibull", "Gamma", "LogLogistic", "LogGaussian")) {
    if (is.null(Xtemp)) X <- matrix(rep(1, dimat)) else X <- Xtemp
    param <- paramtemp
    if (!is.null(Mtemp)) param$mean <- 0
    covmatrix$namesnuis <- unique(c(meantemp, covmatrix$namesnuis))
  } else {
    ## =======================
    ## FIX 2: rispetta X utente
    ## =======================
    if (is.null(X)) X <- covmatrix$X else X <- as.matrix(X)
  }

  num_betas <- ncol(X); NS <- 0
  if (spacetime || bivariate) { NS <- cumsum(covmatrix$ns); NS <- c(0, NS)[-(length(covmatrix$ns) + 1)] }
  if (is.null(Xloc)) Xloc <- as.matrix(rep(1, dimat2)) else { if (spacetime_dyn) Xloc <- as.matrix(Xloc) }

  nuisance <- param[covmatrix$namesnuis]; nuisance <- Filter(Negate(is.null), nuisance)
  sel <- substr(names(nuisance), 1, 4) == "mean"
  betas <- as.numeric(nuisance[sel])
  if (length(betas) > 1 && is.null(X)) stop("Covariates matrix X is missing\n")
  if (bivariate) {
    sel1 <- substr(names(nuisance), 1, 6) == "mean_1"; betas1 <- as.numeric(nuisance[sel1])
    sel2 <- substr(names(nuisance), 1, 6) == "mean_2"; betas2 <- as.numeric(nuisance[sel2])
  }
  other_nuis <- as.numeric(nuisance[!sel])

  cmodel <- corrmodel; cdistance <- distance
  corrmodel_num <- CkCorrModel(covmatrix$corrmodel)
  distance_num  <- CheckDistance(covmatrix$distance)
  corrparam <- covmatrix$param[covmatrix$namescorr]
  if (bivariate && !(which == 1 || which == 2)) stop("which parameter must be 1 or 2")

  ## Costruzione coordinate stacked ccc
  if (!grid) {
    if (!spacetime && !bivariate) {
      if (is.null(covmatrix$coordz))  ccc <- cbind(covmatrix$coordx, covmatrix$coordy, 0)
      else                            ccc <- cbind(covmatrix$coordx, covmatrix$coordy, covmatrix$coordz)
    }
  }
  if (grid) {
    if (is.null(coordz)) { ccc <- as.matrix(expand.grid(covmatrix$coordx, covmatrix$coordy)); ccc <- cbind(ccc, 0) }
    else                 { ccc <- as.matrix(expand.grid(covmatrix$coordx, covmatrix$coordy, covmatrix$coordz)) }
    grid <- FALSE
  } else {
    if ((spacetime || bivariate) && (!spacetime_dyn)) {
      if (!is.null(covmatrix$coordz)) ccc <- cbind(rep(covmatrix$coordx, covmatrix$numtime), rep(covmatrix$coordy, covmatrix$numtime), rep(covmatrix$coordz, covmatrix$numtime))
      else                            ccc <- cbind(rep(covmatrix$coordx, covmatrix$numtime), rep(covmatrix$coordy, covmatrix$numtime), 0)
    }
    if ((spacetime || bivariate) && (spacetime_dyn)) {
      ccc <- do.call(rbind, args = c(covmatrix$coordx_dyn)); if (ncol(ccc) == 2) ccc <- cbind(ccc, 0)
    }
  }

  ## =======================
  ## FIX 1: anisotropia su ccc DOPO la costruzione
  ## =======================
  if (!is.null(anisopars)) ccc <- GeoAniso(ccc, c(anisopars$angle, anisopars$ratio))

  ############################################################################
  ## CALCOLO CORRELAZIONI E COSTRUZIONE CC - MODELLI CONTINUI
  ############################################################################
  if (covmatrix$model %in% c(1,10,59,12,18,21,26,24,22,27,28,29,38,39,9,34,40,20)) {
    cc_call <- dotCall64::.C64(
      "Corr_c",
      SIGNATURE = c("double","double","double","double","double","integer",
                    "integer","double","double","double","integer","integer",
                    "integer","integer","integer","integer","double","integer",
                    "integer","double","integer","integer","double"),
      corri = dotCall64::vector_dc("double", dimat * dimat2),
      ccc[,1], ccc[,2], ccc[,3], covmatrix$coordt, corrmodel_num, 0,
      locx, locy, locz, covmatrix$numcoord, numloc, tloc,
      covmatrix$ns, NS, covmatrix$numtime, corrparam,
      covmatrix$spacetime, covmatrix$bivariate, time, distance_num, which-1, covmatrix$radius,
      INTENT = c("w", rep("r", 22)), NAOK = TRUE, PACKAGE = "GeoModels"
    )

    ## --- trasformazioni modello-specifiche su corri + calcolo vvar ---
    if (bivariate) {
      corri <- cc_call$corri
      vvar <- if (which == 1) covmatrix$param["sill_1"] + covmatrix$param["nugget_1"] else covmatrix$param["sill_2"] + covmatrix$param["nugget_2"]
      M <- 0
    } else {
      rho <- cc_call$corri
      cc  <- (1 - as.numeric(covmatrix$param["nugget"])) * rho
      ## GAUSSIAN
      if (covmatrix$model == 1) {
        vv <- as.numeric(covmatrix$param['sill'])
        if (is.null(covmatrix$copula)) corri <- cc else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov(cc, covmatrix$model, nuisance)
        vvar <- vv; M <- 0
      }
      ## ... (resto del tuo codice invariato) ...
    }

    ## --- COSTRUISCO CC per i pesi ---
    if (is.null(covmatrix$copula)) CC <- matrix(corri * vvar, nrow = dimat, ncol = dimat2) else CC <- matrix(corri, nrow = dimat, ncol = dimat2)
  }

  ############################################################################
  ## CALCOLO CORRELAZIONI E COSTRUZIONE CC - MODELLI DISCRETI
  ############################################################################
  if (covmatrix$model %in% c(2,11,14,19,30,36,16,43,44,45,46,47,57,58)) {
    if (!bivariate) {
      if (is.null(MM)) { mu <- X %*% betas; mu0 <- Xloc %*% betas } else { mu <- MM; mu0 <- Mloc }
    } else {
      X11 <- X[1:covmatrix$ns[1], ]; X22 <- X[(covmatrix$ns[1]+1):(covmatrix$ns[1]+covmatrix$ns[2]), ]
      mu <- c(X11 %*% betas1, X22 %*% betas2)
    }
    KK <- 0
    if (covmatrix$model %in% c(2, 11)) {
      if (is.null(nloc)) { KK <- rep(round(mean(n)), dimat2) }
      else if (is.numeric(nloc)) { KK <- if (length(nloc) == 1) rep(nloc, dimat2) else nloc }
      if (length(n) == 1) n <- rep(n, dimat)
      if (length(KK) != dimat2) stop("dimension of nloc is wrong\n")
    }
    if (covmatrix$model == 19) KK <- min(nloc)
    if (covmatrix$model %in% c(16,45)) KK <- n

    cop <- 0
    if (!is.null(covmatrix$copula) && covmatrix$copula == "Gaussian") cop <- 1

    ccorr <- dotCall64::.C64(
      'Corr_c_bin',
      SIGNATURE = c(rep("double",5),"integer","integer","double","double","double",
                    rep("integer",9),"double","double","double","integer","integer",
                    "double","integer","integer","double","integer"),
      corri = dotCall64::vector_dc("double", dimat*dimat2),
      ccc[,1], ccc[,2], ccc[,3], covmatrix$coordt, corrmodel_num, as.integer(0),
      locx, locy, locz, covmatrix$numcoord, numloc, covmatrix$model, tloc,
      KK, n, covmatrix$ns, NS, covmatrix$numtime,
      rep(c(mu), dimat2), other_nuis, corrparam, covmatrix$spacetime,
      bivariate, time, distance_num, as.integer(which-1), covmatrix$radius, cop,
      INTENT = c("w", rep("r", 28)), PACKAGE = 'GeoModels', VERBOSE = 0, NAOK = TRUE
    )

    if (cop == 1) {
      nuisance$n <- n[1]
      corri <- gaussian_copula_cov(ccorr$corri, covmatrix$model, nuisance)
    } else corri <- ccorr$corri

    CC <- matrix(corri, nrow = dimat, ncol = dimat2)   # per discreti non moltiplichiamo per vvar qui
  }

  ############################################################################
  ## CALCOLO PESI FINALI
  ############################################################################
  weights_result <- getInvC(covmatrix, CC, mse = FALSE)
  krig_weights <- weights_result$a

  ############################################################################
  ## RETURN
  ############################################################################
  result <- list(
    weights = krig_weights,
    model = model,
    corrmodel = cmodel,
    coordx = covmatrix$coordx,
    coordy = covmatrix$coordy,
    coordz = covmatrix$coordz,
    bivariate = bivariate,
    spacetime = covmatrix$spacetime,
    param = covmatrix$param,
    covmatrix = covmatrix$covmatrix,
    xloc = Xloc,
    tloc = tloc,
    CC = CC
  )
  
  return(result)
}
