
####################################################

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
    if (!is.null(X)) X <- covmatrix$X else X <- matrix(1, nrow = dimat, ncol = 1)
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
  if (!is.null(anisopars)) ccc <- GeoAniso(ccc, c(anisopars$angle, anisopars$ratio))
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
        if (is.null(covmatrix$copula)) corri <- cc else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- vv; M <- 0
      }
      ## SKEW GAUSSIAN
      if (covmatrix$model == 10) {
        vv <- as.numeric(covmatrix$param['sill']); sk <- as.numeric(covmatrix$param['skew']); sk2 <- sk^2
        if (is.null(covmatrix$copula)) {
          corr2 <- rho^2
          corri <- (2*sk2)*(sqrt(1-corr2) + rho*asin(rho)-1)/(pi*vv+sk2*(pi-2)) + (cc*vv)/(vv+sk2*(1-2/pi))
        } else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- vv + sk^2 * (1 - 2/pi); M <- sk*sqrt(2/pi)
      }
      ## SKEW LAPLACE
      if (covmatrix$model == 59) {
        vv <- as.numeric(covmatrix$param['sill']); sk <- as.numeric(covmatrix$param['skew'])
        if (is.null(covmatrix$copula)) {
          corri <- rho^2
        } else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- vv*(1/(sk^2) + 1/(1-sk)^2); M <- 1/sk + 1/(1-sk)
      }
      ## STUDENT T
      if (covmatrix$model == 12) {
        vv <- as.numeric(covmatrix$param['sill']); nu <- 1/as.numeric(covmatrix$param['df'])
        if (is.null(covmatrix$copula)) {
          cc2 <- cc^2; gamma_num <- gamma((nu - 1)/2)^2; gamma_den <- gamma(nu/2)^2; Fvals <- hypergeo::hypergeo(0.5, 0.5, nu/2, cc2)
          corri <- ((nu - 2) * gamma_num * Re(Fvals) * cc) / (2 * gamma_den)
        } else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- vv*nu/(nu-2); M <- 0
      }
      ## TUKey h
      if (covmatrix$model == 34) {
        vv <- as.numeric(covmatrix$param['sill']); h <- as.numeric(covmatrix$param['tail'])
        if (h > 0) {
          if (is.null(covmatrix$copula)) { corri <- (cc*(1-2*h)^(1.5))/((1-h)^2-(h*cc)^2)^(1.5) } else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        } else { if (is.null(covmatrix$copula)) corri <- cc else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance) }
        vvar <- vv*(1-2*h)^(-1.5); M <- 0
      }
      ## TUKey h2
      if (covmatrix$model == 40) {
        vv <- as.numeric(covmatrix$param['sill']); tail1 <- as.numeric(covmatrix$param['tail1']); tail2 <- as.numeric(covmatrix$param['tail2'])
        if (is.null(covmatrix$copula)) {
          hr <- tail1; hl <- tail2; y <- cc
          x1 <- 1 - (1 - y^2)*hr; y1 <- 1 - (1 - y^2)*hl
          x2 <- (1-hr)^2 - (y*hr)^2; y2 <- (1-hl)^2 - (y*hl)^2; g <- 1 - hl - hr + (1 - y^2)*hl*hr
          h1 <- sqrt(1 - y^2/x1^2) + (y/x1)*asin(y/x1)
          h2 <- sqrt(1 - y^2/y1^2) + (y/y1)*asin(y/y1)
          h3 <- sqrt(1 - y^2/(x1*y1)) + sqrt(y^2/(x1*y1))*asin(sqrt(y^2/(x1*y1)))
          p1 <- x1*h1/(2*pi*(x2)^(3/2)) + y/(4*(x2)^(3/2))
          p2 <- y1*h2/(2*pi*(y2)^(3/2)) + y/(4*(y2)^(3/2))
          p3 <- -(x1*y1)^(1/2)*h3/(2*pi*(g)^(3/2)) + y/(4*(g)^(3/2))
          mm <- (hr - hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
          vv1 <- 0.5*(1-2*hl)^(-3/2) + 0.5*(1-2*hr)^(-3/2) - (mm)^2
          corri <- (p1 + p2 + 2*p3 - mm^2)/vv1
          vvar <- vv*vv1; M <- sqrt(vv)*mm
        } else if (covmatrix$copula == "Gaussian") { corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance); vvar <- vv; M <- 0 }
      }
      ## SAS
      if (covmatrix$model == 20) {
        vv <- as.numeric(covmatrix$param['sill']); tail <- as.numeric(covmatrix$param['tail']); skew <- as.numeric(covmatrix$param['skew'])
        if (is.null(covmatrix$copula)) {
          d <- tail; e <- skew
          mm <- sinh(e/d)*exp(0.25)*(besselK(.25,(d+1)/(2*d))+besselK(.25,(1-d)/(2*d)))/(sqrt(8*pi))
          vv1 <- cosh(2*e/d)*exp(0.25)*(besselK(.25,(d+2)/(2*d))+besselK(0.25,(2-d)/(2*d)))/(sqrt(32*pi)) - 0.5 - mm^2
          corri <- corrsas(cc, e, d)
          vvar <- vv*vv1; M <- sqrt(vv)*mm
        } else { corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance); vvar <- vv; M <- 0 }
      }
      ## SKEW STUDENT T
      if (covmatrix$model == 18) {
        vv <- as.numeric(covmatrix$param['sill']); nu <- 1/as.numeric(covmatrix$param['df']); sk <- as.numeric(covmatrix$param['skew']); sk2 <- sk^2
        if (is.null(covmatrix$copula)) {
          w <- sqrt(1-sk2); y <- rho
          CorSkew <- (2*sk2/(pi*w*w+sk2*(pi-2)))*(sqrt(1-y*y)+y*asin(y)-1) + w*w*cc/(w*w+sk2*(1-2/pi))
          corri <- (pi*(nu-2)*gamma((nu-1)/2)^2/(2*(pi*gamma(nu/2)^2-sk2*(nu-2)*gamma((nu-1)/2)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,nu/2,y*y))*((1-2*sk2/pi)*CorSkew+2*sk2/pi)-2*sk2/pi)
        } else corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        D1 <- (nu-1)*0.5; D2 <- nu*0.5; mm <- sqrt(nu)*gamma(D1)*sk/(sqrt(pi)*gamma(D2))
        vvar <- vv*(nu/(nu-2)-mm^2); M <- sqrt(vv)*mm
      }
      ## TWO-PIECE t
      if (covmatrix$model == 27) {
        nu <- as.numeric(1/covmatrix$param['df']); sk <- as.numeric(covmatrix$param['skew']); sk2 <- sk^2; vv <- as.numeric(covmatrix$param['sill'])
        if (is.null(covmatrix$copula)) {
          corr2 <- cc^2
          a1 <- Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
          a2 <- cc*asin(cc) + (1-corr2)^(0.5)
          ll <- qnorm((1-sk)/2)
          p11 <- pbivnorm::pbivnorm(ll, ll, rho = rho, recycle = TRUE)
          a3 <- 3*sk2 + 2*sk + 4*p11 - 1
          KK <- (nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2)
          corri <- KK*(a1*a2*a3-4*sk2)
        } else corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        ttemp <- gamma(0.5*(nu-1))/gamma(0.5*nu)
        vvar <- vv*((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*ttemp^2); M <- - sqrt(vv)*2*sk*sqrt(nu/pi)*ttemp
      }
      ## TWO-PIECE Tukey h
      if (covmatrix$model == 38) {
        tail <- as.numeric(covmatrix$param['tail']); sk <- as.numeric(covmatrix$param['skew']); sk2 <- sk^2; vv <- as.numeric(covmatrix$param['sill'])
        if (is.null(covmatrix$copula)) {
          corr2 <- cc^2; gg2 <- (1-(1-corr2)*tail)^2; xx <- corr2/gg2
          A <- (asin(sqrt(xx))*sqrt(xx) + sqrt(1-xx))/(1-xx)^(1.5)
          ll <- qnorm((1-sk)/2)
          p11 <- pbivnorm::pbivnorm(ll, ll, rho = rho, recycle = TRUE)
          a3 <- 3*sk2 + 2*sk + 4*p11 - 1
          mm <- 8*sk2/(pi*(1-tail)^2)
          ff <- (1+3*sk2)/(1-2*tail)^(1.5)
          M  <- (2*(1-corr2)^(3/2))/(pi*gg2)
          corri <- (M*A*a3 - mm)/(ff - mm)
        } else corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- vv*(ff - mm); M <- -sqrt(vv)*2*sk*sqrt(2/pi)/(1-tail)
      }
      ## TWO-PIECE Gaussian
      if (covmatrix$model == 29) {
        vv <- as.numeric(covmatrix$param['sill']); sk <- as.numeric(nuisance['skew']); sk2 <- sk^2
        if (is.null(covmatrix$copula)) {
          corr2 <- sqrt(1-cc^2)
          ll <- qnorm((1-sk)/2)
          p11 <- pbivnorm::pbivnorm(ll, ll, rho = rho, recycle = TRUE)
          KK <- 3*sk2 + 2*sk + 4*p11 - 1
          corri <- (2*((corr2 + cc*asin(cc))*KK) - 8*sk2)/(3*pi*sk2 - 8*sk2 + pi)
        } else corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- vv*((1+3*sk2) - 8*sk2/pi); M <- -sqrt(vv)*2*sk*sqrt(2/pi)
      }
      ## BETA
      if (covmatrix$model == 28) {
        shape1 <- as.numeric(covmatrix$param['shape1']); shape2 <- as.numeric(covmatrix$param['shape2'])
        corr2 <- cc^2
        cc1 <- 0.5*(shape1 + shape2)
        vv  <- shape1*shape2/((cc1+1)*(shape1+shape2)^2)
        idx <- which(abs(corr2) > 1e-10); corr22 <- corr2[idx]; nu2 <- shape1/2; alpha2 <- shape2/2
        res <- 0; A <- 0; k <- 0
        while (k <= 100) {
          p1 <- 2*(lgamma(cc1+k)-lgamma(cc1)+lgamma(nu2+1+k)-lgamma(nu2+1))
          p2 <- lgamma(k+1)+(lgamma(nu2+k)-lgamma(nu2))+2*(lgamma(cc1+1+k)-lgamma(cc1+1))
          b1 <- p1 - p2; b2 <- log(hypergeo::genhypergeo(U=c(cc1+k,cc1+k,alpha2), L=c(cc1+k+1,cc1+k+1), polynomial=TRUE,maxiter=1000, z=corr22)); b3 <- k*log(corr22)
          sum <- exp(b1 + b2 + b3); res <- res + sum
          if (all(sum < 1e-6)) { break } else { A <- res }
          k <- k + 1
        }
        cc[idx] <- A; corri <- shape1*(cc1 + 1 ) * ((1-corr2)^(cc1) * cc - 1)/shape2; corri[-idx] <- 0
        ssup <- as.numeric((covmatrix$param["max"] - covmatrix$param["min"]))
        vvar <- (ssup^2)*vv
        M <- ssup*shape1/(shape1+shape2) + as.numeric(covmatrix$param["min"])
      }
      ## GAMMA
      if (covmatrix$model == 21) {
        sh <- as.numeric(covmatrix$param["shape"])
        if (is.null(covmatrix$copula)) { corri <- cc^2 } else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- 2/sh; M <- 0
      }
      ## WEIBULL
      if (covmatrix$model == 26) {
        sh <- as.numeric(covmatrix$param['shape'])
        if (is.null(covmatrix$copula)) {
          cc2 <- cc^2; c2 <- gamma(1 + 1/sh)^2; denom <- gamma(1 + 2/sh) - c2; Fz <- hypergeo::hypergeo(-1/sh, -1/sh, 1, cc2)
          corri <- (c2 / denom) * (Re(Fz) - 1)
        } else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- gamma(1+2/sh)/gamma(1+1/sh)^2 - 1; M <- 0
      }
      ## LOGLOGISTIC
      if (covmatrix$model == 24) {
        sh <- as.numeric(covmatrix$param['shape'])
        if (is.null(covmatrix$copula)) {
          cc2 <- cc^2 ; num <- pi * sin(2 * pi / sh); den <- 2 * sh * (sin(pi / sh))^2 - pi * sin(2 * pi / sh)
          factor <- num / den
          F1 <- hypergeo::hypergeo(-1/sh, -1/sh, 1, cc2)
          F2 <- hypergeo::hypergeo( 1/sh,  1/sh, 1, cc2)
          corri <- factor * (Re(F1) * Re(F2) - 1)
        } else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- 2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh)) - 1; M <- 0
      }
      ## LOGGAUSSIAN
      if (covmatrix$model == 22) {
        ss <- as.numeric(covmatrix$param['sill'])
        if (is.null(covmatrix$copula)) { corri <- (exp(ss*cc) - 1)/(exp(ss) - 1) } else if (covmatrix$copula == "Gaussian") corri <- gaussian_copula_cov_fast(cc, covmatrix$model, nuisance)
        vvar <- exp(ss) - 1; M <- 0
      }
    } # end !bivariate

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
      corri <- gaussian_copula_cov_fast(ccorr$corri, covmatrix$model, nuisance)
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


