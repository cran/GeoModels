
####################################################
### File name: GeoLik_original_style_fixed.R
####################################################

### Optim call for log-likelihood maximization
Lik <- function(copula,bivariate,coordx,coordy,coordz,coordt,coordx_dyn,corrmodel,data,fixed,flagcor,flagnuis,grid,lower,
                mdecomp,model,namescorr,namesnuis,namesparam,numcoord,numpairs,numparamcor,numtime,
                optimizer,onlyvar,param,radius,setup,spacetime,sparse,varest,taper,type,upper,ns,X,neighb,MM,aniso,score)
{
  #########################################################################
  ## Utilities
  #########################################################################

  combine_param <- function(param, fixed, namesparam) {
    names(param) <- namesparam
    if (!is.null(fixed) && length(fixed) > 0) fixed <- unlist(fixed, use.names = TRUE)
    out <- c(param, fixed)
    out <- as.numeric(out)
    names(out) <- c(namesparam, names(fixed))
    out
  }

  get_mean <- function(X, nuisance, MM = NULL) {
    if (!is.null(MM)) return(as.numeric(MM))
    sel <- substr(names(nuisance), 1, 4) == "mean"
    mm <- as.numeric(nuisance[sel])
    as.vector(X %*% mm)
  }

  build_correlation_matrix <- function(corr_vec, ident) {
    corrmat_full <- ident
    lt <- lower.tri(ident)
    corrmat_full[lt] <- corr_vec
    corrmat_full <- t(corrmat_full)
    corrmat_full[lt] <- corr_vec
    corrmat_full
  }

  get_rhobeg <- function(param, lower, upper) {
    width <- upper - lower
    finite_width <- width[is.finite(width) & width > 0]
    if (length(finite_width) > 0) {
      min(0.1 * min(finite_width), 0.2 * max(abs(param), 1))
    } else {
      0.2 * max(abs(param), 1)
    }
  }

  make_optimize_bounds <- function(param, lower, upper, namesparam) {
    p0 <- as.numeric(param)[1]
    pname <- namesparam[1]

    ll <- if (is.null(lower)) -Inf else as.numeric(lower)[1]
    uu <- if (is.null(upper))  Inf else as.numeric(upper)[1]

    span <- max(10 * (abs(p0) + 1), 10)
    if (startsWith(pname, "mean")) span <- max(100 * (abs(p0) + 1), 100)

    positive_par <- pname %in% c(
      "scale", "scale_s", "scale_t", "sill", "nugget", "power", "power2",
      "shape", "shape1", "shape2", "df", "tail", "tail1", "tail2",
      "min", "max", "rate", "ratio"
    )

    if (!is.finite(ll)) {
      if (positive_par) ll <- .Machine$double.eps else ll <- p0 - span
    }

    if (!is.finite(uu)) {
      if (pname %in% c("scale", "scale_s", "scale_t")) {
        uu <- 4
      } else if (pname == "nugget") {
        uu <- 1 - .Machine$double.eps
      } else if (pname == "angle") {
        uu <- pi
      } else {
        uu <- p0 + span
      }
    }

    if (p0 <= ll) ll <- min(.Machine$double.eps, p0 - span)
    if (p0 >= uu) uu <- p0 + span

    if (!is.finite(ll) || !is.finite(uu) || ll >= uu) {
      stop("invalid lower/upper bounds for optimize", call. = FALSE)
    }

    list(lower = ll, upper = uu)
  }

  make_bobyqa_setup <- function(param, lower, upper, namesparam) {
    p0 <- as.numeric(param)
    ll <- as.numeric(lower)
    uu <- as.numeric(upper)

    if (length(ll) != length(p0) || length(uu) != length(p0)) {
      stop("lower and upper must have the same length as param for bobyqa", call. = FALSE)
    }

    names(ll) <- namesparam
    names(uu) <- namesparam

    span <- pmax(10 * (abs(p0) + 1), 10)
    mean_par <- startsWith(namesparam, "mean")
    span[mean_par] <- pmax(100 * (abs(p0[mean_par]) + 1), 100)

    idx <- is.infinite(ll) & ll < 0
    ll[idx] <- p0[idx] - span[idx]

    idx <- is.infinite(uu) & uu > 0
    uu[idx] <- p0[idx] + span[idx]

    ## standard finite fallbacks for common constrained parameters
    idx <- !is.finite(ll) & namesparam %in% c("scale", "scale_s", "scale_t", "sill", "power", "power2", "shape", "tail", "tail1", "tail2", "ratio")
    ll[idx] <- .Machine$double.eps

    idx <- !is.finite(uu) & namesparam %in% c("scale", "scale_s", "scale_t")
    uu[idx] <- 4
    idx <- !is.finite(uu) & namesparam == "nugget"
    uu[idx] <- 1 - .Machine$double.eps
    idx <- !is.finite(uu) & namesparam == "angle"
    uu[idx] <- pi

    if (any(!is.finite(ll)) || any(!is.finite(uu))) {
      stop("bobyqa requires finite lower and upper bounds", call. = FALSE)
    }

    idx <- p0 <= ll
    if (any(idx)) ll[idx] <- p0[idx] - span[idx]
    idx <- p0 >= uu
    if (any(idx)) uu[idx] <- p0[idx] + span[idx]

    if (any(uu <= ll)) {
      stop("bobyqa requires upper > lower for all parameters", call. = FALSE)
    }

    width <- uu - ll
    finite_width <- width[is.finite(width) & width > 0]
    if (length(finite_width) > 0) {
      rhobeg <- min(0.1 * min(finite_width), 0.2 * max(abs(p0), 1))
    } else {
      rhobeg <- 0.2 * max(abs(p0), 1)
    }
    rhobeg <- max(rhobeg, .Machine$double.eps^0.25)
    rhoend <- rhobeg * 1e-6
    npar <- length(p0)

    list(
      lower = ll,
      upper = uu,
      control = list(
        maxfun = 100000L,
        rhobeg = rhobeg,
        rhoend = rhoend,
        npt = min(2 * npar + 1, (npar + 1) * (npar + 2) / 2),
        iprint = 0
      )
    )
  }

  #########################################################################
  ## C calls for correlations
  #########################################################################

  matr <- function(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
  {
    cc <- dotCall64::.C64(
      as.character(corrmat),
      SIGNATURE = c(rep("double",5),"integer","double","double","double","integer","integer"),
      cr = dotCall64::vector_dc("double", length(corr)),
      coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, radius, ns, NS,
      INTENT = c("w", rep("r", 10)),
      PACKAGE = 'GeoModels', VERBOSE = 0, NAOK = TRUE
    )$cr
    cc
  }

  matr2 <- function(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius,model,mu)
  {
    cc <- dotCall64::.C64(
      as.character(corrmat),
      SIGNATURE = c(rep("double",5),"integer","double","integer","double","double","double","integer","integer","integer"),
      cr = dotCall64::vector_dc("double", length(corr)),
      coordx, coordy, coordz, coordt, corrmodel, c(mu), ns, nuisance, paramcorr, radius, ns, NS, model,
      INTENT = c("rw", rep("r", 13)),
      PACKAGE = 'GeoModels', VERBOSE = 0, NAOK = TRUE
    )$cr
    cc
  }

  #########################################################################
  ## Densities / criteria
  #########################################################################

  LogNormDenRestr <- function(const,cova,dimat,mdecomp,stdata,...) {
    llik <- 1.0e8
    decompvarcov <- MatDecomp(cova, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)
    ivarcov <- MatInv(cova)
    sumvarcov <- sum(ivarcov)
    if (!is.finite(sumvarcov) || sumvarcov <= 0) return(llik)
    p <- ivarcov - array(rowSums(ivarcov), c(dimat, 1)) %*% colSums(ivarcov) / sumvarcov
    llik <- 0.5 * (const + logdetvarcov + log(sumvarcov) + crossprod(t(crossprod(stdata, p)), stdata))
    as.numeric(llik)
  }

  LogNormDenRestr_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata) {
    llik <- 1.0e8
    ident[lower.tri(ident, diag = TRUE)] <- cova
    ident <- t(ident)
    ident[lower.tri(ident, diag = TRUE)] <- cova
    decompvarcov <- MatDecomp(ident, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)
    ivarcov <- MatInv(ident)
    sumvarcov <- sum(ivarcov)
    if (!is.finite(sumvarcov) || sumvarcov <= 0) return(llik)
    p <- ivarcov - array(rowSums(ivarcov), c(dimat, 1)) %*% colSums(ivarcov) / sumvarcov
    llik <- 0.5 * (const + logdetvarcov + log(sumvarcov) + crossprod(t(crossprod(stdata, p)), stdata))
    as.numeric(llik)
  }

  LogNormDenStand <- function(const, cova, dimat, mdecomp, stdata, ...) {
    llik <- 1.0e8
    decompvarcov <- MatDecomp(cova, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)
    y <- try(forwardsolve(decompvarcov, stdata, transpose = FALSE), silent = TRUE)
    if (inherits(y, "try-error")) return(llik)
    quad <- sum(y^2)
    if (!is.finite(quad)) return(llik)
    0.5 * (const + logdetvarcov + quad)
  }

  LogNormDenStand_spam <- function(const, cova, dimat, mdecomp, stdata, ...) {
    llik <- 1.0e8
    mcov <- try(spam::as.spam(cova), silent = TRUE)
    if (inherits(mcov, "try-error")) return(llik)
    cholS <- try(spam::chol.spam(mcov), silent = TRUE)
    if (inherits(cholS, "try-error")) return(llik)
    llik <- 0.5 * (const + 2 * c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) +
                     sum(stdata * spam::solve.spam(cholS, stdata)))
    as.numeric(llik)
  }

  LogNormDenStand22 <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata) {
    llik <- 1.0e8
    decompvarcov <- MatDecomp(cova, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)
    y <- try(forwardsolve(decompvarcov, stdata, transpose = FALSE), silent = TRUE)
    if (inherits(y, "try-error")) return(llik)
    llik <- 0.5 * (const + logdetvarcov + sum(y^2))
    as.numeric(llik)
  }

  LogNormDenStand_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata) {
    llik <- 1.0e8
    ident[lower.tri(ident, diag = TRUE)] <- cova
    ident <- t(ident)
    ident[lower.tri(ident, diag = TRUE)] <- cova
    decompvarcov <- MatDecomp(ident, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)
    y <- try(forwardsolve(decompvarcov, stdata, transpose = FALSE), silent = TRUE)
    if (inherits(y, "try-error")) return(llik)
    llik <- 0.5 * (const + logdetvarcov + sum(y^2))
    as.numeric(llik)
  }

  LogNormDenStand_biv_spam <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata) {
    llik <- 1.0e8
    ident[lower.tri(ident, diag = TRUE)] <- cova
    ident <- t(ident)
    ident[lower.tri(ident, diag = TRUE)] <- cova
    mcov <- try(spam::as.spam(ident), silent = TRUE)
    if (inherits(mcov, "try-error")) return(llik)
    cholS <- try(chol(mcov), silent = TRUE)
    if (inherits(cholS, "try-error")) return(llik)
    llik <- 0.5 * (const + 2 * c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) +
                     sum(stdata * spam::solve.spam(cholS, stdata)))
    as.numeric(llik)
  }

  CVV <- function(const,cova,dimat,mdecomp,stdata,...) {
    llik <- 1.0e8
    decompvarcov <- MatDecomp(cova, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    invarcov <- MatInv(cova)
    D <- diag(1 / diag(invarcov))
    M <- crossprod(invarcov, D)
    C <- tcrossprod(M, M)
    as.numeric(mean(crossprod(t(crossprod(stdata, C)), stdata)))
  }

  CVV_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata) {
    llik <- 1.0e8
    ident[lower.tri(ident, diag = TRUE)] <- cova
    ident <- t(ident)
    ident[lower.tri(ident, diag = TRUE)] <- cova
    decompvarcov <- MatDecomp(ident, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    invarcov <- MatInv(ident)
    D <- diag(1 / diag(invarcov))
    M <- crossprod(invarcov, D)
    C <- tcrossprod(M, M)
    as.numeric(mean(crossprod(t(crossprod(stdata, C)), stdata)))
  }

  ## Tapering functions are kept for compatibility with type 5/6.
  LogNormDenTap <- function(const,cova,dimat,mdecomp,stdata,setup,...) {
    lliktap <- 1.0e8
    varcovtap <- try(new("spam", entries = cova * setup$taps, colindices = setup$ja,
                         rowpointers = setup$ia, dimension = as.integer(rep(dimat, 2))), silent = TRUE)
    if (inherits(varcovtap, "try-error")) return(lliktap)
    cholvctap <- try(spam::update.spam.chol.NgPeyton(setup$struct, varcovtap), silent = TRUE)
    if (inherits(cholvctap, "try-error")) return(lliktap)
    logdet <- c(spam::determinant(cholvctap)$modulus)
    inv <- spam::solve.spam(cholvctap)
    slot(varcovtap, "entries") <- inv[setup$idx] * setup$taps
    as.numeric(0.5 * (const + 2 * logdet + drop(t(stdata) %*% varcovtap %*% stdata)))
  }

  LogNormDenTap1 <- function(const,cova,dimat,mdecomp,stdata,setup,...) {
    lliktap <- 1.0e8
    varcovtap <- try(new("spam", entries = cova * setup$taps, colindices = setup$ja,
                         rowpointers = setup$ia, dimension = as.integer(rep(dimat, 2))), silent = TRUE)
    if (inherits(varcovtap, "try-error")) return(lliktap)
    cholvctap <- try(spam::update.spam.chol.NgPeyton(setup$struct, varcovtap), silent = TRUE)
    if (inherits(cholvctap, "try-error")) return(lliktap)
    logdet <- c(spam::determinant(cholvctap)$modulus)
    as.numeric(0.5 * (const + 2 * logdet + sum(stdata * spam::solve.spam(cholvctap, stdata))))
  }

  LogNormDenTap_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata) {
    lliktap <- 1.0e8
    varcovtap <- try(new("spam", entries = cova * setup$taps, colindices = setup$ja,
                         rowpointers = setup$ia, dimension = as.integer(rep(dimat, 2))), silent = TRUE)
    if (inherits(varcovtap, "try-error")) return(lliktap)
    cholvctap <- try(spam::chol.spam(varcovtap), silent = TRUE)
    if (inherits(cholvctap, "try-error")) return(lliktap)
    logdet <- c(spam::determinant(cholvctap)$modulus)
    inv <- spam::solve.spam(cholvctap)
    slot(varcovtap, "entries") <- inv[setup$idx] * setup$taps
    as.numeric(0.5 * (const + 2 * logdet + drop(t(stdata) %*% varcovtap %*% stdata)))
  }

  LogNormDenTap1_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata) {
    lliktap <- 1.0e8
    cova[cova == nuisance['sill']] <- nuisance['sill'] + nuisance['nugget']
    varcovtap <- try(new("spam", entries = cova * setup$taps, colindices = setup$ja,
                         rowpointers = setup$ia, dimension = as.integer(rep(dimat, 2))), silent = TRUE)
    if (inherits(varcovtap, "try-error")) return(lliktap)
    cholvctap <- try(spam::update.spam.chol.NgPeyton(setup$struct, varcovtap), silent = TRUE)
    if (inherits(cholvctap, "try-error")) return(lliktap)
    logdet <- c(spam::determinant.spam.chol.NgPeyton(cholvctap)$modulus)
    as.numeric(0.5 * (const + 2 * logdet + sum(stdata * spam::solve.spam(cholvctap, stdata))))
  }

  LogshDenTap1 <- function(const,cova,dimat,mdecomp,nuisance,setup,sill,stdata,...) {
    lliktap <- 1.0e8
    varcovtap <- try(new("spam", entries = setup$taps * cova, colindices = setup$ja,
                         rowpointers = setup$ia, dimension = as.integer(rep(dimat, 2))), silent = TRUE)
    if (inherits(varcovtap, "try-error")) return(lliktap)
    cholvctap <- try(spam::update.spam.chol.NgPeyton(setup$struct, varcovtap), silent = TRUE)
    if (inherits(cholvctap, "try-error")) return(lliktap)
    logdet <- 2 * c(spam::determinant(cholvctap)$modulus)
    skew <- as.numeric(nuisance["skew"])
    delta <- as.numeric(nuisance["tail"])
    Z <- sinh(delta * asinh(stdata) - skew)
    C <- delta * sqrt((1 + Z^2) / (stdata^2 + 1))
    as.numeric(0.5 * (const + const * log(sill) / log(2*pi) + logdet - 2 * sum(log(C)) +
                        sum(Z * spam::solve.spam(cholvctap, Z))))
  }

  #########################################################################
  ## Objective functions
  #########################################################################

  loglik <- function(param, const, coordx, coordy, coordz, coordt, corr, corrmat, corrmodel, data, dimat, fixed, fname,
                     grid, ident, mdecomp, model, namescorr, namesnuis, namesparam, radius, setup, X, ns, NS, MM, aniso, namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance <- pram[namesnuis]
    Mean <- get_mean(X, nuisance, MM)

    if (aniso) {
      anisopar <- pram[namesaniso]
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(anisopar))
      coordx <- coords1[, 1]
      coordy <- coords1[, 2]
      coordz <- if (ncol(coords1) == 3) coords1[, 3] else NULL
    }

    sill <- as.numeric(nuisance['sill'])
    nugget <- as.numeric(nuisance['nugget'])
    if (!is.finite(sill) || !is.finite(nugget) || sill <= 0 || nugget < 0 || nugget >= 1) return(llik)

    if (fname %in% c("LogNormDenTap", "LogNormDenTap1")) {
      corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, ns, NS, radius)
      if (length(corrv) == 0 || any(!is.finite(corrv))) return(llik)
      cova_vec <- corrv * sill * (1 - nugget)
      return(do.call(what = fname, args = list(stdata = data - c(Mean), const = const, cova = cova_vec,
                                               dimat = dimat, mdecomp = mdecomp, setup = setup,
                                               nuisance = nuisance, sill = sill)))
    }

    corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, ns, NS, radius)
    if (length(corrv) == 0 || any(!is.finite(corrv))) return(llik)
    corrmat_full <- build_correlation_matrix(corrv, ident)
    cova <- sill * ((1 - nugget) * corrmat_full + nugget * ident)
    stdata <- data - c(Mean)

    do.call(what = fname, args = list(stdata = stdata, const = const, cova = cova,
                                      dimat = dimat, mdecomp = mdecomp, setup = setup))
  }

  loglik_biv <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                         grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance <- pram[namesnuis]
    sel1 <- substr(names(nuisance), 1, 6) == "mean_1"
    sel2 <- substr(names(nuisance), 1, 6) == "mean_2"
    mm1 <- as.numeric(nuisance[sel1])
    mm2 <- as.numeric(nuisance[sel2])
    X1 <- as.matrix(X[1:ns[1], ])
    X2 <- as.matrix(X[(ns[1] + 1):(ns[2] + ns[1]), ])
    Mean <- as.double(c(X1 %*% mm1, X2 %*% mm2))
    if (!is.null(MM)) Mean <- as.numeric(MM)
    stdata <- data - Mean
    corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, ns, NS, radius)
    do.call(what = fname, args = list(stdata = stdata, const = const, cova = corrv, ident = ident,
                                      dimat = dimat, mdecomp = mdecomp, nuisance = nuisance, setup = setup))
  }

  loglik_miss_T <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                            grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance <- pram[namesnuis]
    Mean <- get_mean(X, nuisance, MM)

    if (aniso) {
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(pram[namesaniso]))
      coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- if (ncol(coords1) == 3) coords1[,3] else NULL
    }

    corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, ns, NS, radius)
    df <- 1 / as.numeric(nuisance['df'])
    sill <- as.numeric(nuisance['sill'])
    nugget <- as.numeric(nuisance['nugget'])
    if (!is.finite(df) || df <= 2 || sill <= 0 || nugget < 0 || nugget >= 1) return(llik)

    corrv <- try(exp(log(df - 2) + 2 * lgamma(0.5 * (df - 1)) -
                       (log(2) + 2 * lgamma(df / 2)) +
                       log(Re(hypergeo::hypergeo(0.5, 0.5, df / 2, corrv^2))) + log(abs(corrv))) * sign(corrv),
                 silent = TRUE)
    if (inherits(corrv, "try-error") || any(!is.finite(corrv))) return(llik)
    corrmat_full <- build_correlation_matrix(corrv, ident)
    cova <- sill * ((1 - nugget) * corrmat_full + nugget * ident)
    LogNormDenStand(stdata = data - c(Mean), const = const, cova = cova, dimat = dimat, mdecomp = mdecomp)
  }

  loglik_miss_skewT <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                                grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance <- pram[namesnuis]
    Mean <- get_mean(X, nuisance, MM)

    if (aniso) {
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(pram[namesaniso]))
      coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- if (ncol(coords1) == 3) coords1[,3] else NULL
    }

    corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, ns, NS, radius)
    nu <- 1 / as.numeric(nuisance['df'])
    skew <- as.numeric(nuisance['skew'])
    sill <- as.numeric(nuisance['sill'])
    nugget <- as.numeric(nuisance['nugget'])
    if (!is.finite(nu) || nu <= 2 || abs(skew) >= 1 || sill <= 0 || nugget < 0 || nugget >= 1) return(llik)

    skew2 <- skew^2
    w <- sqrt(1 - skew2)
    KK <- 2 * skew2 / pi
    D1 <- (nu - 1) / 2
    D2 <- nu / 2
    CorSkew <- (2 * skew2 / (pi * w^2 + skew2 * (pi - 2))) *
      (sqrt(1 - corrv^2) + corrv * asin(corrv) - 1) +
      w^2 * corrv / (w^2 + skew2 * (1 - 2 / pi))
    corr3 <- (pi * (nu - 2) * gamma(D1)^2 /
                (2 * (pi * gamma(D2)^2 - skew2 * (nu - 2) * gamma(D1)^2))) *
      (Re(hypergeo::hypergeo(0.5, 0.5, D2, corrv^2)) * ((1 - KK) * CorSkew + KK) - KK)
    if (any(!is.finite(corr3))) return(llik)
    corrmat_full <- build_correlation_matrix(corr3, ident)
    cova <- sill * ((1 - nugget) * corrmat_full + nugget * ident)
    LogNormDenStand(stdata = data - c(Mean), const = const, cova = cova, dimat = dimat, mdecomp = mdecomp)
  }

  loglik_miss_Poisgamma <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                                    grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance0 <- pram[namesnuis]
    Mean <- get_mean(X, nuisance0, MM)
    sel <- substr(names(nuisance0), 1, 4) == "mean"
    nuisance <- nuisance0[!sel]
    nuisance <- nuisance[order(names(nuisance))]

    if (aniso) {
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(pram[namesaniso]))
      coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- if (ncol(coords1) == 3) coords1[,3] else NULL
    }

    if (as.numeric(nuisance['nugget']) < 0 || as.numeric(nuisance['nugget']) > 1) return(llik)
    mu <- Mean
    corrv <- matr2(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, ns, NS, radius, 46, mu)
    cova <- ident
    cova[lower.tri(cova)] <- corrv
    cova <- t(cova)
    cova[lower.tri(cova)] <- corrv
    ff <- exp(mu)
    cova[cbind(seq_len(dimat), seq_len(dimat))] <- ff * (1 + ff / as.numeric(nuisance["shape"]))
    LogNormDenStand22(stdata = data - c(ff), const = const, cova = cova, dimat = dimat,
                      ident = ident, mdecomp = mdecomp, nuisance = nuisance, setup = setup)
  }

  loglik_miss_Pois <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                               grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance0 <- pram[namesnuis]
    Mean <- get_mean(X, nuisance0, MM)
    sel <- substr(names(nuisance0), 1, 4) == "mean"
    nuisance <- nuisance0[!sel]

    if (aniso) {
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(pram[namesaniso]))
      coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- if (ncol(coords1) == 3) coords1[,3] else NULL
    }

    if (as.numeric(nuisance['nugget']) < 0 || as.numeric(nuisance['nugget']) > 1) return(llik)
    mu <- Mean
    corrv <- matr2(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, ns, NS, radius, 30, mu)
    cova <- ident
    cova[lower.tri(cova)] <- corrv
    cova <- t(cova)
    cova[lower.tri(cova)] <- corrv
    diag(cova) <- exp(mu)
    LogNormDenStand22(stdata = data - c(exp(mu)), const = const, cova = cova, dimat = dimat,
                      ident = ident, mdecomp = mdecomp, nuisance = nuisance, setup = setup)
  }

  LogNormDenStand_Tukey2H <- function(const, cova, dimat, mdecomp, nuisance, sill, setup, stdata) {
    llik <- 1.0e8
    delta1 <- as.numeric(nuisance["tail1"])
    delta2 <- as.numeric(nuisance["tail2"])
    if (delta1 <= 0 || delta1 >= 0.5 || delta2 <= 0 || delta2 >= 0.5 || sill <= 0) return(llik)

    decompvarcov <- MatDecomp(cova, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)

    pos_idx <- which(stdata >= 0)
    neg_idx <- which(stdata < 0)
    tau_inv <- numeric(dimat)
    log_jacobian <- 0

    if (length(pos_idx) > 0) {
      st_pos <- stdata[pos_idx]
      W_pos <- VGAM::lambertW(delta1 * st_pos^2)
      g_pos <- sqrt(W_pos / delta1)
      jac_pos <- g_pos / (st_pos * (1 + W_pos))
      tau_inv[pos_idx] <- g_pos
      log_jacobian <- log_jacobian + sum(log(abs(jac_pos)))
    }
    if (length(neg_idx) > 0) {
      st_neg <- stdata[neg_idx]
      W_neg <- VGAM::lambertW(delta2 * st_neg^2)
      g_neg <- sqrt(W_neg / delta2)
      jac_neg <- -g_neg / (st_neg * (1 + W_neg))
      tau_inv[neg_idx] <- -g_neg
      log_jacobian <- log_jacobian + sum(log(abs(jac_neg)))
    }

    y <- try(forwardsolve(decompvarcov, tau_inv, transpose = FALSE), silent = TRUE)
    if (inherits(y, "try-error")) return(llik)
    quad_form <- sum(y^2)
    if (!is.finite(quad_form) || !is.finite(log_jacobian)) return(llik)
    0.5 * (const + dimat * log(sill) + logdetvarcov + quad_form) - log_jacobian
  }

  loglik_tukey2h <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                             grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance <- pram[namesnuis]
    Mean <- get_mean(X, nuisance, MM)

    if (aniso) {
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(pram[namesaniso]))
      coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- if (ncol(coords1) == 3) coords1[,3] else NULL
    }

    sill <- as.numeric(nuisance['sill'])
    nugget <- as.numeric(nuisance['nugget'])
    tail1 <- as.numeric(nuisance['tail1'])
    tail2 <- as.numeric(nuisance['tail2'])
    if (tail1 <= 0 || tail1 >= 0.5 || tail2 <= 0 || tail2 >= 0.5 || nugget < 0 || nugget >= 1 || sill <= 0) return(llik)

    nuisance_temp <- nuisance
    nuisance_temp['sill'] <- 1
    nuisance_temp['nugget'] <- 0
    corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance_temp, paramcorr, ns, NS, radius)
    if (any(!is.finite(corrv))) return(llik)
    corrmat_full <- build_correlation_matrix(corrv, ident)
    cova <- (1 - nugget) * corrmat_full + nugget * ident
    stdata <- (data - c(Mean)) / sqrt(sill)
    LogNormDenStand_Tukey2H(const = const, cova = cova, dimat = dimat, mdecomp = mdecomp,
                            nuisance = nuisance, sill = sill, setup = setup, stdata = stdata)
  }

  LogNormDenStand_TukeyH <- function(const, cova, ident, dimat, mdecomp, delta, sill, setup, stdata) {
    llik <- 1.0e8
    if (delta <= 0 || sill <= 0) return(llik)
    decompvarcov <- MatDecomp(cova, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)
    Wval <- VGAM::lambertW(delta * stdata^2)
    IL <- sign(stdata) * sqrt(Wval / delta)
    IW <- 1 / (stdata * (1 + Wval))
    y <- try(forwardsolve(decompvarcov, IL, transpose = FALSE), silent = TRUE)
    if (inherits(y, "try-error")) return(llik)
    quad_form <- sum(y^2)
    log_jacobian <- -2 * sum(log(abs(IL * IW)))
    if (!is.finite(quad_form) || !is.finite(log_jacobian)) return(llik)
    0.5 * (const + dimat * log(sill) + logdetvarcov + quad_form + log_jacobian)
  }

  loglik_tukeyh <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                            grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance <- pram[namesnuis]
    Mean <- get_mean(X, nuisance, MM)

    if (aniso) {
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(pram[namesaniso]))
      coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- if (ncol(coords1) == 3) coords1[,3] else NULL
    }

    sill <- as.numeric(nuisance['sill'])
    tail <- as.numeric(nuisance['tail'])
    nugget <- as.numeric(nuisance['nugget'])
    if (tail <= 0 || tail > 0.5 || nugget < 0 || nugget >= 1 || sill <= 0) return(llik)
    nuisance_temp <- nuisance
    nuisance_temp['sill'] <- 1
    nuisance_temp['nugget'] <- 0
    corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance_temp, paramcorr, ns, NS, radius)
    if (any(!is.finite(corrv))) return(llik)
    corrmat_full <- build_correlation_matrix(corrv, ident)
    cova <- (1 - nugget) * corrmat_full + nugget * ident
    LogNormDenStand_TukeyH(stdata = ((data - c(Mean)) / sqrt(sill)), const = const, cova = cova,
                           dimat = dimat, ident = ident, mdecomp = mdecomp, delta = tail, sill = sill, setup = setup)
  }

  LogNormDenStand_LG <- function(const, cova, ident, dimat, mdecomp, nuisance, sill, setup, stdata, V_data, mu_s) {
    llik <- 1.0e8
    decompvarcov <- MatDecomp(cova, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)
    y <- try(forwardsolve(decompvarcov, stdata, transpose = FALSE, upper.tri = FALSE), silent = TRUE)
    if (inherits(y, "try-error")) return(llik)
    quadterm <- sum(y^2)
    sigma <- sqrt(sill)
    jacobian_term <- -dimat * log(sigma) - sum(log(V_data)) - sum(log(mu_s))
    loglik_val <- -0.5 * (const + logdetvarcov + quadterm) + jacobian_term
    -loglik_val
  }

  loglik_loggauss <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                              grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance <- pram[namesnuis]
    Mean <- get_mean(X, nuisance, MM)

    if (aniso) {
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(pram[namesaniso]))
      coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- if (ncol(coords1) == 3) coords1[,3] else NULL
    }

    sill <- as.numeric(nuisance["sill"])
    nugget <- as.numeric(nuisance["nugget"])
    if (sill <= 0 || nugget < 0 || nugget >= 1 || any(data <= 0)) return(llik)
    mu_s <- exp(Mean)
    V_data <- data / mu_s
    if (any(V_data <= 0)) return(llik)
    sigma <- sqrt(sill)
    Z_data <- (log(V_data) + sill / 2) / sigma
    nuisance_temp <- nuisance
    nuisance_temp["sill"] <- 1
    corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance_temp, paramcorr, ns, NS, radius)
    if (any(!is.finite(corrv))) return(llik)
    corrmat_full <- build_correlation_matrix(corrv, ident)
    cova <- (1 - nugget) * corrmat_full + nugget * ident
    LogNormDenStand_LG(stdata = Z_data, const = const, cova = cova, dimat = dimat, ident = ident,
                       mdecomp = mdecomp, nuisance = nuisance, sill = sill, setup = setup,
                       V_data = V_data, mu_s = mu_s)
  }

  LogNormDenStand_SH <- function(const,cova,mdecomp,nuisance,sill,setup,stdata, dimat = NULL, ...) {
    llik <- 1.0e8
    decompvarcov <- MatDecomp(cova, mdecomp)
    if (is.logical(decompvarcov) && !decompvarcov) return(llik)
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if (is.na(logdetvarcov)) return(llik)
    skew <- as.numeric(nuisance["skew"])
    delta <- as.numeric(nuisance["tail"])
    Z <- sinh(delta * asinh(stdata) - skew)
    C <- delta * sqrt((1 + Z^2) / (stdata^2 + 1))
    y <- try(forwardsolve(decompvarcov, Z, transpose = FALSE), silent = TRUE)
    if (inherits(y, "try-error") || any(C <= 0)) return(llik)
    llik <- 0.5 * (const + const * log(sill) / log(2*pi) + logdetvarcov - 2 * sum(log(C)) + sum(y^2))
    as.numeric(llik)
  }

  loglik_sh <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                        grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso) {
    llik <- 1.0e8
    pram <- combine_param(param, fixed, namesparam)
    paramcorr <- as.numeric(pram[namescorr])
    nuisance <- pram[namesnuis]
    Mean <- get_mean(X, nuisance, MM)

    if (aniso) {
      coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(pram[namesaniso]))
      coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- if (ncol(coords1) == 3) coords1[,3] else NULL
    }

    sill <- as.numeric(nuisance['sill'])
    nugget <- as.numeric(nuisance['nugget'])
    tail <- as.numeric(nuisance['tail'])
    if (tail <= 0 || sill <= 0 || nugget < 0 || nugget >= 1) return(llik)
    nuisance_temp <- nuisance
    nuisance_temp['sill'] <- 1
    corrv <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance_temp, paramcorr, ns, NS, radius)
    if (any(!is.finite(corrv))) return(llik)
    corrmat_full <- build_correlation_matrix(corrv, ident)
    cova <- (1 - nugget) * corrmat_full + nugget * ident
    do.call(what = fname, args = list(stdata = ((data - c(Mean)) / sqrt(sill)), const = const, cova = cova,
                                      mdecomp = mdecomp, nuisance = nuisance, sill = sill, setup = setup,
                                      dimat = dimat))
  }

  #########################################################################
  ## Main setup
  #########################################################################

  spacetime_dyn <- FALSE
  NS <- 0
  fname <- NULL
  hessian <- FALSE
  namesaniso <- c("angle", "ratio")

  if (!is.null(coordx_dyn)) spacetime_dyn <- TRUE
  if (grid) {
    a <- expand.grid(coordx, coordy)
    coordx <- a[, 1]
    coordy <- a[, 2]
  }

  if (!spacetime_dyn) dimat <- numcoord * numtime
  if (spacetime_dyn)  dimat <- sum(ns)

  if (is.null(dim(X))) {
    if (!bivariate) X <- as.matrix(rep(1, dimat))
    if ( bivariate) X <- as.matrix(rep(1, ns[1] + ns[2]))
  } else {
    if (!bivariate) num_betas <- ncol(X)
    if ( bivariate) num_betas <- c(ncol(X), ncol(X))
  }

  corrmat <- "CorrelationMat"
  if (model == 36 || model == 47) corrmat <- "CorrelationMat_dis"
  if (spacetime) {
    corrmat <- "CorrelationMat_st_dyn"
    if (model == 36 || model == 47) corrmat <- "CorrelationMat_st_dyn_dis"
  }
  if (bivariate) corrmat <- "CorrelationMat_biv_dyn"

  if (spacetime || bivariate) {
    NS <- cumsum(ns)
    if (spacetime_dyn) {
      data <- t(unlist(data))
      NS <- c(0, NS)[-(length(ns) + 1)]
    } else {
      data <- matrix(t(data), nrow = 1)
      NS <- rep(0, numtime)
    }
  }

  dd <- 0
  if (bivariate) dd <- dimat
  numpairstot <- dimat * (dimat - 1) / 2 + dd
  const <- dimat * log(2*pi)

  if (is.null(neighb)) {
    corr <- double(numpairstot)
    ident <- diag(dimat)
  } else {
    corr <- double(numpairstot)
    ident <- diag(dimat)
  }

  #########################################################################
  ## Select objective and density function
  #########################################################################

  if (model == 1 || model == 20) {
    lname <- 'loglik'
    if (bivariate) lname <- 'loglik_biv'

    if (type == 3) {
      if (bivariate) fname <- "LogNormDenRestr_biv" else fname <- "LogNormDenRestr"
      const <- const - log(2*pi)
    }
    if (type == 4) {
      if (bivariate) fname <- 'LogNormDenStand_biv' else fname <- 'LogNormDenStand'
    }
    if (type == 5 || type == 6) {
      corrmat <- "CorrelationMat_tap"
      if (spacetime) corrmat <- "CorrelationMat_st_tap"
      if (bivariate) {
        if (type == 5) fname <- 'LogNormDenTap_biv'
        if (type == 6) fname <- 'LogNormDenTap1_biv'
        corrmat <- "CorrelationMat_biv_tap"
      } else {
        if (type == 5) fname <- 'LogNormDenTap'
        if (type == 6) fname <- 'LogNormDenTap1'
      }
      corr <- double(numpairs)
      tapcorr <- double(numpairs)
      tcor <- matr(corrmat, tapcorr, coordx, coordy, coordz, coordt, setup$tapmodel, c(0,0,1), 1, ns, NS, radius)
      tape <- new("spam", entries = tcor, colindices = setup$ja, rowpointers = setup$ia,
                  dimension = as.integer(rep(dimat, 2)))
      setup$struct <- try(spam::chol.spam(tape, silent = TRUE))
      setup$taps <- tcor
    }
    if (type == 8) {
      if (bivariate) fname <- "CVV_biv" else fname <- "CVV"
    }
    if (sparse) fname <- paste(fname, "_spam", sep = "")
  }

  if (model == 20) {
    lname <- 'loglik_sh'
    if (bivariate) lname <- 'loglik_biv_sh'
    fname <- 'LogNormDenStand_SH'
    if (type == 6) fname <- 'LogshDenTap1'
  }
  if (model == 34) {
    lname <- 'loglik_tukeyh'
    if (bivariate) lname <- 'loglik_biv_tukeyh'
  }
  if (model == 40) {
    lname <- 'loglik_tukey2h'
    if (bivariate) lname <- 'loglik_biv_tukey2h'
  }
  if (model == 47) {
    lname <- 'loglik_miss_Poisgamma'
    if (bivariate) lname <- 'loglik_biv_miss_Poisgamma'
  }
  if (model == 35) {
    lname <- 'loglik_miss_T'
    if (bivariate) lname <- 'loglik_biv_miss_T'
  }
  if (model == 36) {
    lname <- 'loglik_miss_Pois'
    if (bivariate) lname <- 'loglik_biv_miss_Pois'
  }
  if (model == 37) {
    lname <- 'loglik_miss_skewT'
    if (bivariate) lname <- 'loglik_biv_miss_skewT'
  }
  if (model == 22) {
    lname <- 'loglik_loggauss'
    if (bivariate) lname <- 'loglik_biv_loggauss'
  }

  if (type != 5 && type != 6) corrmat <- paste(corrmat, "2", sep = "")
  if (is.null(coordz)) coordz <- double(length(coordx))

  if (!exists(lname, mode = "function")) {
    stop("Objective function not implemented in GeoLik.R: ", lname)
  }
  obj_fun <- get(lname, mode = "function")
  data_opt <- as.numeric(t(data))

  common_args <- list(
    const = const,
    coordx = coordx,
    coordy = coordy,
    coordz = coordz,
    coordt = coordt,
    corr = corr,
    corrmat = corrmat,
    corrmodel = corrmodel,
    data = data_opt,
    dimat = dimat,
    fixed = fixed,
    fname = fname,
    grid = grid,
    ident = ident,
    mdecomp = mdecomp,
    model = model,
    namescorr = namescorr,
    namesnuis = namesnuis,
    namesparam = namesparam,
    radius = radius,
    setup = setup,
    X = X,
    ns = ns,
    NS = NS,
    MM = MM,
    aniso = aniso,
    namesaniso = namesaniso
  )

  #########################################################################
  ## Optimization
  #########################################################################

  if (!onlyvar) {
    ## Safe, conservative controls for full-likelihood optimization.
    ## These match the same spirit used in CompLik2: tighter than the
    ## previous quick defaults, but without changing the likelihood itself.
    opt_reltol <- 1e-12
    opt_pgtol  <- 1e-10
    opt_factr  <- 1e4
    opt_maxit  <- 100000L

    if (length(param) == 1) {
      optimizer <- "optimize"
      opt_bounds <- make_optimize_bounds(param, lower, upper, namesparam)
      Likelihood <- do.call(stats::optimize,
                            c(list(f = obj_fun, lower = opt_bounds$lower, upper = opt_bounds$upper, maximum = FALSE), common_args))
    }

    if (length(param) > 1) {
      if (optimizer == 'L-BFGS-B') {
        Likelihood <- do.call(stats::optim,
                              c(list(par = param, fn = obj_fun, method = 'L-BFGS-B', lower = lower, upper = upper,
                                     hessian = hessian,
                                     control = list(factr = opt_factr, pgtol = opt_pgtol, maxit = opt_maxit, parscale = abs(param) + 1)),
                                common_args))
      }

      if (optimizer %in% c('Nelder-Mead', 'BFGS', 'SANN')) {
        opt_control <- list(reltol = opt_reltol, maxit = opt_maxit, parscale = abs(param) + 1)
        if (optimizer == 'SANN') opt_control <- list(maxit = opt_maxit, parscale = abs(param) + 1)
        Likelihood <- do.call(stats::optim,
                              c(list(par = param, fn = obj_fun, method = optimizer, hessian = hessian,
                                     control = opt_control), common_args))
      }

      if (optimizer == 'nlm') {
        Likelihood <- do.call(stats::nlm,
                              c(list(f = obj_fun, p = param, hessian = hessian, iterlim = opt_maxit, steptol = 1e-6, gradtol = 1e-8),
                                common_args))
      }

      if (optimizer == 'bobyqa') {
        bobyqa_args <- make_bobyqa_setup(param, lower, upper, namesparam)
        bobyqa_args$control$maxfun <- opt_maxit
        Likelihood <- do.call(minqa::bobyqa,
                              c(list(par = param, fn = obj_fun,
                                     lower = bobyqa_args$lower,
                                     upper = bobyqa_args$upper,
                                     control = bobyqa_args$control),
                                common_args))
      }

      if (optimizer == 'nlminb') {
        Likelihood <- do.call(stats::nlminb,
                              c(list(start = param, objective = obj_fun, lower = lower, upper = upper,
                                     control = list(iter.max = opt_maxit, eval.max = opt_maxit,
                                                    rel.tol = opt_reltol, x.tol = 1e-8)),
                                common_args))
      }
    }

    if (!exists("Likelihood")) stop("Optimizer not implemented: ", optimizer)

    #######################################################################
    ## Post-processing
    #######################################################################

    if (optimizer %in% c('Nelder-Mead', 'L-BFGS-B', 'BFGS', 'SANN')) {
      conv_code <- Likelihood$convergence
      conv_msg <- if (!is.null(Likelihood$message)) Likelihood$message else NA_character_
      names(Likelihood$par) <- namesparam
      param <- Likelihood$par
      maxfun <- -Likelihood$value
      Likelihood$value <- maxfun
      Likelihood$convergence_code <- conv_code
      Likelihood$opt_code <- conv_code
      Likelihood$optim_message <- conv_msg
      Likelihood$opt_message <- conv_msg
      Likelihood$counts <- as.numeric(Likelihood$counts[1])

      if (conv_code == 0) {
        Likelihood$convergence <- 'Successful'
      } else if (conv_code == 1) {
        Likelihood$convergence <- paste0(optimizer, ': iteration limit reached')
      } else if (optimizer == 'Nelder-Mead' && conv_code == 10) {
        Likelihood$convergence <- 'Nelder-Mead warning: simplex degenerated; check solution'
      } else {
        Likelihood$convergence <- paste0(optimizer, ' warning: code=', conv_code)
      }
    }

    if (optimizer == 'bobyqa') {
      names(Likelihood$par) <- namesparam
      param <- Likelihood$par
      maxfun <- -Likelihood$fval
      Likelihood$value <- maxfun
      Likelihood$convergence_code <- Likelihood$ierr
      Likelihood$opt_code <- Likelihood$ierr
      Likelihood$optim_message <- Likelihood$msg
      Likelihood$opt_message <- Likelihood$msg
      Likelihood$counts <- if (!is.null(Likelihood$feval)) as.numeric(Likelihood$feval) else NA_real_

      if (Likelihood$ierr == 0) {
        Likelihood$convergence <- 'Successful'
      } else if (Likelihood$ierr == 1) {
        Likelihood$convergence <- 'bobyqa warning: maximum number of function evaluations reached'
      } else if (Likelihood$ierr == 3) {
        Likelihood$convergence <- 'bobyqa warning: trust-region step failed to reduce quadratic model'
      } else if (Likelihood$ierr == 4) {
        Likelihood$convergence <- 'bobyqa warning: bounds too tight relative to rhobeg'
      } else if (Likelihood$ierr == 5) {
        Likelihood$convergence <- 'bobyqa warning: numerical cancellation problem'
      } else {
        Likelihood$convergence <- paste0('bobyqa warning: ierr=', Likelihood$ierr)
      }
    }

    if (optimizer %in% c('nlminb', 'multinlminb')) {
      conv_code <- Likelihood$convergence
      conv_msg <- if (!is.null(Likelihood$message)) Likelihood$message else NA_character_
      names(Likelihood$par) <- namesparam
      param <- Likelihood$par
      maxfun <- -Likelihood$objective
      Likelihood$value <- maxfun
      Likelihood$convergence_code <- conv_code
      Likelihood$opt_code <- conv_code
      Likelihood$optim_message <- conv_msg
      Likelihood$opt_message <- conv_msg
      Likelihood$counts <- as.numeric(Likelihood$iterations)
      if (conv_code == 0) {
        Likelihood$convergence <- 'Successful'
      } else {
        Likelihood$convergence <- paste0(
          'nlminb warning: code=', conv_code,
          if (!is.na(conv_msg)) paste0('; message=', conv_msg) else ''
        )
      }
    }

    if (optimizer == 'nlm') {
      conv_code <- Likelihood$code
      names(Likelihood$estimate) <- namesparam
      param <- Likelihood$estimate
      maxfun <- -as.numeric(Likelihood$minimum)
      Likelihood$value <- maxfun
      Likelihood$par <- param
      Likelihood$param <- param
      Likelihood$convergence_code <- conv_code
      Likelihood$opt_code <- conv_code
      Likelihood$optim_message <- NA_character_
      Likelihood$opt_message <- NA_character_
      Likelihood$counts <- as.numeric(Likelihood$iterations)
      if (conv_code == 1 || conv_code == 2) {
        Likelihood$convergence <- 'Successful'
      } else if (conv_code == 4 || conv_code == 5) {
        Likelihood$convergence <- paste0('nlm warning: code=', conv_code, '; iteration/step limit reached')
      } else {
        Likelihood$convergence <- paste0('nlm warning: code=', conv_code)
      }
    }

    if (optimizer == 'optimize') {
      param <- Likelihood$minimum
      names(param) <- namesparam
      maxfun <- -Likelihood$objective
      Likelihood$value <- maxfun
      Likelihood$par <- param
      Likelihood$param <- param
      Likelihood$convergence <- 'Successful'
      Likelihood$convergence_code <- 0
      Likelihood$optim_message <- NA_character_
      Likelihood$counts <- NA_real_
    }

    serious_failure <- is.null(maxfun) || !is.finite(maxfun) || maxfun <= -1.0e14
    if (serious_failure) {
      Likelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
      warning('Optimization may have failed: try with other starting values', call. = FALSE)
    }

    if (is.null(Likelihood$par)) Likelihood$par <- param
    names(Likelihood$par) <- namesparam

    numparam <- length(param)
    Likelihood$claic <- NULL
    Likelihood$clbic <- NULL
    if (type == 3 || type == 4) {
      Likelihood$claic <- -2 * (maxfun - numparam)
      Likelihood$clbic <- -2 * maxfun + numparam * log(dimat)
    }
  } else {
    Likelihood <- list(value = NA_real_)
    Likelihood$par <- param
    names(Likelihood$par) <- namesparam
    Likelihood$param <- param
    Likelihood$claic <- NULL
    Likelihood$clbic <- NULL
    Likelihood$convergence <- 'None'
    serious_failure <- FALSE
    maxfun <- Likelihood$value
    numparam <- length(param)
  }

  if (!exists("serious_failure", inherits = FALSE)) serious_failure <- FALSE

  #########################################################################
  ## Hessian, score, asymptotic covariance
  #########################################################################

  if (varest) {
    Likelihood$hessian <- numDeriv::hessian(
      func = obj_fun,
      x = as.numeric(Likelihood$par),
      method = "Richardson",
      const = const,
      coordx = coordx,
      coordy = coordy,
      coordz = coordz,
      coordt = coordt,
      corr = corr,
      corrmat = corrmat,
      corrmodel = corrmodel,
      data = data_opt,
      dimat = dimat,
      fixed = fixed,
      fname = fname,
      grid = grid,
      ident = ident,
      mdecomp = mdecomp,
      model = model,
      namescorr = namescorr,
      namesnuis = namesnuis,
      namesparam = namesparam,
      radius = radius,
      setup = setup,
      X = X,
      ns = ns,
      NS = NS,
      MM = MM,
      aniso = aniso,
      namesaniso = namesaniso
    )
    rownames(Likelihood$hessian) <- namesparam
    colnames(Likelihood$hessian) <- namesparam
  }

  Likelihood$score <- NULL
  if (score) {
    Likelihood$score <- -numDeriv::grad(
      func = obj_fun,
      x = as.numeric(Likelihood$par),
      method = "Richardson",
      const = const,
      coordx = coordx,
      coordy = coordy,
      coordz = coordz,
      coordt = coordt,
      corr = corr,
      corrmat = corrmat,
      corrmodel = corrmodel,
      data = data_opt,
      dimat = dimat,
      fixed = fixed,
      fname = fname,
      grid = grid,
      ident = ident,
      mdecomp = mdecomp,
      model = model,
      namescorr = namescorr,
      namesnuis = namesnuis,
      namesparam = namesparam,
      radius = radius,
      setup = setup,
      X = X,
      ns = ns,
      NS = NS,
      MM = MM,
      aniso = aniso,
      namesaniso = namesaniso
    )
    names(Likelihood$score) <- namesparam
  }

  if (!isTRUE(serious_failure) || identical(Likelihood$convergence, 'None')) {
    if (varest) {
      if ((model == 20 || model == 22 || model == 1 || model == 34 || model == 35 || model == 37) && !(type == 5 || type == 6)) {
        aa <- try(abs(det(Likelihood$hessian)), silent = TRUE)
        if (inherits(aa, "try-error") || !is.finite(aa) || aa <= 1e-12) {
          warning("Asymptotic information matrix is singular", immediate. = TRUE)
          Likelihood$varcov <- NULL
        } else {
          ii <- try(solve(Likelihood$hessian), silent = TRUE)
          if (inherits(ii, "try-error")) {
            warning("Asymptotic information matrix is singular", immediate. = TRUE)
            Likelihood$varcov <- NULL
          } else {
            eig <- eigen(ii, symmetric = TRUE, only.values = TRUE)$values
            if (any(!is.finite(eig)) || min(eig) <= 0) {
              warning("Asymptotic information matrix is not positive-definite")
              Likelihood$varcov <- NULL
            } else {
              Likelihood$varcov <- ii
              dimnames(Likelihood$varcov) <- list(namesparam, namesparam)
              Likelihood$stderr <- sqrt(diag(Likelihood$varcov))
              names(Likelihood$stderr) <- namesparam
            }
          }
        }
      } else if (!is.null(Likelihood$varcov)) {
        eig <- eigen(Likelihood$varcov, symmetric = TRUE, only.values = TRUE)$values
        if (any(!is.finite(eig)) || min(eig) < 0) {
          warning("Asymptotic information matrix is not positive-definite")
          Likelihood$varcov <- 'none'
          Likelihood$stderr <- 'none'
        } else {
          dimnames(Likelihood$varcov) <- list(namesparam, namesparam)
          Likelihood$stderr <- diag(Likelihood$varcov)
          if (any(Likelihood$stderr < 0)) {
            Likelihood$stderr <- 'none'
          } else {
            Likelihood$stderr <- sqrt(Likelihood$stderr)
            names(Likelihood$stderr) <- namesparam
          }
        }
      }
    }
  }

  if (varest && is.null(Likelihood$varcov)) {
    Likelihood$varcov <- 'none'
    Likelihood$stderr <- 'none'
  }

  if (varest) Likelihood$sensmat <- Likelihood$hessian

  return(Likelihood)
}
