####################################################
### File name: GeoQQ.r (cleaned; comments in English)
####################################################
GeoQQ <- function(fit, type="Q", add=FALSE, ylim=c(0,1), xlim=NULL, breaks=10, ...)
{
  ###### utilities ##########
  # Log-logistic quantile with shape (>0) and scale (>0)
  qllogis1 <- function(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE) {
    if (missing(shape)) stop("argument 'shape' is missing, with no default")
    if (shape <= 0) stop("'shape' must be positive")
    if (scale <= 0) stop("'scale' must be positive")
    if (log.p) p <- exp(p)
    if (any(p < 0 | p > 1)) stop("probabilities must be between 0 and 1")
    if (!lower.tail) p <- 1 - p
    result <- numeric(length(p)); result[p == 0] <- 0; result[p == 1] <- Inf
    valid <- p > 0 & p < 1
    if (any(valid)) {
      pv <- p[valid]
      result[valid] <- scale * (pv / (1 - pv))^(1/shape)
    }
    result
  }

  # Log-logistic density (fix: density, not quantile)
  dllogis1 <- function(x, shape, rate = 1, scale = 1/rate) {
    if (missing(shape)) stop("argument 'shape' is missing, with no default")
    if (shape <= 0) stop("'shape' must be positive")
    if (scale <= 0) stop("'scale' must be positive")
    out <- numeric(length(x)); out[] <- 0
    pos <- x > 0 & is.finite(x)
    z <- x[pos]/scale
    out[pos] <- (shape/scale) * z^(shape - 1) / (1 + z^shape)^2
    out
  }
  ##################################

  if(!inherits(fit,"GeoFit"))  stop("A GeoFit object is needed as input\n")
  if(type!="Q" && type!="D") stop("Type can be Q or D \n")

  model <- fit$model
  copula <- if(!is.null(fit$copula)) fit$copula else NULL

  fit$param <- unlist(fit$param); fit$fixed <- unlist(fit$fixed)
  pp <- c(fit$param, fit$fixed)

  # Safe extraction: single-bracket returns NA when name is absent (no error)
  MM <- as.numeric(pp["mean"])
  VV <- as.numeric(pp["sill"])   # may be NA for models that don't need it
  if (is.na(MM)) stop("'mean' parameter is required")

  # Which models actually require 'sill' (marginal variance)?
  .requires_VV <- function(m) {
    m %in% c(
      "Gaussian","Gaussian_misp_Binomial","Gaussian_misp_Poisson","Gaussian_misp_BinomialNeg",
      "SkewGaussian","SkewLaplace",
      "StudentT","Gaussian_misp_StudentT",
      "SkewStudentT","Gaussian_misp_SkewStudentT",
      "LogGaussian","Logistic",
      "SinhAsinh","Tukeyh","Tukeyh2",
      "Tukeygh","Gaussian_misp_Tukeygh",
      "TwoPieceGaussian","TwoPieceStudentT","TwoPieceTukeyh"
    )
  }

  # Save and restore par settings
  opar <- par(no.readonly = TRUE); on.exit(par(opar))

  # xlim helpers
  has_xlim <- !is.null(xlim)
  .XL <- function(def) { if (has_xlim) xlim else def }            # choose axis limits
  .GRID <- function(def_range, n = 512) {                         # grid across visible x-range
    xr <- if (has_xlim) xlim else def_range
    seq(xr[1], xr[2], length.out = n)
  }

  # Library checks for skew/t families when needed (clear error if missing)
  need_sn   <- model %in% c("SkewGaussian","SkewStudentT","Gaussian_misp_SkewStudentT")
  need_VGAM <- model %in% c("Tukeyh","Tukeyh2","TwoPieceTukeyh","Tukeygh","Gaussian_misp_Tukeygh")
  if (need_sn && !requireNamespace("sn", quietly = TRUE))
    stop("Package 'sn' is required for skew-* models.")
  if (need_VGAM && !requireNamespace("VGAM", quietly = TRUE))
    stop("Package 'VGAM' is required for Tukey-h/gh models.")

  ########################################
  #### QQ plot
  ########################################
  if(type=="Q") {
    xlab <- "Theoretical Quantiles"
    ylab <- "Sample Quantiles"

    if(!fit$bivariate){
      # Unpack data (finite only)
      dd <- if (is.list(fit$coordx_dyn)) unlist(fit$data) else c(t(fit$data))
      dd <- dd[is.finite(dd)]
      if (!length(dd)) stop("empty data")

      N <- length(dd)
      probabilities  <- (1:N)/(N+1)
      probabilities1 <- c(0.25,0.75)
      q_e  <- quantile(dd, probabilities)
      q_e1 <- quantile(dd, probabilities1)

      # If the model needs VV, enforce presence once (univariate)
      if (.requires_VV(model) && is.na(VV)) {
        stop("'sill' parameter is required for model '", model,
             "' but is missing in 'fit'. Provide 'sill' in fit$param/fixed.")
      }

      # ---- Models ----
      if(model %in% c("Beta2")){
        mm <- 1/(1+exp(-MM)); sh <- as.numeric(pp["shape"])
        pmin <- as.numeric(pp["min"]); pmax <- as.numeric(pp["max"])
        q_t  <- pmin + (pmax-pmin)*qbeta(probabilities,  shape1=mm*sh, shape2=(1-mm)*sh)
        q_t1 <- pmin + (pmax-pmin)*qbeta(probabilities1, shape1=mm*sh, shape2=(1-mm)*sh)
        plot(q_t,q_e, main="Beta qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Kumaraswamy2")){
        mm <- 1/(1+exp(-MM)); sh <- as.numeric(pp["shape"])
        pmin <- as.numeric(pp["min"]); pmax <- as.numeric(pp["max"])
        shape1 <- log(0.5)/log1p(-(mm^sh))
        q_t  <- pmin + (pmax-pmin)*(1 - (1 - (probabilities)^sh)^shape1)
        q_t1 <- pmin + (pmax-pmin)*(1 - (1 - (probabilities1)^sh)^shape1)
        plot(q_t,q_e, main="Kumaraswamy qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Gaussian","Gaussian_misp_Binomial","Gaussian_misp_Poisson","Gaussian_misp_BinomialNeg")) {
        q_t  <- qnorm(probabilities,  mean=MM, sd=sqrt(VV))
        q_t1 <- qnorm(probabilities1, mean=MM, sd=sqrt(VV))
        plot(q_t,q_e, main="Gaussian qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Binomial")) {
        q_t  <- qbinom(probabilities,  size=fit$n, prob=pnorm(pp["mean"]))
        q_t1 <- qbinom(probabilities1, size=fit$n, prob=pnorm(pp["mean"]))
        plot(q_t,q_e, main="Binomial qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("BinomialLogistic")) {
        q_t  <- qbinom(probabilities,  size=fit$n, prob=plogis(MM))
        q_t1 <- qbinom(probabilities1, size=fit$n, prob=plogis(MM))
        plot(q_t,q_e, main="Binomial-Logistic qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("BinomialNeg")) {
        q_t  <- qnbinom(probabilities,  size=fit$n, prob=pnorm(pp["mean"]))
        q_t1 <- qnbinom(probabilities1, size=fit$n, prob=pnorm(pp["mean"]))
        plot(q_t,q_e, main="Binomial Neg qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Poisson")) {
        q_t  <- qpois(probabilities,  lambda=exp(pp["mean"]))
        q_t1 <- qpois(probabilities1, lambda=exp(pp["mean"]))
        plot(q_t,q_e, main="Poisson qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("PoissonGamma")) {
        ff <- exp(pp["mean"])
        sh <- as.numeric(pp["shape"])
        q_t  <- qnbinom(probabilities,  size=sh, prob=sh/(ff+sh))
        q_t1 <- qnbinom(probabilities1, size=sh, prob=sh/(ff+sh))
        plot(q_t,q_e, main="Poisson qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("SkewGaussian")) {
        omega <- as.numeric(sqrt((pp["skew"]^2 + pp["sill"])/pp["sill"]))
        alpha <- as.numeric(pp["skew"]/sqrt(pp["sill"]))
        q_t  <- MM + sqrt(VV)*sn::qsn(probabilities,  xi=0, omega=omega, alpha=alpha)
        q_t1 <- MM + sqrt(VV)*sn::qsn(probabilities1, xi=0, omega=omega, alpha=alpha)
        plot(q_t,q_e, main="Skew Gaussian qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model%in%c("StudentT","Gaussian_misp_StudentT")){
        df <- as.numeric(round(1/pp["df"]))
        q_t  <- MM + sqrt(VV)*qt(probabilities,  df=df)
        q_t1 <- MM + sqrt(VV)*qt(probabilities1, df=df)
        plot(q_t,q_e, main="t qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model%in%c("SkewStudentT","Gaussian_misp_SkewStudentT")){
        alpha <- as.numeric(pp["skew"])
        nu    <- as.numeric(round(1/pp["df"]))
        q_t  <- MM + sqrt(VV)*sn::qst(probabilities,  xi=0, omega=1, alpha=alpha, nu=nu)
        q_t1 <- MM + sqrt(VV)*sn::qst(probabilities1, xi=0, omega=1, alpha=alpha, nu=nu)
        plot(q_t,q_e, main="skewt qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Weibull")){
        shape <- pp["shape"]
        q_t  <- exp(MM)*qweibull(probabilities,  shape=shape, scale=1/(gamma(1+1/shape)))
        q_t1 <- exp(MM)*qweibull(probabilities1, shape=shape, scale=1/(gamma(1+1/shape)))
        plot(q_t,q_e, main="Weibull qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Gamma")){
        shape <- pp["shape"]
        q_t  <- exp(MM)*qgamma(probabilities,  shape=shape/2, rate=shape/2)
        q_t1 <- exp(MM)*qgamma(probabilities1, shape=shape/2, rate=shape/2)
        plot(q_t,q_e, main="Gamma qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("LogGaussian")){
        # Fix: meanlog must be MM - VV/2 to ensure E[Y]=exp(MM)
        q_t  <- qlnorm(probabilities,  meanlog = MM - VV/2, sdlog = sqrt(VV))
        q_t1 <- qlnorm(probabilities1, meanlog = MM - VV/2, sdlog = sqrt(VV))
        plot(q_t,q_e, main="LogGaussian qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("LogLogistic")){
        shape <- pp["shape"]
        cc <- gamma(1+1/shape)*gamma(1-1/shape)
        q_t  <- qllogis1(probabilities,  shape=shape, scale=exp(MM)/cc)
        q_t1 <- qllogis1(probabilities1, shape=shape, scale=exp(MM)/cc)
        plot(q_t,q_e, main="LogLogistic qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Logistic")){
        q_t  <- qlogis(probabilities,  location=MM, scale=sqrt(VV))
        q_t1 <- qlogis(probabilities1, location=MM, scale=sqrt(VV))
        plot(q_t,q_e, main="Logistic qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("SinhAsinh")){
        tail <- as.numeric(pp["tail"]); skew <- as.numeric(pp["skew"])
        q_t  <- MM + sqrt(VV)*sinh( (asinh(qnorm(probabilities))/tail) + skew/tail )
        q_t1 <- MM + sqrt(VV)*sinh( (asinh(qnorm(probabilities1))/tail) + skew/tail )
        plot(q_t,q_e, main="Sas qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Tukeyh")){
        tail <- as.numeric(pp["tail"])
        uu <- qnorm(probabilities); uu1 <- qnorm(probabilities1)
        q_t  <- MM + sqrt(VV)*uu * exp(0.5*tail*uu^2)
        q_t1 <- MM + sqrt(VV)*uu1* exp(0.5*tail*uu1^2)
        plot(q_t,q_e, main="Tukey-h qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Tukeyh2")){
        qtpTukeyh221 <- function(x,tail1,tail2){
          sel1 <- which(x>=0); sel2 <- which(x<0)
          x1 <- x[sel1]; x2 <- x[sel2]
          uu1 <- qnorm(x1,0,1); uu2 <- qnorm(x2,0,1)
          qq1 <- uu1*exp(0.5*tail1*uu1^2)
          qq2 <- uu2*exp(0.5*tail2*uu2^2)
          c(qq2,qq1)
        }
        tail1 <- as.numeric(pp["tail1"]); tail2 <- as.numeric(pp["tail2"])
        q_t  <- MM + sqrt(VV)*qtpTukeyh221(probabilities,  tail1,tail2)
        q_t1 <- MM + sqrt(VV)*qtpTukeyh221(probabilities1, tail1,tail2)
        plot(q_t,q_e, main="Tukey-hh qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }


      if(model %in% c("SkewLaplace")){
      qtpSkewLaplace22 <- function(u, sk){
      mm=0
      vv=1
      s  <- sqrt(vv)
      res <- rep(NA_real_, length(u))
      hi <- u >= sk
      lo <- u <  sk
      res[hi] <- mm - (s / sk) * (log1p(-u[hi]) - log1p(-sk))
      res[lo] <- mm + (s / (1 - sk)) * (log(u[lo]) - log(sk))
      res
      }
        skew <- as.numeric(pp["skew"])
        q_t  <- MM + sqrt(VV)*qtpSkewLaplace22(probabilities,skew)
        q_t1 <- MM + sqrt(VV)*qtpSkewLaplace22(probabilities1,skew)
        plot(q_t,q_e, main="Skew Laplace qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("Tukeygh","Gaussian_misp_Tukeygh")){
        tail <- as.numeric(pp["tail"]); skew <- as.numeric(pp["skew"])
        uu <- qnorm(probabilities); uu1 <- qnorm(probabilities1)
        if (abs(skew) < 1e-8) {
          q_t  <- MM + sqrt(VV) * (uu  * exp(0.5*tail*uu^2))
          q_t1 <- MM + sqrt(VV) * (uu1 * exp(0.5*tail*uu1^2))
        } else {
          q_t  <- MM + sqrt(VV) * ((exp(skew*uu)  - 1) * exp(0.5*tail*uu^2)  / skew)
          q_t1 <- MM + sqrt(VV) * ((exp(skew*uu1) - 1) * exp(0.5*tail*uu1^2) / skew)
        }
        plot(q_t,q_e, main="Tukey-gh qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("TwoPieceGaussian")){
        qtpGaussian <- function(x,skew){
          sel1 <- which(x>0 & x<0.5*(1+skew))
          sel2 <- which(x>=0.5*(1+skew) & x<=1)
          x1 <- x[sel1]; x2 <- x[sel2]
          qq1 <- (1+skew)*qnorm(x1/(1+skew))
          qq2 <- (1-skew)*qnorm((x2-skew)/(1-skew))
          c(qq1,qq2)
        }
        skew <- as.numeric(pp["skew"])
        q_t  <- MM + sqrt(VV)*qtpGaussian(probabilities,  skew)
        q_t1 <- MM + sqrt(VV)*qtpGaussian(probabilities1, skew)
        plot(q_t,q_e, main="Two-Piece Gaussian qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("TwoPieceStudentT")){
        qtpt1 <- function(x,skew,df){
          sel1 <- which(x>0 & x<0.5*(1+skew))
          sel2 <- which(x>=0.5*(1+skew) & x<=1)
          x1 <- x[sel1]; x2 <- x[sel2]
          qq1 <- (1+skew)*qt(x1/(1+skew), df=df)
          qq2 <- (1-skew)*qt((x2-skew)/(1-skew), df=df)
          c(qq1,qq2)
        }
        skew <- as.numeric(pp["skew"]); df <- 1/as.numeric(pp["df"])
        q_t  <- MM + sqrt(VV)*qtpt1(probabilities,  skew,df)
        q_t1 <- MM + sqrt(VV)*qtpt1(probabilities1, skew,df)
        plot(q_t,q_e, main="Two-Piece Student qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      if(model %in% c("TwoPieceTukeyh")){
        qtukh <- function(xx,tail){
          uu <- qnorm(xx); uu*exp(0.5*tail*uu^2)
        }
        qtptukey <- function(x,skew,tail){
          sel1 <- which(x>0 & x<0.5*(1+skew))
          sel2 <- which(x>=0.5*(1+skew) & x<=1)
          x1 <- x[sel1]; x2 <- x[sel2]
          qq1 <- (1+skew)*qtukh(x1/(1+skew), tail=tail)
          qq2 <- (1-skew)*qtukh((x2-skew)/(1-skew), tail=tail)
          c(qq1,qq2)
        }
        skew <- as.numeric(pp["skew"]); tail <- as.numeric(pp["tail"])
        q_t  <- MM + sqrt(VV)*qtptukey(probabilities,  skew,tail)
        q_t1 <- MM + sqrt(VV)*qtptukey(probabilities1, skew,tail)
        plot(q_t,q_e, main="Two-Piece Tukey-h qq-plot", xlab=xlab, ylab=ylab,
             xlim=.XL(range(q_t, finite=TRUE)), ...)
      }

      # Robust reference line from 25%–75%
      x <- q_t1; y <- q_e1
      slope <- diff(y)/diff(x); int <- y[1L] - slope * x[1L]
      abline(int, slope, pch=20)
    }

    if(fit$bivariate){
      par(mfrow=c(1,2))
      if(is.list(fit$coordx_dyn)){ dd1 <- fit$data[[1]]; dd2 <- fit$data[[2]]}
      else  {dd1 <- fit$data[1,]; dd2 <- fit$data[2,];}
      dd1 <- dd1[is.finite(dd1)]; dd2 <- dd2[is.finite(dd2)]

      if(model %in% c("Gaussian")) {
        if (has_xlim) {
          qqnorm(dd1,main="First Gaussian qq-plot", xlim=xlim); abline(0,1)
          qqnorm(dd2,main="Second Gaussian qq-plot", xlim=xlim); abline(0,1)
        } else {
          qqnorm(dd1,main="First Gaussian qq-plot"); abline(0,1)
          qqnorm(dd2,main="Second Gaussian qq-plot"); abline(0,1)
        }
      }

      if(model %in% c("SkewGaussian")){
        omega1 <- sqrt((pp["skew_1"]^2 + pp["sill_1"])/pp["sill_1"])
        alpha1 <- pp["skew_1"]/sqrt(pp["sill_1"])
        probabilities <- (1:length(dd1))/(length(dd1)+1)
        q_t1 <- sn::qsn(probabilities, xi=0, omega=as.numeric(omega1), alpha=as.numeric(alpha1))
        q_e1 <- quantile(dd1, probabilities)
        plot(q_t1,q_e1, main="First Skew Gaussian qq-plot",
             xlim=.XL(range(q_t1, finite=TRUE)), ...)
        aa <- lm(q_e1~1+q_t1); abline(as.numeric(aa$coefficients[1]),as.numeric(aa$coefficients[2]))

        omega2 <- sqrt((pp["skew_2"]^2 + pp["sill_2"])/pp["sill_2"])
        alpha2 <- pp["skew_2"]/sqrt(pp["sill_2"])
        probabilities <- (1:length(dd2))/(length(dd2)+1)
        q_t2 <- sn::qsn(probabilities, xi=0, omega=as.numeric(omega2), alpha=as.numeric(alpha2))
        q_e2 <- quantile(dd2, probabilities)
        plot(q_t2,q_e2, main="Second Skew Gaussian qq-plot",
             xlim=.XL(range(q_t2, finite=TRUE)), ...)
        aa <- lm(q_e2~1+q_t2); abline(as.numeric(aa$coefficients[1]),as.numeric(aa$coefficients[2]))
      }
    }
  } # end QQ


  ##########################################################
  #### Density / histogram
  ##########################################################
  if(type=="D") {
    if(!fit$bivariate){

      dd <- if (is.list(fit$coordx_dyn)) unlist(fit$data) else c(t(fit$data))
      dd <- dd[is.finite(dd)]
      if (!length(dd)) stop("empty data")

      # If the model needs VV, enforce presence once (univariate)
      if (.requires_VV(model) && is.na(VV)) {
        stop("'sill' parameter is required for model '", model,
             "' but is missing in 'fit'. Provide 'sill' in fit$param/fixed.")
      }

      # Helper: draw histogram once; overlay theoretical curve with lines()
      .hist_once <- function(main_title, xr_default) {
        if(!add) {
          hist(dd, freq=FALSE, xlim=.XL(xr_default), ylab="Density", xlab="",
               main=main_title, ylim=ylim, breaks=breaks, ...)
        }
      }

      # ------------ Continuous models ------------
      if(model %in% c("Beta2")){
        mm <- 1/(1+exp(-MM)); sh <- as.numeric(pp["shape"])
        pmin <- as.numeric(pp["min"]); pmax <- as.numeric(pp["max"])
        xr <- .XL(c(pmin, pmax))
        ll <- .GRID(xr, n = 800)
        ds <- ifelse(ll>=pmin & ll<=pmax,
                     dbeta((ll-pmin)/(pmax-pmin), shape1=mm*sh, shape2=(1-mm)*sh)/(pmax-pmin), 0)
        .hist_once("Beta Histogram", xr_default = c(pmin,pmax))
        lines(ll, ds, ...)
      }

      if(model %in% c("Kumaraswamy2")){
        sh <- as.numeric(pp["shape"])
        pmin <- as.numeric(pp["min"]); pmax <- as.numeric(pp["max"])
        xr <- .XL(c(pmin,pmax))
        ll <- .GRID(xr, n = 800)
        dkuma <- function(x,MM,sh,pmin,pmax){
          out <- numeric(length(x))
          idx <- which(x>=pmin & x<=pmax)
          if (length(idx)) {
            q <- (x[idx]-pmin)/(pmax-pmin)
            m1 <- 1/(1+exp(-MM))
            shapei <- log(0.5)/log1p(-m1^sh)
            out[idx] <- shapei*sh*q^(sh-1)*(1-q^sh)^(shapei-1)/(pmax-pmin)
          }
          out
        }
        ds <- dkuma(ll, MM, sh, pmin, pmax)
        .hist_once("Kumaraswamy Histogram", xr_default = c(pmin,pmax))
        lines(ll, ds, ...)
      }

      if(model%in%c("Gaussian")){
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        .hist_once("Gaussian Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, dnorm(ll, mean=MM, sd=sqrt(VV)), ...)
      }

      if(model%in%c("StudentT","Gaussian_misp_StudentT")){
        df <- as.numeric(round(1/pp["df"]))
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        .hist_once("Student T Histogram", xr_default = range(dd, finite=TRUE))
        # Fix: divide by sqrt(VV) inside lines()
        lines(ll, dt((ll - MM)/sqrt(VV), df = df) / sqrt(VV), ...)
      }

      if(model%in%c("SkewStudentT","Gaussian_misp_SkewStudentT")){
        alpha <- as.numeric(pp["skew"]); nu <- as.numeric(round(1/pp["df"]))
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        d_st <- sn::dst((ll-MM)/sqrt(VV), xi=0, omega=1, alpha=alpha, nu=nu)/sqrt(VV)
        .hist_once("Skew-T Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, d_st, ...)
      }

      if(model %in% c("SkewGaussian")){
        skew <- as.numeric(pp["skew"]); sill <- as.numeric(pp["sill"])
        omega <- sqrt((skew^2 + sill)/sill); alpha <- skew/sqrt(sill)
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        d_sn <- sn::dsn((ll-MM)/sqrt(VV), xi=0, omega=omega, alpha=alpha)/sqrt(VV)
        .hist_once("Skew Gaussian Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, d_sn, ...)
      }

      if(model %in% c("Weibull")){
        shape <- pp["shape"]
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        d_w <- ifelse(ll>0, dweibull(ll, shape=shape, scale=exp(MM)/(gamma(1+1/shape))), 0)
        .hist_once("Weibull Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, d_w, ...)
      }

      if(model %in% c("Gamma")){
        shape <- pp["shape"]
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        d_g <- ifelse(ll>0, dgamma(ll, shape=shape/2, rate=shape/(2*exp(MM))), 0)
        .hist_once("Gamma Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, d_g, ...)
      }

      if(model %in% c("LogGaussian")){
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        d_l <- ifelse(ll>0, dlnorm(ll, meanlog = MM - VV/2, sdlog = sqrt(VV)), 0)
        .hist_once("LogGaussian Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, d_l, ...)
      }

      if(model %in% c("Logistic")){
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        d_l <- dlogis(ll, location = MM, scale = sqrt(VV))
        .hist_once("Logistic Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, d_l, ...)
      }

      if(model %in% c("LogLogistic")){
        shape <- pp["shape"]
        cc <- gamma(1+1/shape)*gamma(1-1/shape)
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        d_l <- dllogis1(ll, shape = shape, scale = exp(MM)/cc)  # fix
        .hist_once("LogLogistic Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, d_l, ...)
      }

      if(model %in% c("SinhAsinh")){
        tail <- as.numeric(pp["tail"]); skew <- as.numeric(pp["skew"])
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        dsas <- function(z, skew, tail){
          s <- sinh(tail*asinh(z) - skew)
          a <- (2*pi*(1+z^2))^(-0.5)
          c0 <- sqrt(1+s^2)
          d0 <- exp(-0.5*s^2)
          tail*c0*d0*a
        }
        ds <- dsas((ll-MM)/sqrt(VV), skew, tail)/sqrt(VV)
        .hist_once("SAS Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, ds, ...)
      }



    if(model %in% c("SkewLaplace")){
       skew <- as.numeric(pp["skew"])
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        dslaplace <- function(x, sk) {
        if (sk <= 0 || sk >= 1) stop("sk debe estar en (0,1)")
        f <- numeric(length(x))
        # Parte izquierda: x < 0
        sel1 <- x < 0
        f[sel1] <- sk * (1 - sk) * exp((1 - sk) * x[sel1])
        # Parte derecha: x >= 0
        sel2 <- !sel1
        f[sel2] <- sk * (1 - sk) * exp(-sk * x[sel2])
        return(f)
        }


        #dslaplace <- function(z, skew){
        #aa=1:length(x)
        #sel1=I(x>=0)*aa
        #sel2=I(x<0)*aa
        #x1=x[sel1]        
        #x2=x[sel2]
        #ds1=skew*(1-skew)*exp(-skew*x1)
        #ds2=skew*(1-skew)*exp((1-skew)*x2)
        #return(c(ds2,ds1))
        #}
        ds <- dslaplace((ll-MM)/sqrt(VV), skew)/sqrt(VV)
        .hist_once("Skew Laplace Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, ds, ...)
      }



      if(model %in% c("Tukeyh")){
        inverse_lamb <- function(x,tail) {
          val <- sqrt(VGAM::lambertW(tail*x*x)/tail)
          sign(x)*val
        }
        dTukey <- function(z,tail){
          a <- z*(1 + VGAM::lambertW(tail*z*z))
          b <- inverse_lamb(z,tail)
          c0 <- dnorm(inverse_lamb(z,tail),0,1)
          b*c0/a
        }
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        ds <- dTukey((ll-MM)/sqrt(VV), tail=as.numeric(pp["tail"]))/sqrt(VV)
        .hist_once("Tukey-h Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, ds, ...)
      }

      if(model %in% c("Tukeyh2")){
        tail1 <- as.numeric(pp["tail1"]); tail2 <- as.numeric(pp["tail2"])
        inverse_lamb <- function(x,tail){
          val <- sqrt(VGAM::lambertW(tail*x*x)/tail)
          sign(x)*val
        }
        d_piece <- function(z, tail){
          a <- z*(1 + VGAM::lambertW(tail*z*z))
          b <- inverse_lamb(z,tail)
          c0 <- dnorm(inverse_lamb(z,tail),0,1)
          b*c0/a
        }
        xr <- .XL(range(dd, finite=TRUE))
        ll <- .GRID(xr)
        zs <- (ll-MM)/sqrt(VV)
        ds <- numeric(length(zs))
        idx1 <- which(zs>=0); idx2 <- which(zs<0)
        ds[idx1] <- d_piece(zs[idx1], tail1)
        ds[idx2] <- d_piece(zs[idx2], tail2)
        ds <- ds / sqrt(VV)
        .hist_once("Tukey-hh Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, ds, ...)
      }

      if(model %in% c("TwoPieceGaussian")){
        skew <- as.numeric(pp["skew"])
        xr <- .XL(range(dd, finite=TRUE)); ll <- .GRID(xr)
        z <- (ll-MM)/sqrt(VV)
        # Fix: include Jacobian 1/(1±skew)
        ds <- numeric(length(z))
        idx1 <- which(z>=0); idx2 <- which(z<0)
        ds[idx1] <- dnorm(z[idx1]/(1-skew),0,1) / (sqrt(VV)*(1-skew))
        ds[idx2] <- dnorm(z[idx2]/(1+skew),0,1) / (sqrt(VV)*(1+skew))
        .hist_once("Two-Piece Gaussian Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, ds, ...)
      }

      if(model %in% c("TwoPieceStudentT")){
        skew <- as.numeric(pp["skew"]); df <- 1/as.numeric(pp["df"])
        xr <- .XL(range(dd, finite=TRUE)); ll <- .GRID(xr)
        z <- (ll-MM)/sqrt(VV)
        # Fix: include Jacobian 1/(1±skew)
        ds <- numeric(length(z))
        idx1 <- which(z>=0); idx2 <- which(z<0)
        ds[idx1] <- dt(z[idx1]/(1-skew), df=df) / (sqrt(VV)*(1-skew))
        ds[idx2] <- dt(z[idx2]/(1+skew), df=df) / (sqrt(VV)*(1+skew))
        .hist_once("Two-Piece Student Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, ds, ...)
      }

      if(model %in% c("TwoPieceTukeyh")){
        skew <- as.numeric(pp["skew"]); tail <- as.numeric(pp["tail"])
        inverse_lamb <- function(x,tail){
          val <- sqrt(VGAM::lambertW(tail*x*x)/tail)
          sign(x)*val
        }
        dTukeyh <- function(z,tail){
          a <- z*(1 + VGAM::lambertW(tail*z*z))
          b <- inverse_lamb(z,tail)
          c0 <- dnorm(inverse_lamb(z,tail),0,1)
          b*c0/a
        }
        xr <- .XL(range(dd, finite=TRUE)); ll <- .GRID(xr)
        z <- (ll-MM)/sqrt(VV)
        # Fix: include Jacobian 1/(1±skew)
        ds <- numeric(length(z))
        idx1 <- which(z>=0); idx2 <- which(z<0)
        ds[idx1] <- dTukeyh(z[idx1]/(1-skew), tail) / (sqrt(VV)*(1-skew))
        ds[idx2] <- dTukeyh(z[idx2]/(1+skew), tail) / (sqrt(VV)*(1+skew))
        .hist_once("Two-Piece Tukey-h Histogram", xr_default = range(dd, finite=TRUE))
        lines(ll, ds, ...)
      }

      # ------------ Discrete models ------------
      .disc_support <- function() {
        if (has_xlim) {
          ys <- seq(ceiling(xlim[1]), floor(xlim[2]), by=1)
          if (length(ys) == 0) ys <- sort(unique(as.numeric(dd)))
        } else {
          ys <- sort(unique(as.numeric(dd)))
        }
        ys
      }

      if(model %in% c("BinomialNeg")) {
        count <- as.numeric(table(dd)); ll_emp <- count/sum(count)
        y_emp <- sort(unique(as.numeric(dd)))
        ys <- .disc_support()
        ds <- dnbinom(ys, size=fit$n, prob=pnorm(MM))
        if(!add) plot(y_emp, ll_emp, type="h", col="blue", main="Binomial Negative Histogram",
                      ylab="Density", xlab="", xlim=.XL(range(c(y_emp, ys))), ylim=ylim, lwd=2, ...)
        points(ys, ds, type="p", lwd=3)
        lines(ys, ds)
      }

      if(model %in% c("Binomial")) {
        count <- as.numeric(table(dd)); ll_emp <- count/sum(count)
        y_emp <- sort(unique(as.numeric(dd)))
        ys <- .disc_support()
        ds <- dbinom(ys, size=fit$n, prob=pnorm(MM))
        if(!add) plot(y_emp, ll_emp, type="h", col="blue", main="Binomial Histogram",
                      ylab="Density", xlab="", xlim=.XL(range(c(y_emp, ys))), ylim=ylim, lwd=2, ...)
        points(ys, ds, type="p", lwd=3)
        lines(ys, ds)
      }

      if(model %in% c("BinomialLogistic")) {
        count <- as.numeric(table(dd)); ll_emp <- count/sum(count)
        y_emp <- sort(unique(as.numeric(dd)))
        ys <- .disc_support()
        ds <- dbinom(ys, size=fit$n, prob=plogis(MM))
        if(!add) plot(y_emp, ll_emp, type="h", col="blue", main="Binomial-Logistic Histogram",
                      ylab="Density", xlab="", xlim=.XL(range(c(y_emp, ys))), ylim=ylim, lwd=2, ...)
        points(ys, ds, type="p", lwd=3)
        lines(ys, ds)
      }

      if(model %in% c("Poisson")) {
        count <- as.numeric(table(dd)); ll_emp <- count/sum(count)
        y_emp <- sort(unique(as.numeric(dd)))
        ys <- .disc_support()
        ds <- dpois(ys, lambda=exp(MM))
        if(!add) plot(y_emp, ll_emp, type="h", col="blue", main="Poisson Histogram",
                      ylab="Density", xlab="", xlim=.XL(range(c(y_emp, ys))), ylim=ylim, lwd=2, ...)
        points(ys, ds, type="p", lwd=3)
        lines(ys, ds)
      }

      if(model %in% c("PoissonGamma")) {
        count <- as.numeric(table(dd)); ll_emp <- count/sum(count)
        y_emp <- sort(unique(as.numeric(dd)))
        ys <- .disc_support()
        ff <- exp(MM); sh <- as.numeric(pp["shape"])
        ds <- dnbinom(ys, size=sh, prob=sh/(ff+sh))
        if(!add) plot(y_emp, ll_emp, type="h", col="blue", main="Poisson Gamma (NegBin) Histogram",
                      ylab="Density", xlab="", xlim=.XL(range(c(y_emp, ys))), ylim=ylim, lwd=2, ...)
        points(ys, ds, type="p", lwd=3)
        lines(ys, ds)
      }
    } # end univariate
  } # end D

  invisible()
}
