## Regression tests for the replacement fixed-margin models 85 and 89.
## Run AFTER installing this patched GeoModels:
## Rscript -e 'source(system.file("validation", "Test_FixedMarginST.R", package="GeoModels"))'
## The native-kernel tests were run independently in the development environment.
## This R integration suite requires R and the package dependencies.
local({
  stopifnot(requireNamespace("GeoModels", quietly=TRUE))
  library(GeoModels)
  ns <- asNamespace("GeoModels")
  internal <- function(name) get(name, envir=ns, inherits=FALSE)
  equal <- function(x, y, tolerance=2e-8) {
    stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y),
                              tolerance=tolerance, check.attributes=FALSE)))
  }
  must_error <- function(expr) {
    stopifnot(inherits(tryCatch(force(expr), error=identity), "error"))
  }
  code <- internal("CkCorrModel")
  stopifnot(code("Genwend_Genwend_nosep")==85,
            code("GenWend_GenWend_nosep")==85,
            code("Matern_Matern_nosep")==89)
  stopifnot(identical(internal("CorrelationPar")(85),
            c("scale_s","scale_t","power2_s","power2_t","smooth_s","smooth_t","sep")))
  eval_native <- function(model, h, u, p, d=2L) {
    p <- unlist(p[internal("CorrelationPar")(code(model))], use.names=FALSE)
    .Call("GeoFixedMarginSTEvaluate", as.integer(code(model)),
          as.double(h), as.double(u), as.double(p), as.integer(d), PACKAGE="GeoModels")
  }
  coefficients <- function(p, d=2L, matern=FALSE) {
    if (matern) {
      A <- max(d, 4*p$smooth_s*p$smooth_t)
      B <- max(2*d*p$smooth_t, 2*p$smooth_s)
    } else {
      Ls <- d/(p$power2_s+2*p$smooth_s+d+1)
      Lt <- 1/(p$power2_t+2*p$smooth_t+2)
      Us <- (2*p$smooth_s+1)/p$power2_s
      Ut <- (2*p$smooth_t+1)/p$power2_t
      A <- max(Ls*Lt,Us*Ut); B <- max(Ls*Ut,Us*Lt)
    }
    2*p$sep/(A*(1-p$sep)+B*(1+p$sep))
  }
  parts <- function(x,mu,k) {
    v <- pmax(1-x,0)
    if (k==0) return(list(w=v^mu, D=-x*v^mu))
    list(w=v^(mu+1)*(1+(mu+1)*x), D=-(mu+2)*x^2*v^(mu+1))
  }
  gw_ref <- function(h,u,p,d=2L) {
    a <- parts(h/p$scale_s,p$power2_s,p$smooth_s)
    b <- parts(u/p$scale_t,p$power2_t,p$smooth_t)
    a$w*b$w+coefficients(p,d)*a$D*b$D
  }
  pg <- list(mean=0,sill=1,nugget=0,scale_s=1,scale_t=1,
             power2_s=6,power2_t=6,smooth_s=0,smooth_t=1,sep=-0.7)
  pm <- list(mean=0,sill=1,nugget=0,scale_s=0.7,scale_t=1.5,
             smooth_s=1.5,smooth_t=0.5,sep=0.7)
  models <- c("Genwend_Genwend_nosep","Matern_Matern_nosep")
  ps <- list(pg,pm)
  h <- c(0,1e-6,0.1,0.5,0.99,1,1.3)
  u <- c(0,0.2,0.6,1.2)
  hu <- expand.grid(h=h,u=u)
  for (d in 1:3) for (s in c(-1,-0.5,0,0.5,1)) {
    p <- pg;p$sep <- s
    equal(eval_native(models[1],hu$h,hu$u,p,d),gw_ref(hu$h,hu$u,p,d))
  }
  ## Same margins for every sep, with exact zeros at the compact cutoffs.
  for (m in seq_along(models)) {
    p <- ps[[m]]
    for (s in c(-1,0,1)) {
      q <- p;q$sep <- s
      equal(eval_native(models[m],h,rep(0,length(h)),q),
            eval_native(models[m],h,rep(0,length(h)),p))
      equal(eval_native(models[m],rep(0,length(u)),u,q),
            eval_native(models[m],rep(0,length(u)),u,p))
      z <- GeoCorrFct(h,t=u,corrmodel=models[m],param=q)
      equal(z$corr,eval_native(models[m],hu$h,hu$u,q))
    }
    bad <- p;bad$sep <- 1.01
    must_error(GeoCorrFct(h,t=u,corrmodel=models[m],param=bad))
    must_error(GeoCorrFct(h,t=u,corrmodel=models[m],param=p,distance="Geod"))
    ## A single positive temporal lag must not receive an extra diagonal nugget.
    q <- p;q$nugget <- 0.3;q$sill <- 2
    z <- GeoCorrFct(c(0,0.2),t=0.25,corrmodel=models[m],param=q,covariance=TRUE)
    equal(z$corr,2*0.7*eval_native(models[m],c(0,0.2),c(0.25,0.25),p))
  }
  bad <- pg;bad$smooth_s <- -0.5
  must_error(GeoCorrFct(h,t=u,corrmodel=models[1],param=bad))
  bad <- pg;bad$power2_s <- 3.5;bad$smooth_s <- 0
  must_error(GeoCorrFct(h,t=u,corrmodel=models[1],param=bad))
  ## Explicit CheckCor constraints (also defended inside the native evaluator).
  for (model in models) {
    good <- if (identical(model, models[1])) pg else pm
    for (bad_value in c(-1.01, 1.01, Inf, NaN)) {
      bad <- good; bad$sep <- bad_value
      must_error(eval_native(model, 0.15, 0.2, bad))
    }
    for (which_scale in c("scale_s", "scale_t")) {
      bad <- good; bad[[which_scale]] <- 0
      must_error(eval_native(model, 0.15, 0.2, bad))
    }
  }
  for (which_smooth in c("smooth_s", "smooth_t")) {
    bad <- pm; bad[[which_smooth]] <- 0
    must_error(eval_native(models[2], 0.15, 0.2, bad))
  }
  bad <- pg; bad$smooth_t <- -0.5
  must_error(eval_native(models[1], 0.15, 0.2, bad))
  bad <- pg; bad$power2_t <- pg$smooth_t + 3
  must_error(eval_native(models[1], 0.15, 0.2, bad))
  ## Public interfaces are identical to the original supplied package.
  stopifnot(identical(names(formals(GeoCorrFct)),
    c("x", "t", "corrmodel", "model", "distance", "param", "radius", "n",
      "covariance", "variogram")))
  stopifnot(identical(names(formals(GeoCorrFct_Cop)),
    c("x", "t", "corrmodel", "model", "copula", "distance", "param", "radius",
      "n", "covariance", "variogram")))

  ## Dimension inference and restoration, using only existing coordinates/metric.
  equal(internal(".GeoSTCurrentDimension")(),2)
  xy_dim <- matrix(c(0,0, 0.3,0.2), ncol=2, byrow=TRUE)
  xyz_dim <- cbind(xy_dim,c(0.1,0.4))
  equal(internal(".GeoSTDimension")(coordx=xy_dim),2)
  equal(internal(".GeoSTDimension")(coordx=xyz_dim),3)
  equal(internal(".GeoSTDimension")(coordx=xyz_dim[,1],
    coordy=xyz_dim[,2],coordz=xyz_dim[,3]),3)
  equal(internal(".GeoSTDimension")(coordx_dyn=list(xyz_dim,xyz_dim)),3)
  equal(internal(".GeoSTDimension")(distance="Chor"),3)

  with_coordinates <- function(model, coordinates, expr) {
    old <- internal(".GeoSTEnter")(model,coordx=coordinates)
    on.exit(internal(".GeoSTRestore")(old),add=TRUE)
    force(expr)
  }
  for (m in seq_along(models)) {
    p <- ps[[m]]
    ## Existing chordal metric selects dimension three without an extra argument.
    z <- GeoCorrFct(h,t=u,corrmodel=models[m],param=p,distance="Chor")
    equal(z$corr,eval_native(models[m],hu$h,hu$u,p,3))
    equal(internal(".GeoSTCurrentDimension")(),2)
    zc <- GeoCorrFct_Cop(h,t=u,corrmodel=models[m],param=p)
    equal(zc$corr,eval_native(models[m],hu$h,hu$u,p))
    zc <- GeoCorrFct_Cop(h,t=u,corrmodel=models[m],param=p,distance="Chor")
    equal(zc$corr,eval_native(models[m],hu$h,hu$u,p,3))
    equal(internal(".GeoSTCurrentDimension")(),2)

    ## A curve inside a coordinate-based operation inherits the 3D context.
    with_coordinates(models[m],xyz_dim,{
      equal(internal(".GeoSTCurrentDimension")(),3)
      z <- GeoCorrFct(h,t=u,corrmodel=models[m],param=p)
      equal(z$corr,eval_native(models[m],hu$h,hu$u,p,3))
      zc <- GeoCorrFct_Cop(h,t=u,corrmodel=models[m],param=p)
      equal(zc$corr,eval_native(models[m],hu$h,hu$u,p,3))
      ## A nested 2D coordinate operation must not erase the enclosing 3D scope.
      with_coordinates(models[m],xy_dim,{
        z <- GeoCorrFct(h,t=u,corrmodel=models[m],param=p)
        equal(z$corr,eval_native(models[m],hu$h,hu$u,p,2))
      })
      equal(internal(".GeoSTCurrentDimension")(),3)
      bad <- p;bad$sep <- 2
      must_error(GeoCorrFct(h,t=u,corrmodel=models[m],param=bad))
      must_error(GeoCorrFct_Cop(h,t=u,corrmodel=models[m],param=bad))
      equal(internal(".GeoSTCurrentDimension")(),3)
    })
    equal(internal(".GeoSTCurrentDimension")(),2)
    bad <- p;bad$sep <- 2
    must_error(with_coordinates(models[m],xyz_dim,
      GeoCorrFct(h,t=u,corrmodel=models[m],param=bad)))
    equal(internal(".GeoSTCurrentDimension")(),2)
    ## No dimension setting leaks into later standalone Euclidean evaluations.
    z <- GeoCorrFct(h,t=u,corrmodel=models[m],param=p)
    equal(z$corr,eval_native(models[m],hu$h,hu$u,p,2))
  }

  set.seed(27102026)
  xy <- cbind(runif(6),runif(6))
  times <- c(0,0.25,0.7)
  n <- nrow(xy);nt <- length(times)
  coords <- xy[rep(seq_len(n),nt),,drop=FALSE]
  ts <- rep(times,each=n)
  hs <- as.matrix(dist(coords));us <- abs(outer(ts,ts,"-"))
  for (m in seq_along(models)) {
    p <- ps[[m]]
    cv <- GeoCovmatrix(coordx=xy,coordt=times,corrmodel=models[m],param=p)
    ref <- matrix(eval_native(models[m],c(hs),c(us),p),n*nt)
    equal(as.matrix(cv$covmatrix),ref)
    if (m==1L) {
      cvs <- GeoCovmatrix(coordx=xy,coordt=times,corrmodel=models[m],param=p,sparse=TRUE)
      equal(as.matrix(cvs$covmatrix),ref)
    }
    ## Three-dimensional parameter normalization and restoration to 2D.
    xyz <- cbind(xy,runif(n));c3 <- xyz[rep(seq_len(n),nt),,drop=FALSE]
    h3 <- as.matrix(dist(c3))
    cv3 <- GeoCovmatrix(coordx=xyz,coordt=times,corrmodel=models[m],param=p)
    equal(as.matrix(cv3$covmatrix),matrix(eval_native(models[m],c(h3),c(us),p,3),n*nt))
    equal(internal(".GeoSTCurrentDimension")(),2)
    ## Irregular/dynamic design is a submatrix of the complete design.
    subsets <- list(1:6,c(1,3,5),c(2,4,6))
    dyn <- lapply(subsets,function(ii) xy[ii,,drop=FALSE])
    keep <- unlist(Map(function(ii,j) ii+(j-1L)*n,subsets,seq_len(nt)))
    cvd <- GeoCovmatrix(coordx_dyn=dyn,coordt=times,corrmodel=models[m],param=p)
    equal(as.matrix(cvd$covmatrix),ref[keep,keep])
    ## Simulate a Gaussian vector, then verify public kriging against dense algebra.
    y <- drop(t(chol(ref)) %*% rnorm(n*nt))
    data <- matrix(y,nrow=nt,byrow=TRUE)
    loc <- matrix(c(0.35,0.45),1,2);t0 <- 0.4
    k <- eval_native(models[m],sqrt(rowSums(sweep(coords,2,loc[1,])^2)),abs(ts-t0),p)
    kr <- GeoKrig(data=data,coordx=xy,coordt=times,corrmodel=models[m],
                  param=p,loc=loc,time=t0,mse=TRUE)
    equal(kr$pred,drop(crossprod(k,solve(ref,y))))
    equal(kr$mse,1-drop(crossprod(k,solve(ref,k))))
    ## Full likelihood, interaction-only optimization. No assertion of statistical
    ## recovery from this tiny sample: check finite output and the parameter domain.
    fixed <- p;fixed$sep <- NULL
    fit <- GeoFit(data=data,coordx=xy,coordt=times,corrmodel=models[m],
                  fixed=fixed,start=list(sep=0),lower=list(sep=-1),upper=list(sep=1),
                  likelihood="Full",type="Standard",optimizer="optimize",varest=FALSE)
    stopifnot(is.finite(fit$logCompLik),is.finite(fit$param$sep),abs(fit$param$sep)<=1)
    if (identical(Sys.getenv("GEOMODELS_TEST_CL"),"true")) {
      fitcl <- GeoFit(data=data,coordx=xy,coordt=times,corrmodel=models[m],
         fixed=fixed,start=list(sep=0),lower=list(sep=-1),upper=list(sep=1),
         likelihood="Marginal",type="Pairwise",optimizer="optimize",neighb=3,
         maxtime=0.7,score=TRUE,varest=FALSE)
      stopifnot(is.finite(fitcl$logCompLik),abs(fitcl$param$sep)<=1)
    }
  }
  cat("Fixed-margin space-time R regression checks passed.\n")
})
