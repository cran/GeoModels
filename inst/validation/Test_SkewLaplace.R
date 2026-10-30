# Post-installation checks for GeoModels 2.3.0.1.
# These are reproducible R integration tests supplied for the user's environment.
# They were NOT run as an R script in the build environment (R unavailable).
# Native algorithm tests run during preparation are described in the report.
library(GeoModels)
stopifnot(utils::packageVersion("GeoModels") >= package_version("2.3.0.1"))
cat("Loaded GeoModels", as.character(utils::packageVersion("GeoModels")),
    "from", find.package("GeoModels"), "\n")
D <- getFromNamespace(".GeoSkewLaplaceLogDensity", "GeoModels")

# 1. Exact independence and exact density at zero (formerly broken).
logmarg <- function(x, p) log(p) + log1p(-p) + ifelse(x >= 0, -p*x, (1-p)*x)
x <- c(-5, -.2, 0, .1, 1, 10)
y <- rev(x)
for (p in c(.1, .2, .5, .75, .9)) {
  stopifnot(max(abs(D(x, y, 0, p) - logmarg(x,p) - logmarg(y,p))) < 1e-11)
  for (r in c(.1, .7, .95, .99)) {
    f00 <- p^2*(1-p)^2 / sqrt((1-r^2)*(1-r^2*(2*p-1)^2))
    stopifnot(abs(D(0,0,r,p)-log(f00)) < 1e-10,
              max(abs(D(x,y,r,p)-D(y,x,r,p))) < 1e-10,
              max(abs(D(x,y,r,p)-D(-x,-y,r,1-p))) < 1e-9,
              max(abs(D(x,y,r,p)-D(x,y,-r,p))) < 1e-10)
  }
}
stopifnot(abs(D(1,2,.99,.2,log=FALSE) / .0237389248007547 - 1) < 1e-9)
stopifnot(abs(D(2,4,.7,.2,sill=4)-D(1,2,.7,.2)+log(4)) < 1e-10)

# 2. Integration over one component recovers the correct marginal density.
for (r in c(.2, .9, .99)) {
  integrand <- function(z) D(rep(1,length(z)), z, r, .2, log=FALSE)
  value <- integrate(integrand, -Inf, 0, rel.tol=1e-7)$value +
           integrate(integrand, 0, Inf, rel.tol=1e-7)$value
  expected <- exp(logmarg(1,.2))
  stopifnot(abs(value/expected-1) < 2e-6)
}
cat("Bivariate density: independence, zero, symmetry, reflection, scale, high correlation and marginals passed.\n")

# 3. TWO correlated observations: posterior radial and angular-product moments.
# References are independent 128/256-node Kibble quadrature, not the Gibbs code.
reference <- list(
  list(p=0.75,r=0.6,y=c(1,-0.7),target=c(0.8847424204898078,1.073269716319773,1.325981518937683,0.1845742808057538)),
  list(p=0.2,r=0.85,y=c(1,2),target=c(0.8895962651025641,0.8429300333437326,0.7256373590617335,1.104338062580795)),
  list(p=0.5,r=0.95,y=c(-1,-2),target=c(1.064802693352979,0.8543856203504906,0.8653019490091481,2.35760876413698))
)
set.seed(20261027)
for (case in reference) {
  nrep <- 30000L
  Q <- solve(matrix(c(1,case$r,case$r,1),2,2))
  z <- .Call("geo_skewlaplace_latent_gibbs", as.double(case$y), Q,
             as.double(case$p), 3000L, nrep, 3L, PACKAGE="GeoModels")
  A <- (z[,1,]^2 + z[,2,]^2)/2
  B <- (z[,3,]^2 + z[,4,]^2)/2
  residual <- sweep(A/case$p-B/(1-case$p),1,case$y,"-")
  stopifnot(max(abs(residual)) < 1e-10)
  radial <- sweep(A/case$p,1,pmax(case$y,0),"-")
  values <- cbind(t(radial),
                  z[1,1,]*z[2,1,]+z[1,2,]*z[2,2,],
                  z[1,3,]*z[2,3,]+z[1,4,]*z[2,4,])
  nb <- 100L; block <- nrep %/% nb
  means <- vapply(seq_len(nb),function(k)
    colMeans(values[(k-1L)*block+seq_len(block),,drop=FALSE]), numeric(4))
  se <- apply(means,1,stats::sd)/sqrt(nb)
  estimate <- colMeans(values)
  print(data.frame(p=case$p,r=case$r,moment=c("T1","T2","plus product","minus product"),
                   estimate=estimate,reference=case$target,MCSE=se))
  stopifnot(all(abs(estimate-case$target) < 6*se + .001))
}
cat("Multi-observation latent Gibbs tests passed.\n")

# 4. Full R -> C -> R conditional paths, with covariates, exact matches,
# duplicate prediction sites, and a prediction site outside compact support.
set.seed(920)
n <- 40L
coords <- cbind(runif(n),runif(n))
X <- cbind(1,runif(n,-1,1))
pars <- list(mean=.1,mean1=-.25,nugget=0,sill=1,skew=.75,
             scale=.3,power2=4,smooth=0)
y_direct <- as.numeric(GeoSim(coordx=coords,corrmodel="GenWend",
                    model="SkewLaplace",param=pars,X=X,progress=FALSE)$data)
y_copula <- as.numeric(GeoSimCopula(coordx=coords,corrmodel="GenWend",
                    model="SkewLaplace",copula="Gaussian",param=pars,X=X,
                    progress=FALSE)$data)
loc <- rbind(coords[1:2,,drop=FALSE], c(.33,.52), c(.33,.52), c(4,4))
Xloc <- rbind(X[1:2,,drop=FALSE], c(1,.2), c(1,.2), c(1,.3))
base <- list(coordx=coords,loc=loc,X=X,Xloc=Xloc,corrmodel="GenWend",
             model="SkewLaplace",param=pars,method="Cholesky",nrep=2000L,
             local=FALSE,progress=FALSE)
direct <- do.call(GeoSimcond,c(base,list(data=y_direct,copula=NULL,
                                      n_iter=1500L,mcmc_thin=3L)))
copula <- do.call(GeoSimcond,c(base,list(data=y_copula,copula="Gaussian")))
for (tag in c("direct","copula")) {
  obj <- get(tag); observations <- if(tag=="direct") y_direct else y_copula
  draws <- do.call(rbind,obj$condsim)
  stopifnot(identical(dim(draws),c(2000L,5L)),all(is.finite(draws)),
            max(abs(draws[,1]-observations[1])) < 1e-8,
            max(abs(draws[,2]-observations[2])) < 1e-8,
            max(abs(draws[,3]-draws[,4])) < 1e-10,
            obj$simulation_details$observation_factorizations == 1L,
            obj$simulation_details$prediction_factorizations == 1L)
  eta <- drop(Xloc[5,,drop=FALSE] %*% c(pars$mean,pars$mean1))
  mu <- eta + sqrt(pars$sill)*(1/pars$skew-1/(1-pars$skew))
  se <- stats::sd(draws[,5])/sqrt(nrow(draws))
  stopifnot(abs(mean(draws[,5])-mu) < 6*se)
  prob <- mean(draws[,5] <= eta)
  stopifnot(abs(prob-pars$skew) < .06)
  cat(tag,": exact matches, repeated sites, far-site marginal, cached factors passed.\n")
}

# 5. Independence likelihood agrees in check loss with rq(tau=skew).
if (requireNamespace("quantreg",quietly=TRUE)) {
  fit <- GeoFit(data=y_direct,coordx=coords,corrmodel="GenWend",X=X,
                model="SkewLaplace",likelihood="Marginal",type="Independence",
                start=list(mean=.1,mean1=-.25),
                fixed=pars[setdiff(names(pars),c("mean","mean1"))])
  beta <- unlist(fit$param[c("mean","mean1")],use.names=FALSE)
  qr <- quantreg::rq(y_direct ~ X[,2],tau=pars$skew)
  loss <- function(b) {u <- y_direct-drop(X%*%b);sum(u*(pars$skew-(u<0)))}
  print(rbind(GeoFit=beta,rq=unname(coef(qr))))
  stopifnot(abs(loss(beta)-loss(coef(qr))) < 1e-4*(1+loss(coef(qr))))
  cat("Independence likelihood / quantile regression check passed.\n")
} else cat("quantreg not installed: optional rq comparison skipped.\n")

cat("ALL REQUESTED POST-INSTALLATION CHECKS PASSED.\n")
