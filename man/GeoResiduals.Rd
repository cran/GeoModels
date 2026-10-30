\name{GeoResiduals}
\alias{GeoResiduals}
\encoding{UTF-8}
\title{Compute Residuals from a Fitted GeoModels Model}
\description{
 The procedure return a GeoFit object associated to the estimated residuals. For a random field Y defined on the real line (Gaussian, Skew Gaussian, Tukeyh etcc) they are computed as (Y-m)/sqrt(v) where m and v are the estimated mean and variance respectively.
 For a random field Y defined on the positive real line (Gamma, Weibull, Log-Gaussian) they are computed as Y/m where m is the estimated mean.
 In the first case residuals have zero mean and unut variance with a specific distribution defined on the real line.
 In the second case residuals have unit mean with a specific distribution defined on the positive real line.
 For the discrete Bernoulli/Binary, Binomial, Geometric/negative-binomial,
 Poisson, and Poisson--Gamma models, the function returns Pearson residuals
 computed from the fitted marginal mean and variance. The paper specifies these marginal moments but does
 not prescribe a residual definition; Pearson residuals are an explicit package
 convention.
For continuous models, the returned object can be coupled with
\code{\link{GeoQQ}} and \code{\link{GeoCovariogram}}. For discrete Pearson
residuals, use the returned object for residual maps and variograms, and apply
\code{GeoQQ()} or \code{GeoPit()} to the original fitted object for marginal
assessment.
}
\usage{GeoResiduals(fit)}
\arguments{
 \item{fit}{A fitted object obtained from the
 \code{\link{GeoFit}}.}
}

\details{

The fitted location is evaluated from the full mean specification: a
site-specific fixed mean vector when present, otherwise \eqn{X\beta} with
coefficients \code{mean}, \code{mean1}, and so on matched in order to the
columns of \code{X}.  The returned residual object is converted to an
intercept-only, standardized representation for subsequent diagnostics.  All
other fields of the original \code{GeoFit} object, including neighbourhood,
thinning, anisotropy, convergence, and pair-selection metadata, are preserved.

For \code{Beta2} and \code{Kumaraswamy2}, the residual transformation uses
the fitted support \code{min}/\code{max} and the corresponding fitted marginal
CDF, so the standardized residual scale is the unit interval.

For the discrete models described in the package paper, Pearson residuals are
used: \eqn{(Y-E(Y))/\sqrt{Var(Y)}}.  Thus \eqn{p=\Phi(\mu)} for
Bernoulli, Binomial and negative-binomial fields, while
\eqn{\lambda=\exp(\mu)} for Poisson fields.  For Poisson--Gamma fields,
with fitted mean \eqn{\lambda} and shape \eqn{\kappa}, the fitted variance is
\eqn{\lambda(1+\lambda/\kappa)}.  These deterministic residuals
are intended for residual maps and variograms; they are not randomized quantile
residuals, and \code{GeoQQ()} intentionally rejects the Pearson-residual object.  Bivariate residuals are currently supported only
for Gaussian models.
}

\value{
 Returns an (updated) object of class \code{GeoFit}.  The component
 \code{residual_type} records whether the data are Pearson or model-standardized residuals.
}


\seealso{\code{\link{GeoFit}}.}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}



\examples{
library(GeoModels)



###########################
###Example 1: Residuals using a Gaussian RF
###########################
set.seed(211)
model="Gaussian";
N=150 # number of location sites
# Set the coordinates of the points:
x = runif(N, 0, 1)
y = runif(N, 0, 1)
coords=cbind(x,y)

# regression parameters
mean = 5
mean1=0.8

X=cbind(rep(1,N),runif(N))
# correlation parameters:
corrmodel = "Wend0"
sill = 1
nugget = 0
scale = 0.3
power2=4

param=list(mean=mean,mean1=mean1, sill=sill, nugget=nugget, 
 scale=scale,power2=power2)
# Simulation of the Gaussian RF:
data = GeoSim(coordx=coords, corrmodel=corrmodel, X=X,model=model,param=param)$data

start=list(mean=mean,mean1=mean1, scale=scale,sill=sill)
fixed=list(nugget=nugget,power2=power2)
# Maximum composite-likelihood fitting 
fit = GeoFit(data,coordx=coords, corrmodel=corrmodel,model=model,X=X,
 likelihood="Conditional",type='Pairwise',start=start,
 fixed=fixed,neighb=3)

res=GeoResiduals(fit)
mean(res$data) # should be approx 0
var(res$data) # should be approx 1
# checking goodness of fit marginal model
GeoQQ(res);GeoQQ(res,type="D",col="red",ylim=c(0,0.5),breaks=20);
# Empirical estimation of the variogram for the residuals:
vario = GeoVariogram(res$data,coordx=coords,maxdist=0.5)
# Comparison between empirical amd estimated semivariogram for the residuals
GeoCovariogram(res, show.vario=TRUE, vario=vario,pch=20)





###########################
###Example 2: Residuals using a Weibull RF
###########################
model="Weibull";shape=4
N=150 # number of location sites
# Set the coordinates of the points:
x = runif(N, 0, 1)
y = runif(N, 0, 1)
coords=cbind(x,y)


# regression parameters
mean = 5
mean1=0.8

X=cbind(rep(1,N),runif(N))
# correlation parameters:
corrmodel = "Wend0"
sill = 1
nugget = 0
scale = 0.3
power2=4

param=list(mean=mean,mean1=mean1, sill=sill, nugget=nugget, 
	 scale=scale,shape=shape,power2=power2)
# Simulation of the Gaussian RF:
data = GeoSim(coordx=coords, corrmodel=corrmodel, X=X,model=model,param=param)$data

I=Inf
start=list(mean=mean,mean1=mean1, scale=scale,shape=shape)
lower=list(mean=-I,mean1=-I, scale=0,shape=0)
upper=list(mean= I,mean1= I, scale=I,shape=I)
fixed=list(nugget=nugget,sill=sill,power2=power2)
# Maximum composite-likelihood fitting 
fit = GeoFit(data,coordx=coords, corrmodel=corrmodel,model=model,X=X,
 likelihood="Conditional",type='Pairwise',start=start,
 optimizer="nlminb", lower=lower,upper=upper,
 fixed=fixed,neighb=3)


res=GeoResiduals(fit)
mean(res$data) # should be approx 1
# checking goodness of fit marginal model
GeoQQ(res);GeoQQ(res,type="D",lwd=2,ylim=c(0,1.7),breaks=20);
# Empirical estimation of the variogram for the residuals:
vario = GeoVariogram(res$data,coordx=coords,maxdist=0.5)
# Comparison between empirical amd estimated semivariogram for the residuals
GeoCovariogram(res, show.vario=TRUE, vario=vario,pch=20)

}

\keyword{Composite}
