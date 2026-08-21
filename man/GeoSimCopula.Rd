\name{GeoSimCopula}
\alias{GeoSimCopula}
\alias{print.GeoSimCopula}
\encoding{UTF-8}
\title{Simulation of Gaussian and non-Gaussian random fields using copula.}
\description{
 Simulation of Gaussian and some non-Gaussian univariate spatial and
 spatio-temporal random fields using Gaussian, skew-Gaussian, Clayton-like or
 AMH constructions. Bivariate correlation models are not supported.
 The function returns a realization of a random field for a given covariance model and covariance parameters.
 Exact simulation is available through Cholesky or SVD decomposition; for supported purely spatial correlation models, turning-bands simulation is available through \code{method="TB"}.
}
\usage{
GeoSimCopula(coordx=NULL, coordy=NULL,coordz=NULL, coordt=NULL, 
coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, 
 method="cholesky", model='Gaussian', n=1, param,
 anisopars=NULL,radius=1, sparse=FALSE,
 copula="Gaussian",X=NULL,spobj=NULL,nrep=1,progress=FALSE,check.duplicates=FALSE,
 L=10000,parallel=FALSE,ncores=6)
}
\arguments{
 \item{coordx}{Optional when \code{coordx_dyn} or \code{spobj} is supplied. A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix 
 Coordinates on a sphere for a fixed radius \code{radius} 
 are passed in lon/lat format expressed in decimal degrees.}
 \item{coordy}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordz}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}

 \item{coordt}{A numeric vector giving the temporal coordinates at which the field is simulated. Optional argument; the default is \code{NULL}, in which case a spatial random field is expected. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}
 \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
 containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
 }
 \item{corrmodel}{String; the name of a correlation model, for the
 see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}
 \item{distance}{String; the name of the spatial distance. The default
 is \code{Eucl}, the Euclidean distance. See \code{\link{GeoFit}} for details.}
 \item{grid}{Logical; if \code{FALSE} (the default) the data
 are interpreted as spatial or spatio-temporal realisations on a set
 of non-equispaced spatial sites (irregular grid).}
 \item{method}{String; simulation engine for the latent fields. Use \code{"cholesky"} (default) or \code{"svd"} for exact simulation through \code{GeoSim}, and \code{"TB"} for turning-bands simulation through \code{GeoSimapprox}.} 
 \item{model}{String; the type of RF and therefore the densities associated to the likelihood
 objects. \code{Gaussian} is the default, see the Section
 Details.}
 \item{n}{Numeric; the number of trials for binomial random fields. The number of successes in the negative Binomial random fields. Default is \eqn{1}.} 
 \item{param}{A list of parameter values required in the
 simulation procedure of random fields, see \bold{Examples}.}
 \item{anisopars}{ A list of two elements "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
 \item{radius}{Numeric; a value indicating the radius of the sphere when using the great 
 circle distance. Default value is 1.} 
 \item{sparse}{Logical; if \code{TRUE} then cholesky decomposition is performed
 using sparse matrices algorithms (spam packake). It should be used
 with compactly supported covariance models.\code{FALSE} is the default.}
 \item{copula}{String; one of \code{"Gaussian"}, \code{"SkewGaussian"}, \code{"Clayton"}, or \code{"AMH"}.}
 \item{X}{Numeric; Matrix of space-time covariates.}
 \item{spobj}{An object of class \code{sp} or \code{spacetime}. For space-time objects, the current \code{sp2Geo()} conversion uses sequential temporal indices and does not preserve irregular original time spacing; to retain irregular temporal distances, use the explicit-coordinate interface with numeric \code{coordt}.}
 \item{nrep}{Numeric; Numbers of indipendent replicates.}
 \item{progress}{Logic; If TRUE then a progress bar is shown.}
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
 \item{L}{Positive integer; number of turning-band lines used when \code{method="TB"}. Ignored by exact methods.}
 \item{parallel}{Logical; default \code{FALSE}. For \code{method="TB"}, set \code{TRUE} to allow \code{GeoSimapprox} to select a parallel turning-bands execution path.}
 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE} and \code{method="TB"}, an explicit integer requests that many workers for the TB backend, subject to detected-core and job-count limits. Set \code{ncores=NULL} for automatic selection.}
}


\details{

Only univariate correlation models are accepted. Dynamic coordinates are
supported and the returned data preserve the list structure and row order of
\code{coordx_dyn}.

With \code{method="TB"}, \code{GeoSimCopula} keeps the copula and marginal
construction unchanged but generates each required latent field through
\code{GeoSimapprox(method="TB")}. Consequently, TB is available only for
purely spatial models supported by the turning-bands backend, requires
\code{distance="Eucl"}, and inherits the TB restrictions on correlation
families and coordinate dimension. The \code{sparse} argument is ignored by
the TB backend. Exact \code{"cholesky"} and \code{"svd"} simulations continue
to use \code{GeoSim}.

For the constructive Clayton-like simulator, \code{param$nu} must be a
positive integer and is never rounded silently. For \code{copula="SkewGaussian"},
\code{param$nu} is the bounded reflection-asymmetry parameter
\eqn{\eta\in(-1,1)} and is internally mapped to
\eqn{\gamma_\eta=\eta/\sqrt{1-\eta^2}}.  The continuous margins audited
for copula simulation and pairwise copula fitting with the Gaussian,
Clayton-like and skew--Gaussian copulas are
\code{Gaussian}, \code{StudentT}, \code{LogGaussian}, \code{Gamma},
\code{Weibull}, \code{Beta}, \code{Beta2}, \code{Kumaraswamy},
\code{Kumaraswamy2}, \code{Logistic}, and \code{SkewLaplace}.
The Gaussian copula additionally has audited simulation and pairwise fitting for
\code{Poisson}, \code{Binomial}, and \code{BinomialNeg}. Exact conditional
simulation for these discrete copula margins is not currently implemented.
For \code{model="StudentT"},
\code{param$df} is the reciprocal degrees-of-freedom parameter: the Student-t
degrees of freedom are exactly \code{1/param$df}, without integer rounding.
The same no-rounding convention is used for the Student-t degrees of freedom
of a \code{"SkewStudentT"} marginal. Because copula margins are obtained through
inverse CDFs, the integer latent-field restrictions used by direct \code{GeoSim} /
\code{GeoSimapprox} Gamma, Beta, and Student-t constructions do not apply here.
Gamma margins require a finite
positive \code{shape} and use shape \code{shape/2} and scale
\code{2/shape}, hence the multiplicative factor has unit mean. With dynamic
coordinates, \code{X} may be a row-stacked matrix or a list of matrices aligned
with \code{coordx_dyn}.  For \code{SkewStudentT}, \code{skew} is the native
\eqn{\delta\in(-1,1)} parameter and is converted to the skew-t shape
\eqn{\alpha=\delta/\sqrt{1-\delta^2}}.  The Tukey-g-and-h implementation
uses its continuous limit when the skew parameter is zero.  Kumaraswamy margins
use the inverse of the CDF employed by \code{GeoPit}.

Aliases accepted by the main modelling API, including \code{Gauss},
\code{SkewGauss}, \code{LogGauss}, and \code{TwoPieceGauss}, are canonicalized
before the marginal transform. \code{nrep} must be a positive integer. For the AMH
construction, \code{nu} is checked to be a finite scalar; no additional theoretical
parameter range is imposed by this simulator. The returned \code{param} component is
the parameter list supplied by the user.

\code{Binary} and \code{Bernoulli} are aliases of \code{Binomial} with
\code{n=1}; \code{Geom} and \code{Geometric} are aliases of
\code{BinomialNeg} with \code{n=1}.

For univariate simulation, the marginal location is
\eqn{\mu=X\beta}, with coefficients named \code{mean}, \code{mean1}, and so
on in the order of the columns of \code{X}.  If \code{X=NULL}, the model is
intercept-only.  A vector \code{param$mean} supplies a known location value at
each observation and cannot be combined with \code{X}.
}

\value{
 Returns an object of class \code{GeoSimCopula}.
 An object of class \code{GeoSimCopula} is a list containing
 at most the following components:
 
 \item{bivariate}{Always \code{FALSE}; bivariate correlation models are rejected.}
 \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
 \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
 \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
 \item{coordx_dyn}{A list of dynamical (in time) spatial coordinates;}
 \item{corrmodel}{The correlation model; see \code{\link{GeoCovmatrix}}.}
 \item{data}{The simulated data. Dynamic space-time output is a list with one vector per time, aligned with \code{coordx_dyn}; fixed-location space-time output preserves the matrix layout of \code{GeoSim}.}
 \item{distance}{The type of spatial distance;}
 \item{method}{The method of simulation}
 \item{model}{The type of RF, see \code{\link{GeoFit}}.}
 \item{n}{The number of trial for Binomial random fields;the number of successes in a negative Binomial random fields;} 
 \item{numcoord}{The number of spatial coordinates;}
 \item{numtime}{The number the temporal realisations of the RF;}
 \item{param}{A list of the parameters}
 \item{radius}{The radius of the sphere if coordinates are passed in lon/lat format;}
 \item{randseed}{The seed used for the random simulation;}
 \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if
 spatial RF;}
 \item{copula}{The type of copula}
}



\references{

 Bevilacqua M., Alvarado E., Caamano C. (2024) A flexible Clayton-like spatial copula with application to bounded support data. \emph{Journal of Multivariate Analysis} \bold{201}


}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{
library(GeoModels)

################################################################
###
### Example: Simulation of a reparametrized Beta RF
### for beta regression
### with Gaussian and Clayton Copula 
### with underlying Wendland correlation.
###
###############################################################
set.seed(261)
NN=1400
x <- runif(NN);y <- runif(NN)
coords=cbind(x,y)

corrmodel="GenWend"
X=cbind(rep(1,NN),runif(NN))

NuisParam("Beta2",num_betas=2,copula="Gaussian")
CorrParam("GenWend")
#### Gaussian copula
param=list(smooth=0,power2=4, min=0,max=1,
 mean=0.1,mean1=0.1,scale=0.3,nugget=0,shape=5)

data <- GeoSimCopula(coordx=coords, corrmodel=corrmodel, model="Beta2",param=param,
 copula="Gaussian",sparse=TRUE,X=X)$data

if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(coords,data)


#### Clayton copula
NuisParam("Beta2",num_betas=2,copula="Clayton")
CorrParam("GenWend")
param=list(smooth=0,power2=4, min=0,max=1,
 mean=0.2,mean1=0.1,scale=0.3,nugget=0,shape=6,nu=4)
data1 <- GeoSimCopula(coordx=coords, corrmodel=corrmodel, model="Beta2",param=param,
 copula="Clayton",sparse=TRUE,X=X)$data

hist(data1,freq=FALSE)
if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(coords,data1)

}
\keyword{Simulation}
