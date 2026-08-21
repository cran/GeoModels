\name{GeoSimapprox}
\alias{GeoSimapprox}
\encoding{UTF-8}
\title{Fast simulation of Gaussian and non-Gaussian random fields.}
\description{
 Simulation of Gaussian and some non-Gaussian spatial, spatio-temporal and spatial bivariate 
 random fields using two approximate methods of simulation: circulant embeeding and spectral turning band. (see Examples). 
}
\usage{
GeoSimapprox(coordx=NULL, coordy=NULL, coordz=NULL,coordt=NULL, 
coordx_dyn=NULL,corrmodel, distance="Eucl",
grid=FALSE, max.ext=1,
method="TB", L=10000,model='Gaussian',parallel=FALSE,ncores=6,
n=1,param,anisopars=NULL, radius=6371,X=NULL,spobj=NULL,
nrep=1,progress=TRUE,check.duplicates=FALSE)
}
\arguments{
 \item{coordx}{May be \code{NULL} when \code{coordx_dyn} or \code{spobj} is supplied. Otherwise, a numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix 
 Coordinates on a sphere for a fixed radius \code{radius} 
 are passed in lon/lat format expressed in decimal degrees.}
 \item{coordy}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordz}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}

 \item{coordt}{A numeric vector giving the temporal coordinates. Optional argument; the default is \code{NULL}, in which case a spatial random field is expected. For space-time \code{method = "CE"}, \code{coordt} must be finite, strictly increasing, and equally spaced; space-time turning-bands simulation is not implemented.}
 \item{coordx_dyn}{For dynamic simulation sites, a list with one two- or three-column coordinate matrix per element of \code{coordt}.  The output at time \eqn{t} follows the row order of \code{coordx_dyn[[t]]}.  See \code{\link{GeoModels-spacetime-ordering}}.}
 \item{corrmodel}{String; the name of a correlation model, for the
 see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}
 \item{parallel}{Logical; default \code{FALSE}. For \code{method = "TB"}, \code{TRUE} enables an
 adaptive scheduler that chooses among serial execution, parallelization over
 spatial chunks within each turning-bands simulation, and parallelization over
 independent replicates. Small workloads can remain serial when worker overhead
 is expected to dominate.}
 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE}, an explicit integer requests that many workers, capped by detected cores and the available TB jobs. Set \code{ncores=NULL} for automatic selection, capped at six workers and normally leaving one detected core free.}
 \item{distance}{String; the name of the spatial distance. The default
 is \code{Eucl}, the Euclidean distance. Turning-bands simulation requires
 \code{distance="Eucl"}. See \code{\link{GeoFit}} for details.}
 \item{grid}{Logical; if \code{FALSE} (the default) the data
 are interpreted as spatial or spatio-temporal realisations on a set
 of non-equispaced spatial sites (irregular grid).}
 \item{max.ext}{Positive integer; maximum number of successive doubling attempts used by the spatial or temporal CE embedding.}
 \item{method}{String; the approximation method. The default is \code{TB}
 (turning bands). The alternative \code{CE} uses circulant embedding. Spatial
 CE requires an exact regular two-dimensional grid, Euclidean distance and no
 \code{anisopars}. Spatial TB accepts explicit irregular two- or three-dimensional
 Euclidean coordinates for univariate fields. Spatio-temporal CE requires fixed
 spatial locations, an equally spaced time vector and a supported separable
 correlation model.} 

 \item{L}{Numeric; the number of lines in the turning band method.} 
 \item{model}{String; the type of RF and therefore the densities associated to the likelihood
 objects. \code{Gaussian} is the default, see the Section
 Details.}
 \item{n}{Positive integer size parameter. For Binomial it may be scalar or site-specific; for direct Negative Binomial it is the single common number \eqn{r} of successes, with Geometric corresponding to \eqn{r=1}.} 
 \item{param}{A list of parameter values required in the
 simulation procedure of random fields, see \bold{Examples}.}
 \item{anisopars}{ A list of two elements "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
 \item{radius}{Numeric; a value indicating the radius of the sphere when using the great 
 circle distance. Default value is the radius of the earth in Km (i.e. 6371)} 
 \item{X}{Numeric design matrix for the mean.  For fixed locations, rows are ordered by time blocks, with all sites at the first time followed by all sites at the second time.  For dynamic locations, supply either a stacked matrix in temporal-block order or a list with \code{X[[t]]} aligned with \code{coordx_dyn[[t]]}.}
 \item{spobj}{An object of class \code{sp} or \code{spacetime}. For space-time objects, the current \code{sp2Geo()} conversion uses sequential temporal indices and does not preserve irregular original time spacing; to retain irregular temporal distances, use the explicit-coordinate interface with numeric \code{coordt}.}
 \item{nrep}{Numeric; Numbers of indipendent replicates.}
 \item{progress}{Logic; If TRUE then a progress bar is shown.}
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
}


\details{

For \code{method = "CE"}, spatial simulations require \code{grid = TRUE},
separate equally spaced axes in \code{coordx} and \code{coordy},
\code{distance = "Eucl"}, and \code{anisopars = NULL}. The simulated grid is
exactly the grid supplied by the user. The \code{sill} and \code{nugget}
parameters are applied once to the latent Gaussian field.

Spatio-temporal CE is a separable hybrid method: Cholesky decomposition is
used for the spatial correlation matrix and circulant embedding/FFT for the
temporal correlation. It requires fixed spatial locations and a finite,
strictly increasing, equally spaced \code{coordt}. Dynamic coordinates supplied
through \code{coordx_dyn} are not supported by the approximate methods; use
\code{GeoSim(..., method = "cholesky")} instead. A single temporal instant is
handled as a spatial draw at the supplied locations. Bivariate approximate
simulation is not available with \code{method = "CE"}; use \code{method = "TB"}
or \code{GeoSim()}.  Turning bands are implemented only for purely spatial
models; \code{method = "TB"} is rejected for all spatio-temporal correlation
models and requires \code{distance="Eucl"}. Purely spatial univariate TB is restricted to the Matern, generalized-Wendland /
hypergeometric, and Kummer families implemented by the spectral sampler; unsupported
correlation families are rejected before entering the TB kernel. \code{L} must be a
positive integer.

Bivariate approximate simulation is restricted to \code{model="Gaussian"},
\code{corrmodel="Bi_matern"}, common spatial support for the two variables, and
\code{method="TB"}. The two marginal sills and nugget effects are applied to the
standardized bivariate TB output before the component means are added. Non-Gaussian
bivariate approximate simulation is rejected explicitly.

\code{Binary} and \code{Bernoulli} are aliases of a binomial field with
\code{n = 1}; \code{Geom} and \code{Geometric} are aliases of a
negative-binomial field with \code{n = 1}. Direct simulation of
\code{model = "Beta2"} is not implemented by \code{GeoSimapprox}; use
\code{\link{GeoSimCopula}} with an explicit copula instead. The function
stops explicitly rather than returning a latent Gaussian draw. Inferential
\code{Gaussian_misp_*} model names and any other model without an implemented direct
simulator are also rejected explicitly.


For the Tukey transformed-Gaussian models, \code{Tukeyh} requires
\code{0 <= tail < 0.5}, while \code{Tukeyh2} requires both
\code{0 <= tail1 < 0.5} and \code{0 <= tail2 < 0.5}.  In \code{Tukeyh2},
\code{tail1} is the right-tail parameter and \code{tail2} is the left-tail
parameter. Zero is an allowed boundary and recovers the Gaussian transformation
on the corresponding side. \code{SinhAsinh} requires a strictly positive
\code{tail}.

The integer latent-field restrictions of \code{GeoSim} apply unchanged:
\code{Gamma} requires integer \code{shape}; \code{Beta} requires integer
\code{shape1} and \code{shape2}; \code{1/df} must be an integer at least 3 for
the direct Student-t constructions; and \code{2*shape} must be integer for the direct
Poisson-Gamma constructions. No such parameter is rounded silently. Binomial and
negative-binomial \code{n} values must be positive integers. The negative-binomial
simulator uses per-location counters rather than storing the complete Bernoulli history,
and count-process simulators accumulate counts without a growing event-indicator matrix.
\code{nrep} must be a positive integer. The returned \code{param} component preserves
the parameter list supplied by the user.

For univariate simulation, the location parameter is
\eqn{\mu=X\beta}, with coefficients named \code{mean}, \code{mean1}, and so
on in the order of the columns of \code{X}.  If \code{X=NULL}, the simulation
is intercept-only and uses scalar \code{mean}.  Alternatively,
\code{param$mean} may be a vector with one value per simulated observation;
this external mean is mutually exclusive with \code{X}.
}

\section{Spatio-temporal ordering}{
With fixed locations, a simulated space-time realization is returned as a
\eqn{T \times N}{T x N} matrix: rows correspond to \code{coordt} and columns
to rows of \code{coordx}.  The corresponding internal and \code{X} row order
is \code{c(t(data))}, i.e. time then site.

With dynamic locations, the simulated realization is a list of length
\eqn{T}.  Element \code{data[[t]]} has one value per row of
\code{coordx_dyn[[t]]}, in the same row order.  Replicates, when requested,
contain objects with this same fixed or dynamic layout.  See
\code{\link{GeoModels-spacetime-ordering}}.
}

\value{
 Returns an object of class \code{GeoSim}.
 An object of class \code{GeoSim} is a list containing
 at most the following components:
 
 \item{bivariate}{Logical:\code{TRUE} if the Gaussian RF is bivariate, otherwise \code{FALSE};}
 \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
 \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
 \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
 \item{coordx_dyn}{A list of dynamical (in time) spatial coordinates;}
 \item{corrmodel}{The correlation model; see \code{\link{GeoCovmatrix}}.}
 \item{data}{The simulated data.  For fixed-location space-time simulation this is a matrix with times in rows and sites in columns; for dynamic locations it is a list with one vector per time, aligned with \code{coordx_dyn}.  See \code{\link{GeoModels-spacetime-ordering}}.}
 \item{distance}{The type of spatial distance;}
 \item{method}{The method of simulation}
 \item{model}{The type of RF, see \code{\link{GeoFit}}.}
 \item{n}{The Binomial number of trials; for direct Negative Binomial, the common number \eqn{r} of successes.} 
 \item{numcoord}{The number of spatial coordinates;}
 \item{numtime}{The number the temporal realisations of the RF;}
 \item{param}{The parameter list supplied to the simulation call;}
 \item{radius}{The radius of the sphere if coordinates are passed in lon/lat format;}
 \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if
 spatial RF;}
 \item{nrep}{The number of indipendent replicates;}
}

\references{

T. Gneiting, H. Sevcikova, D. B. Percival, M. Schlather and Y. Jiang (2006)
Fast and Exact Simulation of Large Gaussian Lattice Systems in R2: Exploring the Limits
\emph{Journal of Computational and Graphical Statistics} 15 (3)

M. Bevilacqua, X. Emery, F. Cuevas Pacheco (2025) Fast simulation of Gaussian random fields with
flexible correlation models in Euclidean spaces \emph{arxiv} 


}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}



\section{Three-dimensional coordinates}{
For a purely spatial univariate field, \code{method = "TB"} supports explicit
irregular \eqn{N \times 3} Euclidean coordinates with \code{grid = FALSE} and
\code{anisopars = NULL}. Three-dimensional TB is not implemented for
spatio-temporal or bivariate models. CE remains restricted to regular
two-dimensional grids. For unsupported three-dimensional cases, use
\code{GeoSim(..., method = "cholesky")} instead.
}
\examples{
library(GeoModels)


################################################################
###
### Example 1. Simulation of a large spatial Gaussian RF 
### with Matern covariance model
### using circulant embeeding method
### It works only for regular grid
###############################################################
set.seed(68)
x = seq(0,1,0.005)
y = seq(0,1,0.005)
param=list(smooth=1.5,mean=0,sill=1,scale=0.2/3,nugget=0)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSimapprox(coordx=x,coordy=y, grid=TRUE,corrmodel="Matern", model="Gaussian",
 method="CE",param=param)$data
if (requireNamespace("fields", quietly = TRUE)) {
  fields::image.plot(matrix(data1, length(x), length(y), byrow = TRUE))
}

################################################################
###
### Example 2. Simulation of a large spatial Tukey-h RF 
### with Matern covariance model
### using spectral Turning band method
### It works for (ir)regular grid
###############################################################
set.seed(68)
x = runif(50000)
y = runif(50000)
coords=cbind(x,y)
param=list(smooth=0.5,mean=0,sill=1,scale=0.06,nugget=0,tail=0.15)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSimapprox(coords, corrmodel="Matern", model="Tukeyh",
                      method="TB",L=1000,param=param)$data
if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(coords,data1)


################################################################
###
### Example 3. Simulation of a large spacetime Gaussian RF 
### with separable matern covariance model
### using Circular embeeding method
### It works for (large) regular time grid
###############################################################
set.seed(68)
coordt <- (0:100)
coords <- cbind( runif(100, 0 ,1), runif(100, 0 ,1))
param <- list(mean = 0, sill = 1, nugget = 0.25,
 scale_s = 0.05, scale_t = 2, 
 smooth_s = 0.5, smooth_t = 0.5)
# Simulation of a spatial Gaussian RF with Matern correlation function
param<-list(nugget=0,mean=0,scale_s=0.2/3,scale_t=2/3,sill=1,smooth_s=0.5,smooth_t=0.5)

data <- GeoSimapprox(coordx=coords, coordt=coordt, corrmodel="Matern_Matern",
 model="Gaussian",method="CE",param=param)$data
dim(data)

################################################################
###
### Example 4. Simulation of a large spacetime Gaussian RF 
### with separable GenWend covariance model
### using Circular embeeding method in time
###############################################################
set.seed(68)
# Simulation of a spatial Gaussian RF with Matern correlation function
param<-list(nugget=0,mean=0,scale_s=0.2,scale_t=3,sill=1,
 smooth_s=0,smooth_t=0, power2_s=4,power2_t=4)

data <- GeoSimapprox(coordx=coords, coordt=coordt, corrmodel="GenWend_GenWend",
 model="Gaussian",method="CE",param=param)$data
dim(data)


################################################################
###
### Example 6. Simulation of a large bivariate Gaussian RF
### with bivariate Matern correlation model 
### using spectral Turning band method
###############################################################

# Define the spatial-coordinates of the points:
#x <- runif(20000, 0, 2)
#y <- runif(20000, 0, 2)
#coords <- cbind(x,y)

# Simulation of a bivariate spatial Gaussian RF:
# with a Bivariate Matern
#set.seed(12)
#param=list(mean_1=4,mean_2=2,smooth_1=0.5,smooth_2=0.5,smooth_12=0.5,
# scale_1=0.12,scale_2=0.1,scale_12=0.15,
# sill_1=1,sill_2=1,nugget_1=0,nugget_2=0,pcol=0.5)
#data <- GeoSimapprox(coordx=coords,corrmodel="Bi_matern",
# param=param,method="TB",L=1000)$data
#opar=par(no.readonly = TRUE)
#par(mfrow=c(1,2))
#fields::quilt.plot(coords,data[1,],col=terrain.colors(100),main="1",xlab="",ylab="")
#fields::quilt.plot(coords,data[2,],col=terrain.colors(100),main="2",xlab="",ylab="")

}
\keyword{Simulation}
