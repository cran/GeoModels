\name{GeoSim}
\alias{GeoSim}
\alias{print.GeoSim}
\encoding{UTF-8}
\title{Simulation of Gaussian and non-Gaussian random fields.}
\description{
  Simulates a realization of a Gaussian or non-Gaussian spatial, spatio-temporal,
  or spherical random field, and Gaussian bivariate random fields, for a specified
  covariance or correlation model.
  The covariance parameters are supplied through \code{param}; available correlation
  models are documented in \code{\link{GeoCovmatrix}}.
}
\usage{
GeoSim(coordx=NULL, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, 
 distance="Eucl", grid=FALSE, method="cholesky", 
 model='Gaussian', n=1, param,anisopars=NULL,radius=1, 
 sparse=FALSE,X=NULL,spobj=NULL,nrep=1,progress=TRUE,check.duplicates=FALSE)
}
\arguments{
 \item{coordx}{May be \code{NULL} when \code{coordx_dyn} or \code{spobj} is supplied. Otherwise, a numeric (\eqn{d \times 2}{d x 2})-matrix or 
 (\eqn{d \times 3}{d x 3})-matrix. Coordinates on a sphere for a fixed radius \code{radius} 
 are passed in lon/lat format expressed in decimal degrees.}
 \item{coordy}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordz}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordt}{A numeric vector giving the temporal coordinates at which the field is simulated. Optional argument; the default is \code{NULL}, in which case a spatial random field is expected. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}
 \item{coordx_dyn}{For dynamic simulation sites, a list with one two- or three-column coordinate matrix per element of \code{coordt}.  The output at time \eqn{t} follows the row order of \code{coordx_dyn[[t]]}.  See \code{\link{GeoModels-spacetime-ordering}}.}
 \item{corrmodel}{String; the name of a correlation model, for the
 see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}
 \item{distance}{String; the name of the spatial distance. The default
 is \code{Eucl}, the Euclidean distance. See \code{\link{GeoFit}} for details.}
 \item{grid}{Logical; if \code{FALSE} (the default) the data
 are interpreted as spatial or spatio-temporal realisations on a set
 of non-equispaced spatial sites (irregular grid).}
 \item{method}{String; the type of matrix decomposition used in the simulation. 
 Default is cholesky. The other possible choices is \code{svd}.} 
 \item{model}{String; the type of RF and therefore the densities associated to the likelihood
 objects. \code{Gaussian} is the default, see the Section
 Details.}
 \item{n}{Positive integer size parameter. For Binomial it may be scalar or site-specific; for direct Negative Binomial it is the single common number \eqn{r} of successes, with Geometric corresponding to \eqn{r=1}.} 
 \item{param}{A list of parameter values required in the
 simulation procedure of random fields, see \bold{Examples}.}
 \item{anisopars}{ A list of two elements "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
 \item{radius}{Numeric; a value indicating the radius of the sphere when using the great 
 circle distance. Default value is 1.} 
 \item{sparse}{Logical; if \code{TRUE} then cholesky decomposition is performed
 using sparse matrices algorithms (spam packake). It should be used
 with compactly supported covariance models.\code{FALSE} is the default.}
 \item{X}{Numeric design matrix for the mean.  For fixed locations, rows are ordered by time blocks, with all sites at the first time followed by all sites at the second time.  For dynamic locations, supply either a stacked matrix in temporal-block order or a list with \code{X[[t]]} aligned with \code{coordx_dyn[[t]]}.}
 \item{spobj}{An object of class \code{sp} or \code{spacetime}. For space-time objects, the current \code{sp2Geo()} conversion uses sequential temporal indices and does not preserve irregular original time spacing; to retain irregular temporal distances, use the explicit-coordinate interface with numeric \code{coordt}.}
 \item{nrep}{Numeric; Numbers of indipendent replicates.}
 \item{progress}{Logic; If TRUE then a progress bar is shown.}
 
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
}


\details{

Bivariate simulation is currently implemented only for \code{model = "Gaussian"};
other bivariate marginal models are rejected explicitly.  \code{Binary} and
\code{Bernoulli} are aliases of a binomial field with \code{n = 1}; \code{Geom}
and \code{Geometric} are aliases of a negative-binomial field with
\code{n = 1}. Direct simulation of \code{model = "Beta2"} is not implemented
by \code{GeoSim}; use \code{\link{GeoSimCopula}} with an explicit copula instead.
The function stops explicitly rather than returning a latent Gaussian draw.
Models whose names contain \code{Gaussian_misp_} are inferential working models and
are rejected as data-generating mechanisms. Other model names that are valid elsewhere
in GeoModels but do not have a direct simulator are also rejected explicitly.

Several direct non-Gaussian constructions use an integer number of independent latent
Gaussian fields. Consequently, \code{Gamma} requires integer \code{shape};
\code{Beta} requires integer \code{shape1} and \code{shape2}; and for
\code{StudentT}, \code{SkewStudentT}, and \code{TwoPieceStudentT},
\code{1/df} must be an integer at least 3. For \code{PoissonGamma} and
\code{PoissonGammaZIP}, \code{2*shape} must be a positive integer. These
parameters are never rounded silently. These restrictions apply to the direct latent-field
constructions in \code{GeoSim}; they do not apply to quantile-transformed copula margins
in \code{GeoSimCopula}.

\code{nrep} must be a positive integer. For binomial and negative-binomial fields,
\code{n} must contain positive integers and may be scalar or have one value per
simulated observation. The returned \code{param} component preserves the parameter
list supplied by the user rather than the standardized latent-Gaussian working parameters.

For univariate simulation, the location parameter is
\eqn{\mu=X\beta}, with coefficients named \code{mean}, \code{mean1}, and so
on in the order of the columns of \code{X}.  If \code{X=NULL}, the simulation
is intercept-only and uses scalar \code{mean}.  Alternatively,
\code{param$mean} may be a vector with one value per simulated observation;
this external mean is mutually exclusive with \code{X}.

For the Tukey transformed-Gaussian models, \code{Tukeyh} requires
\code{0 <= tail < 0.5}, while \code{Tukeyh2} requires both
\code{0 <= tail1 < 0.5} and \code{0 <= tail2 < 0.5}.  In \code{Tukeyh2},
\code{tail1} is the right-tail parameter and \code{tail2} is the left-tail
parameter. Zero is an allowed boundary and recovers the Gaussian transformation
on the corresponding side. \code{SinhAsinh} requires a strictly positive
\code{tail}.
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
 \item{coordz}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
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



\seealso{
\code{\link{GeoCovmatrix}} for covariance matrix construction,
  \code{\link{GeoFit}} for parameter estimation,
  \code{\link{GeoSimcond}} and \code{\link{GeoSimapprox}} for conditional and approximate simulation.
}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{
library(GeoModels)


################################################################
###
### Example 1. Simulation of a spatial Gaussian RF 
### with Matern and Generalized Wendland correlations
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500);y <- runif(500)
coords=cbind(x,y)
set.seed(261)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(coordx=coords, corrmodel="Matern", param=list(smooth=0.5,
 mean=0,sill=1,scale=0.4/3,nugget=0))$data

set.seed(261)
data2 <- GeoSim(coordx=coords, corrmodel="GenWend", param=list(smooth=0,
 power2=4,mean=0,sill=1,scale=0.4,nugget=0))$data
opar=par(no.readonly = TRUE)
par(mfrow=c(1,2))
if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(coords, data1, main = "Matern", xlab = "", ylab = "")
}
if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(coords, data2, main = "Wendland", xlab = "", ylab = "")
}
par(opar)
 

################################################################
###
### Example 2. Simulation of a spatial geometric RF 
### with underlying Wend0 correlation
###
################################################################

# Define the spatial-coordinates of the points:
x <- runif(800);y <- runif(800)
coords <- cbind(x,y)
set.seed(251)
# Simulation of a spatial Binomial RF:
sim <- GeoSim(coordx=coords, corrmodel="Wend0",
 model="BinomialNeg",n=1,sparse=TRUE,
 param=list(nugget=0,mean=0,scale=.2,power2=4))

if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(
  coords, sim$data, nlevel = max(sim$data),
  col = terrain.colors(max(sim$data + 1))
 )
}

################################################################
###
### Example 3. Simulation of a spatial Weibull RF
### with underlying Matern correlation on a regular grid
###
###############################################################
# Define the spatial-coordinates of the points:
x <- seq(0,1,0.032)
y <- seq(0,1,0.032)
set.seed(261)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(x,y,grid=TRUE, corrmodel="Matern",model="Weibull", 
 param=list(shape=1.2,mean=0,scale=0.3/3,nugget=0,smooth=0.5))$data
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(x, y, data1, main = "Weibull RF", xlab = "", ylab = "")
}

################################################################
###
### Example 4. Simulation of a spatial t RF
### with with underlying Generalized Wendland correlation 
###
###############################################################
# Define the spatial-coordinates of the points:
x <- seq(0,1,0.03)
y <- seq(0,1,0.03)
set.seed(268)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(x,y,grid=TRUE, corrmodel="GenWend",model="StudentT", sparse=TRUE,
 param=list(df=1/4,mean=0,sill=1,scale=0.3,nugget=0,smooth=1,power2=5))$data
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(
  x, y, data1, col = terrain.colors(100), main = "Student-t RF",
  xlab = "", ylab = ""
 )
}


################################################################
###
### Example 5. Simulation of a sinhasinh RF
### with underlying Wend0 correlation.
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)
coords <- cbind(x,y)
set.seed(261)
corrmodel="Wend0"
# Simulation of a spatial Gaussian RF:
param=list(power2=4,skew=0,tail=1,
 mean=0,sill=1,scale=0.2,nugget=0) ## gaussian case
data0 <- GeoSim(coordx=coords, corrmodel=corrmodel,
 model="SinhAsinh", param=param,sparse=TRUE)$data
plot(density(data0),xlim=c(-7,7))

param=list(power2=4,skew=0,tail=0.7,
 mean=0,sill=1,scale=0.2,nugget=0) ## heavy tails
data1 <- GeoSim(coordx=coords, corrmodel=corrmodel,
 model="SinhAsinh", param=param,sparse=TRUE)$data
lines(density(data1),lty=2)

param=list(power2=4,skew=0.5,tail=1,
 mean=0,sill=1,scale=0.2,nugget=0) ## asymmetry
data2 <- GeoSim(coordx=coords, corrmodel=corrmodel,
 model="SinhAsinh", param=param,sparse=TRUE)$data
lines(density(data2),lty=3)

################################################################
###
### Example 6. Simulation of a bivariate Gaussian RF
### with bivariate Matern correlation model
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)
coords <- cbind(x,y)

# Simulation of a bivariate spatial Gaussian RF:
# with a separable Bivariate Matern
param=list(mean_1=4,mean_2=2,smooth_1=0.5,smooth_2=0.5,smooth_12=0.5,
 scale_1=0.12,scale_2=0.1,scale_12=0.15,
 sill_1=1,sill_2=1,nugget_1=0,nugget_2=0,pcol=0.5)
data <- GeoSim(coordx=coords,corrmodel="Bi_matern",
 param=param)$data
opar=par(no.readonly = TRUE)
par(mfrow=c(1,2))
if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(
  coords, data[1, ], col = terrain.colors(100), main = "1",
  xlab = "", ylab = ""
 )
}
if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(
  coords, data[2, ], col = terrain.colors(100), main = "2",
  xlab = "", ylab = ""
 )
}
par(opar)


################################################################
###
### Example 7. Simulation of a spatio temporal Gaussian random field.
### observed on fixed location sites with double Matern correlation 
###
###############################################################



coordt=1:5

# Define the spatial-coordinates of the points:
x <- runif(50, 0, 2)
y <- runif(50, 0, 2)
coords <- cbind(x,y)

param<-list(nugget=0,mean=0,scale_s=0.2/3,scale_t=2/3,sill=1,smooth_s=0.5,smooth_t=0.5)
data <- GeoSim(coordx=coords, coordt=coordt, corrmodel="Matern_Matern",
 param=param)$data
dim(data)

################################################################
###
### Example 8. Simulation of a spatio temporal Gaussian random field.
### observed on dynamical location sites with double Matern correlation 
###
###############################################################

# Define the dynamical spatial-coordinates of the points:

coordt=1:5
coordx_dyn=list()
maxN=30
set.seed(8)
for(k in 1:length(coordt))
{
NN=sample(1:maxN,size=1)
x <- runif(NN, 0, 1)
y <- runif(NN, 0, 1)
coordx_dyn[[k]]=cbind(x,y)
}
coordx_dyn

param<-list(nugget=0,mean=0,scale_s=0.2/3,scale_t=2/3,sill=1,smooth_s=0.5,smooth_t=0.5)
data <- GeoSim(coordx_dyn=coordx_dyn, coordt=coordt, corrmodel="Matern_Matern",
 param=param)$data
## spatial realization at first temporal instants
data[[1]]
## spatial realization at third temporal instants
data[[3]]




################################################################
###
### Example 9. Simulation of a Gaussian RF 
### with a Wend0 correlation in the north emisphere of the planet earth
### using geodesic distance
###############################################################
distance="Geod";radius=6371

NN=3000 ## total point on the sphere on lon/lat format
set.seed(80)
coords=cbind(runif(NN,-180,180),runif(NN,0,90))
## Set the wendland parameters
corrmodel <- "Wend0"
param<-list(mean=0,sill=1,nugget=0,scale=1000,power2=3)
# Simulation of a spatial Gaussian RF on the sphere
#set.seed(2)
data <- GeoSim(coordx=coords,corrmodel=corrmodel,sparse=TRUE,
 distance=distance,radius=radius,param=param)$data
#require(globe)
#globe::globeearth(eye=place("newyorkcity"))
#globe::globepoints(loc=coords,pch=20,col = cm.colors(length(data),alpha=0.4)[rank(data)])




}
\keyword{Simulation}
