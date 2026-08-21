\name{GeoNA}
\alias{GeoNA}
\encoding{UTF-8}
\title{Deleting NA values (missing values) from a spatial or spatio-temporal dataset.}
\description{The function deletes NA values from a spatial or spatio-temporal dataset }
\usage{GeoNA(data, coordx, coordy=NULL,coordz=NULL, coordt=NULL,
coordx_dyn=NULL, grid=FALSE, X=NULL, setting="spatial")}
\arguments{
 \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid) or a
 (\eqn{t \times d}{t x d})-matrix (a single spatio-temporal realisation) or an (\eqn{d \times d \times t \times n }{d x d x t})-array
 (a single spatio-temporal realisation on regular grid) giving the data.}
 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix 
 Coordinates on a sphere for a fixed radius \code{radius} 
 are passed in lon/lat format expressed in decimal degrees.}
 \item{coordy}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordz}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordt}{A numeric vector giving one dimension of temporal coordinates; the default is \code{NULL}, in which case a spatial random field is expected. Temporal coordinates may be irregularly spaced and their supplied numeric values are preserved.}
 \item{coordx_dyn}{For dynamic locations, a list with one two- or three-column coordinate matrix per time.  It must be aligned with the list \code{data}; missing responses remove the corresponding rows from both objects.}
 \item{grid}{Logical; if \code{FALSE} (the default) the data 
 are interpreted as spatial or spatio-temporal realisations on a set
 of non-equispaced spatial sites (irregular grid).}
 \item{X}{Covariates aligned with \code{data}.  For fixed-coordinate space-time data, use an \eqn{TN \times p}{TN x p} matrix ordered as \code{c(t(data))}.  For dynamic data, \code{X} must be a list with \code{X[[t]]} having one row per element of \code{data[[t]]} and per row of \code{coordx_dyn[[t]]}.}
 \item{setting}{String; are data spatial, spatio-temporal or spatial bivariate
 (respectively \code{spatial}, \code{spacetime}, \code{bivariate}) }

}

\section{Spatio-temporal ordering}{
For fixed sites, \code{data} has times in rows and sites in columns; rows of
\code{X} follow time-major order.  Removing a site removes the corresponding
column of \code{data}, coordinate row, and its row from every temporal block
of \code{X}.  For dynamic sites, \code{data}, \code{coordx_dyn}, and
\code{X} are aligned lists and missing values are removed independently
within each temporal block.  See \code{\link{GeoModels-spacetime-ordering}}.
}

\value{
 Returns a list containing the following components:
 \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
 \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
 \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
 \item{data}{The data without NAvalues }
 \item{grid}{\code{TRUE} if the spatial data are observed in a regular grid, otherwise \code{FALSE};}
 \item{perc}{The percentage of NA values.}
 \item{setting}{Are data of spatial or spatio-temporal or spatial bivariate type}
 \item{X}{The covariate matrix after deleting the rows corresponding to removed observations. For fixed-coordinate spatio-temporal data, deleting a site removes its row from every temporal block.} 
}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

library(GeoModels)

# Define the spatial-coordinates of the points:
set.seed(79)
x = runif(200, 0, 1)
y = runif(200, 0, 1)
coords=cbind(x,y)
# Set the exponential cov parameters:
corrmodel = "Matern"
mean=0
sill=1
nugget=0
scale=0.3/3
smooth=0.5
param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,
 param=param)$data


data[1:100]=NA
# removing NA
a=GeoNA(data,coordx=coords)
a$perc # percentage of NA values 
#a$coordx# spatial coordinates without missing values
#a$data # data without missinng values


}
\keyword{Composite}
