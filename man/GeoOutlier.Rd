\name{GeoOutlier}
\alias{GeoOutlier}
\encoding{UTF-8}
\title{Spatial outlier detection}
\description{
Detects outliers in a univariate purely spatial data set using the selected local procedure. Space-time and bivariate outlier detection are not currently implemented.
}
\usage{
GeoOutlier(data, coordx, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL, 
 distance="Eucl", grid=FALSE, neighb=10,alpha=0.001,
 method="Z-Median", radius=1, bivariate=FALSE,X=NULL)
}
\arguments{
 \item{data}{A numeric vector containing a single spatial realization, or a matrix representing a single spatial realization on a regular grid. Space-time data are not currently supported.}
 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix 
 Coordinates on a sphere for a fixed radius \code{radius} 
 are passed in lon/lat format expressed in decimal degrees.}
 \item{coordy}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordz}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}

 \item{coordt}{Retained for API compatibility. \code{GeoOutlier} is currently implemented only for purely spatial data and rejects space-time input.}

 \item{coordx_dyn}{Retained for API compatibility. Dynamic space-time coordinates are not currently implemented.}

 \item{distance}{String; the name of the spatial distance. The default
 is \code{Eucl}, the Euclidean distance. See \code{\link{GeoFit}} for details.}

 \item{grid}{Logical; if \code{FALSE} (the default), the data are interpreted as a spatial realization on explicit, possibly non-equispaced, spatial sites.}

 \item{neighb}{Numeric; an optional positive integer indicating the 
 order of neighborhoodused for Z-Median algorithm.}

 \item{alpha}{Numeric; a numeric value between 0 and 1 used for Z-Median algorithm.}

 \item{method}{String; The name of the algorithm for detecting spatial ouliers. Default is Z-median proposed in Chen et al. (2008)}
\item{radius}{Numeric; a value indicating the radius of the sphere when using the great 
 circle distance. Default value is 1.} 

 \item{bivariate}{Retained for API compatibility. Bivariate outlier detection is not currently implemented and \code{TRUE} is rejected.}

 \item{X}{Numeric; an optional matrix of spatial covariates with one row per observation.}
}


\value{
 Returns a matrix containing the detected spatial outliers
}

\references{
 Chen D, Lu C, Kou Y, Chen F (2008) On detecting spatial outliers. Geoinformatica 12:455–475
 
 Bevilacqua M., Caamaño C., Arellano-Valle R. B., Camilo Gomez C. (2022) A class of random fields with two-piece marginal distributions for modeling point-referenced data with spatial outliers. \emph{Test} 10.1007/s11749-021-00797-5
}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\section{Three-dimensional coordinates}{
Purely spatial irregular coordinates may be supplied as an
\eqn{N \times 3} matrix or as separate \code{coordx}, \code{coordy}, and
\code{coordz} vectors. Three-dimensional coordinates require
\code{distance = "Eucl"}. Space-time and bivariate outlier detection are
not implemented.
}
