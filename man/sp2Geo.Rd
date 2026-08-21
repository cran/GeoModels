\name{sp2Geo}
\alias{sp2Geo}
\encoding{UTF-8}
\title{Extracting information from an sp or spacetime object}
\description{Extracting information from an sp or spacetime object}
\usage{sp2Geo(spobj,spdata = NULL)}
\arguments{
 \item{spobj}{An object of class sp or spacetime}
 \item{spdata}{Character: The name of data in the sp or spacetime object}

}
\value{A list with spatio-temporal information}


\details{
The function accepts an \pkg{sp} or \pkg{spacetime} object and returns information used by GeoModels functions. Spatial objects should have a coordinate reference system (CRS) defined; if the CRS is missing, \code{sp2Geo} issues a warning because projected/geographic status cannot be determined reliably. For \code{STFDF} and \code{STIDF} objects, the current conversion represents the temporal positions by the sequential indices \code{1, 2, ..., T}; the original spacing of an irregular time index is therefore not preserved. To retain irregular temporal distances, use the explicit-coordinate interfaces of the relevant GeoModels functions and supply the numeric \code{coordt} values directly rather than passing \code{spobj}.
}
\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}
\examples{
if (requireNamespace("sp", quietly = TRUE)) {
# Define the spatial-coordinates of the points:
set.seed(3)
N <- 30 # number of location sites
x <- runif(N, 0, 1)
set.seed(6)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

# Define spatial matrix covariates and regression parameters
X <- cbind(rep(1,N),runif(N))
# Define spatial matrix dependent variable
Y <- rnorm(nrow(X))

crs <- sp::CRS("EPSG:4326")
obj1 <- sp::SpatialPoints(coords, proj4string = crs)
obj2 <- sp::SpatialPointsDataFrame(
  coords, data = data.frame(X, Y), proj4string = crs
)

# sp2Geo info extraction
b <- sp2Geo(obj2,spdata = "Y")
class(b)
b
}
}

\keyword{Composite}
