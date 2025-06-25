\name{CheckDistance}
\alias{CheckDistance}
\encoding{UTF-8}
\title{Checking Distance Type}
\description{
Checks the validity and type of the specified distance.
}
\usage{
CheckDistance(distance)
}
\arguments{
  \item{distance}{A character string indicating the type of distance. 
  Available options are: \code{"Eucl"} (Euclidean), \code{"Geod"} (Geodesic), and \code{"Chor"} (Chordal).
  See also \code{\link{GeoCovmatrix}}.}
}
\details{
This function checks whether the specified distance type is valid.
}
\value{
An integer:
\itemize{
  \item \code{0} for Euclidean distance
  \item \code{1} for Geodesic distance
  \item \code{2} for Chordal distance
}
If the input is not recognized, the function returns \code{NULL}.
}
\author{
Moreno Bevilacqua \email{moreno.bevilacqua89@gmail.com} \cr
\url{https://sites.google.com/view/moreno-bevilacqua/home} \cr
Víctor Morales Oñate \email{victor.morales@uv.cl} \cr
\url{https://sites.google.com/site/moralesonatevictor/} \cr
Christian Caamaño-Carrillo \email{chcaaman@ubiobio.cl} \cr
\url{https://www.researchgate.net/profile/Christian-Caamano}
}
\keyword{composite}
