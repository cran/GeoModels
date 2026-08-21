\name{CkType}
\alias{CkType}
\encoding{UTF-8}
\title{Checking Likelihood Objects}
\description{
Subroutine called by \code{InitParam}. \\
The procedure controls the type of likelihood objects inserted by the users.
}
\usage{
CkType(type)
}
\arguments{
 \item{type}{String; the type of likelihood object. Supported values are
 \code{Pairwise}, \code{Standard}, \code{Tapering}, \code{Tapering1},
 \code{Tapering2}, \code{GeoWLS}, \code{CV}, and \code{Independence}.}
}
\value{
The function returns a numeric positive integer, or \code{NULL} if the type of likelihood is invalid.
}
\details{
The procedure checks the likelihood object.
}
\seealso{\code{\link{GeoFit}}}
\author{
Moreno Bevilacqua \email{moreno.bevilacqua89@gmail.com} \url{https://sites.google.com/view/moreno-bevilacqua/home} \\
Víctor Morales Oñate \email{victor.morales@uv.cl} \url{https://sites.google.com/site/moralesonatevictor/} \\
Christian Caamaño-Carrillo \email{chcaaman@ubiobio.cl} \url{https://www.researchgate.net/profile/Christian-Caamano}
}
\keyword{Composite}
