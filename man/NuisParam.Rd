\name{NuisParam}
\alias{NuisParam}
\encoding{UTF-8}
\title{Lists the Nuisance Parameters of a random field}
\description{
Internal function handling Nuisance Parameters of a random field.
}
\usage{
NuisParam(model, bivariate=FALSE,num_betas=c(1,1),copula=NULL)
}
\arguments{
 \item{model}{String; the name of a random field.}
 \item{bivariate}{Logical; if \code{FALSE} (the default) the correlation model is univariate spatial or spatio-temporal. 
 Otherwise is bivariate.}
 \item{num_betas}{Numerical; the nunber of mean parameters in the linear specification (default is 1) }
 \item{copula}{The type of copula.} 
}

\value{Return a vector string of nuisance parameters.}

\details{The function returns a list with the nuisance parameters of a given
 random field model.}
\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Composite}
