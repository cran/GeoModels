\name{MargParam}
\alias{MargParam}
\encoding{UTF-8}
\title{Lists the Nuisance Parameters of a Random Field}
\description{
The procedure returns a list with the nuisance parameters of a given
  random field model.
}
\usage{
MargParam(model, bivariate=FALSE,num_betas=c(1,1),copula=NULL)
}
\arguments{
  \item{model}{String; the name of a random field.}
  \item{bivariate}{Logical; if \code{FALSE} (the default) the correlation  model is univariate spatial or spatial-temporal.  
       Otherwise is bivariate.}
  \item{num_betas}{Numerical; the nunber of mean parameters in the linear specification (default is 1)     }
  \item{copula}{The type of copula.} 
}

\value{Return a vector string of nuisance parameters.}

\details{The function returns a list with the nuisance parameters of a given
  random field model.}
\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{
library(GeoModels)

MargParam("Gaussian")

MargParam("Binomial")

MargParam("Weibull",num_betas=2)

MargParam("SkewGaussian", num_betas=3)

MargParam("SinhAsinh")

MargParam("Beta2",copula="Clayton")

MargParam("StudentT")
## note that in the bivariate case sill_1 e sill_2 are considered as correlation parameteres
MargParam("Gaussian", bivariate=TRUE)

}
\keyword{Composite}
