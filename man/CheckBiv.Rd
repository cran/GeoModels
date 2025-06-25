\name{CheckBiv}
\alias{CheckBiv}
\encoding{UTF-8}
\title{Checking Bivariate Covariance Models}
\description{
Checks whether the correlation model is bivariate.
}
\usage{
CheckBiv(numbermodel)
}
\arguments{
  \item{numbermodel}{A numeric value; the number associated with a given correlation model.}
}
\details{ 
This function checks whether the correlation model is bivariate.
}
\value{
A logical value: \code{TRUE} if the correlation model is bivariate, and \code{FALSE} otherwise.
} 
\author{
Moreno Bevilacqua \email{moreno.bevilacqua89@gmail.com} \cr
\url{https://sites.google.com/view/moreno-bevilacqua/home} \cr
Víctor Morales Oñate \email{victor.morales@uv.cl} \cr
\url{https://sites.google.com/site/moralesonatevictor/} \cr
Christian Caamaño-Carrillo \email{chcaaman@ubiobio.cl} \cr
\url{https://www.researchgate.net/profile/Christian-Caamano}
}
\examples{
library(GeoModels)
CheckBiv(CkCorrModel("Bi_matern_sep"))
}
\keyword{composite}
