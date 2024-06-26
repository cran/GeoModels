\name{GeoDoScores}
\alias{GeoDoScores}
\encoding{UTF-8}
\title{
Computation of drop-one predictive scores}
\description{
  The function computes RMSE, MAE, LSCORE, CRPS predictive scores based on drop-one 
prediction for a spatial, spatiotemporal and bivariate Gaussian RFs}
\usage{
GeoDoScores(data,  method="cholesky", matrix)
}

\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) 
   or a a(\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation).
    or a a(\eqn{2 \times d}{2 x d})-matrix (a single bivariate realisation).}
  \item{method}{String; the type of matrix decomposition used in the computation of the predictive scores. Default is \code{cholesky}. The other possible choices is \code{svd}.} 
   \item{matrix}{An object of class GeoCovmatrix. See the Section \bold{Details}.}
}

\details{
  For a given covariance matrix object (\code{\link{GeoCovmatrix}})
  and a given spatial, spatiotemporal or bivariare realization
  from a Gaussian random field,
  the function computes four predictive scores based on drop-one prediction.}

\value{
  Returns a list  containing the following informations:
  \item{RMSE}{Root-mean-square error predictive score}
   \item{MAE}{Mean absolute error predictive score}
  \item{LSCORE}{Logarithmic predictive score}
  \item{CRPS}{Continuous ranked probability predictive score}
}


\references{
 Zhang H. and Wang Y. (2010). \emph{Kriging
   and cross-validation for massive spatial data}.
  Environmetrics, \bold{21}, 290--304.
Gneiting T. and Raftery A. \emph{Strictly Proper Scoring Rules, 
   Prediction, and Estimation}.
  Journal of the American Statistical Association, \bold{102}
}

\seealso{\code{\link{GeoCovmatrix}}}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

library(GeoModels)

################################################################
######### Examples of predictive score computation  ############
################################################################
set.seed(8)
  # Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)
coords=cbind(x,y)
matrix1 <- GeoCovmatrix(coordx=coords, corrmodel="Matern", param=list(smooth=0.5,
                      sill=1,scale=0.2,nugget=0))
 
data <- GeoSim(coordx=coords, corrmodel="Matern", param=list(mean=0,smooth=0.5,
                      sill=1,scale=0.2,nugget=0))$data

Pr_scores <- GeoDoScores(data,matrix=matrix1)

Pr_scores

}
\keyword{Predictive scores}
