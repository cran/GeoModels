\name{GeoDoScores}
\alias{GeoDoScores}
\encoding{UTF-8}
\title{
Computation of drop-one predictive scores}
\description{
 The function computes RMSE, MAE, MAD, logarithmic score and CRPS from exact drop-one linear-prediction identities based on a \code{GeoCovmatrix} object.}
\usage{
GeoDoScores(data, method="cholesky", matrix)
}

\arguments{
 \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) 
 or a a(\eqn{t \times d}{t x d})-matrix (a single spatio-temporal realisation).
 or a a(\eqn{2 \times d}{2 x d})-matrix (a single bivariate realisation).}
 \item{method}{String; matrix decomposition used for dense covariance matrices. Possible values are \code{"cholesky"} (default) and \code{"svd"}. Sparse covariance matrices use the sparse Cholesky path and therefore require \code{"cholesky"}.} 
 \item{matrix}{An object returned by \code{\link{GeoCovmatrix}}.}
}

\details{
 Let \eqn{Q=\Sigma^{-1}} and let \eqn{r} denote the data after subtraction of the model mean. The drop-one residual and conditional variance are computed without repeatedly refitting the model:
 \deqn{e_i^{(-i)} = (Qr)_i/Q_{ii}, \qquad v_i^{(-i)}=1/Q_{ii}.}
 The standardized residual is \eqn{z_i=(Qr)_i/\sqrt{Q_{ii}}}. The logarithmic score uses \eqn{\frac{1}{2}\{\log(2\pi v_i^{(-i)})+z_i^2\}}, and the Gaussian CRPS uses its standard closed-form expression with the normal density \eqn{\phi} and distribution function \eqn{\Phi}. The SVD path solves the linear system through the SVD itself rather than passing an SVD object to triangular solvers.}

\value{
 Returns a list containing the following information:
 \item{RMSE}{Root-mean-square error predictive score}
 \item{MAE}{Mean absolute drop-one prediction error.}
 \item{MAD}{Median absolute drop-one prediction error.}
 \item{LSCORE}{Mean Gaussian negative log predictive density for the drop-one predictions.}
 \item{CRPS}{Mean Gaussian continuous ranked probability score for the drop-one predictions.}
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
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

library(GeoModels)

################################################################
######### Examples of predictive score computation ############
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
