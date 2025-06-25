\name{corrsas}
\alias{corrsas}
\encoding{UTF-8}
\title{Correlation Function for Sinh-Arcsinh Random Fields}
\description{
  Computes the correlations f for a random field transformed via the sinh-arcsinh (SAS) distribution. This transformation introduces flexible skewness and tail behavior to an underlying Gaussian field. The resulting correlation is derived via an infinite Hermite expansion, as described in Equation (16) of Blasi et al. (2022).
}
\usage{
corrsas(corr, skew, tail, max_coeff = NULL)
}
\arguments{
  \item{corr}{A numeric vector of correlation values of the underlying standard Gaussian random field.}
  \item{skew}{A numeric value representing the skewness parameter \eqn{\alpha} of the sinh-arcsinh transformation. Positive values induce right-skewness, negative values left-skewness.}
  \item{tail}{A positive numeric value representing the tailweight parameter \eqn{\kappa}. Values less than 1 yield heavier tails than Gaussian, while values greater than 1 produce lighter tails.}
  \item{max_coeff}{Optional integer. The maximum number of Hermite coefficients used in the infinite series expansion. If \code{NULL}, a default truncation value is used internally.}
}
\details{
  The correlation of the sinh-arcsinh transformed field is computed as:
  \deqn{
    \rho_{SAS}(h) = \sum_{j=1}^{\infty} \frac{\xi_j^2(\alpha, \kappa)}{j!} \rho(h)^j
  }
  where \eqn{\rho(h)} is the correlation function of the underlying Gaussian field and \eqn{\xi_j(\alpha, \kappa)} are Hermite coefficients depending on the skewness and tail parameters. This series is truncated at \code{max_coeff} terms for computational feasibility.

  See Equation (16) in Blasi et al. (2022) for the full derivation.
}
\value{
  A numeric vector of adjusted correlation values corresponding to the SAS-transformed process.
}
\references{
  Blasi, F., Caamaño-Carrillo, C., Bevilacqua, M., Furrer, R. (2022). A selective view of climatological data and likelihood estimation. \emph{Spatial Statistics}, 50, 100596. \doi{10.1016/j.spasta.2022.100596}
}

\examples{
# Example usage:
rho <- seq(0, 1, length.out = 50)
rho_sas <- corrsas(rho, skew = 0.5, tail = 0.8, max_coeff = 20)
plot(rho, rho_sas, type = "l", main = "SAS Correlation", 
         xlab = "Original Correlation", ylab = "Transformed Correlation")
}
