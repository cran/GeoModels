\name{GeoOMP}
\alias{GeoOMP}
\title{Configure Native OpenMP Acceleration}
\usage{
GeoOMP(threads = NULL, min_pairs = NULL)
}
\arguments{
  \item{threads}{Optional non-negative integer. One selects serial execution,
  values greater than one request that many OpenMP threads, and zero requests
  the OpenMP runtime default.}
  \item{min_pairs}{Optional positive integer giving the minimum number of
  retained pairs required before the OpenMP path is used.}
}
\description{
Configures or inspects optional native OpenMP acceleration. The parallel path
covers univariate Gaussian spatial pairwise marginal and pairwise conditional
likelihoods with the Matern correlation model for any positive smoothness.
Exact half-integer smoothness uses elementary formulas; general smoothness uses
an internal allocation-free Bessel-K kernel that does not call the R API from
OpenMP workers. Other combinations retain the serial implementation.
}
\details{
The default is one thread. This prevents accidental oversubscription in
workflows that already use process-level or future-based parallelism. When
using multiple outer workers, keep \code{GeoOMP(threads = 1L)}. For a single
large \code{GeoFit} call, native OpenMP can instead be enabled explicitly.
}
\value{
Returns a list containing \code{compiled}, \code{threads},
\code{max_threads}, and \code{min_pairs}.
}
\examples{
GeoOMP()
GeoOMP(threads = 1L)
\dontrun{GeoOMP(threads = 8L)}
}
