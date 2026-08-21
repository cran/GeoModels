\name{GeoModels-spacetime-ordering}
\alias{GeoModels-spacetime-ordering}
\alias{spacetime-ordering}
\alias{spacetime-data-order}
\encoding{UTF-8}
\title{Ordering Conventions for Spatio-temporal Data, Covariates and Predictions}
\description{
Defines the ordering used by \pkg{GeoModels} for fixed-location and
dynamic-location spatio-temporal data, observation-level means and
covariates, and kriging prediction tasks.
}
\details{
The conventions below apply consistently to data, coordinates, design
matrices, known mean vectors, and prediction outputs.
}
\section{Irregular temporal coordinates}{
For the explicit-coordinate space-time interfaces, \code{coordt} is a numeric
coordinate vector and does not need to be equally spaced. Covariance,
likelihood, pair-selection, variogram, neighborhood and kriging calculations
use temporal lags derived from the supplied values, typically
\eqn{|t_i-t_j|}. Consequently, \code{maxtime} is a temporal-distance threshold
expressed in the same units as \code{coordt}, not an index lag or an order of
neighboring time points.

Important exceptions are explicit. Space-time
\code{GeoSimapprox(..., method = "CE")} requires finite, strictly increasing,
equally spaced temporal coordinates, and space-time turning-bands simulation
is not implemented. \code{GeoSimcond} and \code{GeoOutlier} currently reject
space-time input.

The \code{spobj} route has a separate legacy limitation: the current
\code{sp2Geo()} conversion of \code{STFDF}/\code{STIDF} objects replaces their
time positions by sequential indices \code{1, 2, ..., T}. Therefore an
irregular original time spacing is not preserved through \code{spobj}. When
irregular temporal distances matter, pass the data and coordinates through the
explicit interface and supply the numeric \code{coordt} values directly.
}

\section{Fixed spatial locations}{
Let \eqn{T} be \code{length(coordt)} and let \eqn{N} be the number of rows of
\code{coordx}. The data must be a \eqn{T \times N}{T x N} matrix. Row
\eqn{t} corresponds to \code{coordt[t]}, and column \eqn{i} corresponds to
row \eqn{i} of \code{coordx}.

The internal observation order is time-major and is equivalent to
\code{c(t(data))}:
\preformatted{
time 1: site 1, site 2, ..., site N
time 2: site 1, site 2, ..., site N
...
time T: site 1, site 2, ..., site N
}
Any observation-level vector, including a known mean vector, must use this
order. A design matrix \code{X} must have \eqn{NT} rows in the same order.
For example, if \code{X_time} has one row per time,
\code{X_time[rep(seq_len(T), each = N), , drop = FALSE]} creates the
corresponding observation-level design.
}

\section{Dynamic spatial locations}{
For dynamic locations, \code{coordx_dyn} must be a list of length \eqn{T}.
The element \code{coordx_dyn[[t]]} is an \eqn{N_t \times 2}{N_t x 2} or
\eqn{N_t \times 3}{N_t x 3} matrix containing the locations observed at
\code{coordt[t]}.

The data must be a list of length \eqn{T}; \code{data[[t]]} must contain
\eqn{N_t} values, and its \eqn{i}-th value corresponds to row \eqn{i} of
\code{coordx_dyn[[t]]}. The internal order is the time-wise concatenation
\code{unlist(data, use.names = FALSE)}, matched by
\code{do.call(rbind, coordx_dyn)}:
\preformatted{
all observations at time 1,
then all observations at time 2,
...,
then all observations at time T.
}
A design matrix \code{X} can be supplied either as a
\eqn{(\sum_t N_t) \times p}{sum(N_t) x p} matrix in this concatenated order,
or as a list of length \eqn{T} with \code{X[[t]]} having \eqn{N_t} rows and
the same row order as \code{coordx_dyn[[t]]}. A known observation-level mean
uses the same time-wise concatenation.
}

\section{Mean parameters}{
For a univariate linear mean,
\deqn{\mu = X\beta.}
The regression vector \eqn{\beta} is represented by the parameters
\code{mean}, \code{mean1}, \code{mean2}, and so on. The columns of \code{X}
correspond to these coefficients in exactly this order. An intercept-only
model is equivalent to a one-column matrix of ones and the single coefficient
\code{mean}. A supplied observation-level mean vector is a known mean and is
not interpreted as regression coefficients.
}

\section{Kriging prediction order}{
Let \eqn{N_{loc}} be the number of rows of \code{loc} and let
\eqn{T_{loc}} be \code{length(time)}. Kriging prediction tasks are ordered
location-major:
\preformatted{
location 1: time 1, time 2, ..., time Tloc
location 2: time 1, time 2, ..., time Tloc
...
location Nloc: time 1, time 2, ..., time Tloc
}
Therefore the rows of \code{Xloc} and the elements of \code{Mloc} must follow
this order. If \code{X_time_pred} has one row per prediction time, use
\code{X_time_pred[rep(seq_len(Tloc), times = Nloc), , drop = FALSE]}.
Equivalently, if \code{mu_time} contains one prediction mean per time, use
\code{rep(mu_time, times = Nloc)} for \code{Mloc}.

Dynamic coordinates describe the observed sites only and do not change the
prediction-task order. Spatio-temporal kriging output is a
\eqn{T_{loc} \times N_{loc}}{Tloc x Nloc} matrix: rows correspond to
\code{time} and columns correspond to rows of \code{loc}.
}
\seealso{
\code{\link{GeoFit}}, \code{\link{GeoSim}}, \code{\link{GeoVariogram}},
\code{\link{GeoKrig}}, and \code{\link{GeoKrigloc}}.
}
\examples{
T <- 3
N <- 2
data <- matrix(seq_len(T * N), nrow = T, ncol = N)
c(t(data))

X_time <- cbind(1, seq_len(T))
X <- X_time[rep(seq_len(T), each = N), , drop = FALSE]

Nloc <- 2
Xloc <- X_time[rep(seq_len(T), times = Nloc), , drop = FALSE]
Mloc <- rep(seq_len(T), times = Nloc)

## Irregular observation times are allowed in the explicit interface.
coordt_irregular <- c(0, 0.25, 1.1, 2.8)
## Temporal lags are based on differences between these values; for example,
abs(outer(coordt_irregular, coordt_irregular, "-"))
}
\keyword{spatial}
