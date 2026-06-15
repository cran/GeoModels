\name{GeoNeighIndex}
\alias{GeoNeighIndex}
\encoding{UTF-8}
\title{Spatial, spatio-temporal, or bivariate nearest-neighbour indices}
\description{
  The function returns nearest-neighbour pair indices for spatial, spatio-temporal,
  or bivariate data. Optionally, a stochastic thinning mechanism can be applied
  to retain only a subset of the candidate nearest-neighbour pairs.
}
\usage{
GeoNeighIndex(coordx, coordy=NULL, coordz=NULL, coordt=NULL,
              coordx_dyn=NULL, distance="Eucl", neighb=4,
              maxdist=NULL, maxtime=1, radius=1,
              bivariate=FALSE, p_neighb=1,
              thin_method="bernoulli")
}
\arguments{
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or
    (\eqn{d \times 3}{d x 3})-matrix. Coordinates on a sphere for a fixed radius
    \code{radius} are passed in longitude/latitude format expressed in decimal degrees.}

  \item{coordy}{A numeric vector giving one dimension of spatial coordinates; optional argument,
    default is \code{NULL}.}

  \item{coordz}{A numeric vector giving one dimension of spatial coordinates; optional argument,
    default is \code{NULL}.}

  \item{coordt}{A numeric vector giving the temporal coordinates. Optional argument,
    default is \code{NULL}; if \code{NULL}, a purely spatial random field is expected.}

  \item{coordx_dyn}{A list of numeric coordinate matrices containing spatial coordinates that may vary
    over time. For spatio-temporal data, the list length must equal the number of time points. For
    bivariate data with different spatial supports, the list must have length two, with one coordinate
    matrix for each variable. Optional argument, default is \code{NULL}.}

  \item{distance}{String; the name of the spatial distance. Default is \code{"Eucl"}
    (Euclidean distance). See the Section \bold{Details} of \code{\link{GeoFit}}.}

  \item{neighb}{Numeric; a positive integer indicating the nearest-neighbour order. In the bivariate
    case, it may also be a vector of length three, corresponding to within-variable 1, cross-variable,
    and within-variable 2 neighbourhood sizes.}

  \item{maxdist}{A numeric value denoting the maximum spatial distance; see \bold{Details}. In the
    bivariate case, it may also be a vector of length three, corresponding to within-variable 1,
    cross-variable, and within-variable 2 distance thresholds.}

  \item{maxtime}{A numeric value denoting the maximum temporal distance; see \bold{Details}.}

  \item{radius}{Numeric; a value indicating the radius of the sphere when using great-circle distances.
    Default value is \code{1}.}

  \item{bivariate}{Logical; if \code{FALSE} (default), the data are interpreted as univariate spatial or
    spatio-temporal realisations. If \code{TRUE}, the data are interpreted as a realization from a
    bivariate field.}

  \item{p_neighb}{Numeric; a value in \eqn{(0,1]} controlling stochastic thinning. Its interpretation
    depends on \code{thin_method}. If \code{thin_method="bernoulli"}, \code{p_neighb} controls the
    expected retained fraction of candidate pairs through calibrated independent Bernoulli inclusion
    probabilities. If \code{thin_method="TargetBalanced"}, \code{p_neighb} is interpreted as a nominal target
    fraction of candidate pairs for the hard-core greedy TargetBalanceding; the final number of retained pairs
    is capped by the endpoint-disjoint TargetBalanceding constraint and may be smaller than the target.}

  \item{thin_method}{String; stochastic thinning scheme. Available options are \code{"bernoulli"}
    and \code{"TargetBalanced"}. The default is \code{"bernoulli"}. With \code{"bernoulli"}, the function uses
    independent Bernoulli thinning, possibly with pair-specific probabilities depending on spatial or
    temporal lags. With \code{"TargetBalanced"}, the function uses hard-core greedy TargetBalanceding and retains only
    endpoint-disjoint pairs.}
}
\details{
  The function first builds a candidate set of directed nearest-neighbour pairs. For purely spatial data,
  the candidate set contains spatial nearest-neighbour pairs. For spatio-temporal data, the function includes
  within-time spatial pairs, pure temporal same-site pairs, and cross-time spatio-temporal pairs up to
  \code{maxtime}. For bivariate data, the function includes within-variable and cross-variable pairs.

  If \code{thin_method="bernoulli"} and \code{p_neighb<1}, candidate pairs are retained independently with
  calibrated Bernoulli probabilities. These probabilities may depend on pair features, such as spatial or
  temporal lag, but they are calibrated so that the expected number of retained pairs is approximately
  \eqn{p_neighb} times the number of candidate pairs.

  If \code{thin_method="TargetBalanced"}, the function applies a hard-core greedy TargetBalanceding procedure. A random
  permutation of the candidate pairs is scanned, and a pair is retained only if neither endpoint has already
  been used by a previously retained pair. Therefore no observation index is used in more than one retained
  pair. In this case \code{p_neighb} is not a marginal inclusion probability. It defines the nominal target
  \eqn{round(p_neighb d)}, where \eqn{d} is the number of candidate pairs, but the final number of retained
  pairs is bounded above by \eqn{\lfloor n/2 \rfloor}, where \eqn{n} is the number of observation indices,
  and may be smaller due to TargetBalanceding feasibility.

  If \code{thin_method="bernoulli"} and \code{p_neighb=1}, no thinning is applied. If
  \code{thin_method="TargetBalanced"} and \code{p_neighb=1}, the function attempts to retain as many pairs as allowed
  by the hard-core TargetBalanceding constraint; this is not equivalent to no thinning.
}
\value{
  Returns a list containing some of the following components:
  \item{colidx}{Vector of neighbour indices.}
  \item{rowidx}{Vector of target indices.}
  \item{lags}{Vector of spatial distances.}
  \item{lagt}{Vector of temporal distances, returned for spatio-temporal data.}
  \item{first}{Variable indicator for the first component of a bivariate pair, returned for bivariate data.}
  \item{second}{Variable indicator for the second component of a bivariate pair, returned for bivariate data.}
  \item{maxdist}{Maximum spatial distance used to construct the candidate pairs, when available.}
  \item{neighb}{Nearest-neighbour order used to construct the candidate pairs, when available.}
  \item{n_candidates}{Number of candidate pairs before thinning.}
  \item{n_retained}{Number of pairs retained after thinning or TargetBalanceding.}
  \item{target_retained}{Target number of retained pairs. For Bernoulli thinning this is the expected
    retained count. For hard-core TargetBalanceding this is the capped target count.}
  \item{target_retained_raw}{Uncapped target number of retained pairs, returned for hard-core TargetBalanceding.}
  \item{TargetBalanceding_cap}{Maximum number of endpoint-disjoint pairs, returned for hard-core TargetBalanceding.}
  \item{effective_fraction}{Observed retained fraction, \code{n_retained/n_candidates}.}
  \item{expected_retained}{Expected retained count under calibrated Bernoulli thinning.}
  \item{thin_method}{Thinning method used.}
  \item{p_neighb_interpretation}{Text description of how \code{p_neighb} was interpreted.}
}
\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com}, \url{https://sites.google.com/view/moreno-bevilacqua/home},
Victor Morales Onate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamano-Carrillo, \email{chcaaman@ubiobio.cl}, \url{https://www.researchgate.net/profile/Christian-Caamano}
}
\examples{
require(GeoModels)
NN <- 400
coords <- cbind(runif(NN), runif(NN))

corrmodel <- "Matern"
scale <- 0.5/3
param <- list(mean=0, sill=1, nugget=0, scale=scale, smooth=0.5)

set.seed(951)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,
               model="Gaussian", param=param)$data

sel <- GeoNeighIndex(coordx=coords, neighb=5)

data1 <- data[sel$colidx]
data2 <- data[sel$rowidx]

## plotting pairs that are neighbours of order 5
plot(data1, data2, xlab="", ylab="",
     main="h-scatterplot, neighb=5")

## Bernoulli thinning: p_neighb controls the expected retained fraction
sel_ber <- GeoNeighIndex(coordx=coords, neighb=5,
                         p_neighb=0.2,
                         thin_method="bernoulli")

data1 <- data[sel_ber$colidx]
data2 <- data[sel_ber$rowidx]

## plotting a random fraction of pairs that are neighbours of order 5 
plot(data1, data2, xlab="", ylab="",
     main="h-scatterplot, neighb=5")
}
\keyword{Scatterplot}
