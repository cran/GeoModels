\name{GeoNeighIndex}
\alias{GeoNeighIndex}
\encoding{UTF-8}
\title{Spatial, spatio-temporal, or bivariate nearest-neighbour indices}
\description{
Builds directed nearest-neighbour candidate pairs for spatial, spatio-temporal,
or bivariate data. Candidate pairs can optionally be thinned using independent
Bernoulli sampling or an exact fixed-budget, stratum-wise design.
}
\usage{
GeoNeighIndex(coordx=NULL, coordy=NULL, coordz=NULL, coordt=NULL,
 coordx_dyn=NULL, distance="Eucl", neighb=4,
 maxdist=NULL, maxtime=1, radius=1,
 bivariate=FALSE, p_neighb=1,
 thin_method="bernoulli", check.duplicates=FALSE)
}
\arguments{
 \item{coordx}{A numeric \eqn{N \times 2}{N x 2} or \eqn{N \times 3}{N x 3}
 coordinate matrix. Longitude/latitude coordinates can be supplied for spherical
 distances. It may be omitted when \code{coordx_dyn} is supplied.}

 \item{coordy}{Optional numeric vector giving a second spatial coordinate when
 coordinates are supplied separately.}

 \item{coordz}{Optional numeric vector giving a third spatial coordinate when
 coordinates are supplied separately.}

 \item{coordt}{A numeric vector giving the temporal coordinates. Optional
 argument, default is \code{NULL}; if \code{NULL}, a purely spatial random
 field is expected. Temporal coordinates may be irregularly spaced; temporal
 lags are computed from the supplied coordinate values.}

 \item{coordx_dyn}{A list of numeric coordinate matrices. For spatio-temporal
 data, the list must contain one coordinate matrix per element of \code{coordt}.
 For bivariate data with different spatial supports, it must contain the two
 variable-specific coordinate matrices. Optional argument, default is
 \code{NULL}.}

 \item{distance}{Character string specifying the spatial distance. Default is
 \code{"Eucl"}. See \code{\link{GeoFit}}.}

 \item{neighb}{A positive integer giving the nearest-neighbour candidate size.
 In the bivariate case, a length-three positive integer vector can be supplied
 for within-variable 1, cross-variable, and within-variable 2 pairs.}

 \item{maxdist}{Optional non-negative maximum spatial distance. In the bivariate
 case, a length-three vector can be supplied for the three pair types.}

 \item{maxtime}{A non-negative numeric value giving the maximum temporal
 distance in the same units as \code{coordt}. Spatio-temporal candidate pairs
 satisfy the corresponding temporal-distance restriction; see Details.}

 \item{radius}{Numeric radius used for spherical distances. Default is
 \code{1}.}

 \item{bivariate}{Logical; if \code{FALSE} (default), construct univariate
 spatial or spatio-temporal pairs. If \code{TRUE}, construct bivariate spatial
 pairs.}

 \item{p_neighb}{Numeric scalar in \eqn{(0,1]}. If \code{1}, the candidate set
 is returned without thinning. For \code{thin_method="bernoulli"}, it is the
 marginal inclusion probability in the constant-probability design and hence
 controls the retained size in expectation. For
 \code{thin_method="FixedBudget"}, it defines the exact global retained-pair
 budget \eqn{K=\mathrm{round}(p_{neighb} d)}, where \eqn{d} is the candidate
 pair count.}

 \item{thin_method}{Character string selecting the thinning design. Accepted
 values are \code{"bernoulli"} and \code{"FixedBudget"}, case-insensitively.
 The legacy name \code{"TargetBalanced"} is accepted as a backward-compatible
 alias for \code{"FixedBudget"}.}
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
}
\details{
The unthinned output is a directed candidate graph: \code{rowidx} identifies
candidate targets and \code{colidx} their neighbours. For a purely spatial
configuration with \eqn{N} locations and candidate size \eqn{m}, the regular
case contains approximately \eqn{Nm} directed candidate pairs.

For chordal and geodesic distances, nearest-neighbour ordering is built in
three-dimensional unit-sphere coordinates. Because chordal distance is monotone
in great-circle distance, this gives the correct global nearest-neighbour order,
including near the dateline and poles. Reported lags are then evaluated in the
requested spherical metric and scaled by \code{radius}.

For spatio-temporal data, temporal coordinates are treated as numeric
coordinates, not as equally spaced indices. Candidate temporal separations are
computed as
\deqn{|t_i-t_j|,}
and \code{maxtime} is interpreted as a temporal-distance threshold in the same
units as \code{coordt}. Thus irregularly spaced temporal coordinates are
supported directly.

With \code{thin_method="bernoulli"} and \code{p_neighb < 1}, candidate pairs
are retained independently. In the current constant-weight implementation the
inclusion probability is \code{p_neighb}, so the retained count is random with
expectation approximately \eqn{p_{neighb}d}.

With \code{thin_method="FixedBudget"} and \code{p_neighb < 1}, exactly
\eqn{K=\mathrm{round}(p_{neighb}d)} candidate pairs are retained, subject only
to truncation to the available candidate count. The budget is allocated across
natural candidate strata using proportional quotas with randomized residual
allocation, followed by sampling without replacement within each stratum. For
purely spatial data, strata are the target-specific nearest-neighbour lists.
For spatio-temporal data, target and temporal-lag descriptors define the
strata; for bivariate data, the variable-pair descriptors are also included.
The resulting global retained size is exact.

For \code{thin_method="FixedBudget"}, the exact target is \eqn{round(p d)}. Hence very small \code{p_neighb} values may legitimately return zero retained pairs. This is a valid graph-selection result; \code{GeoFit} rejects such a graph because a pairwise likelihood cannot be estimated without pair contributions. Exact duplicate observation locations are rejected before nearest-neighbor construction.

When \code{p_neighb = 1}, no thinning is applied, independently of
\code{thin_method}.
}
\value{
A list containing the pair indices and distances required by downstream
GeoModels procedures. Depending on the setting, components include:
\item{rowidx}{Directed target indices.}
\item{colidx}{Directed neighbour indices.}
\item{lags}{Spatial distances for the retained pairs.}
\item{lagt}{Temporal distances for spatio-temporal pairs.}
\item{first}{First variable indicator for bivariate pairs.}
\item{second}{Second variable indicator for bivariate pairs.}
\item{maxdist}{Spatial-distance threshold, when stored for the selected case.}
\item{neighb}{Nearest-neighbour candidate size, when stored for the selected case.}

The function returns only the retained pair structure; it does not append
summary fields such as candidate or retained counts. These can be obtained from
the lengths of the returned index vectors.
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

## Full directed 5-NN candidate graph.
sel <- GeoNeighIndex(coordx=coords, neighb=5)
length(sel$rowidx)

## Bernoulli thinning: retained size is random with expectation 20 percent
## of the candidate graph.
set.seed(1)
sel_ber <- GeoNeighIndex(coordx=coords, neighb=5,
 p_neighb=0.2, thin_method="bernoulli")
length(sel_ber$rowidx)

## Fixed-budget thinning: the global retained size is exact.
set.seed(1)
sel_fix <- GeoNeighIndex(coordx=coords, neighb=5,
 p_neighb=0.2, thin_method="FixedBudget")
length(sel_fix$rowidx)

## Irregular temporal coordinates: maxtime is a distance threshold.
times <- c(0, 0.25, 1.1, 2.8)
sel_st <- GeoNeighIndex(coordx=coords[1:30, ], coordt=times,
 neighb=3, maxtime=1)
head(sel_st$lagt)
}
\keyword{nearest neighbour}
\keyword{spatio-temporal}
