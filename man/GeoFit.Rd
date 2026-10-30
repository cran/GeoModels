\name{GeoFit}
\alias{GeoFit}
\alias{print.GeoFit}
\encoding{UTF-8}
\title{Maximum-Likelihood-Based Fitting of Gaussian and non-Gaussian random fields.}
\description{
 Maximum weighted (stochastic) composite-likelihood fitting for Gaussian and some non-Gaussian univariate spatial, spatio-temporal and bivariate spatial random fields. The function allows fixing any of the parameters and setting upper/lower bounds in the optimization.
 Different optimization methods can be used.
}
\usage{
GeoFit(data, coordx, coordy=NULL, coordz=NULL, coordt=NULL,
 coordx_dyn=NULL, copula=NULL, corrmodel=NULL, distance="Eucl",
 fixed=NULL, anisopars=NULL, est.aniso=c(FALSE,FALSE),
 grid=FALSE, likelihood="Marginal", lower=NULL, maxdist=Inf,
 neighb=NULL, p_neighb=1, maxtime=Inf, memdist=TRUE,
 method="cholesky", model="Gaussian", n=1, onlyvar=FALSE,
 optimizer="Nelder-Mead", radius=1, score=FALSE,
 sensitivity=FALSE, sparse=FALSE, start=NULL,
 thin_method="bernoulli", type="Pairwise", upper=NULL, 
 varest=FALSE, weighted=FALSE,
 X=NULL, spobj=NULL, spdata=NULL, check.duplicates=FALSE)
}

\arguments{
 \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
 or a (\eqn{t \times d}{t x d})-matrix (a single spatio-temporal realisation) or a (\eqn{d \times d \times t \times n}{d x d x t x n})-array
 (a single spatio-temporal realisation on regular grid).
 See Details for the accepted data layouts.}

 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix.
 Coordinates on a sphere for a fixed radius \code{radius} are passed in lon/lat format expressed in decimal degrees.}

 \item{coordy}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, default is \code{NULL}.}

 \item{coordz}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, default is \code{NULL}.}

 \item{coordt}{A numeric vector assigning one dimension of the observation-time coordinates. Optional argument, default is \code{NULL}; if \code{NULL}, a spatial random field is expected. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}

 \item{coordx_dyn}{For dynamic locations, a list with one numeric coordinate matrix per temporal instant.  If \eqn{T} is \code{length(coordt)}, the list must have length \eqn{T}; \code{coordx_dyn[[t]]} must have \eqn{N_t} rows and two or three columns.  Its row order must match \code{data[[t]]} and, when supplied as a list, \code{X[[t]]}.  See \code{\link{GeoModels-spacetime-ordering}}.}

 \item{copula}{String; the copula used by pairwise copula likelihoods. Supported values are \code{"Gaussian"}, \code{"Clayton"} (the constructive Clayton-like spatial copula), and \code{"SkewGaussian"}.}

 \item{corrmodel}{String; the name of a correlation model; see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}

 \item{distance}{String; the name of the spatial distance. Default is \code{"Eucl"} (Euclidean distance). See Details for the accepted options.}

 \item{fixed}{An optional named list giving the values of the parameters that will be considered as known values. Parameter names must match exactly the nuisance/marginal or correlation parameters supported by the selected model and correlation model; invalid names are reported explicitly.
 The listed parameters for a given model/correlation function will not be estimated.}

 \item{anisopars}{A list of two elements: \code{"angle"} and \code{"ratio"}, i.e. the anisotropy angle and the anisotropy ratio, respectively. Geometric anisotropy is available only with \code{distance="Eucl"}.}

 \item{est.aniso}{A bivariate logical vector providing which anisotropy parameters must be estimated.}

 \item{grid}{Logical; if \code{FALSE} (default) the data are interpreted as spatial or spatio-temporal realisations on a set of non-equispaced spatial sites. If \code{TRUE}, \code{coordx}, \code{coordy} (and optionally \code{coordz}) are interpreted as grid axes and the data dimensions must match those axes; grid inputs are internally converted to the canonical explicit-coordinate layout before fitting.}

 \item{likelihood}{String; the configuration of the composite likelihood. \code{"Marginal"} is the default; see Details for the accepted options.}

 \item{lower}{An optional named list giving lower bounds for parameters when the optimizer is \code{L-BFGS-B}, \code{nlminb}, \code{bobyqa} or \code{optimize}.
 Names must match parameters that are being estimated.}

 \item{maxdist}{Numeric; an optional positive value indicating the maximum spatial distance considered in the composite computation.
 See Details for more information.}

 \item{neighb}{Numeric; an optional positive integer indicating the order of neighborhood in the composite likelihood computation.
 See Details for more information.}

 \item{p_neighb}{Numeric scalar in \eqn{(0,1]}. If \code{1} (default), no thinning is applied. For \code{thin_method="bernoulli"}, it is the candidate-pair inclusion probability in the constant-weight design and controls the retained size in expectation. For \code{thin_method="FixedBudget"}, it defines the exact global retained-pair budget \eqn{K=\mathrm{round}(p_{neighb}d)}, where \eqn{d} is the candidate-pair count.}

 \item{maxtime}{Numeric; an optional non-negative maximum temporal-distance threshold, expressed in the same units as \code{coordt}, used to select pairwise composite-likelihood contributions.}

 \item{memdist}{Deprecated logical argument retained for backward compatibility. The selected pair structure is always precomputed and reused during composite-likelihood optimization. Supplying \code{FALSE} produces a warning and is treated as \code{TRUE}.}

 \item{method}{String; the type of matrix decomposition/linear algebra backend used in likelihood computations.
 Default is \code{"cholesky"}. Another possible choice is \code{"svd"} for dense
 likelihood computations. Sparse Gaussian and tapered likelihoods require
 \code{"cholesky"}; requesting \code{"svd"} for these paths raises an error.}

 \item{model}{String; the type of random field (and associated density) used in the likelihood objects. Default is \code{"Gaussian"}; see Details for the accepted options.}

 \item{n}{Positive integer size parameter. For direct Binomial models it may be scalar or contain one value per observation. For direct Negative-Binomial models it is the common number \eqn{r} of successes and must be scalar; Geometric corresponds to \eqn{r=1}. Non-integer or non-positive values are rejected.}

 \item{onlyvar}{Logical; if \code{TRUE} (and \code{varest=TRUE}) only the variance-covariance matrix is computed without optimizing.
 Default is \code{FALSE}.}

 \item{optimizer}{String; the optimization algorithm (see \code{\link{optim}} for details). Default is \code{"Nelder-Mead"}.
 Other possible choices are \code{"nlm"}, \code{"BFGS"}, \code{"SANN"}, \code{"L-BFGS-B"}, \code{"nlminb"}, \code{"bobyqa"}.
 For \code{"L-BFGS-B"}, \code{"nlminb"} and \code{"bobyqa"} bounds can be passed via \code{lower} and \code{upper}.
 In the one-dimensional case, \code{optimize} is used.}

 \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates. Default is \code{1}.}

 \item{score}{Logical; if \code{TRUE} the score function is computed. Default is \code{FALSE}.}

 \item{sensitivity}{Logical; if \code{TRUE} the sensitivity matrix is computed. Default is \code{FALSE}.}

 \item{sparse}{Logical; if \code{TRUE} then maximum likelihood / composite likelihood may exploit sparse-matrix algorithms (e.g., \pkg{spam}).
 Typically used with compactly supported covariance models. Default is \code{FALSE}.}

 \item{start}{An optional named list with initial values for parameters to be estimated. Default is \code{NULL}. Parameter names must match exactly the nuisance/marginal or correlation parameters supported by the selected model and correlation model; invalid names are reported explicitly. Parameters omitted from \code{start} receive conservative internal initial values; explicitly supplied values are preserved (see Details).}

 \item{thin_method}{String; thinning scheme in stochastic weighted pairwise
 likelihood (used when \code{p_neighb < 1}). Default is \code{"bernoulli"}
 (independent Bernoulli thinning). An alternative thinning scheme is
 \code{"FixedBudget"}. The legacy name \code{"TargetBalanced"} is accepted as an alias for \code{"FixedBudget"}.}

 \item{type}{String; the type of likelihood objects. If \code{"Pairwise"} (default) the composite likelihood is formed by pairwise components (see Details).}

 \item{upper}{An optional named list giving upper bounds for parameters when the optimizer is \code{L-BFGS-B}, \code{nlminb}, \code{bobyqa} or \code{optimize}.
 Names must match parameters that are being estimated.}

 \item{varest}{Logical; if \code{TRUE} the estimates variances and standard errors are returned. For composite likelihood estimation it is deprecated.
 Use \code{sensitivity=TRUE} and update the object using \code{GeoVarest}. Default is \code{FALSE}.}

 \item{weighted}{Logical; if \code{TRUE} the likelihood objects are weighted; see Details for the accepted options. Default is \code{FALSE}.}

 \item{X}{Numeric design matrix for the linear mean \eqn{X\beta}.  Its columns correspond in order to \code{mean}, \code{mean1}, \ldots.  For fixed-location space-time data, rows follow \code{c(t(data))}: all sites at the first time, then all sites at the second time, and so on.  For dynamic locations, use either a stacked matrix in time-wise order or a list with \code{X[[t]]} aligned with \code{coordx_dyn[[t]]}.  See \code{\link{GeoModels-spacetime-ordering}}.}

 \item{spobj}{An object of class \code{sp} or \code{spacetime}. For space-time objects, the current \code{sp2Geo()} conversion uses sequential temporal indices and does not preserve irregular original time spacing; to retain irregular temporal distances, use the explicit-coordinate interface with numeric \code{coordt}.}

 \item{spdata}{Character; the name of the data component in the \code{sp} or \code{spacetime} object.}
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
}
\details{
For the fixed-margin models \code{Genwend_Genwend_nosep} (85) and
\code{Matern_Matern_nosep} (89), \code{sep} is the normalized interaction
on \eqn{[-1,1]}; \code{sep=0} is separable and is not a boundary value.
Automatic starts use zero for this parameter and automatic bounds admit
both signs. The GW shape constraints are coupled to the smoothness and
spatial dimension. Invalid trial parameter vectors are rejected by the
native covariance evaluator; user-specified starts must be valid.
See \code{\link{GeoFixedMarginST}} for parameter meanings and examples.
These definitions replace the previous models at codes 85 and 89, so old
fits using those codes must be refitted.



For a univariate model the mean is specified as
\deqn{\mu = X\beta.}
The entries of \eqn{\beta} are represented by the parameters \code{mean},
\code{mean1}, and so on.  The columns of \code{X} correspond to these
coefficients in exactly this order.  When \code{X=NULL}, an intercept-only model is used, which is
equivalent to a one-column matrix of ones and the single coefficient
\code{mean}.  Missing starting values for regression coefficients and other
estimated parameters are initialized internally. Explicit entries in \code{start}
always define starting values for parameters that remain estimated, including
\code{likelihood="Full", type="Standard"}; they are never converted to fixed
mean coefficients. A coefficient is fixed only when it is supplied in \code{fixed}.
Supplied coefficient names must be compatible with \code{ncol(X)}. For difficult
non-Gaussian models, user-supplied starting values can still improve numerical
optimization.

For \code{model="Wrapped"}, omitted automatic starting values use circular
statistics. If \eqn{\bar C=n^{-1}\sum_i\cos(\theta_i)},
\eqn{\bar S=n^{-1}\sum_i\sin(\theta_i)}, and
\eqn{R=(\bar C^2+\bar S^2)^{1/2}}, the intercept is initialized from the
sample mean direction through the model link
\eqn{\mu=2\arctan(\eta)+\pi}, while \code{sill} is initialized as
\eqn{-2\log R}. Additional automatically initialized mean coefficients are
set to zero. Explicit entries in \code{start} or \code{fixed} retain
precedence. The wrapped pairwise and conditional densities use winding
indices from \eqn{-3} through \eqn{3}.

A site-specific known mean can instead be supplied as a vector in
\code{fixed$mean}, with one value per observation.  This external mean is
mutually exclusive with \code{X} and with estimated mean coefficients in
\code{start}.  A scalar \code{fixed$mean} is not an external vector: it is the
fixed intercept coefficient.

 \code{GeoFit} provides weighted and stochastic weighted composite-likelihood
 estimation based on pairs for Gaussian and non-Gaussian random fields,
 including nearest-neighbor and stochastic nearest-neighbor pairwise likelihoods;
 see Caamaño-Carrillo et al. (2024) and Bevilacqua et al. (2026).
 It also provides independence composite-likelihood estimation. The accepted
likelihood/type combinations are checked explicitly: \code{likelihood="Full"}
uses the full-likelihood types (including \code{"Standard"}).  Full-likelihood
objectives are currently implemented for Gaussian, SinhAsinh, LogGaussian,
Tukeyh, Tukeyh2, and the internally supported misspecified-Gaussian full
likelihoods. Other margins, including Gamma and Weibull, are rejected before
optimization with an explicit capability message;
\code{likelihood="Marginal"} is used with \code{"Pairwise"} or
\code{"Independence"}; and \code{likelihood="Conditional"} is used with
\code{"Pairwise"}. The historical Difference composite likelihood is no longer
supported. For space-time pairwise fitting, GeoFit checks that a registered
space-time native kernel exists for the requested marginal model and stops
before optimization when it does not. Bivariate pairwise fitting is currently
implemented only for the Gaussian model with marginal pairwise likelihood.

For pairwise copula fitting, the current implementation is univariate and purely
spatial.  Continuous margins supported with the Gaussian, Clayton-like, and
skew-Gaussian copulas are \code{"Gaussian"}, \code{"StudentT"},
\code{"LogGaussian"}, \code{"Gamma"}, \code{"Weibull"},
\code{"Beta"}, \code{"Beta2"}, \code{"Kumaraswamy"},
\code{"Kumaraswamy2"}, \code{"Logistic"}, and \code{"SkewLaplace"}.
The Gaussian and skew-Gaussian copulas additionally support the discrete margins
\code{"Poisson"}, \code{"Binomial"}, and \code{"BinomialNeg"}. For the
skew-Gaussian copula these probabilities are evaluated as copula-rectangle
probabilities; when \code{nu=0} the calculation reduces exactly to the Gaussian
copula likelihood.
For \code{copula="Clayton"}, \code{nu} is the positive-integer parameter of
the constructive Clayton-like random field and must be supplied in
\code{fixed}; it is not a continuously estimable copula parameter.  For
\code{copula="SkewGaussian"}, \code{nu} is the bounded asymmetry parameter
\eqn{\eta\in(-1,1)} and invalid optimizer proposals are rejected by the
objective function.

 The optimization method is specified using \code{optimizer}. The default method is \code{Nelder-Mead}; other available methods are
 \code{nlm}, \code{BFGS}, \code{SANN}, \code{L-BFGS-B}, \code{bobyqa}, and \code{nlminb}. In the last three cases,
 bounds can be specified using \code{lower} and \code{upper}.

 Depending on the dimension of \code{data} and on the name of the correlation model,
 the observations are assumed to be a realization of a spatial, spatio-temporal or bivariate random field.
 Specifically, with \code{data}, \code{coordx}, \code{coordy}, \code{coordt}:
 \itemize{
 \item If \code{data} is a numeric \eqn{d}{d}-dimensional vector and \code{coordx}, \code{coordy} are two numeric \eqn{d}{d}-dimensional vectors
 (or \code{coordx} is a (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}), then the data are interpreted as a single spatial realisation observed on \eqn{d}{d} spatial sites;
 \item If \code{data} is a numeric (\eqn{t \times d}{t x d})-matrix and \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, then the data are interpreted as a single spatio-temporal realisation observed on \code{d} sites and \code{t} times;
 \item If \code{data} is a numeric (\eqn{2 \times d}{2 x d})-matrix, then the data are interpreted as a single bivariate spatial realisation observed on \code{d} spatial sites;
 \item If \code{data} is a list, \code{coordx_dyn} is a list and \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, then the data are interpreted as a spatio-temporal realisation observed on dynamical spatial sites (different locations for each time) and for \code{t} times.
 }

 It is also possible to specify a matrix of covariates using \code{X}. Specifically:
 \itemize{
 \item In the spatial case, \code{X} must be a (\eqn{d \times k}{d x k}) matrix associated to \code{data} (a \eqn{d}{d}-vector);
 \item In the spatio-temporal case, \code{X} must be a (\eqn{N \times k}{N x k}) matrix associated to \code{data} (a \eqn{t \times d}{t x d}-matrix), where \eqn{N=t\times d};
 \item In the bivariate case, \code{X} must be a (\eqn{N \times k}{N x k}) matrix associated to \code{data} (a \eqn{2 \times d}{2 x d}-matrix), where \eqn{N=2\times d}.
 }

 The \code{distance} parameter allows different kinds of spatial distances:
 \enumerate{
 \item \code{Eucl}, Euclidean distance (default);
 \item \code{Chor}, chordal distance;
 \item \code{Geod}, geodesic distance.
 }

 The \code{likelihood} parameter represents the composite-likelihood configuration:
 \enumerate{
 \item \code{Conditional}, composite likelihood formed by conditionals;
 \item \code{Marginal}, composite likelihood formed by marginals (default);
 \item \code{Full}, standard likelihood.
 }

 It must be coupled with \code{type}:
 \enumerate{
 \item \code{Pairwise}, composite likelihood based on pairs;
 \item \code{Independence}, composite likelihood based on independence;
 \item \code{Standard}, standard likelihood.
 }

 Observation coordinates are checked for exact duplicates before fitting. For space-time data, the full space-time point must be unique; the same spatial site observed at different times is allowed. The check is hash-based and is performed once at the user-facing fit, not during bootstrap/refit iterations that reuse an already validated design. If stochastic thinning retains no pairwise contributions, fitting stops with an explicit error rather than optimizing an empty criterion.

For \code{model="PoissonGamma"} and \code{model="PoissonGammaZIP"}, pairwise fitting accepts every finite \code{shape>0}. This is the continuous Kibble--Gamma pairwise extension used by the marginal moments, correlation function, and bivariate probabilities. The constructive random-field simulation based on a finite sum of squared Gaussian fields is more restrictive and requires \eqn{2\,shape} to be a positive integer; \code{GeoSim} enforces that simulation constraint explicitly.

Stochastic thinning of nearest-neighbor pairs can be enabled via \code{p_neighb<1}. The argument \code{thin_method}
 controls the thinning scheme (default \code{"bernoulli"}).
}
\section{Spatio-temporal ordering}{
For fixed spatial locations, \code{data} is a \eqn{T \times N}{T x N}
matrix: row \eqn{t} corresponds to \code{coordt[t]} and column \eqn{i} to
row \eqn{i} of \code{coordx}.  Observation-level quantities are ordered as
\code{c(t(data))}, hence \code{X} and a vector \code{fixed$mean} use the
order time then site.

For dynamic locations, \code{coordx_dyn}, \code{data}, and optionally a
list-valued \code{X} have one aligned element per time.  Element \eqn{t}
contains the coordinates, responses, and covariate rows observed at
\code{coordt[t]}; internal concatenation is by increasing list index.
See \code{\link{GeoModels-spacetime-ordering}} for the complete convention.
}

\value{
 Returns an object of class \code{GeoFit}.
 An object of class \code{GeoFit} is a list containing at most the following components:
 \item{bivariate}{Logical: \code{TRUE} if the random field is bivariate, otherwise \code{FALSE}.}
 \item{clic}{The composite information criterion after a \code{GeoVarest} call; if the full likelihood is considered then it coincides with AIC.}
 \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates.}
 \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates.}
 \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates.}
 \item{coordx_dyn}{A list of dynamical (in time) spatial coordinates.}
 \item{conf.int}{Confidence intervals for standard maximum likelihood estimation.}
 \item{convergence}{A string that denotes if convergence is reached.}
 \item{copula}{The type of copula.}
 \item{corrmodel}{The correlation model.}
 \item{data}{The vector/matrix/array (or list) of data.}
 \item{distance}{The type of spatial distance.}
 \item{fixed}{A list of fixed parameters.}
 \item{iterations}{The number of iterations used by the numerical routine.}
 \item{likelihood}{The configuration of the composite likelihood.}
 \item{logCompLik}{The value of the log composite-likelihood at the maximum.}
 \item{maxdist}{The maximum spatial distance used in the weighted composite likelihood (or \code{NULL}).}
 \item{maxtime}{The maximum temporal-distance threshold used in the composite likelihood, expressed in the same units as \code{coordt}.}
 \item{message}{Extra message passed from the numerical routines.}
 \item{model}{The density associated to the likelihood objects.}
 \item{estimation_model}{The inferential model/working likelihood supplied to
 \code{GeoFit}. For misspecified fits this differs from \code{model}, which
 records the data-generating response model used by simulation, diagnostics,
 and prediction routines.}
 \item{missp}{\code{TRUE} if a misspecified Gaussian model is used in the composite likelihood.}
 \item{n}{The Binomial number of trials, or the common Negative-Binomial number \eqn{r} of successes.}
 \item{neighb}{The order of spatial neighborhood in the composite likelihood computation.}
 \item{ns}{The number of (different) location sites in the bivariate case.}
 \item{numcoord}{The number of spatial coordinates.}
 \item{numtime}{The number of temporal realisations.}
 \item{param}{A list of parameter estimates.}
 \item{radius}{The radius of the sphere in the case of great-circle distance.}
 \item{stderr}{Standard errors for standard maximum likelihood estimation.}
 \item{sensmat}{The sensitivity matrix.}
 \item{varcov}{The variance-covariance matrix of the estimates.}
 \item{type}{The type of likelihood objects.}
 \item{X}{The matrix of covariates.}
}

\references{

 \bold{General Composite-likelihood}:

 Varin, C., Reid, N. and Firth, D. (2011). An overview of composite
 likelihood methods. \emph{Statistica Sinica}, \bold{21}, 5--42.

 Varin, C. and Vidoni, P. (2005). A note on composite likelihood inference
 and model selection. \emph{Biometrika}, \bold{92}, 519--528.

 \bold{Non-Gaussian random fields}:

 Alegría, A., Caro, S., Bevilacqua, M., Porcu, E. and Clarke, J. (2017).
 Estimating covariance functions of multivariate skew-Gaussian random fields
 on the sphere. \emph{Spatial Statistics}, \bold{22}, 388--402.
 \doi{10.1016/j.spasta.2017.05.003}

 Alegría, A., Bevilacqua, M. and Porcu, E. (2016). Likelihood-based inference
 for multivariate space-time wrapped-Gaussian fields.
 \emph{Journal of Statistical Computation and Simulation},
 \bold{86}(13), 2583--2597.

 Bevilacqua, M., Caamaño-Carrillo, C. and Gaetan, C. (2020). On modelling
 positive continuous data with spatio-temporal dependence.
 \emph{Environmetrics}, \bold{31}(7), e2628. \doi{10.1002/env.2628}

 Bevilacqua, M., Caamaño-Carrillo, C., Arellano-Valle, R. B. and
 Morales-Oñate, V. (2021). Non-Gaussian geostatistical modeling using
 (skew) t processes. \emph{Scandinavian Journal of Statistics},
 \bold{48}(1), 212--245. \doi{10.1111/sjos.12447}

 Blasi, F., Caamaño-Carrillo, C., Bevilacqua, M. and Furrer, R. (2022).
 A selective view of climatological data and likelihood estimation.
 \emph{Spatial Statistics}, \bold{50}, 100596.
 \doi{10.1016/j.spasta.2022.100596}

 Bevilacqua, M., Caamaño-Carrillo, C., Arellano-Valle, R. B. and
 Gómez, C. (2022). A class of random fields with two-piece marginal
 distributions for modeling point-referenced data with spatial outliers.
 \emph{TEST}, \bold{31}(3), 644--674. \doi{10.1007/s11749-021-00797-5}

 Morales-Navarrete, D., Bevilacqua, M., Caamaño-Carrillo, C. and
 Castro, L. M. (2024). Modelling point referenced spatial count data:
 A Poisson process approach. \emph{Journal of the American Statistical
 Association}, \bold{119}(545), 664--677.
 \doi{10.1080/01621459.2022.2140053}

 Bevilacqua, M., Alvarado, E. and Caamaño-Carrillo, C. (2024).
 A flexible Clayton-like spatial copula with application to bounded support
 data. \emph{Journal of Multivariate Analysis}, \bold{201}, 105277.
 \doi{10.1016/j.jmva.2023.105277}

 \bold{Weighted composite-likelihood for (non-)Gaussian random fields}:

 Bevilacqua, M., Gaetan, C., Mateu, J. and Porcu, E. (2012).
 Estimating space and space-time covariance functions for large data sets:
 a weighted composite likelihood approach.
 \emph{Journal of the American Statistical Association, Theory and Methods},
 \bold{107}, 268--280. \doi{10.1080/01621459.2011.646928}

 Bevilacqua, M. and Gaetan, C. (2015). Comparing composite likelihood methods
 based on pairs for spatial Gaussian random fields.
 \emph{Statistics and Computing}, \bold{25}(5), 877--892.
 \doi{10.1007/s11222-014-9471-3}

 Caamaño-Carrillo, C., Bevilacqua, M., López, C. and Morales-Oñate, V. (2024).
 Nearest neighbours weighted composite likelihood based on pairs for
 (non-)Gaussian massive spatial data with an application to Tukey-hh random
 fields estimation. \emph{Computational Statistics and Data Analysis},
 \bold{191}, 107887. \doi{10.1016/j.csda.2023.107887}

 Bevilacqua, M., Cuevas-Pacheco, F. and Caamaño-Carrillo, C. (2026).
Fast stochastic nearest neighbor pairwise composite likelihood for massive
spatial datasets. \emph{arXiv preprint}, arXiv:2607.06142.
\doi{10.48550/arXiv.2607.06142}

}


\seealso{
\code{\link{GeoCovmatrix}} for covariance matrix construction,
  \code{\link{GeoSim}} for simulation,
  \code{\link{GeoKrig}} for prediction,
  \code{\link{GeoVariogram}} and \code{\link{GeoWLS}} for variogram-based tools,
  \code{\link{GeoVarest}} for bootstrap variance estimation.
}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}
 
\examples{
library(GeoModels)

###############################################################
############ Examples of spatial Gaussian random fields ################
###############################################################


# Define the spatial-coordinates of the points:
set.seed(3)
N=100 # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

# Define spatial matrix covariates and regression parameters
X=cbind(rep(1,N),runif(N))
mean <- 0.2
mean1 <- -0.5

# Set the covariance model's parameters:
corrmodel <- "Matern"
sill <- 1
nugget <- 0
scale <- 0.2/3
smooth=0.5


param<-list(mean=mean,mean1=mean1,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

# Simulation of the spatial Gaussian random field:
data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param,X=X)$data



################################################################
###
### Example 0. Maximum independence composite likelihood fitting of
### a Gaussian random field (no dependence parameters)
### 
###############################################################
# setting starting parameters to be estimated
start<-list(mean=mean,mean1=mean1,sill=sill)

fit1 <- GeoFit(data=data,coordx=coords,likelihood="Marginal",
 type="Independence", start=start,X=X)
print(fit1)


################################################################
###
### Example 1. Maximum conditional pairwise likelihood fitting of
### a Gaussian random field using Nelder-Mead
### 
###############################################################
# setting fixed and starting parameters to be estimated
fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill)

fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, 
 neighb=3,likelihood="Conditional",optimizer="Nelder-Mead",
 type="Pairwise", start=start,fixed=fixed,X=X)
print(fit1)

################################################################
###
### Example 2. Maximum stochastic marginal pairwise likelihood 
### fitting of a Gaussian random field using Nelder-Mead
### 
###############################################################

#N=100000 # number of location sites
#x <- runif(N, 0, 1)
#y <- runif(N, 0, 1)
#coords <- cbind(x,y)
#X=cbind(rep(1,N),runif(N))
#data <- GeoSimapprox(coordx=coords,method="TB",L=20000,parallel=TRUE,
#	  corrmodel=corrmodel, param=param,X=X)$data
#fixed<-list(nugget=nugget)
#start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill,smooth=smooth)

#fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
# neighb=3,likelihood="Marginal",optimizer="Nelder-Mead",
# p_neighb=0.2, thin_method="bernoulli",
# type="Pairwise", start=start,fixed=fixed,X=X)
#print(fit2)

################################################################
###
### Example 3. Standard Maximum likelihood fitting of
### a Gaussian random field using nlminb
###
###############################################################
# Define the spatial-coordinates of the points:
set.seed(3)
N=100 # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

param<-list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param)$data

# setting fixed and parameters to be estimated
fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,scale=scale,sill=sill)

I=Inf
lower<-list(mean=-I,scale=0,sill=0)
upper<-list(mean=I,scale=I,sill=I)
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
 optimizer="nlminb",upper=upper,lower=lower,
 likelihood="Full",type="Standard", 
 start=start,fixed=fixed)
print(fit2)


###############################################################
############ Examples of spatial non-Gaussian random fields #############
###############################################################


################################################################
###
### Example 4. Maximum pairwise likelihood fitting of a Weibull random field 
### with Generalized Wendland correlation with Nelder-Mead
### 
###############################################################
set.seed(524)
# Define the spatial-coordinates of the points:
N=100
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
X=cbind(rep(1,N),runif(N))
mean=1; mean1=2 # regression parameters
nugget=0
shape=2
scale=0.2
smooth=0

model="Weibull"
corrmodel="GenWend"
param=list(mean=mean,mean1=mean1,scale=scale,
 shape=shape,nugget=nugget,power2=4,smooth=smooth)
# Simulation of a non-stationary Weibull random field:
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model,X=X,
 param=param)$data


fixed<-list(nugget=nugget,power2=4,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,shape=shape)

# Maximum independence likelihood:
fit <- GeoFit(data=data, coordx=coords, X=X,
 likelihood="Marginal", type="Independence", corrmodel=corrmodel,
 model=model, start=start, fixed=fixed)
print(unlist(fit$param))

## estimating dependence parameter fixing vector mean parameter
Xb <- as.numeric(X \%*\% unlist(fit$param)[1:2])
fixed<-list(nugget=nugget,power2=4,smooth=smooth,mean=Xb)
start<-list(scale=scale,shape=shape)

# Maximum conditional composite-likelihood fitting of the random fields:
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 neighb=3,likelihood="Conditional",type="Pairwise",
 optimizer="Nelder-Mead",
 start=start,fixed=fixed)
print(unlist(fit1$param))



### joint estimation of the dependence parameter and mean parameters
fixed<-list(nugget=nugget,power2=4,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,shape=shape)
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 neighb=3,likelihood="Conditional",type="Pairwise",X=X,
 optimizer="Nelder-Mead",
 start=start,fixed=fixed)
print(unlist(fit2$param))



################################################################
###
### Example 5. Maximum pairwise likelihood fitting of
### a Skew-Gaussian spatial random fields with Wendland correlation
###
###############################################################
set.seed(261)
model="SkewGaussian"
# Define the spatial-coordinates of the points:
x <- runif(100, 0, 1)
y <- runif(100, 0, 1)
coords <- cbind(x,y)

corrmodel="Wend0"
mean=0;nugget=0
sill=1
skew=-4.5
power2=4
c_supp=0.2

# model parameters
param=list(power2=power2,skew=skew,
 mean=mean,sill=sill,scale=c_supp,nugget=nugget)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model, param=param)$data

plot(density(data))
fixed=list(power2=power2,nugget=nugget)
start=list(scale=c_supp,skew=skew,mean=mean,sill=sill)
lower=list(scale=0,skew=-I,mean=-I,sill=0)
upper=list(scale=I,skew=I,mean=I,sill=I)
# Maximum marginal pairwise likelihood:
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 neighb=3,likelihood="Marginal",type="Pairwise",
 optimizer="bobyqa",lower=lower,upper=upper,
 start=start,fixed=fixed)
print(unlist(fit1$param))


################################################################
###
### Example 6. Maximum pairwise likelihood fitting of 
### a Bernoulli random field with exponential correlation 
###
###############################################################

set.seed(422)
N=100
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
mean=0.1; mean1=0.8; mean2=-0.5 # regression parameters
X=cbind(rep(1,N),runif(N),runif(N)) # matrix covariates
corrmodel <- "Wend0"
param=list(mean=mean,mean1=mean1,mean2=mean2,nugget=0,scale=0.2,power2=4)
# Simulation of the spatial Binomial-Gaussian random field:
data <- GeoSim(coordx=coords, corrmodel=corrmodel, model="Binomial", n=1,X=X,
 param=param)$data


## estimating the marginal parameters using independence cl
fixed <- list(power2=4,scale=0.2,nugget=0)
start <- list(mean=mean,mean1=mean1,mean2=mean2)

# Maximum independence likelihood:
fit <- GeoFit(data=data, coordx=coords, n=1, X=X,
 likelihood="Marginal", type="Independence", corrmodel=corrmodel,
 model="Binomial", start=start, fixed=fixed)
 
print(fit)


## estimating dependence parameter fixing vector mean parameter
Xb <- as.numeric(X \%*\% unlist(fit$param))
fixed <- list(nugget=0,power2=4,mean=Xb)
start <- list(scale=0.2)
lower <- list(scale=0)
upper <- list(scale=2)
# Maximum Marginal pairwise likelihood:
fit1 <- GeoFit(data=data, coordx=coords, corrmodel=corrmodel, n=1,
 likelihood="Marginal", type="Pairwise", neighb=3,
 model="Binomial", start=start, fixed=fixed,
 lower=list(scale=0.001), upper=list(scale=1))
 
print(fit1)


## estimating jointly marginal and dependence parameters
fixed <- list(nugget=0,power2=4)
start <- list(mean=mean,mean1=mean1,mean2=mean2,scale=0.2)

# Maximum conditional pairwise likelihood:
fit2 <- GeoFit(data=data, coordx=coords, corrmodel=corrmodel, n=1, X=X,
 likelihood="Marginal", type="Pairwise", neighb=3,
 model="Binomial", start=start, fixed=fixed)
 
print(fit2)


###############################################################
######### Examples of Gaussian spatio-temporal random fields ###########
###############################################################
set.seed(52)
# Define the temporal sequence:
time <- seq(1, 5, 1)

# Define the spatial-coordinates of the points:
x <- runif(20, 0, 1)
y <- runif(20, 0, 1)
coords=cbind(x,y)

# Set the covariance model's parameters:
scale_s=0.2/3;scale_t=1
smooth_s=0.5;smooth_t=0.5
sill=1
nugget=0
mean=0

param<-list(mean=0,scale_s=scale_s,scale_t=scale_t,
 smooth_t=smooth_t, smooth_s=smooth_s ,sill=sill,nugget=nugget)

# Simulation of the spatio-temporal Gaussian random field:
data <- GeoSim(coordx=coords,coordt=time,corrmodel="Matern_Matern",
 param=param)$data

################################################################
###
### Example 7. Maximum pairwise likelihood fitting of a
### space time Gaussian random fields with double-exponential correlation
###
###############################################################
# Fixed parameters
fixed<-list(nugget=nugget,smooth_s=smooth_s,smooth_t=smooth_t)
# Starting value for the estimated parameters
start<-list(mean=mean,scale_s=scale_s,scale_t=scale_t,sill=sill)

# Maximum composite-likelihood fitting of the random fields:
fit <- GeoFit(data=data,coordx=coords,coordt=time,
 corrmodel="Matern_Matern",maxtime=1,neighb=3,
 likelihood="Marginal",type="Pairwise",
 start=start,fixed=fixed)
print(fit)



###############################################################
######### Examples of a bivariate Gaussian random field ###########
###############################################################

################################################################
### Example 8. Maximum pairwise likelihood fitting of a
### bivariate Gaussian random fields with separable Bivariate matern 
### (cross) correlation model 
###############################################################

# Define the spatial-coordinates of the points:
set.seed(89)
x <- runif(100, 0, 1)
y <- runif(100, 0, 1)
coords=cbind(x,y)
# parameters
param=list(mean_1=0,mean_2=0,scale=0.1,smooth=0.5,sill_1=1,sill_2=1,
 nugget_1=0,nugget_2=0,pcol=0.2)

# Simulation of a spatial bivariate Gaussian random field:
data <- GeoSim(coordx=coords, corrmodel="Bi_Matern_sep", 
 param=param)$data

# selecting fixed and estimated parameters
fixed=list(mean_1=0,mean_2=0,nugget_1=0,nugget_2=0,smooth=0.5)
start=list(sill_1=var(data[1,]),sill_2=var(data[2,]),
 scale=0.1,pcol=cor(data[1,],data[2,]))


# Maximum marginal pairwise likelihood
fitcl<- GeoFit(data=data, coordx=coords, corrmodel="Bi_Matern_sep",
 likelihood="Marginal",type="Pairwise",
 start=start,fixed=fixed,
 neighb=3)
print(fitcl)

}
\keyword{Composite}
