\name{GeoCovmatrix}
\alias{GeoCovmatrix}
\encoding{UTF-8}
\title{Spatial and Spatio-temporal Covariance Matrix of (non-)Gaussian random fields}
\description{
 The function computes the covariance matrix associated to
 a spatial or spatio-temporal Gaussian or supported non-Gaussian random field
 with given underlying covariance model and a set of spatial location sites (and temporal instants).
}
\usage{
GeoCovmatrix(estobj = NULL, coordx = NULL, coordy = NULL, coordz = NULL, coordt = NULL,
 coordx_dyn = NULL, corrmodel, distance = "Eucl", grid = FALSE,
 model = "Gaussian", n = 1, param, anisopars = NULL, radius = 1,
 sparse = FALSE, copula = NULL, X = NULL, spobj = NULL,
 check.duplicates = FALSE)
}
\arguments{
 \item{estobj}{An object of class \code{Geofit} that includes information about data, model and estimates.}
 \item{coordx}{A numeric (\eqn{d \times 2})-matrix or (\eqn{d \times 3})-matrix.
 Coordinates on a sphere for a fixed radius \code{radius}
 are passed in lon/lat format expressed in decimal degrees. The default is \code{NULL};
 this is useful when coordinates are supplied through \code{coordx_dyn}, \code{estobj}, or \code{spobj}.}
 \item{coordy}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, the default is \code{NULL}.}
 \item{coordz}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, the default is \code{NULL}.}
 \item{coordt}{A numeric vector giving one dimension of temporal coordinates. Optional argument; the default is \code{NULL}, in which case a spatial random field is expected. At the moment implemented only for the Gaussian case. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}
 \item{coordx_dyn}{For dynamic locations, a list of length \eqn{T}, where \eqn{T} is \code{length(coordt)}.  Element \code{coordx_dyn[[t]]} contains the two- or three-dimensional locations at \code{coordt[t]}; covariance-matrix rows and columns concatenate these blocks by time.}
 \item{corrmodel}{String; the name of a correlation model; see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}
 \item{distance}{String; the name of the spatial distance. The default is \code{"Eucl"}, the Euclidean distance. See \code{\link{GeoFit}}.}
 \item{grid}{Logical; if \code{FALSE} (the default) the data
 are interpreted as spatial or spatio-temporal realisations on a set
 of non-equispaced spatial sites (irregular grid). See \code{\link{GeoFit}}.}
 \item{n}{Positive integer size parameter. For a direct Binomial random field it may be a scalar or one value per observation. For a direct Negative-Binomial random field it is the common number \eqn{r} of successes and must be scalar; the Geometric model is the special case \eqn{r=1}. Default is \eqn{1}.}
 \item{model}{String; the type of random field. See \code{\link{GeoFit}}.}
 \item{param}{A list of parameter values required for the covariance model.}
 \item{anisopars}{A list of two elements \code{"angle"} and \code{"ratio"}, i.e. the anisotropy angle and the anisotropy ratio, respectively.}
 \item{radius}{Numeric; a value indicating the radius of the sphere
 when using covariance models valid for the great circle distance. Default value is \code{1}.}
 \item{sparse}{Logical; if \code{TRUE} the function returns an object of class \code{spam}.
 This option should be used when a parametric compactly supported covariance is
 used. Default is \code{FALSE}.}
 \item{copula}{String; optional copula specification. Observed-scale covariance-matrix construction is implemented for \code{"Gaussian"}, \code{"Clayton"}, and \code{"SkewGaussian"} for the supported continuous copula margins.}
 \item{X}{Numeric design matrix used to evaluate the mean when required by the marginal covariance.  For fixed sites, rows follow time-major order; for dynamic sites, supply a stacked matrix in temporal-block order or a list aligned with \code{coordx_dyn}.}
 \item{spobj}{An object of class \code{sp} or \code{spacetime}. For space-time objects, the current \code{sp2Geo()} conversion uses sequential temporal indices and does not preserve irregular original time spacing; to retain irregular temporal distances, use the explicit-coordinate interface with numeric \code{coordt}.}
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
}
\section{Spatio-temporal ordering}{
For fixed locations, covariance rows and columns are ordered
\preformatted{
(s[1], t[1]), ..., (s[N], t[1]),
(s[1], t[2]), ..., (s[N], t[2]), ...
}
For dynamic locations, all rows of \code{coordx_dyn[[1]]} come first, then all
rows of \code{coordx_dyn[[2]]}, and so on.  The rows of \code{X} and an
external \code{param$mean} vector use this same order.  See
\code{\link{GeoModels-spacetime-ordering}}.
}

\value{
 Returns an object of class \code{GeoCovmatrix}.
 An object of class \code{GeoCovmatrix} is a list containing
 at most the following components:
 \item{bivariate}{Logical: \code{TRUE} if the Gaussian random field is bivariate, otherwise \code{FALSE}.}
 \item{coordx}{A \eqn{d}-dimensional vector of spatial coordinates.}
 \item{coordy}{A \eqn{d}-dimensional vector of spatial coordinates.}
 \item{coordt}{A \eqn{t}-dimensional vector of temporal coordinates.}
 \item{coordx_dyn}{A list of \eqn{t} matrices of spatial coordinates.}
 \item{covmatrix}{The covariance matrix if \code{type} is \code{Standard}. An object of
 class \code{spam} if \code{type} is \code{Tapering} or \code{Standard} and
 \code{sparse} is \code{TRUE}. For \code{model="Wrapped"}, this is the
 covariance matrix of \eqn{\sin(\Theta-\mu)}, not an ordinary covariance
 matrix of numerically represented angles.}
 \item{wrapped_cos_covmatrix}{For \code{model="Wrapped"}, the covariance
 matrix of \eqn{\cos(\Theta-\mu)}. This component is present only for \code{model="Wrapped"}.}
 \item{wrapped_circular_correlation}{For \code{model="Wrapped"}, the
 Jammalamadaka--Sarma circular correlation matrix. This component is present only for \code{model="Wrapped"}.}
 \item{corrmodel}{String: the correlation model.}
 \item{distance}{String: the type of spatial distance.}
 \item{grid}{Logical: the value supplied by the user, indicating whether the original spatial coordinates described a regular grid.}
 \item{nozero}{In the case of tapered matrix the percentage of non-zero values in the
 covariance matrix; otherwise \code{NULL}.}
 \item{n}{The Binomial number of trials, or the common Negative-Binomial number \eqn{r} of successes.}
 \item{namescorr}{String: the names of the correlation parameters.}
 \item{numcoord}{Numeric: the number of spatial coordinates.}
 \item{numtime}{Numeric: the number of temporal coordinates.}
 \item{model}{The type of RF, see \code{\link{GeoFit}}.}
 \item{param}{Numeric: the scalar covariance and marginal parameters.}
 \item{mean_external}{For a univariate model, the optional known mean vector supplied through \code{param$mean}. It is stored separately from \code{param}, in observation order, so that covariance and prediction routines do not interpret its entries as regression coefficients.}
 \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if spatial covariance model.}
 \item{sparse}{Logical: is the returned object of class \code{spam}?}
 \item{call}{The matched call used to create the object.}
}
\details{


For plain (non-copula) bivariate correlation models, observed-scale covariance is currently validated only for Gaussian margins; unsupported non-Gaussian bivariate combinations are rejected explicitly. For direct Negative-Binomial fields, the public argument \code{n} denotes the single common integer \eqn{r} in the waiting-time construction.
For direct \code{StudentT} and \code{SkewStudentT} fields, the common scale-mixture
denominator uses the nugget-free latent correlation. In \code{SkewStudentT}, the
absolute-Gaussian component also uses that correlation; the additive Gaussian
component carries the nugget. This agrees with \code{\link{GeoSim}} and the
cross-covariances used by \code{\link{GeoKrig}}. Setting \code{skew=0} recovers
the direct Student-t covariance, including when \code{nugget>0}.
For \code{copula="Gaussian"}, \code{copula="Clayton"}, and \code{copula="SkewGaussian"}, continuous-margin covariance matrices are constructed directly on the observed scale. The implementation supports \code{"Gaussian"}, \code{"StudentT"}, \code{"LogGaussian"}, \code{"Gamma"}, \code{"Weibull"}, \code{"Beta"}, \code{"Beta2"}, \code{"Kumaraswamy"}, \code{"Kumaraswamy2"}, \code{"Logistic"}, \code{"SkewLaplace"}, \code{"Tukeyh"}, \code{"Tukeyh2"}, and \code{"SinhAsinh"}. Under the Gaussian copula, covariance for the last three is evaluated with the exact covariance of the corresponding monotone transformed-Gaussian field. For location-dependent copula margins, pairwise covariances use the two site-specific marginal predictors supplied through \code{X} or an external mean vector. Clayton-like covariance uses cached deterministic copula quadrature; skew-Gaussian covariance uses a cached bivariate Hermite expansion of the latent representation, avoiding pair-by-pair adaptive two-dimensional integration. The count margins \code{"Poisson"}, \code{"Binomial"}, and \code{"BinomialNeg"} are also available for Gaussian, Clayton-like, and skew-Gaussian copula covariance in the global prediction path; Binomial currently requires a common trial count \code{n}.
The same observed-scale covariance engines are used by \code{\link{GeoCovariogram}} and \code{\link{GeoKrig}}, so covariance matrices, fitted semivariogram curves, and linear prediction are numerically consistent.

For \code{model="Wrapped"}, let \eqn{\Theta=(\mu+Z) \bmod 2\pi}, with
latent Gaussian variance \eqn{\sigma^2} given by \code{sill} and off-diagonal latent
correlation \eqn{r_{ij}=(1-\tau)\rho_{ij}}, where \eqn{\tau} is \code{nugget}. Ordinary covariance
of the numeric angles is not used because it depends on the arbitrary branch
cut. Instead, \code{covmatrix} contains
\deqn{\operatorname{Cov}[\sin(\Theta_i-\mu_i),\sin(\Theta_j-\mu_j)]
= e^{-\sigma^2}\sinh(\sigma^2 r_{ij}),}
with the exact marginal variance on the diagonal. The returned
\code{wrapped_cos_covmatrix} contains
\deqn{e^{-\sigma^2}
[\cosh(\sigma^2 r_{ij})-1],}
and \code{wrapped_circular_correlation} contains
\deqn{\frac{\sinh(\sigma^2 r_{ij})}{\sinh(\sigma^2)}.}
These matrices are evaluated with numerically stable exponential formulas.


For univariate models whose covariance depends on the marginal location, the
mean is evaluated as \eqn{X\beta}, where coefficients \code{mean},
\code{mean1}, and so on correspond in order to the columns of \code{X}.  With
\code{X=NULL}, an intercept-only design is used.  A vector
\code{param$mean} may instead supply one known mean per observation and is
mutually exclusive with \code{X}.

 In the spatial case, the covariance matrix of the
 random vector \deqn{[Z(s_1),\ldots,Z(s_n)]^T}
 with a specific spatial covariance model is computed. Here \eqn{n} is the number of spatial location sites.

 In the space-time case, the covariance matrix of the
 random vector \deqn{[Z(s_1,t_1),Z(s_2,t_1),\ldots,Z(s_n,t_1),\ldots,Z(s_n,t_m)]^T}
 with a specific space-time covariance model is computed. Here \eqn{m} is the number of temporal instants.

 In the bivariate case, rows and columns are in variable-block order.
 For common locations the covariance matrix is associated with
 \deqn{[Z_1(s_1),\ldots,Z_1(s_n),Z_2(s_1),\ldots,Z_2(s_n)]^T.}
 For heterotopic locations supplied by \code{coordx_dyn=list(coords1,coords2)},
 the first block contains all observations of variable 1 at \code{coords1}
 and the second contains all observations of variable 2 at \code{coords2}.

 The location site \eqn{s_i} can be a point in the \eqn{d}-dimensional Euclidean space with \eqn{d=2} or \eqn{d=3}
 or a point (given in lon/lat degree format) on a sphere of arbitrary radius.

 A list with all implemented spatial, space-time and bivariate
 correlation models is given below.
 The argument \code{param} is a list including all the parameters of a given
 correlation model specified by the argument \code{corrmodel}.
 For each correlation model one can check the associated parameters' names using \code{\link{CorrParam}}.
 In what follows
 \eqn{\kappa>0}, \eqn{\beta>0}, \eqn{\alpha, \alpha_s, \alpha_t \in (0,2]} and \eqn{\gamma \in [0,1]}.
 The associated parameters in the argument \code{param} are
 \code{smooth}, \code{power2}, \code{power}, \code{power_s}, \code{power_t} and \code{sep} respectively.
 Moreover let \eqn{1(A)=1} when \eqn{A} is true and \eqn{0} otherwise.

 \itemize{
 \item Spatial correlation models:
 \enumerate{
 \item \code{GenCauchy} (generalised Cauchy in Gneiting and Schlather 2004) defined as:
 \deqn{R(h) = ( 1+h^{\alpha} )^{-\beta / \alpha}}
 If \eqn{h} is the geodesic distance then \eqn{\alpha \in (0,1]}.

 \item \code{Matern} defined as:
 \deqn{R(h) = 2^{1-\kappa} \Gamma(\kappa)^{-1} h^\kappa K_\kappa(h)}
 If \eqn{h} is the geodesic distance then \eqn{\kappa \in (0,0.5]}.

 \item \code{Kummer} (Kummer hypergeometric in Ma and Bhadra 2022) defined as:
 \deqn{R(h) = \Gamma(\kappa+\alpha) U(\alpha,1-\kappa,0.5 h^2 ) / \Gamma(\kappa+\alpha)}
 where \eqn{U(.,.,.)} is the Kummer hypergeometric function.
 If \eqn{h} is the geodesic distance then \eqn{\kappa \in (0,0.5]}.

 \item \code{Kummer_Matern} It is a rescaled version of the \code{Kummer} model,
 i.e. \eqn{h} must be divided by \eqn{(2(1+\alpha))^{0.5}}. When \eqn{\alpha} goes to infinity it is the Matern model.

 \item \code{Wave} defined as:
 \deqn{R(h)=\sin(h)/h}
 This model is valid only for dimensions less than or equal to 3.

 \item \code{GenWend} (Generalized Wendland in Bevilacqua et al. 2019) defined as:
 \deqn{R(h) = A (1-h^2)^{\beta+\kappa} F(\beta/2,(\beta+1)/2,2\beta+\kappa+1,1-h^2) 1(h \in [0,1])}
 where \eqn{\mu \ge 0.5(d+1)+\kappa}
 and \eqn{A=(\Gamma(\kappa)\Gamma(2\kappa+\beta+1))/(\Gamma(2\kappa)\Gamma(\beta+1-\kappa)2^{\beta+1})} and \eqn{F(.,.,.)}
 is the Gaussian hypergeometric function.
 The cases \eqn{\kappa=0,1,2} correspond to the \code{Wend0}, \code{Wend1} and \code{Wend2} models respectively.

 \item \code{GenWend_Matern} (Generalized Wendland Matern in Bevilacqua et al. 2022).
 It is defined as a rescaled version of the Generalized Wendland, i.e.
 \eqn{h} must be divided by \eqn{(\Gamma(\beta+2\kappa+1)/\Gamma(\beta))^{1/(1+2\kappa)}}.
 When \eqn{\beta} goes to infinity it is the Matern model.

 \item \code{GenWend_Matern2} (Generalized Wendland Matern second parametrisation).
 It is defined as a rescaled version of the Generalized Wendland, i.e.
 \eqn{h} must be multiplied by \eqn{\beta} and the smoothness parameter is \eqn{\kappa-0.5}.
 When \eqn{\beta} goes to infinity it is the Matern model.

 \item \code{Hypergeometric} (Hypergeometric model in Bevilacqua et al. 2025).

 \item \code{Hypergeometric_Matern} (Hypergeometric model first parametrisation).

 \item \code{Hypergeometric_Matern2} (Hypergeometric model second parametrisation).

 \item \code{Multiquadric} defined as:
 \deqn{R(h) = (1-\alpha 0.5)^{2\beta}/(1+(\alpha 0.5)^2-\alpha \cos(h))^{\beta}, \quad h \in [0,\pi]}
 This model is valid on the unit sphere and \eqn{h} is the geodesic distance.

 \item \code{Sinpower} defined as:
 \deqn{R(h) = 1-(\sin(h/2))^{\alpha},\quad h \in [0,\pi]}
 This model is valid on the unit sphere and \eqn{h} is the geodesic distance.

 \item \code{F_Sphere} (F family in Alegria et al. 2021) defined as:
 \deqn{R(h) = K F(1/\alpha,1/\alpha+0.5,2/\alpha+0.5+\kappa),\quad h \in [0,\pi]}
 where \eqn{K =(\Gamma(a)\Gamma(i))/(\Gamma(i)\Gamma(o))}.
 This model is valid on the unit sphere and \eqn{h} is the geodesic distance.
 }

 \item Spatio-temporal correlation models:
 \itemize{
 \item Non-separable models:
 \enumerate{
 \item \code{Gneiting} defined as:
 \deqn{R(h, u) = \exp(-h^{\alpha_s}/((1+u^{\alpha_t})^{0.5 \gamma \alpha_s}))/(1+u^{\alpha_t})}

 \item \code{Gneiting_GC} defined as:
 \deqn{R(h, u) = \exp(-u^{\alpha_t}/((1+h^{\alpha_s})^{0.5 \gamma \alpha_t}))/(1+h^{\alpha_s})}
 where \eqn{h} can be either Euclidean or geodesic distance.

 \item \code{Iacocesare} defined as:
 \deqn{R(h, u) = (1+h^{\alpha_s}+u^{\alpha_t})^{-\beta}}

 \item \code{Porcu} defined as:
 \deqn{R(h, u) = (0.5 (1+h^{\alpha_s})^\gamma + 0.5 (1+u^{\alpha_t})^\gamma)^{-\gamma^{-1}}}

 \item \code{Porcu1} defined as:
 \deqn{R(h, u) = \exp(-h^{\alpha_s}(1+u^{\alpha_t})^{0.5 \gamma \alpha_s})/((1+u^{\alpha_t})^{1.5})}

 \item \code{Stein} defined as:
 \deqn{R(h, u) = (h^{\psi(u)}K_{\psi(u)}(h))/(2^{\psi(u)}\Gamma(\psi(u)+1))}
 where \eqn{\psi(u)=\nu+u^{0.5\alpha_t}}.

 \item \code{Gneiting_mat_S} defined as:
 \deqn{R(h, u) = \phi(u)^{\tau_t} \mathrm{Mat}(h \phi(u)^{-\beta},\nu_s)}
 where \eqn{\phi(u)=(1+u^{0.5\alpha_t})}, \eqn{\tau_t \ge 3.5+\nu_s}, \eqn{\beta \in [0,1]}.

 \item \code{Gneiting_mat_T} defined by interchanging \eqn{h} with \eqn{u} in \code{Gneiting_mat_S}.

 \item \code{Gneiting_wen_S} defined as:
 \deqn{R(h, u) = \phi(u)^{\tau_t} \mathrm{GenWend}(h \phi(u)^{\beta},\nu_s,\mu_s)}
 where \eqn{\phi(u)=(1+u^{0.5\alpha_t})}, \eqn{\tau_t \ge 2.5+2\nu_s}, \eqn{\beta \in [0,1]}.

 \item \code{Gneiting_wen_T} defined by interchanging \eqn{h} with \eqn{u} in \code{Gneiting_wen_S}.

 \item \code{Matern_Matern_nosep} (code 89) is the fixed-margin model
 \deqn{R(h,u)=M_s(h)M_t(u)+\zeta Q_s(h)Q_t(u),}
 where \eqn{Q_j(v)=2\nu_j\{M_{\nu_j,\beta_j}(v)-M_{\nu_j+1,\beta_j}(v)\}}.
 The existing parameters \code{scale_s}, \code{scale_t}, \code{smooth_s},
 \code{smooth_t}, and \code{sep} represent \eqn{\beta_s,\beta_t,\nu_s,\nu_t,\rho}.
 All scales and smoothness parameters are positive. The coefficient
 \eqn{\zeta} is obtained from the dimension-dependent normalization of
 \eqn{\rho=\code{sep}\in[-1,1]}; zero gives the separable product.
 See \code{\link{GeoFixedMarginST}} for the formula and valid metrics.

 \item \code{Genwend_Genwend_nosep} (code 85; alias
 \code{GenWend_GenWend_nosep}) is the fixed-margin compactly supported model
 \deqn{R(h,u)=W_s(h)W_t(u)+\eta D_s(h)D_t(u),}
 where \eqn{W_j(v)=\mathcal{GW}_{\kappa_j,\mu_j}(v/a_j)} and
 \eqn{D_j(v)=\mathcal{GW}_{\kappa_j,\mu_j+1}(v/a_j)-W_j(v)}.
 The two terms in each difference have the same support radius.
 The existing parameters \code{scale_s}, \code{scale_t} are the actual support
 radii \eqn{a_s,a_t}; \code{power2_s}, \code{power2_t} are \eqn{\mu_s,\mu_t};
 and \code{smooth_s}, \code{smooth_t} are \eqn{\kappa_s,\kappa_t}.
 Require \eqn{\kappa_j>-1/2},
 \eqn{\mu_s>(d+1)/2+\kappa_s+2}, and \eqn{\mu_t>\kappa_t+3}.
 The parameter \code{sep} is the normalized interaction \eqn{\rho\in[-1,1]},
 not the coefficient \eqn{\eta}; zero gives separability.
 The correlation is exactly zero if \eqn{h\ge a_s} or \eqn{u\ge a_t}.
 See \code{\link{GeoFixedMarginST}} for the normalization.

 \item \code{Multiquadric_st} defined as:
 \deqn{R(h, u)= ((1-0.5\alpha_s)^2/(1+(0.5\alpha_s)^2-\alpha_s \psi(u) \cos(h)))^{a_s},\quad h \in [0,\pi]}
 where \eqn{\psi(u)=(1+(u/a_t)^{\alpha_t})^{-1}}.
 This model is valid on the unit sphere and \eqn{h} is the geodesic distance.

 \item \code{Sinpower_st} defined as:
 \deqn{R(h, u)=(\exp(\alpha_s \cos(h) \psi(u)/a_s)(1+\alpha_s \cos(h) \psi(u)/a_s))/k}
 where \eqn{\psi(u)=(1+(u/a_t)^{\alpha_t})^{-1}} and \eqn{k=(1+\alpha_s/a_s)\exp(\alpha_s/a_s)}, \eqn{h \in [0,\pi]}.
 This model is valid on the unit sphere and \eqn{h} is the geodesic distance.
 }

 \item Separable models:

 Space-time separable correlation models are easily obtained as the product of a spatial and a temporal correlation model, that is
 \deqn{R(h,u)=R(h) R(u)}
 Several combinations are possible:
 \enumerate{
 \item \code{Exp_Exp} defined as:
 \deqn{R(h, u) = \exp(-h)\exp(-u)}

 \item \code{Matern_Matern} defined as:
 \deqn{R(h, u) = \mathrm{Matern}(h;\kappa_s)\mathrm{Matern}(u;\kappa_t)}

 \item \code{GenWend_GenWend} defined as:
 \deqn{R(h, u) = \mathrm{GenWend}(h;\kappa_s,\mu_s)\mathrm{GenWend}(u;\kappa_t,\mu_t)}

 \item \code{Stable_Stable} defined as:
 \deqn{R(h, u) = \exp(-h^{\alpha_s})\exp(-u^{\alpha_t})}
 }
 Note that some models are nested (e.g. \code{Exp_Exp} within \code{Matern_Matern}).
 }

 \item Spatial bivariate correlation models:
 \enumerate{
 \item \code{Bi_Matern} (Bivariate full Matern model)
 \item \code{Bi_Matern_contr} (Bivariate Matern model with constraints)
 \item \code{Bi_Matern_sep} (Bivariate separable Matern model)
 \item \code{Bi_LMC} (Bivariate linear model of coregionalization)
 \item \code{Bi_LMC_contr} (Bivariate LMC with constraints)
 \item \code{Bi_Wendx} (Bivariate full Wendland model)
 \item \code{Bi_Wendx_contr} (Bivariate Wendland model with constraints)
 \item \code{Bi_Wendx_sep} (Bivariate separable Wendland model)
 \item \code{Bi_F_Sphere} (Bivariate full F model on the unit sphere)
 }
 }

 \bold{Remarks}:\cr
 In what follows we assume
 \eqn{\sigma^2,\sigma_1^2,\sigma_2^2,\tau^2,\tau_1^2,\tau_2^2,a,a_s,a_t,a_{11},a_{22},a_{12},\kappa_{11},\kappa_{22},\kappa_{12},f_{11},f_{12},f_{21},f_{22}}
 positive.

 The associated names of the parameters in \code{param} are
 \code{sill}, \code{sill_1}, \code{sill_2},
 \code{nugget}, \code{nugget_1}, \code{nugget_2},
 \code{scale}, \code{scale_s}, \code{scale_t}, \code{scale_1}, \code{scale_2}, \code{scale_12},
 \code{smooth_1}, \code{smooth_2}, \code{smooth_12}, \code{a_1}, \code{a_12}, \code{a_21}, \code{a_2}
 respectively.

 Let \eqn{R(h)} be a spatial correlation model given in standard notation.
 Then the covariance model applied with arbitrary variance, nugget and scale equals to \eqn{\sigma^2} if \eqn{h=0} and
 \deqn{C(h)=\sigma^2(1-\tau^2)R(h/a,\ldots), \quad h>0}
 with nugget parameter \eqn{\tau^2} between 0 and 1.

 Similarly, if \eqn{R(h,u)} is a spatio-temporal correlation model given in standard notation,
 then the covariance model is \eqn{\sigma^2} if \eqn{h=0} and \eqn{u=0} and
 \deqn{C(h,u)=\sigma^2(1-\tau^2)R(h/a_s,u/a_t,\ldots), \quad h>0, u>0}

 Here \sQuote{\ldots} stands for additional parameters.

 The bivariate models implemented are the following:
 \enumerate{
 \item \code{Bi_Matern} defined as:
 \deqn{C_{ij}(h)=\rho_{ij}(\sigma_i\sigma_j+\tau_i^2 1(i=j,h=0))\mathrm{Matern}(h/a_{ij},\kappa_{ij}), \quad i,j=1,2,\; h\ge 0}
 where \eqn{\rho=\rho_{12}=\rho_{21}} is the colocated correlation parameter and \eqn{\rho_{ii}=1}.
 The model \code{Bi_Matern_sep} (separable Matern) is a special case when
 \eqn{a=a_{11}=a_{12}=a_{22}} and \eqn{\kappa=\kappa_{11}=\kappa_{12}=\kappa_{22}}.
 The model \code{Bi_Matern_contr} (constrained Matern)
 is a special case when \eqn{a_{12}=0.5(a_{11}+a_{22})} and \eqn{\kappa_{12}=0.5(\kappa_{11}+\kappa_{22})}.

 \item \code{Bi_GenWend} defined as:
 \deqn{C_{ij}(h)=\rho_{ij}(\sigma_i\sigma_j+\tau_i^2 1(i=j,h=0))\mathrm{GenWend}(h/a_{ij},\nu_{ij},\kappa_{ij}), \quad i,j=1,2,\; h\ge 0}
 where \eqn{\rho=\rho_{12}=\rho_{21}} is the colocated correlation parameter and \eqn{\rho_{ii}=1}.
 The model \code{Bi_GenWend_sep} (separable GenWendland) is a special case when
 \eqn{a=a_{11}=a_{12}=a_{22}} and \eqn{\mu=\mu_{11}=\mu_{12}=\mu_{22}}.
 The model \code{Bi_GenWend_contr} (constrained GenWendland)
 is a special case when \eqn{a_{12}=0.5(a_{11}+a_{22})} and \eqn{\mu_{12}=0.5(\mu_{11}+\mu_{22})}.

 \item \code{Bi_LMC} defined as:
 \deqn{C_{ij}(h)=\sum_{k=1}^{2}(f_{ik}f_{jk}+\tau_i^2 1(i=j,h=0))R(h/a_k)}
 where \eqn{R(h)} is a correlation model.
 The model \code{Bi_LMC_contr}
 is a special case when \eqn{f=f_{12}=f_{21}}.
 Bivariate LMC models, in the current version of the package, are obtained with \eqn{R(h)} equal to the exponential correlation model.
 }
}
\references{
 Alegria, A., Cuevas-Pacheco, F., Diggle, P. and Porcu, E. (2021).
 The F-family of covariance functions: A Matérn analogue for modeling random fields on spheres.
 \emph{Spatial Statistics} 43, 100512.

 Bevilacqua, M., Faouzi, T., Furrer, R. and Porcu, E. (2019).
 Estimation and prediction using generalized Wendland functions under fixed domain asymptotics.
 \emph{Annals of Statistics} 47(2), 828--856.

 Bevilacqua, M., Caamaño-Carrillo, C. and Porcu, E. (2022).
 Unifying compactly supported and Matérn covariance functions in spatial statistics.
 \emph{Journal of Multivariate Analysis} 189, 104949.

 Daley, D. J., Porcu, E. and Bevilacqua, M. (2015).
 Classes of compactly supported covariance functions for multivariate random fields.
 \emph{Stochastic Environmental Research and Risk Assessment} 29(4), 1249--1263.

 Emery, X. and Alegria, A. (2022).
 The Gauss hypergeometric covariance kernel for modeling second-order stationary random fields in Euclidean spaces: its compact support, properties and spectral representation.
 \emph{Stochastic Environmental Research and Risk Assessment} 36, 2819--2834.

 Gneiting, T. (2002).
 Nonseparable, stationary covariance functions for space-time data.
 \emph{Journal of the American Statistical Association} 97, 590--600.

 Gneiting, T., Kleiber, W. and Schlather, M. (2010).
 Matérn cross-covariance functions for multivariate random fields.
 \emph{Journal of the American Statistical Association} 105, 1167--1177.

 Ma, P. and Bhadra, A. (2022).
 Beyond Matérn: on a class of interpretable confluent hypergeometric covariance functions.
 \emph{Journal of the American Statistical Association}, 1--14.

 Porcu, E., Bevilacqua, M. and Genton, M. (2015).
 Spatio-temporal covariance and cross-covariance functions of the great circle distance on a sphere.
 \emph{Journal of the American Statistical Association}.
 DOI: 10.1080/01621459.2015.1072541.

 Gneiting, T. and Schlather, M. (2004).
 Stochastic models that separate fractal dimension and the Hurst effect.
 \emph{SIAM Review} 46, 269--282.
}
\seealso{
\code{\link{GeoFit}} for parameter estimation,
  \code{\link{GeoSim}} for simulation,
  \code{\link{GeoKrig}} for prediction,
  \code{\link{CorrParam}} for correlation-parameter names.
}
\author{
 Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com}, \url{https://sites.google.com/view/moreno-bevilacqua/home}

 Víctor Morales-Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}

 Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}, \url{https://www.researchgate.net/profile/Christian-Caamano}
}
\examples{
library(GeoModels)
################################################################
###
### Example 1. Estimated spatial covariance matrix associated to
### a WCL estimates using the Matern correlation model
###
###############################################################

set.seed(3)
N <- 100 # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x, y)
# Set the covariance model's parameters:
corrmodel <- "Matern"
mean <- 0.5
sill <- 1
nugget <- 0
scale <- 0.2/3
smooth <- 0.5

param <- list(mean = mean, sill = sill, nugget = nugget, scale = scale, smooth = smooth)
data <- GeoSim(coordx = coords, corrmodel = corrmodel, param = param)$data
fixed <- list(nugget = nugget, smooth = smooth)
start <- list(mean = mean, scale = scale, sill = sill)

fit0 <- GeoFit(data = data, coordx = coords, corrmodel = corrmodel,
 neighb = 3, likelihood = "Conditional", optimizer = "BFGS",
 type = "Pairwise", start = start, fixed = fixed)
print(fit0)

# estimated covariance matrix using Geofit object
mm <- GeoCovmatrix(fit0)$covmatrix
# estimated covariance matrix
mm1 <- GeoCovmatrix(coordx = coords, corrmodel = corrmodel,
 param = c(fit0$param, fit0$fixed))$covmatrix
sum(mm - mm1)

################################################################
###
### Example 2. Spatial covariance matrix associated to
### the Generalized Wendland-Matern correlation model
###
###############################################################

# Correlation Parameters for Gen Wendland model
CorrParam("GenWend_Matern")
# Gen Wendland Parameters
param <- list(sill = 1, scale = 0.04, nugget = 0, smooth = 0, power2 = 1/1.5)

matrix2 <- GeoCovmatrix(coordx = coords, corrmodel = "GenWend_Matern",
 param = param, sparse = TRUE)

# Percentage of non-zero values
matrix2$nozero

################################################################
###
### Example 3. Spatial covariance matrix associated to
### the Kummer correlation model
###
###############################################################

# Correlation Parameters for Kummer model
CorrParam("Kummer")
param <- list(sill = 1, scale = 0.2, nugget = 0, smooth = 0.5, power2 = 1)

matrix3 <- GeoCovmatrix(coordx = coords, corrmodel = "Kummer", param = param)

matrix3$covmatrix[1:4, 1:4]

################################################################
###
### Example 4. Covariance matrix associated to
### the space-time double Matern correlation model
###
###############################################################

# Define the temporal coordinates:
times <- seq(1, 4, 1)

# Correlation Parameters for double Matern model
CorrParam("Matern_Matern")

# Define covariance parameters
param <- list(scale_s = 0.3, scale_t = 0.5, sill = 1, smooth_s = 0.5, smooth_t = 0.5)

# Simulation of a spatial Gaussian random field:
matrix4 <- GeoCovmatrix(coordx = coords, coordt = times,
 corrmodel = "Matern_Matern", param = param)

dim(matrix4$covmatrix)

################################################################
###
### Example 5. Spatial covariance matrix associated to
### a skew Gaussian RF with Matern correlation model
###
###############################################################

param <- list(sill = 1, scale = 0.3/3, nugget = 0, skew = 4, smooth = 0.5)
# Simulation of a spatial Gaussian random field:
matrix5 <- GeoCovmatrix(coordx = coords, corrmodel = "Matern", param = param,
 model = "SkewGaussian")

# covariance matrix
matrix5$covmatrix[1:4, 1:4]

################################################################
###
### Example 6. Spatial covariance matrix associated to
### a Weibull RF with GenWend correlation model
###
###############################################################

param <- list(scale = 0.3, nugget = 0, shape = 4, mean = 0, smooth = 1, power2 = 5)
# Simulation of a spatial Gaussian random field:
matrix6 <- GeoCovmatrix(coordx = coords, corrmodel = "GenWend", param = param,
 sparse = TRUE, model = "Weibull")

# Percentage of non-zero values
matrix6$nozero

################################################################
###
### Example 7. Spatial covariance matrix associated to
### a binomial Gaussian RF with Generalized Wendland correlation model
###
###############################################################

param <- list(mean = 0.2, scale = 0.2, nugget = 0, power2 = 4, smooth = 0)
# Simulation of a spatial Gaussian random field:
matrix7 <- GeoCovmatrix(coordx = coords, corrmodel = "GenWend", param = param,
 n = 5, sparse = TRUE, model = "Binomial")

as.matrix(matrix7$covmatrix)[1:4, 1:4]

################################################################
###
### Example 8. Covariance matrix associated to
### a bivariate Matern exponential correlation model
###
###############################################################

set.seed(8)
# Define the spatial coordinates of the points:
x <- runif(4, 0, 1)
y <- runif(4, 0, 1)
coords <- cbind(x, y)

# Parameters
param <- list(mean_1 = 0, mean_2 = 0, sill_1 = 1, sill_2 = 2,
 scale_1 = 0.1, scale_2 = 0.1, scale_12 = 0.1,
 smooth_1 = 0.5, smooth_2 = 0.5, smooth_12 = 0.5,
 nugget_1 = 0, nugget_2 = 0, pcol = -0.25)

# Covariance matrix
matrix8 <- GeoCovmatrix(coordx = coords, corrmodel = "Bi_matern", param = param)$covmatrix

matrix8
}
