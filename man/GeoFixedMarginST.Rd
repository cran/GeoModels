\name{GeoFixedMarginST}
\alias{GeoFixedMarginST}
\title{Fixed-margin nonseparable generalized Wendland and Matern space-time models}
\description{
Covariance definitions for correlation-model codes 85 and 89. The two models
have a normalized signed interaction parameter \code{sep} in \eqn{[-1,1]}.
Varying it preserves both marginal correlations. These definitions replace
the earlier models at the same codes; stored fits must be refitted.
}
\details{
For \code{Genwend_Genwend_nosep} (alias \code{GenWend_GenWend_nosep}), use
\code{scale_s}, \code{scale_t}, \code{power2_s}, \code{power2_t},
\code{smooth_s}, \code{smooth_t}, \code{sep}, in their existing order.
They represent \eqn{a_s,a_t,\mu_s,\mu_t,\kappa_s,\kappa_t,\rho}.
The scales are the actual compact support radii and the shape parameters
are not reciprocals. With \eqn{d_s=d}, \eqn{d_t=1}, require
\eqn{a_j>0}, \eqn{\kappa_j>-1/2}, and
\eqn{\mu_j>(d_j+1)/2+\kappa_j+2}. Define
\deqn{W_j(v)=\mathcal{GW}_{\kappa_j,\mu_j}(v/a_j),\quad
D_j(v)=\mathcal{GW}_{\kappa_j,\mu_j+1}(v/a_j)-W_j(v),}
\deqn{L_j=\frac{d_j}{\mu_j+2\kappa_j+d_j+1},\quad
U_j=\frac{2\kappa_j+1}{\mu_j},}
\deqn{A=\max(L_sL_t,U_sU_t),\quad B=\max(L_sU_t,U_sL_t),\quad
\eta(\rho)=\frac{2\rho}{A(1-\rho)+B(1+\rho)}.}
The correlation is
\deqn{R_\rho(h,u)=W_s(h)W_t(u)+\eta(\rho)D_s(h)D_t(u).}
It vanishes exactly outside the fixed support rectangle; internal zeros
and negative values are allowed. Both terms in a shape difference use
the same radius. For \eqn{\kappa=0},
\eqn{\mathcal{GW}_{0,\mu}(x)=(1-x)_+^\mu}; for \eqn{\kappa=1},
\eqn{\mathcal{GW}_{1,\mu}(x)=(1-x)_+^{\mu+1}\{1+(\mu+1)x\}}.

For \code{Matern_Matern_nosep}, use \code{scale_s}, \code{scale_t},
\code{smooth_s}, \code{smooth_t}, \code{sep}, representing
\eqn{\beta_s,\beta_t,\nu_s,\nu_t,\rho}. Require positive scales and
\eqn{\nu_j>0}. There is no \eqn{\sqrt{2\nu}} factor in the Matern argument:
\deqn{M_j(v)=\frac{2^{1-\nu_j}}{\Gamma(\nu_j)}(v/\beta_j)^{\nu_j}
K_{\nu_j}(v/\beta_j),\quad M_j(0)=1.}
Set \eqn{Q_j(v)=2\nu_j\{M_{\nu_j,\beta_j}(v)-M_{\nu_j+1,\beta_j}(v)\}},
\eqn{Q_j(0)=0}, and
\deqn{A_M=\max(d,4\nu_s\nu_t),\quad B_M=\max(2d\nu_t,2\nu_s),\quad
\zeta(\rho)=\frac{2\rho}{A_M(1-\rho)+B_M(1+\rho)}.}
Then \eqn{R_\rho(h,u)=M_s(h)M_t(u)+\zeta(\rho)Q_s(h)Q_t(u)}.
This model is globally supported. In both families \eqn{\rho=0} means
separability and \eqn{\rho=\pm1} reaches the admissible interaction endpoints.
The normalized parameter is not itself a correlation coefficient.

Coordinate-based functions infer \eqn{d} from the existing coordinate
arguments (two or three Euclidean coordinates). Chordal distance uses the
three-dimensional Euclidean embedding. Geodesic distance is rejected for
these two families because its positive definiteness is not established here.
No public argument, option, or covariance parameter is added for dimension.
Standalone Euclidean \code{GeoCorrFct} and \code{GeoCorrFct_Cop} calls use
\eqn{d=2}; internal lag-only calls inherit the enclosing model dimension.
The dimension context is restored on exit, including after errors.
For a three-dimensional Euclidean fitted model, \code{GeoCovariogram}
uses the dimension of its stored coordinates. A standalone lag vector cannot
identify that dimension and should not be used as a three-dimensional check.

The Gaussian nuisance parameters retain their existing meaning. With
\code{sill} \eqn{=\sigma^2} and \code{nugget} \eqn{=g}, distinct observed
site-time points have covariance matrix
\eqn{\sigma^2\{(1-g)R+gI\}}. Adding a nugget does not preserve the special
exact two-product diagonalization of the noise-free model. This update
uses the package's existing likelihood, simulation and kriging algorithms;
it does not introduce a new Kronecker likelihood engine.
}
\seealso{\code{\link{GeoCorrFct}}, \code{\link{GeoCovmatrix}},
\code{\link{GeoFit}}, \code{\link{GeoKrig}}}
\examples{
p <- list(mean=0, sill=1, nugget=0, scale_s=1, scale_t=1,
          power2_s=4, power2_t=4, smooth_s=0, smooth_t=0, sep=0.7)
z <- GeoCorrFct(x=c(0,0.2,0.6,1.1), t=c(0,0.3,1.2),
                corrmodel="Genwend_Genwend_nosep", param=p)
stopifnot(all(is.finite(z$corr)))

pm <- list(mean=0, sill=1, nugget=0, scale_s=1, scale_t=1,
           smooth_s=1.5, smooth_t=1.5, sep=-0.5)
zm <- GeoCorrFct(x=c(0,0.2,0.6), t=c(0,0.3),
                 corrmodel="Matern_Matern_nosep", param=pm)
stopifnot(all(is.finite(zm$corr)))
}
