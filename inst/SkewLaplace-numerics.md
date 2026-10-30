# SkewLaplace numerical update: GeoModels 2.3.0.1

This is a development source update of the supplied DUAL_COND archive, not a
claim of a CRAN release or a completed R CMD check. Both dependence models are
retained: `copula=NULL` is the four-Gaussian/two-Gamma random field;
`copula="Gaussian"` is the Gaussian copula with asymmetric Laplace margins.
The unconditional models, p convention, sill convention and correlation models
have not been replaced.

## 1. Bivariate density of the direct field

Let q=1-p, a=1-r^2, and let x,y be location-centered residuals divided by
sigma=sqrt(sill). Define positive alpha,beta by

    alpha-beta = p-q,   alpha*beta = p*q/a,
    D = sqrt((p-q)^2 + 4*p*q/a),
    alpha=(D+p-q)/2, beta=(D-p+q)/2.

Rationalized formulas are used in C to avoid cancellation when p is near an
endpoint. Let C=p*q*alpha*beta/(alpha+beta).

For x>=0, y<=0 the exact standardized joint density is

    f(x,y)=C*exp(-alpha*x + beta*y).

Exchangeability covers the other opposite-sign quadrant. This includes the
axes and the origin. In particular,

    f(0,0)=p^2*q^2/sqrt((1-r^2)*(1-r^2*(2*p-1)^2)).

For x,y>0 define B=p/a+beta and A=p/a-alpha. The code computes A from
A*B=p^2*r^2/a^2 rather than by subtracting nearby quantities. The exact density is

    f(x,y)/C = exp(beta*x-alpha*y) * P(Pois(A*y)>=Pois(B*x))
               + exp(-alpha*x+beta*y) * P(Pois(A*x)>Pois(B*y)),

where the two Poisson variables within each probability are independent.
The negative quadrant follows by replacing (x,y,p) by (-x,-y,1-p).
At r=0 the product of the two correctly parameterized marginal densities is
used exactly. The final log-density subtracts log(sill), as required for two
location-scale variables.

### Derivation (equivalent to the manuscript's Kibble double series)

At one observed site with response x>=0, write
A1=p(T+x), B1=q*T. Conditioning independent Exp(1) variables A1,B1 on their
scaled difference shows that T is exactly Exp(1). Conditional on the first
site's latent Gaussian values, the two squared-radius fields at a second
site are independent scaled noncentral chi-squares with two degrees of freedom.
Integrating T gives the conditional characteristic function

    exp(i*u*r^2*x/(1-i*a*u/p))
    / [1 + i*a*u*(p-q)/(p*q) + a*u^2/(p*q)].

The denominator is the characteristic function of E_alpha-E_beta, where the
exponentials have rates alpha,beta. The numerator is the characteristic function
of an independent compound-Poisson sum J of positive Exp(c) jumps, with

    c=p/a,   lambda=r^2*p*x/a.

Consequently Y2 | Y1=x is distributed as E_alpha-E_beta+J. Convolution with the
asymmetric exponential density splits at J=y. Exponential tilting changes the
compound-Poisson rate and intensity. Since 2*c*J has a noncentral chi-square law
with zero degrees of freedom, including its atom at zero, these tilted
probabilities reduce to noncentral chi-square probabilities, equivalently the
Poisson-difference probabilities above. Multiplication by the marginal density
p*q*exp(-p*x) yields the displayed joint density. Thus no new dependence model
or approximation has been introduced.

### Numerical implementation

The old two-dimensional truncation K<=50 and its repeated R_alloc buffers have
been removed. Opposite quadrants require only elementary functions. Same-sign
quadrants require two one-dimensional positive sums, not a double series.
For lambda<=mu and d=0 or 1, the evaluated probability is

    T_d = sum(k>=d) dpois(k,lambda) * ppois(k-d,mu).

Summation starts near k=sqrt(lambda*mu). Poisson masses and cumulative
probabilities are updated by recurrences, periodically re-anchored to R's dpois
and log upper-gamma pgamma evaluations. These summands are log-concave, so a
decreasing adjacent-term ratio bounds the omitted tail by a geometric sum.
Both sides must meet a relative truncation bound of 2e-13 each. This is a
truncation bound in exact arithmetic, not an interval-arithmetic guarantee on
all floating-point or special-function error. Probabilities near one use the
complement of the SMALL directly computed opposite probability, using log1mexp.
No noncentral chi-square upper-tail subtraction is used.

All working state is scalar: no per-pair heap allocations and no static mutable
cache. If the iteration safeguard (100000 terms per side) or numerical range
checks fail, the routine returns NaN, NEVER the last non-converged finite sum.
The internal vector diagnostic reports this as an error. Existing composite
likelihood wrappers reject non-finite evaluations. Extremely nearly singular
correlations with enormous residuals can still exceed the numerical safeguards.
The |r|=1 singular joint law is not an ordinary two-dimensional density.

## 2. Direct conditional sampler

The C radial-angular sampler now evaluates

    log(I0(x)) = x + log(bessel_i(x,0,2))

using the scaled R C API (expo=2). The old expo=1 plus x was incorrect.
The separate hand-truncated large-x log-Bessel approximation has also been
removed. The randomized slice bracket computes its right endpoint BEFORE
clipping a negative left endpoint at zero. The full conditional radial target
and angular von Mises updates otherwise retain the direct field construction
and its correct conditioning Jacobian.

These changes do not turn finite MCMC output into independent exact posterior
draws. Burn-in, mixing, Monte Carlo error and effective sample size still need
to be assessed for each dataset, particularly with strongly correlated sites.

## 3. Reused Gaussian conditional factors

For global Cholesky propagation, one observed covariance factor is used for the
linear solve and the precision. One conditional prediction covariance factor
is then reused over all four latent fields and all retained draws:

    W=R_*o R_oo^{-1}, C=R_**-W R_o*, C=L L',
    Z_k^*=W Z_k + L*epsilon_k.

Only bounded matrix batches of independent Gaussian innovations are generated.
The same prepared conditional Gaussian engine is used for global Cholesky
simulation of the SkewLaplace Gaussian copula (with one latent Gaussian field
and no Gibbs). Duplicate prediction locations share their latent random value;
locations identical to observations are mapped back to the conditioned value.
Only roundoff-scale negative eigenvalues can be clipped in a semidefinite
fallback; no nugget is added to make a materially invalid covariance pass.

The direct Gibbs still uses a dense n x n precision and O(n^2) work per sweep.
It is not parallelized by parallel=TRUE. Dense conditioning and the returned
latent array remain memory limitations for large n. The TB route retains its
previous approximate Gaussian propagation. Local conditioning of the direct
SkewLaplace field and space-time conditional simulation remain unsupported.

## 4. Validation status

The numerical C cores were compiled and exercised in a standalone C++ harness.
The harness uses Boost/long-double gamma and Poisson primitives in place of R's
pgamma/dpois; a checked scaled I0 implementation in place of R's bessel_i; and
an independent uniform RNG. The production C algorithm source is shared with
these tests. This verifies the numerical algorithms, but NOT the full R API,
R package loading, full R dependency integration, or Mac compilation.
R is not installed in the execution environment; no complete GeoModels
installation or R CMD check was run. Post-installation R tests are supplied in
inst/validation/Test_SkewLaplace.R and as a separate attachment.

Validation performed:
- 1287 density cases over p=.001,...,.999, |r|=0,...,.9999, zero, opposite-sign,
  same-sign and negative residual pairs; agreement with an independent SciPy
  noncentral-chi-square expression to maximum absolute log error 5.4e-13.
- Exact exchangeability, reflection and r-sign checks on the same grid.
- Six numerical marginal integrations, maximum relative discrepancy 1.1e-13.
- Independent positive one-dimensional compound-Poisson convolution at
  p=.2,r=.99,(x,y)=(1,2): density 0.0237389248007547, relative discrepancy
  about 4.1e-14. This is more accurate than the older two-dimensional quadrature
  reported during the preceding audit.
- Correct zero-residual density at p=.2,r=.7: 0.0394999647709364.
- Three TWO-observation Gibbs tests at (p,r)=(.75,.6),(.2,.85),(.5,.95), each
  with burn-in 5000, thinning 3 and 160000 retained draws. Radial means and
  within-field cross-site products agree with independent 128/256-node
  Kibble quadrature within 1.1 estimated Monte Carlo standard errors.
- Predictive means for a third site agree with quadrature-based expectations
  within 2.2 estimated Monte Carlo standard errors. Constraint errors below
  8e-15 on these tests.
- Independent-site radial Exp(1) check (100000 draws), including negligible
  empirical cross-site radial correlation.

The standalone benchmarks are not end-to-end R or Mac timings and are not a
claim that the entire fitting procedure has the same speedup.

API reference for scaled Bessel functions:
https://cran.r-project.org/doc/manuals/r-release/R-exts.html
