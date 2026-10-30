GeoModels 2.3.0 - SkewLaplace dual conditional simulation (development patch)
==========================================================================

What is supported
-----------------
(1) GeoSimcond(..., model="SkewLaplace", copula=NULL)
    Conditions the DIRECT Gamma-difference / four-Gaussian RF of the paper.
    Exact quadratic constraints at all conditioning sites; radial Gibbs-within-
    slice and von Mises angles (C backend); Gaussian conditional substitution
    to the predictions. n_iter is MCMC burn-in, mcmc_thin controls thinning.
    Only global spatial conditioning. method="Cholesky" is exact Gaussian
    propagation; "TB" is an approximate Gaussian propagation.

(2) GeoSimcond(..., model="SkewLaplace", copula="Gaussian")
    Existing Gaussian-copula conditional simulation with SkewLaplace margins;
    unchanged, including local=TRUE support and the marginal PIT transform.

IMPORTANT: these are DIFFERENT joint distributions with the same AL marginals.
The direct model has correlation rho^2 regardless of p. The Gaussian-copula
model uses a Gaussian copula and does not in general have this correlation.

Notes
-----
- Retains marginal independence likelihood skew fix in both R and C from the
  corrected user package.
- The direct sampler is a new radial/von-Mises method, not the geodesic MCMC
  implementation described in the current manuscript. It targets the correct
  latent Gaussian law conditional on the quadratic constraints, including the
  coarea factor, but finite-run MCMC error must be monitored.
- Native routine smoke-tested with standalone, mock-R-API C harness (one and
  two conditioning locations): latent constraints within 7e-15, radial mean
  about 1.004 vs exact 1, angular moments correct. Gaussian single-site
  conditional predictive mean and variance independently checked by Monte Carlo.
- The complete R package could not be installed or tested in this environment
  because R is not installed. Run the accompanying R scripts on your Mac.
- Unrelated issues in the existing direct pairwise density (e.g. high-rho
  truncated series and zero-residual handling) have NOT been changed here.

Examples / tests
----------------
inst/examples/skewlaplace_two_conditionals.R
inst/examples/skewlaplace_conditional_tests.R

Install with R >= 4.4, restart R before and after:
install.packages("~/Downloads/GeoModels_2.3.0_SkewLaplace_DUAL_COND.tar.gz",
                 repos=NULL, type="source")

Then, e.g.:
sys.source(system.file("examples/skewlaplace_conditional_tests.R",
                       package="GeoModels"), envir=new.env(parent=globalenv()))
