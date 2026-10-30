# Targeted fixes to the GeoModels 2.2.9 source snapshot

Source: GeoModels(10).zip. Patch date: 2026-09-27.
The package version and public function signatures are preserved.

## Corrections

1. **Binomial pair trials:** preserve the site-level `n` vector and create
   `n_col`/`n_row` once. Remove the obsolete concatenation followed by a second
   site-index lookup. Applies to Binomial, BinomialLogistic and the
   Gaussian-misspecified Binomial path.
2. **Bivariate Turning Bands:** pass all 24 registered `for_c` arguments,
   including the unused Matérn `params_other` placeholder; align SIGNATURE and
   INTENT. The C entry point and its arity are unchanged.
3. **Tukey-gh at zero skew:** correct the variance to
   `sill / (1 - 2*tail)^1.5`, retaining the existing correlation normalization.
4. **Student covariance at high degrees of freedom:** add a shared native
   correlation evaluator with log-gamma/asymptotic gamma ratios and a stable
   positive hypergeometric series for large parameters. The missing `nu>=170`
   branch is removed. Finite nu still uses the Student formula, not a Gaussian
   approximation.
5. **Student covariance/kriging agreement:** distinguish the nugget-free
   denominator correlation from the signed, nugget-attenuated numerator
   correlation. Use the same helper in GeoCovmatrix, GeoKrig, GeoCorrFct and
   GeoCovariogram. Copula covariance formulas are not changed.
6. **Gaussian-misspecified Student likelihoods:** keep the sign outside the
   logarithm in spatial, space-time and anisotropic pair/conditional kernels.
   Use the same stable Student evaluator. The affected anisotropic paths now
   attenuate the nugget only once; the conditional marginal uses its specified
   mean rather than the observation as its own mean.
7. **GeoDoScores:** center at the observed marginal mean via an internal
   shared helper. Correct inverse-df conversion for TwoPieceStudentT and
   positive, count, transformed and bounded marginal means. Preserve direct
   versus copula parametrizations. Unsupported mean formulas fail explicitly
   instead of silently leaving data uncentered. Wrapped Euclidean scores are
   explicitly rejected: this patch does not invent an angular scoring rule.
8. **Wrapped likelihood:** use direct log-density evaluation and log-sum-exp
   for bivariate and univariate wrapped normals. All existing pair/conditional
   likelihood callers use the log API; the ordinary-density API is retained.
   Alpha remains **3**. Fixed truncation accuracy at arbitrarily large sills
   is not addressed by this numerical-underflow patch. Circular prediction,
   covariance fields and the `wrapped_*` output convention are unchanged.

No public function or legacy model has been removed. Optimizer controls have
not been changed. The ZIP excludes macOS metadata only (`__MACOSX`, `.DS_Store`).

## Verification actually performed

The modified native function bodies were compiled in an isolated C harness.
R mathematical macros/functions were replaced by C99 equivalents in that
harness, with controlled correlation inputs; this is **not** a full package
compilation or an R integration test. The package's own unchanged 2F1
implementation was included in these native tests.

- 330 Student correlation comparisons against high-precision calculations:
  maximum absolute error approximately 4.69e-13; 60 additional large-nu/limit
  checks passed.
- 900 tests across all six modified Student likelihood kernels, with signed
  correlations, positive nugget, several degrees of freedom and weights:
  maximum absolute error approximately 4.89e-12.
- 120 wrapped log-density comparisons: maximum scaled discrepancy
  approximately 1.57e-14. The underflow example gives log-density
  -994.9301217874274, rather than logging a zero density.
- 360 tests across all five modified wrapped likelihood kernels:
  maximum absolute error approximately 1.82e-12.
- No failures in those 1,770 checks. AddressSanitizer/UndefinedBehaviorSanitizer
  smoke tests on the isolated native harness passed.
- Lexical delimiter/quote checks on the R sources and native-call signature
  checks were performed. Lexical checks are not R parsing.

**R is not installed in the patch environment.** Package installation,
`R CMD check`, the full examples, and end-to-end R regression tests were not
executed. No claim of CRAN or cross-platform validation is made.

## Included R regression tests

After installing this source snapshot and its dependencies, run:

```sh
R CMD INSTALL GeoModels
Rscript GeoModels/tests/regression-audit.R
R CMD check --no-manual GeoModels
```

`tests/regression-audit.R` includes public-interface checks for all eight
corrections: exact heterogeneous-binomial pairs and permutations, an actual
small bivariate TB simulation, Tukey-gh/Tukey-h equality, high-nu Student
covariance, Student BLUP versus a joint covariance matrix, marginal centering,
and finite wrapped likelihoods in the density-underflow regime. These tests
are supplied but have not been run in this environment.
