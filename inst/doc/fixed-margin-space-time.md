# Fixed-margin space-time models: local GeoModels patch

This source package is based on the uploaded GeoModels 2.3.0 tree. Its local
version is **2.3.0**. The two model codes and existing parameter names are
retained, but the covariance definitions intentionally replace the old ones.
**Stored fits using codes 85 or 89 must be refitted.** Other model definitions
are not replaced.

## Model 85 — `Genwend_Genwend_nosep`

`GenWend_GenWend_nosep` is also accepted as an alias.

| Parameter | Meaning |
|---|---|
| `scale_s`, `scale_t` | actual support radii a_s, a_t |
| `power2_s`, `power2_t` | shape parameters mu_s, mu_t (not reciprocal shapes) |
| `smooth_s`, `smooth_t` | GW smoothness kappa_s, kappa_t > -1/2 |
| `sep` | normalized interaction rho in [-1,1] |

Let W_j(v)=GW(kappa_j,mu_j,v/a_j),
D_j(v)=GW(kappa_j,mu_j+1,v/a_j)-W_j(v). Both terms in each difference have the
same radius. The correlation is W_s W_t + eta(rho) D_s D_t, where

    L_j = d_j / (mu_j + 2*kappa_j + d_j + 1)
    U_j = (2*kappa_j + 1) / mu_j
    A = max(L_s*L_t, U_s*U_t)
    B = max(L_s*U_t, U_s*L_t)
    eta(rho) = 2*rho / (A*(1-rho) + B*(1+rho))

Here d_s is the spatial Euclidean dimension and d_t=1. The implemented theorem
region is strict:

    mu_s > (d_s+1)/2 + kappa_s + 2
    mu_t > kappa_t + 3

The spatial and temporal scales are positive. Correlation is exactly zero if
h >= scale_s OR u >= scale_t. No arbitrary threshold is used for the compact
support. Negative correlations inside that support are retained.

## Model 89 — `Matern_Matern_nosep`

The parameter order remains `scale_s`, `scale_t`, `smooth_s`, `smooth_t`, `sep`.
The scales beta_s,beta_t and smoothness nu_s,nu_t are positive. The Matérn
argument is lag/beta (no sqrt(2*nu) factor).

    Q_j(v) = 2*nu_j * (M_{nu_j,beta_j}(v) - M_{nu_j+1,beta_j}(v))
    A_M = max(d_s, 4*nu_s*nu_t)
    B_M = max(2*d_s*nu_t, 2*nu_s)
    zeta(rho) = 2*rho / (A_M*(1-rho) + B_M*(1+rho))
    R(h,u) = M_s(h)*M_t(u) + zeta(rho)*Q_s(h)*Q_t(u)

This covariance is globally supported. In both models, sep=0 is the separable
product, and sep=-1/1 are the two sharp interaction endpoints. `sep` is neither
the raw coefficient eta/zeta nor a correlation coefficient.

## Dimension and metrics

Coordinate-based functions infer dimension from their coordinate representation
(two or three Euclidean coordinates). Chordal distance uses Euclidean dimension
three. These two models reject geodesic distance; the proposed Euclidean theorem
does not establish validity for that metric.

All public function signatures retain the arguments, order and defaults of
the original uploaded package. There is no new dimension argument, global
option or additional item in `param`.

The lag-only functions `GeoCorrFct()` and `GeoCorrFct_Cop()` use dimension 2
for standalone Euclidean calls and dimension 3 for chordal distance, selected
through the existing `distance` argument. Internally called curves inherit
the surrounding coordinate-based operation's dimension. The dimension is
set and restored internally, including when an operation raises an error.

A lag vector alone cannot determine the Euclidean spatial dimension. For a
three-dimensional fitted model, `GeoCovariogram()` uses the dimension from
the fitted coordinates; standalone Euclidean lag-only curves keep the 2D
convention and are not direct 3D comparisons. Existing coordinate interfaces
continue to accept two or three columns without any extra argument.

## Fitting, prediction and simulation

The replacement is in the shared native correlation dispatch, so the existing
Gaussian/full/pairwise likelihood, covariance, kriging and exact simulation paths
use these kernels. Automatic interaction starts are zero and bounds admit both
signs. GW automatic starts respect the coupled shape/smoothness inequalities.
Invalid user starts get an explicit error; invalid optimizer trials are rejected.
The existing parameter-only covariance derivatives use boundary-aware finite
differences for codes 85 and 89.

The fixed GW support is registered for `GeoCovmatrix(..., sparse=TRUE)`.
Global and local kriging paths propagate the model dimension, including pair-cache
calls in future workers. The new source does not add a new fast two-Kronecker
likelihood solver: it uses GeoModels' existing matrix and likelihood algorithms.
No approximate space-time turning-bands/CE capability is newly claimed.

Gaussian `sill` and `nugget` retain their current package meaning:
K = sill * [(1-nugget)*R + nugget*I] for distinct observed site-time points.
The Gaussian lag-only curve at h=0 and a positive temporal lag must not receive
an extra nugget; that scalar-temporal-lag branch has been corrected.

## Installation from the supplied source archive

In a fresh R session, with the existing GeoModels dependencies and compiler
available:

```r
install.packages("GeoModels_2.3.0.tar.gz", repos=NULL, type="source")
library(GeoModels)
packageVersion("GeoModels")
?GeoFixedMarginST
```

Alternatively unzip the source tree and run `R CMD INSTALL GeoModels`.
Restart an R session that has the old compiled GeoModels DLL loaded before
using the replacement.


## Compiler diagnostic: `Rconfig.h` not found

`src/header.h` already included `<Rconfig.h>` in the user's original GeoModels
2.3.0 source. This is an R header and must be provided by the R development
installation; it is also included indirectly by `<R.h>`. The R manual
recommends it before BLAS/LAPACK declarations when using `USE_FC_LEN_T`.
Do **not** remove the include to suppress the error.

Compile the *package*, not individual `.c` files with `clang`/`make`, from a
terminal using R's configured compilation system:

```sh
R CMD INSTALL --preclean GeoModels_2.3.0.tar.gz
```

On the user's computer, check the active R installation and its headers:

```r
R.home()
R.home("include")
file.exists(file.path(R.home("include"), "Rconfig.h"))
```

From a terminal, `R CMD config --cppflags` shows R's include flags. If the file
is missing, the local R installation/development headers need repair. If it
exists but the compiler still cannot locate it, inspect the full compilation
command for the `-I` path, and the user's `~/.R/Makevars` overrides. Replacing
this R header with a package-local copy is not a portable solution.

## Validation scope

A previous development build's numerical validation report is retained at
`inst/validation/fixed_margin_native_results.json`. It reports 1,335 checks
of the native kernels with independent special-function substitutes. Those
checks are historical: they are not presented as a fresh full-package test
of this interface-only revision. The covariance formulas and their C
implementations are unchanged; only three dimension-related diagnostic
strings change in the C file.

For this revision, source checks compare exported function argument lists
and defaults against the original uploaded package, verify removal of the
added dimension argument and its calls/documentation, and check the updated
R sources for balanced delimiters. The native numerical source is compared
with the preceding version after normalizing the three diagnostic strings.
See `inst/validation/no_new_arguments_checks.json` for the recorded scope. The version identifier in that historical JSON records the development build originally tested; it is not the package version.

R/Rscript is not available in this environment. The full R package has not
been installed, `R CMD check` has not been run, and the R integration suite
has not been executed here. Passing source checks is not a substitute for
these runtime tests.

## R integration checks to run after installation

```r
source(system.file("validation", "Test_FixedMarginST.R", package="GeoModels"))
# Optional composite-likelihood score test:
Sys.setenv(GEOMODELS_TEST_CL="true")
source(system.file("validation", "Test_FixedMarginST.R", package="GeoModels"))
```

The suite includes checks for the original curve-function signatures,
coordinate/metric-based dimension inference, nested-context restoration,
restoration after errors, both public curve evaluators, and no dimension
leakage after a 3D operation. It also retains the checks for fixed margins,
closed forms, nugget, covariance matrices, sparse support, kriging and
interaction-only fitting. These R tests have been supplied but not run here.

A self-contained example is installed at
`system.file("examples", "fixed_margin_space_time.R", package="GeoModels")`.
