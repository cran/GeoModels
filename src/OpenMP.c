/*
 * Optional OpenMP support.
 *
 * R headers define a macro named `match`, while recent omp.h versions use
 * match(...) in OpenMP declare-variant pragmas.  Therefore, when OpenMP is
 * actually available, undefine that R macro immediately before including
 * omp.h.  If the compiler defines _OPENMP but omp.h is not available, build
 * the serial fallback rather than failing installation.
 */
#include "header.h"
#include "matern_kernel.h"
#include "matern_bessel.h"

#include <math.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>

#if defined(_OPENMP)
# if defined(__has_include)
#  if __has_include(<omp.h>)
#   ifdef match
#    undef match
#   endif
#   include <omp.h>
#   define GEOMODELS_HAVE_OPENMP 1
#  else
#   define GEOMODELS_HAVE_OPENMP 0
#  endif
# else
#  ifdef match
#   undef match
#  endif
#  include <omp.h>
#  define GEOMODELS_HAVE_OPENMP 1
# endif
#else
# define GEOMODELS_HAVE_OPENMP 0
#endif

/*
 * OpenMP support is deliberately opt-in.  The default is one thread so that
 * existing GeoModels workflows, including outer future/parallel loops, keep
 * their previous behaviour and cannot become oversubscribed unexpectedly.
 *
 * The native parallel path covers the univariate Gaussian spatial pairwise
 * marginal and pairwise conditional likelihoods with the Whittle--Matern
 * correlation model (cormod == 14), including general positive smoothness.
 * All other likelihood/correlation combinations fall back to the serial code.
 */
static int geomodels_omp_threads = 1;
static int geomodels_omp_min_pairs = 50000;

int GeoModels_omp_threads_internal(void)
{
#if GEOMODELS_HAVE_OPENMP
    return geomodels_omp_threads > 1 ? geomodels_omp_threads : 1;
#else
    return 1;
#endif
}

int GeoModels_omp_min_pairs_internal(void)
{
    return geomodels_omp_min_pairs > 0 ? geomodels_omp_min_pairs : 1;
}

void GeoModels_set_omp_threads(int *n)
{
    int requested = (n != NULL) ? *n : 1;

#if GEOMODELS_HAVE_OPENMP
    if (requested == 0)
        requested = omp_get_max_threads();
    if (requested < 1)
        requested = 1;
    geomodels_omp_threads = requested;
#else
    (void) requested;
    geomodels_omp_threads = 1;
#endif

    if (n != NULL)
        *n = geomodels_omp_threads;
}

void GeoModels_set_omp_min_pairs(int *n)
{
    int requested = (n != NULL) ? *n : 50000;
    if (requested < 1)
        requested = 1;
    geomodels_omp_min_pairs = requested;
    if (n != NULL)
        *n = geomodels_omp_min_pairs;
}

void GeoModels_openmp_info(int *compiled, int *threads,
                           int *max_threads, int *min_pairs)
{
#if GEOMODELS_HAVE_OPENMP
    if (compiled != NULL)    *compiled = 1;
    if (threads != NULL)     *threads = GeoModels_omp_threads_internal();
    if (max_threads != NULL) *max_threads = omp_get_max_threads();
#else
    if (compiled != NULL)    *compiled = 0;
    if (threads != NULL)     *threads = 1;
    if (max_threads != NULL) *max_threads = 1;
#endif
    if (min_pairs != NULL) *min_pairs = GeoModels_omp_min_pairs_internal();
}

#if GEOMODELS_HAVE_OPENMP

/* Worker-safe Matérn evaluation.  All Rmath quantities that depend only on
 * smoothness are computed by the main thread before entering OpenMP.  General
 * smoothness uses GeoModels' local allocation-free RKbesl kernel; workers do
 * not call bessel_k[_ex], lgammafn, warning(), error(), R_alloc(), or any other
 * R API.
 */
static inline double geomodels_matern_corr_omp(double lag,
                                               double scale,
                                               double smooth,
                                               int half_integer_k,
                                               double log_den,
                                               double small_log_ratio,
                                               double *bessel_work,
                                               int *ok)
{
    if (ok != NULL) *ok = 1;
    if (!(scale > 0.0) || !(smooth > 0.0) || isnan(lag) || lag < 0.0) {
        if (ok != NULL) *ok = 0;
        return NAN;
    }
    if (lag == 0.0) return 1.0;

    const double a = lag / scale;
    if (isinf(a)) return 0.0;

    if (half_integer_k >= 0)
        return geo_matern_half_corr(a, half_integer_k);

    if (a <= 1.0e-10) {
        const double log_a = a > 0.0 ? log(a) : log(lag) - log(scale);
        if (smooth >= 0.99) return 1.0;
        const double corr = -expm1(small_log_ratio +
                                   2.0 * smooth * (log_a - log(2.0)));
        return fmin(1.0, fmax(0.0, corr));
    }

    if (bessel_work == NULL) {
        if (ok != NULL) *ok = 0;
        return NAN;
    }

    double scaled_bessel = 0.0;
    if (!GeoModels_bessel_k_scaled(a, smooth, bessel_work, &scaled_bessel)) {
        if (ok != NULL) *ok = 0;
        return NAN;
    }

    const double log_corr = smooth * log(a) + log(scaled_bessel) - a - log_den;
    return exp(fmin(0.0, log_corr));
}

#endif /* GEOMODELS_HAVE_OPENMP */

/*
 * Evaluate the Gaussian spatial pair objective in parallel.
 *
 * conditional = 0: bivariate marginal pairwise log likelihood.
 * conditional = 1: log f(z1,z2) - log f(z2).
 *
 * Return 1 when the OpenMP path was used, and 0 when the caller should use
 * the original serial implementation (unsupported model, too few pairs,
 * one requested thread, no OpenMP compiler support, or allocation failure).
 */
int GeoModels_gauss_pair_omp(int conditional,
                             int cormod,
                             const double *data1,
                             const double *data2,
                             const double *par,
                             int weighted,
                             const double *mean1,
                             const double *mean2,
                             int n_pairs,
                             double max_dist,
                             double sill,
                             double nugget,
                             double *result)
{
#if !GEOMODELS_HAVE_OPENMP
    (void) conditional; (void) cormod; (void) data1; (void) data2;
    (void) par; (void) weighted; (void) mean1; (void) mean2;
    (void) n_pairs; (void) max_dist; (void) sill; (void) nugget;
    (void) result;
    return 0;
#else
    if (result == NULL || cormod != 14 || n_pairs < GeoModels_omp_min_pairs_internal())
        return 0;

    int nthreads = GeoModels_omp_threads_internal();
    if (nthreads <= 1)
        return 0;
    if (nthreads > n_pairs)
        nthreads = n_pairs;

    const double corr_scale = par[0];
    const double corr_smooth = par[1];
    if (!(corr_scale > 0.0) || !(corr_smooth > 0.0))
        return 0;

    const int half_integer_k = geo_matern_half_order(corr_smooth);
    if (!isfinite(corr_scale) || !isfinite(corr_smooth) ||
        !isfinite(sill) || sill <= 0.0 || !isfinite(nugget) ||
        nugget < 0.0 || nugget > 1.0)
        return 0;

    const double log_den = (corr_smooth - 1.0) * log(2.0) + lgammafn(corr_smooth);
    double small_log_ratio = 0.0;
    if (half_integer_k < 0 && corr_smooth < 0.99) {
        if (corr_smooth < 0.001) {
            const double n2 = corr_smooth * corr_smooth;
            small_log_ratio = 2.0 * corr_smooth *
                (0.57721566490153286061 +
                 n2 * (0.40068563438653142847 +
                 n2 * (0.20738555102867398527 +
                 n2 * 0.14404989676884611812)));
        } else {
            small_log_ratio = lgammafn(1.0 - corr_smooth) -
                              lgammafn(1.0 + corr_smooth);
        }
    }

    int nb = 0;
    if (half_integer_k < 0) {
        if (corr_smooth > (double) (INT_MAX - 2)) return 0;
        nb = 1 + (int) floor(corr_smooth);
        if (nb < 1) return 0;
    }

    double *partial = (double *) calloc((size_t) nthreads, sizeof(double));
    double *bessel_workspace = NULL;
    if (half_integer_k < 0) {
        const size_t nth = (size_t) nthreads;
        const size_t nwork = nth * (size_t) nb;
        if ((size_t) nb != 0U && nwork / (size_t) nb != nth) {
            free(partial);
            return 0;
        }
        bessel_workspace = (double *) calloc(nwork, sizeof(double));
    }
    if (partial == NULL || (half_integer_k < 0 && bessel_workspace == NULL)) {
        free(partial);
        free(bessel_workspace);
        return 0;
    }

    const double correlation_multiplier = 1.0 - nugget;
    const double log_normal_constant = log(2.0 * M_PI) + log(sill);
    int failed = 0;

#pragma omp parallel num_threads(nthreads) reduction(|:failed)
    {
        const int tid = omp_get_thread_num();
        double local_total = 0.0;
        double *thread_bessel = bessel_workspace == NULL ? NULL :
            bessel_workspace + (size_t) tid * (size_t) nb;

#pragma omp for schedule(static)
        for (int i = 0; i < n_pairs; ++i) {
            const double d1 = data1[i];
            const double d2 = data2[i];

            if (!isnan(d1) && !isnan(d2)) {
                const double lag = lags[i];
                int corr_ok = 1;
                const double corr0 = geomodels_matern_corr_omp(
                    lag, corr_scale, corr_smooth, half_integer_k,
                    log_den, small_log_ratio, thread_bessel, &corr_ok);
                if (!corr_ok || !isfinite(corr0)) {
                    failed = 1;
                    continue;
                }
                const double corr = correlation_multiplier * clamp_corr(corr0);
                const double weight = weighted ? CorFunBohman(lag, max_dist) : 1.0;
                const double log_biv = log_biv_Norm(
                    corr, d1, d2, mean1[i], mean2[i], sill, 0.0);
                if (!isfinite(log_biv)) {
                    failed = 1;
                    continue;
                }

                if (conditional) {
                    const double z2 = (d2 - mean2[i]) / sqrt(sill);
                    const double log_marginal_2 =
                        -0.5 * (log_normal_constant + z2 * z2);
                    local_total += (log_biv - log_marginal_2) * weight;
                } else {
                    local_total += log_biv * weight;
                }
            }
        }

        partial[tid] = local_total;
    }

    if (failed) {
        free(partial);
        free(bessel_workspace);
        return 0; /* rare numerical failure: use the original serial path */
    }

    /* Fixed thread-order combination gives repeatable summation for a fixed
       thread count and static schedule. */
    double total = 0.0;
    for (int t = 0; t < nthreads; ++t)
        total += partial[t];

    free(partial);
    free(bessel_workspace);

    *result = total;
    return 1;
#endif
}
