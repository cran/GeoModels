#ifndef GEOMODELS_MATERN_KERNEL_H
#define GEOMODELS_MATERN_KERNEL_H
#include <math.h>

/* Exact half-integers only: nearby smoothness values are different models. */
static inline int geo_matern_half_order(double nu)
{
    if (nu == 0.5) return 0;
    if (nu == 1.5) return 1;
    if (nu == 2.5) return 2;
    if (nu == 3.5) return 3;
    return -1;
}

/* Elementary, allocation-free and independent of the R API. */
static inline double geo_matern_half_corr(double a, int k)
{
    double poly;
    if (isnan(a) || a < 0.0 || k < 0 || k > 3) return NAN;
    if (a > 1.0e4) return 0.0; /* also handles +Inf without Inf * 0 */
    if (k == 0) return exp(-a);
    if (k == 1) poly = 1.0 + a;
    else if (k == 2) poly = 1.0 + a * (1.0 + a / 3.0);
    else poly = 1.0 + a * (1.0 + a * (0.4 + a / 15.0));
    return fmin(1.0, a < 700.0 ? exp(-a) * poly : exp(log(poly) - a));
}
#endif
