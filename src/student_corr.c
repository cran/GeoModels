#include "header.h"

/* Log of A_nu = (nu-2)/2 * [Gamma((nu-1)/2)/Gamma(nu/2)]^2.
   The asymptotic expansion is for the gamma ratio, NOT a replacement of
   the Student model by a Gaussian model.  At x >= 1000 the first omitted
   term, 17/(7168*x^7), is below 2.4e-24. */
double GeoStudentCorrLogFactor(double nu)
{
    if (ISNAN(nu) || nu <= 2.0) return R_NaN;
    if (!R_FINITE(nu)) return 0.0;  /* nu = +Inf */
    const double x = 0.5 * (nu - 2.0);
    if (x >= 1000.0) {
        const double r = 1.0 / x;
        const double r2 = r * r;
        return r * (-0.25 + r2 * (1.0 / 96.0 - r2 / 320.0));
    }
    return log(x) + 2.0 * (lgammafn(x + 0.5) - lgammafn(x + 1.0));
}

/* rho_base controls the common denominator field; rho_num is the signed,
   nugget-attenuated numerator correlation.  Never take log(rho_num). */
double GeoStudentCorrFromLogFactor(double rho_base, double rho_num,
                                   double nu, double log_factor)
{
    if (!R_FINITE(rho_base) || !R_FINITE(rho_num) ||
        ISNAN(nu) || nu <= 2.0 || !R_FINITE(log_factor) ||
        fabs(rho_base) > 1.0 + 64.0 * DBL_EPSILON ||
        fabs(rho_num) > 1.0 + 64.0 * DBL_EPSILON)
        return R_NaN;
    rho_base = fmax(-1.0, fmin(1.0, rho_base));
    rho_num = fmax(-1.0, fmin(1.0, rho_num));
    if (rho_num == 0.0) return 0.0;
    if (!R_FINITE(nu)) return rho_num;
    const double z = rho_base * rho_base;
    /* Gauss's endpoint identity: A_nu * 2F1(1/2,1/2;nu/2;1) = 1. */
    if (z == 1.0) return rho_num;
    const double c = 0.5 * nu;
    double hg;
    if (c >= 50.0) {
        /* Positive series, uniformly safe for z in [0,1].  After term k,
           the remaining sum is <= term_k*(k+1)/(c-2): bound each later
           ratio by j/(c+j-1), then use Gauss's summation formula.
           This avoids large gamma parameters in generic 2F1 backends. */
        double term = 1.0, sum = 1.0;
        int converged = 0;
        for (int k = 1; k <= 10000; ++k) {
            const double a = (double)k - 0.5;
            term *= (a / (double)k) * (a / (c + (double)k - 1.0)) * z;
            sum += term;
            if (term * (((double)k + 1.0) / (c - 2.0)) <=
                4.0 * DBL_EPSILON * sum) {
                converged = 1;
                break;
            }
        }
        if (!converged) return R_NaN;
        hg = sum;
    } else {
        hg = hypergeo(0.5, 0.5, c, z);
    }
    if (!R_FINITE(hg) || hg <= 0.0) return R_NaN;
    return rho_num * exp(log_factor + log(hg));
}

void GeoStudentCorrelationVec(double *rho_base, double *rho_num, double *nu,
                              double *result, int *n)
{
    const double log_factor = GeoStudentCorrLogFactor(*nu);
    for (int i = 0; i < *n; ++i)
        result[i] = GeoStudentCorrFromLogFactor(rho_base[i], rho_num[i],
                                               *nu, log_factor);
}

/* Skew-Student construction: the absolute Gaussian and denominator copies
 * use rho_base; only the signed Gaussian numerator uses rho_num.  Reusing
 * the Student factor also gives the exact nu=+Inf and skew=0 limits. */
double GeoSkewStudentCorrFromLogFactor(double rho_base, double rho_num,
                                      double nu, double skew, double log_factor)
{
    if (!R_FINITE(skew) || fabs(skew) > 1.0 || !R_FINITE(rho_num) ||
        fabs(rho_num) > 1.0 + 64.0 * DBL_EPSILON) return R_NaN;
    if (skew == 0.0)
        return GeoStudentCorrFromLogFactor(rho_base, rho_num, nu, log_factor);
    const double ah = GeoStudentCorrFromLogFactor(rho_base, 1.0, nu, log_factor);
    if (!R_FINITE(ah)) return R_NaN;
    rho_base = fmax(-1.0, fmin(1.0, rho_base));
    rho_num = fmax(-1.0, fmin(1.0, rho_num));
    if (rho_base == 0.0 && rho_num == 0.0) return 0.0;
    if (rho_base == 1.0 && rho_num == 1.0) return 1.0;
    const double a = exp(log_factor), skew2 = skew * skew;
    const double abs_product = sqrt(fmax(0.0, 1.0-rho_base*rho_base)) +
                               rho_base * asin(rho_base);
    const double numerator = ah * ((1.0-skew2)*rho_num +
                               skew2*(2.0/M_PI)*(abs_product-1.0)) +
                             skew2*(2.0/M_PI)*(ah-a);
    return numerator / (1.0-skew2*(2.0/M_PI)*a);
}

double corr_skewt2(double rho_base, double rho_num, double nu, double skew)
{
    return GeoSkewStudentCorrFromLogFactor(rho_base, rho_num, nu, skew,
                                          GeoStudentCorrLogFactor(nu));
}
