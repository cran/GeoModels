/*
 * Direct asymmetric-Laplace random field: allocation-free bivariate density.
 *
 * This is an EXACT reduction of the two-Kibble convolution, not a copula
 * approximation and not the old finite double-series truncation.
 *
 * Let a=1-rho^2, q=1-p, alpha-beta=p-q, alpha*beta=p*q/a.
 * In opposite quadrants f(x,y)=p*q*alpha*beta/(alpha+beta)
 *                              *exp(-alpha*x_positive+beta*x_negative).
 * In the positive quadrant, put B=p/a+beta and A=p/a-alpha. Then
 * f(x,y)/C = exp(beta*x-alpha*y) P(Pois(A*y)>=Pois(B*x))
 *           +exp(-alpha*x+beta*y) P(Pois(A*x)>Pois(B*y)),
 * where the Poisson variables are independent and C=p*q*alpha*beta/(alpha+beta).
 * The negative quadrant follows by (x,y,p)->(-x,-y,1-p).
 *
 * The Poisson-difference tails are noncentral chi-square probabilities.
 * They are evaluated directly as positive, log-concave ONE-dimensional
 * series, in log space, centered near their mode, with geometric remainder
 * bounds. This avoids subtraction of nearly equal noncentral chi-square
 * tails, allocation in a pair loop, and the old Kmax=50 cutoff.
 *
 * Derivation and numerical-validation details: inst/SkewLaplace-numerics.md.
 */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <float.h>

#define GEO_SL_MAX_TAIL_TERMS 100000
#define GEO_SL_LOG_TAIL_TOL (-29.240459028362647) /* log(2e-13) per side */

static double geo_sl_logadd(double x, double y)
{
    if (x == R_NegInf) return y;
    if (y == R_NegInf) return x;
    if (x < y) { double z = x; x = y; y = z; }
    return x + log1p(exp(y - x));
}

/* log(1-exp(x)), x <= 0, including the limiting values. */
static double geo_sl_log1mexp(double x)
{
    if (x == R_NegInf) return 0.0;
    if (x >= 0.0) return R_NegInf;
    return x < -0.69314718055994530942 ? log1p(-exp(x)) : log(-expm1(x));
}

/* log P(Pois(mu) <= j). Gamma upper tails are evaluated directly in log scale. */
static double geo_sl_logpcdf(double j, double mu)
{
    if (j < 0.0) return R_NegInf;
    if (mu == 0.0) return 0.0;
    if (j == 0.0) return -mu;
    const double ans = pgamma(mu, j + 1.0, 1.0, 0, 1);
    return ans > 0.0 ? 0.0 : ans;
}

/*
 * T_d(lambda,mu) = P(Pois(lambda) >= Pois(mu)+d), lambda <= mu, d=0 or 1.
 * Summands w_k = dpois(k,lambda) * ppois(k-d,mu) are log-concave.
 * After the mode, a decreasing adjacent-term ratio bounds the uncomputed
 * remainder by a geometric series. Both tails must meet the relative bound.
 * Failure returns NaN: a non-converged finite sum is NEVER accepted.
 */
static double geo_sl_log_poisson_tail(double lambda, double mu, int d)
{
    if (!R_FINITE(lambda) || !R_FINITE(mu) || lambda < 0.0 || mu < lambda)
        return R_NaN;
    if (lambda == 0.0) return d == 0 ? -mu : R_NegInf;

    double k0 = floor(sqrt(lambda) * sqrt(mu));
    if (k0 < d) k0 = d;
    /* Beyond this range an integer-index recurrence is not a useful kernel. */
    if (!R_FINITE(k0) || k0 > 1e12) return R_NaN;
    const double j0 = k0 - d;
    const double loglambda = log(lambda), logmu = log(mu);
    const double lw0 = dpois(k0, lambda, 1);
    const double lm0 = dpois(j0, mu, 1);
    const double lf0 = geo_sl_logpcdf(j0, mu);
    if (!R_FINITE(lw0) || !R_FINITE(lm0) || !R_FINITE(lf0)) return R_NaN;
    const double lt0 = lw0 + lf0;
    double total = lt0;
    int up_done = 0, down_done = (k0 == d);

    double k = k0, j = j0, lw = lw0, lm = lm0, lf = lf0, previous = lt0;
    for (int it = 1; it <= GEO_SL_MAX_TAIL_TERMS; ++it) {
        k += 1.0; j += 1.0;
        lw += loglambda - log(k);
        lm += logmu - log(j);
        lf = geo_sl_logadd(lf, lm);
        if ((it & 31) == 0) {
            lw = dpois(k, lambda, 1);
            lm = dpois(j, mu, 1);
            lf = geo_sl_logpcdf(j, mu);
        }
        if (lf > 0.0) lf = 0.0;
        const double term = lw + lf;
        if (!R_FINITE(term)) return R_NaN;
        total = geo_sl_logadd(total, term);
        /* Small upward padding protects the geometric bound from roundoff. */
        const double ratio = term - previous + 1e-12;
        if (ratio < 0.0 && term + ratio - geo_sl_log1mexp(ratio)
                              < total + GEO_SL_LOG_TAIL_TOL) {
            up_done = 1; break;
        }
        previous = term;
    }
    if (!up_done) return R_NaN;

    k = k0; j = j0; lw = lw0; lm = lm0; lf = lf0; previous = lt0;
    for (int it = 1; !down_done && it <= GEO_SL_MAX_TAIL_TERMS; ++it) {
        /* First remove the mass at j from F(j). Re-anchor when subtraction
           could amplify rounding, and at least every eight downward steps. */
        double nextlf;
        if ((it & 7) == 0 || lm - lf > -0.69314718055994530942)
            nextlf = geo_sl_logpcdf(j - 1.0, mu);
        else
            nextlf = lf + geo_sl_log1mexp(lm - lf);
        lw += log(k) - loglambda;
        lm += log(j) - logmu;
        k -= 1.0; j -= 1.0;
        lf = nextlf;
        if ((it & 7) == 0) {
            lw = dpois(k, lambda, 1);
            lm = dpois(j, mu, 1);
        }
        const double term = lw + lf;
        if (!R_FINITE(term)) return R_NaN;
        total = geo_sl_logadd(total, term);
        if (k == d) { down_done = 1; break; }
        const double ratio = term - previous + 1e-12;
        if (ratio < 0.0 && term + ratio - geo_sl_log1mexp(ratio)
                              < total + GEO_SL_LOG_TAIL_TOL) {
            down_done = 1; break;
        }
        previous = term;
    }
    if (!down_done || !R_FINITE(total)) return R_NaN;
    return total > 0.0 ? 0.0 : total;
}

static double geo_sl_log_poisson_ge(double lambda, double mu)
{
    if (lambda <= mu) return geo_sl_log_poisson_tail(lambda, mu, 0);
    const double opposite = geo_sl_log_poisson_tail(mu, lambda, 1);
    return geo_sl_log1mexp(opposite);
}

static double geo_sl_standard_logdensity(double x, double y, double p, double r)
{
    const double q = 1.0 - p;
    const double logpq = log(p) + log1p(-p);
    if (r == 0.0) {
        return 2.0 * logpq + (x >= 0.0 ? -p*x : q*x)
                           + (y >= 0.0 ? -p*y : q*y);
    }
    const double a = (1.0 - r) * (1.0 + r);
    const double product = p*q/a;
    const double difference = p-q;
    const double root = hypot(difference, 2.0*sqrt(product));
    double alpha, beta;
    if (difference >= 0.0) {
        alpha = 0.5*(root + difference); beta = product/alpha;
    } else {
        beta = 0.5*(root - difference); alpha = product/beta;
    }
    const double logconstant = 2.0*logpq - log(a) - log(root);
    /* Includes axes and (0,0), with no 0*log(0). */
    if ((x <= 0.0 && y >= 0.0) || (y <= 0.0 && x >= 0.0))
        return logconstant - alpha*fmax(x,y) + beta*fmin(x,y);

    if (x < 0.0) {
        x = -x; y = -y; p = q;
        const double tmp = alpha; alpha = beta; beta = tmp;
    }
    /* Canonical order makes exchangeability exact up to the input conversion. */
    if (x > y) { const double tmp = x; x = y; y = tmp; }
    const double B = p/a + beta;
    /* A*B=p^2*r^2/a^2; use this identity, avoiding cancellation at r=0. */
    const double pr = p*r/a;
    const double A = (pr/B)*pr;
    const double l1 = beta*x - alpha*y + geo_sl_log_poisson_ge(A*y, B*x);
    const double l2 = -alpha*x + beta*y + geo_sl_log_poisson_tail(A*x, B*y, 1);
    if (ISNAN(l1) || ISNAN(l2)) return R_NaN;
    return logconstant + geo_sl_logadd(l1, l2);
}

double log_biv_skewlaplace(double rho, double z1, double z2,
                          double m1, double m2, double sill, double p)
{
    if (!R_FINITE(rho) || fabs(rho) >= 1.0 || !R_FINITE(sill) || sill <= 0.0 ||
        !R_FINITE(p) || p <= 0.0 || p >= 1.0 || !R_FINITE(z1) || !R_FINITE(z2) ||
        !R_FINITE(m1) || !R_FINITE(m2)) return R_NegInf;
    const double sigma = sqrt(sill);
    const double x = (z1-m1)/sigma, y = (z2-m2)/sigma;
    if (!R_FINITE(x) || !R_FINITE(y)) return R_NegInf;
    return geo_sl_standard_logdensity(x, y, p, fabs(rho)) - log(sill);
}

/* Internal, registered vectorized diagnostic interface; no new public API. */
SEXP geo_skewlaplace_logdensity(SEXP sx, SEXP sy, SEXP sr, SEXP sp, SEXP ss)
{
    if (TYPEOF(sx)!=REALSXP || TYPEOF(sy)!=REALSXP || TYPEOF(sr)!=REALSXP ||
        TYPEOF(sp)!=REALSXP || TYPEOF(ss)!=REALSXP || XLENGTH(sp)!=1 || XLENGTH(ss)!=1)
        Rf_error("Invalid argument types for the SkewLaplace density evaluator.");
    const R_xlen_t n = XLENGTH(sx), ny=XLENGTH(sy), nr=XLENGTH(sr);
    if (!((ny==1 || ny==n) && (nr==1 || nr==n)))
        Rf_error("y and rho must have length one or length(x).");
    SEXP out = PROTECT(Rf_allocVector(REALSXP,n));
    for (R_xlen_t i=0; i<n; ++i) {
        REAL(out)[i] = log_biv_skewlaplace(REAL(sr)[nr==1?0:i], REAL(sx)[i],
                            REAL(sy)[ny==1?0:i], 0.0, 0.0, REAL(ss)[0], REAL(sp)[0]);
        if (ISNAN(REAL(out)[i]))
            Rf_error("SkewLaplace density: numerical tail evaluation did not converge.");
        if ((i & 1023)==0) R_CheckUserInterrupt();
    }
    UNPROTECT(1);
    return out;
}
