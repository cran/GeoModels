#include "header.h"

/*
   Balanced implementation of Gauss hypergeometric 2F1(a,b;c;x)
   for the main intended case a,b,c > 0 and |x| <= 1.

   Strategy:
   - direct compensated series away from x = 1;
   - Pfaff transformation for sufficiently negative x;
   - analytic continuation in y = 1-x near x = 1;
   - fallback to the existing robust routine hypergeo() in delicate cases:
       * delta = c-a-b close to 0;
       * delta close to an integer;
       * non-positive parameters;
       * non-finite intermediate quantities.

   The tolerance is intentionally balanced, not ultra-aggressive.
*/

#ifndef MAX_ITERATIONS
#define MAX_ITERATIONS 10000
#endif

#ifndef MACHEP
#define MACHEP DBL_EPSILON
#endif

#define CONV_TOL_BALANCED 1.0e-11
#define DELTA_CRIT_TOL    1.0e-8
#define X_DIRECT_LIMIT    0.75
#define X_NEAR_ONE_LIMIT  0.85
#define X_NEG_PFAFF_LIMIT -0.60

static double hyp2f1_series_balanced(double a, double b, double c, double x);
static double hyp2f1_negative_x_balanced(double a, double b, double c, double x);
static double hyp2f1_near_one_balanced(double a, double b, double c, double x);
static int is_half_integer(double x);

double hypergeom_2f1(double a, double b, double c, double x);
double hypergeo2(double a, double b, double c, double x);

static int is_near_integer(double x, double tol)
{
    double r = round(x);
    return fabs(x - r) < tol;
}

static int is_half_integer(double x)
{
    double twice = 2.0 * x;
    double rounded = round(twice);

    if (fabs(twice - rounded) < 1.0e-8) {
        int int_val = (int) rounded;
        return (int_val % 2 != 0);
    }

    return 0;
}

static double exp_lgamma_ratio4(double u1, double u2, double v1, double v2)
{
    double logv = lgammafn(u1) + lgammafn(u2) - lgammafn(v1) - lgammafn(v2);

    if (!R_FINITE(logv)) return NA_REAL;
    if (logv > 700.0) return R_PosInf;
    if (logv < -745.0) return 0.0;

    return exp(logv);
}

/* Direct 2F1 series with Kahan compensated summation. */
static double hyp2f1_series_balanced(double a, double b, double c, double x)
{
    double sum = 1.0;
    double term = 1.0;
    double comp = 0.0;

    double a_n = a;
    double b_n = b;
    double c_n = c;

    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(c) || !R_FINITE(x)) return NA_REAL;
    if (c == 0.0) return NA_REAL;

    for (int n = 1; n < MAX_ITERATIONS; n++) {
        double denom = c_n * (double)n;
        double ratio;

        if (denom == 0.0) return NA_REAL;

        ratio = (a_n * b_n * x) / denom;
        term *= ratio;

        if (!R_FINITE(term)) return NA_REAL;

        /* Kahan summation */
        double y = term - comp;
        double t = sum + y;
        comp = (t - sum) - y;
        sum = t;

        if (!R_FINITE(sum)) return NA_REAL;

        if (fabs(term) <= CONV_TOL_BALANCED * fmax(1.0, fabs(sum))) {
            break;
        }

        /* Optional geometric remainder estimate. */
        {
            double next_denom = (c_n + 1.0) * ((double)n + 1.0);
            if (next_denom != 0.0) {
                double r_est = fabs(((a_n + 1.0) * (b_n + 1.0) * x) / next_denom);
                if (r_est < 1.0) {
                    double rem = fabs(term) * r_est / (1.0 - r_est);
                    if (rem <= CONV_TOL_BALANCED * fmax(1.0, fabs(sum))) {
                        break;
                    }
                }
            }
        }

        a_n += 1.0;
        b_n += 1.0;
        c_n += 1.0;
    }

    return sum;
}

/* Pfaff transformation for x < 0:
   2F1(a,b;c;x) = (1-x)^(-a) 2F1(a,c-b;c;x/(x-1)). */
static double hyp2f1_negative_x_balanced(double a, double b, double c, double x)
{
    double y = x / (x - 1.0);
    double factor = R_pow(1.0 - x, -a);
    double s;

    if (!R_FINITE(factor)) return hypergeo(a, b, c, x);

    s = hyp2f1_series_balanced(a, c - b, c, y);
    if (!R_FINITE(s)) return hypergeo(a, b, c, x);

    return factor * s;
}

/* Analytic continuation near x = 1.

   For y = 1-x and delta = c-a-b:

   2F1(a,b;c;x) =
       A 2F1(a,b;1-delta;y)
     + B y^delta 2F1(c-a,c-b;1+delta;y),

   where
       A = Gamma(c) Gamma(delta) / [Gamma(c-a) Gamma(c-b)],
       B = Gamma(c) Gamma(-delta) / [Gamma(a) Gamma(b)].

   This formula is fast near x = 1 because y is small. However it is
   numerically delicate when delta is near 0 or near an integer. In those
   cases we use hypergeo() as a robust fallback.
*/
static double hyp2f1_near_one_balanced(double a, double b, double c, double x)
{
    double y = 1.0 - x;
    double delta = c - a - b;
    double A, B, s1, s2, yd, result;

    if (y < 0.0) return NA_REAL;

    if (y == 0.0) {
        if (delta > 0.0) {
            double val = exp_lgamma_ratio4(c, delta, c - a, c - b);
            return val;
        }
        return R_PosInf;
    }

    /* Critical cases: avoid cancellation-prone continuation formula. */
    if (fabs(delta) < DELTA_CRIT_TOL || is_near_integer(delta, DELTA_CRIT_TOL)) {
        return hypergeo(a, b, c, x);
    }

    /* The formula also requires non-pole gamma arguments. For the intended
       a,b,c > 0 case this is usually safe, but c-a or c-b may be non-positive
       integers. If so, use fallback. */
    if (is_near_integer(c - a, DELTA_CRIT_TOL) && (c - a) <= 0.0) return hypergeo(a, b, c, x);
    if (is_near_integer(c - b, DELTA_CRIT_TOL) && (c - b) <= 0.0) return hypergeo(a, b, c, x);

    A = exp_lgamma_ratio4(c, delta, c - a, c - b);
    B = exp_lgamma_ratio4(c, -delta, a, b);

    if (!R_FINITE(A) || !R_FINITE(B)) return hypergeo(a, b, c, x);

    s1 = hyp2f1_series_balanced(a, b, a + b - c + 1.0, y);
    s2 = hyp2f1_series_balanced(c - a, c - b, c - a - b + 1.0, y);

    if (!R_FINITE(s1) || !R_FINITE(s2)) return hypergeo(a, b, c, x);

    yd = R_pow(y, delta);
    if (!R_FINITE(yd)) return hypergeo(a, b, c, x);

    /* Simple compensated two-term sum. */
    {
        double term1 = A * s1;
        double term2 = B * yd * s2;
        double comp = 0.0;
        double tmp;
        double yy;

        if (!R_FINITE(term1) || !R_FINITE(term2)) return hypergeo(a, b, c, x);

        result = term1;
        yy = term2 - comp;
        tmp = result + yy;
        comp = (tmp - result) - yy;
        result = tmp;
    }

    if (!R_FINITE(result)) return hypergeo(a, b, c, x);

    return result;
}

/* Main dispatcher. */
double hypergeom_2f1(double a, double b, double c, double x)
{
    double delta;

    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(c) || !R_FINITE(x)) return NA_REAL;
    if (fabs(x) > 1.0) return NA_REAL;

    if (x == 0.0) return 1.0;
    if (a == 0.0 || b == 0.0) return 1.0;

    /* The optimized path is intended for positive a,b,c. For anything else,
       use the existing robust implementation. */
    if (a <= 0.0 || b <= 0.0 || c <= 0.0) {
        return hypergeo(a, b, c, x);
    }

    delta = c - a - b;

    if (fabs(x - 1.0) < MACHEP) {
        if (delta > 0.0) {
            return exp_lgamma_ratio4(c, delta, c - a, c - b);
        }
        return R_PosInf;
    }

    /* For negative x, Pfaff is helpful only when x is sufficiently negative.
       For small negative x the direct series is faster. */
    if (x < X_NEG_PFAFF_LIMIT) {
        return hyp2f1_negative_x_balanced(a, b, c, x);
    }

    if (fabs(x) <= X_DIRECT_LIMIT) {
        return hyp2f1_series_balanced(a, b, c, x);
    }

    /* Intermediate region. Direct series is still efficient enough. */
    if (x < X_NEAR_ONE_LIMIT) {
        return hyp2f1_series_balanced(a, b, c, x);
    }

    /* Near x=1: use analytic continuation, with robust fallback inside. */
    return hyp2f1_near_one_balanced(a, b, c, x);
}

/* Wrapper keeping your original external behavior:
   for half-integer c, use the existing hypergeo(); otherwise use the balanced path. */
double hypergeo2(double a, double b, double c, double x)
{
    if (is_half_integer(c)) {
        return hypergeo(a, b, c, x);
    }

    return hypergeom_2f1(a, b, c, x);
}
