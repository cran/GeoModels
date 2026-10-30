#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

/* Gaussian-copula covariance polynomial, sum_k coefficient[k] * rho^k.
 * Keep the R helper's Horner order and handling of non-finite inputs, while
 * avoiding one temporary vector per degree.  Inputs are coerced by R.
 */
SEXP GeoGaussianCopulaPolynomial(SEXP rho, SEXP coefficient)
{
    if (TYPEOF(rho) != REALSXP || TYPEOF(coefficient) != REALSXP)
        error("rho and coefficient must be double vectors.");

    const R_xlen_t n = XLENGTH(rho);
    const R_xlen_t degree = XLENGTH(coefficient);
    SEXP result = PROTECT(allocVector(REALSXP, n));
    const double *rr = REAL(rho);
    const double *cc = REAL(coefficient);
    double *out = REAL(result);

    for (R_xlen_t i = 0; i < n; ++i) {
        if (i % 65536 == 0)
            R_CheckUserInterrupt();
        const double r = rr[i];
        if (!R_FINITE(r)) {
            out[i] = NA_REAL;
            continue;
        }
        if (degree == 0) {
            out[i] = 0.0;
            continue;
        }

        double acc = cc[degree - 1];
        for (R_xlen_t k = degree - 1; k > 0; --k)
            acc = cc[k - 1] + r * acc;
        const double value = r * acc;
        out[i] = (fabs(r) < 1e-12 || !R_FINITE(value)) ? 0.0 : value;
    }

    UNPROTECT(1);
    return result;
}

/*
 * Evaluate
 *
 *   sum_{k=1}^degree coefficient[k] * scale[i] * rho[i]^k
 *
 * for every i.  When use_scale == 0 the scale factor is omitted.
 *
 * The loop order intentionally matches the historical R implementation used
 * by the skew-Gaussian covariance engine: powers are accumulated from degree
 * one upward.  This keeps the floating-point evaluation order as close as
 * possible to the previous code while avoiding degree full-vector passes and
 * their temporary allocations.
 */
void GeoCovSeriesForward(double *rho, int *n,
                         double *coefficient, int *degree,
                         double *scale, int *use_scale,
                         double *out)
{
    const int nn = *n;
    const int dd = *degree;
    const int scaled = *use_scale;
    int i, k;

    for (i = 0; i < nn; ++i) {
        const double r = rho[i];
        const double s = scaled ? scale[i] : 1.0;
        double rp = r;
        double acc = 0.0;

        for (k = 0; k < dd; ++k) {
            double ck = coefficient[k];
            if (scaled)
                ck *= s;
            acc += ck * rp;
            rp *= r;
        }

        if (fabs(r) < 1e-15)
            acc = 0.0;
        out[i] = acc;
    }
}


/*
 * Evaluate the location-dependent Gaussian-copula Hermite covariance series
 * without expanding endpoint coefficient matrices to one row per pair.
 *
 * A1 and A2 are nrow1 x degree and nrow2 x degree R matrices in column-major
 * order. idx1/idx2 are one-based row indices for each covariance pair.  The
 * accumulation order in k matches the historical R loop exactly.
 */
void GeoGaussianCovSeriesIndexed(double *rho, int *n,
                                 double *A1, int *nrow1, int *idx1,
                                 double *A2, int *nrow2, int *idx2,
                                 int *degree, double *fac,
                                 double *out)
{
    const int nn = *n;
    const int nr1 = *nrow1;
    const int nr2 = *nrow2;
    const int dd = *degree;
    int i, k;

    for (i = 0; i < nn; ++i) {
        const int r1 = idx1[i] - 1;
        const int r2 = idx2[i] - 1;
        const double r = rho[i];
        double rp = r;
        double acc = 0.0;

        for (k = 0; k < dd; ++k) {
            const double a = A1[r1 + k * nr1];
            const double b = A2[r2 + k * nr2];
            double term = (a * b / fac[k]) * rp;
            if (!isfinite(term))
                term = 0.0;
            acc += term;
            rp *= r;
        }
        out[i] = acc;
    }
}


/*
 * Evaluate a location-dependent covariance series whose coefficient at every
 * degree is obtained by bilinear interpolation on the same eta grid.
 *
 * The gram array contains degree consecutive ng x ng matrices in R column-
 * major order.  i1/i2 are one-based lower-cell indices computed once in R;
 * a1/a2 are the corresponding interpolation fractions.  Keeping these
 * quantities outside the degree loop removes the historical 36 repeated
 * findInterval/interpolation passes for Beta2 and Kumaraswamy2.
 */
void GeoCovSeriesInterp2Forward(double *rho, int *n,
                                double *gram, int *ng, int *degree,
                                int *i1, int *i2,
                                double *a1, double *a2,
                                double *out)
{
    const int nn = *n;
    const int gg = *ng;
    const int dd = *degree;
    const int plane = gg * gg;
    int i, k;

    for (i = 0; i < nn; ++i) {
        const double r = rho[i];
        const int r0 = i1[i] - 1;
        const int c0 = i2[i] - 1;
        const int r1 = r0 + 1;
        const int c1 = c0 + 1;
        const double x = a1[i];
        const double y = a2[i];
        const double w00 = (1.0 - x) * (1.0 - y);
        const double w10 = x * (1.0 - y);
        const double w01 = (1.0 - x) * y;
        const double w11 = x * y;
        double rp = r;
        double acc = 0.0;

        for (k = 0; k < dd; ++k) {
            const int base = k * plane;
            const double c00 = gram[base + r0 + c0 * gg];
            const double c10 = gram[base + r1 + c0 * gg];
            const double c01 = gram[base + r0 + c1 * gg];
            const double c11 = gram[base + r1 + c1 * gg];
            const double ck = w00 * c00 + w10 * c10 +
                              w01 * c01 + w11 * c11;
            acc += ck * rp;
            rp *= r;
        }

        if (fabs(r) < 1e-15)
            acc = 0.0;
        out[i] = acc;
    }
}


/*
 * Return the lower interpolation cell for a monotone grid, matching the
 * clamped findInterval(..., all.inside = TRUE) behavior used by the R
 * reference helpers.  The returned index is zero based and always lies in
 * [0, ng - 2].
 */
static int GeoLowerInterval(const double *grid, const int ng, const double x)
{
    int lo, hi, mid;

    if (x <= grid[0])
        return 0;
    if (x >= grid[ng - 1])
        return ng - 2;

    lo = 0;
    hi = ng - 1;
    while (hi - lo > 1) {
        mid = lo + (hi - lo) / 2;
        if (grid[mid] <= x)
            lo = mid;
        else
            hi = mid;
    }
    if (lo > ng - 2)
        lo = ng - 2;
    return lo;
}

static double GeoClamp(const double x, const double lo, const double hi)
{
    if (x < lo)
        return lo;
    if (x > hi)
        return hi;
    return x;
}


/*
 * Clayton warm-cache one-dimensional lookup.  For multiplicative margins the
 * exp(eta1 + eta2) scaling is fused into the same pass, so the R caller does
 * not allocate a million-element scaling vector before interpolation.
 */
void GeoClaytonInterp1(double *grid, int *ng,
                       double *values, double *x,
                       double *eta1, double *eta2,
                       int *n, int *multiplicative,
                       double *out)
{
    const int gg = *ng;
    const int nn = *n;
    const int mult = *multiplicative;
    int i;

    for (i = 0; i < nn; ++i) {
        const double xx = GeoClamp(x[i], grid[0], grid[gg - 1]);
        const int j = GeoLowerInterval(grid, gg, xx);
        const double d = grid[j + 1] - grid[j];
        const double a = d > 0.0 ? (xx - grid[j]) / d : 0.0;
        double z = (1.0 - a) * values[j] + a * values[j + 1];

        if (mult)
            z *= exp(eta1[i] + eta2[i]);
        out[i] = z;
    }
}


/*
 * Clayton location-dependent trilinear lookup for Beta2/Kumaraswamy2.
 * The table is an ne x ne x nx R array in column-major order.  Cell searches,
 * interpolation weights, the eight table reads, and the interpolation itself
 * are fused in one native pass over covariance pairs.
 */
void GeoClaytonInterp3(double *eta_grid, int *ne,
                       double *x_grid, int *nx,
                       double *tab,
                       double *eta1, double *eta2, double *x,
                       int *n, double *out)
{
    const int ee = *ne;
    const int rr = *nx;
    const int nn = *n;
    const int plane = ee * ee;
    int i;

    for (i = 0; i < nn; ++i) {
        const double z1 = GeoClamp(eta1[i], eta_grid[0], eta_grid[ee - 1]);
        const double z2 = GeoClamp(eta2[i], eta_grid[0], eta_grid[ee - 1]);
        const double xx = GeoClamp(x[i], x_grid[0], x_grid[rr - 1]);
        const int i1 = GeoLowerInterval(eta_grid, ee, z1);
        const int i2 = GeoLowerInterval(eta_grid, ee, z2);
        const int ir = GeoLowerInterval(x_grid, rr, xx);
        const double d1 = eta_grid[i1 + 1] - eta_grid[i1];
        const double d2 = eta_grid[i2 + 1] - eta_grid[i2];
        const double dr = x_grid[ir + 1] - x_grid[ir];
        const double a1 = d1 > 0.0 ? (z1 - eta_grid[i1]) / d1 : 0.0;
        const double a2 = d2 > 0.0 ? (z2 - eta_grid[i2]) / d2 : 0.0;
        const double ar = dr > 0.0 ? (xx - x_grid[ir]) / dr : 0.0;
        const int b0 = ir * plane;
        const int b1 = (ir + 1) * plane;
        const double c000 = tab[b0 + i1 + i2 * ee];
        const double c100 = tab[b0 + (i1 + 1) + i2 * ee];
        const double c010 = tab[b0 + i1 + (i2 + 1) * ee];
        const double c110 = tab[b0 + (i1 + 1) + (i2 + 1) * ee];
        const double c001 = tab[b1 + i1 + i2 * ee];
        const double c101 = tab[b1 + (i1 + 1) + i2 * ee];
        const double c011 = tab[b1 + i1 + (i2 + 1) * ee];
        const double c111 = tab[b1 + (i1 + 1) + (i2 + 1) * ee];
        const double c00 = (1.0 - a1) * c000 + a1 * c100;
        const double c10 = (1.0 - a1) * c010 + a1 * c110;
        const double c01 = (1.0 - a1) * c001 + a1 * c101;
        const double c11 = (1.0 - a1) * c011 + a1 * c111;
        const double c0 = (1.0 - a2) * c00 + a2 * c10;
        const double c1 = (1.0 - a2) * c01 + a2 * c11;

        out[i] = (1.0 - ar) * c0 + ar * c1;
    }
}
