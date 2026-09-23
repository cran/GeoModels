/*
 * Thread-safe modified Bessel K kernel for GeoModels.
 *
 * Numerical core adapted from R Mathlib src/nmath/bessel_k.c:
 * Copyright (C) 1998-2014 Ross Ihaka and the R Core team
 * Copyright (C) 2002-3 The R Foundation
 * R Mathlib is distributed under the GNU General Public License, version 2
 * or later.  GeoModels is distributed under GPL-3-or-later.
 *
 * The recurrence is based on the Campbell/Cody/Stoltz RKbesl algorithm.
 * This adaptation removes all R API calls, allocation and warning side
 * effects.  It is used only with exponential scaling (exp(x) K_nu(x)) by
 * the OpenMP Matérn likelihood kernel.
 */

#include "matern_bessel.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>

#define GEO_MIN0(x,y) (((x) <= (y)) ? (x) : (y))
#define GEO_MAX0(x,y) (((x) <= (y)) ? (y) : (x))
#define GEO_SQXMIN_BESS_K 1.49e-154
#define GEO_SQRT_2DPI 0.79788456080286535587989211986876373695

static void geo_K_bessel(double *x, double *alpha, int *nb,
                         int *ize, double *bk, int *ncalc)
{
    static const double a = .11593151565841244881;
    static const double p[8] = {
        .805629875690432845,20.4045500205365151,
        157.705605106676174,536.671116469207504,
        900.382759291288778,730.923886650660393,
        229.299301509425145,.822467033424113231
    };
    static const double q[7] = {
        29.4601986247850434,277.577868510221208,
        1206.70325591027438,2762.91444159791519,
        3443.74050506564618,2210.63190113378647,
        572.267338359892221
    };
    static const double r[5] = {
        -.48672575865218401848,13.079485869097804016,
        -101.96490580880537526,347.65409106507813131,
        3.495898124521934782e-4
    };
    static const double s[4] = {
        -25.579105509976461286,212.57260432226544008,
        -610.69018684944109624,422.69668805777760407
    };
    static const double t[6] = {
        1.6125990452916363814e-10,2.5051878502858255354e-8,
        2.7557319615147964774e-6,1.9841269840928373686e-4,
        .0083333333333334751799,.16666666666666666446
    };
    static const double estm[6] = {
        52.0583,5.7607,2.7782,14.4303,185.3004,9.3715
    };
    static const double estf[7] = {
        41.8341,7.1075,6.4306,42.511,1.35633,84.5096,20.
    };

    int iend, i, j, k, m, ii, mplus1;
    double x2by4, twox, c, blpha, ratio, wminf;
    double d1, d2, d3, f0, f1, f2, p0, q0, t1, t2, twonu;
    double dm, ex, bk1, bk2, nu;

    ii = 0;
    ex = *x;
    nu = *alpha;
    *ncalc = GEO_MIN0(*nb, 0) - 2;

    if (*nb > 0 && (0. <= nu && nu < 1.) && (1 <= *ize && *ize <= 2)) {
        if (ex <= 0) {
            for (i = 0; i < *nb; ++i) bk[i] = INFINITY;
            *ncalc = *nb;
            return;
        }

        k = 0;
        if (nu < GEO_SQXMIN_BESS_K) {
            nu = 0.;
        } else if (nu > .5) {
            k = 1;
            nu -= 1.;
        }
        twonu = nu + nu;
        iend = *nb + k - 1;
        c = nu * nu;
        d3 = -c;

        if (ex <= 1.) {
            d1 = 0.; d2 = p[0];
            t1 = 1.; t2 = q[0];
            for (i = 2; i <= 7; i += 2) {
                d1 = c * d1 + p[i - 1];
                d2 = c * d2 + p[i];
                t1 = c * t1 + q[i - 1];
                t2 = c * t2 + q[i];
            }
            d1 = nu * d1;
            t1 = nu * t1;
            f1 = log(ex);
            f0 = a + nu * (p[7] - nu * (d1 + d2) / (t1 + t2)) - f1;
            q0 = exp(-nu * (a - nu * (p[7] + nu * (d1 - d2) / (t1 - t2)) - f1));
            f1 = nu * f0;
            p0 = exp(f1);

            d1 = r[4];
            t1 = 1.;
            for (i = 0; i < 4; ++i) {
                d1 = c * d1 + r[i];
                t1 = c * t1 + s[i];
            }
            if (fabs(f1) <= .5) {
                f1 *= f1;
                d2 = 0.;
                for (i = 0; i < 6; ++i) d2 = f1 * d2 + t[i];
                d2 = f0 + f0 * f1 * d2;
            } else {
                d2 = sinh(f1) / nu;
            }
            f0 = d2 - nu * d1 / (t1 * p0);

            if (ex <= 1e-10) {
                bk[0] = f0 + ex * f0;
                if (*ize == 1) bk[0] -= ex * bk[0];
                ratio = p0 / f0;
                c = ex * DBL_MAX;
                if (k != 0) {
                    *ncalc = -1;
                    if (bk[0] >= c / ratio) return;
                    bk[0] = ratio * bk[0] / ex;
                    twonu += 2.;
                    ratio = twonu;
                }
                *ncalc = 1;
                if (*nb == 1) return;
                *ncalc = -1;
                for (i = 1; i < *nb; ++i) {
                    if (ratio >= c) return;
                    bk[i] = ratio / ex;
                    twonu += 2.;
                    ratio = twonu;
                }
                *ncalc = 1;
                goto L420;
            } else {
                c = 1.;
                x2by4 = ex * ex / 4.;
                p0 = .5 * p0;
                q0 = .5 * q0;
                d1 = -1.;
                d2 = 0.;
                bk1 = 0.;
                bk2 = 0.;
                f1 = f0;
                f2 = p0;
                do {
                    d1 += 2.;
                    d2 += 1.;
                    d3 = d1 + d3;
                    c = x2by4 * c / d2;
                    f0 = (d2 * f0 + p0 + q0) / d3;
                    p0 /= d2 - nu;
                    q0 /= d2 + nu;
                    t1 = c * f0;
                    t2 = c * (p0 - d2 * f0);
                    bk1 += t1;
                    bk2 += t2;
                } while (fabs(t1 / (f1 + bk1)) > DBL_EPSILON ||
                         fabs(t2 / (f2 + bk2)) > DBL_EPSILON);
                bk1 = f1 + bk1;
                bk2 = 2. * (f2 + bk2) / ex;
                if (*ize == 2) {
                    d1 = exp(ex);
                    bk1 *= d1;
                    bk2 *= d1;
                }
                wminf = estf[0] * ex + estf[1];
            }
        } else if (DBL_EPSILON * ex > 1.) {
            *ncalc = *nb;
            bk1 = 1. / (GEO_SQRT_2DPI * sqrt(ex));
            for (i = 0; i < *nb; ++i) bk[i] = bk1;
            return;
        } else {
            twox = ex + ex;
            blpha = 0.;
            ratio = 0.;
            if (ex <= 4.) {
                d2 = trunc(estm[0] / ex + estm[1]);
                m = (int) d2;
                d1 = d2 + d2;
                d2 -= .5;
                d2 *= d2;
                for (i = 2; i <= m; ++i) {
                    d1 -= 2.;
                    d2 -= d1;
                    ratio = (d3 + d2) / (twox + d1 - ratio);
                }
                d2 = trunc(estm[2] * ex + estm[3]);
                m = (int) d2;
                c = fabs(nu);
                d3 = c + c;
                d1 = d3 - 1.;
                f1 = DBL_MIN;
                f0 = (2. * (c + d2) / ex + .5 * ex / (c + d2 + 1.)) * DBL_MIN;
                for (i = 3; i <= m; ++i) {
                    d2 -= 1.;
                    f2 = (d3 + d2 + d2) * f0;
                    blpha = (1. + d1 / d2) * (f2 + blpha);
                    f2 = f2 / ex + f1;
                    f1 = f0;
                    f0 = f2;
                }
                f1 = (d3 + 2.) * f0 / ex + f1;
                d1 = 0.;
                t1 = 1.;
                for (i = 1; i <= 7; ++i) {
                    d1 = c * d1 + p[i - 1];
                    t1 = c * t1 + q[i - 1];
                }
                p0 = exp(c * (a + c * (p[7] - c * d1 / t1) - log(ex))) / ex;
                f2 = (c + .5 - ratio) * f1 / ex;
                bk1 = p0 + (d3 * f0 - f2 + f0 + blpha) /
                            (f2 + f1 + f0) * p0;
                if (*ize == 1) bk1 *= exp(-ex);
                wminf = estf[2] * ex + estf[3];
            } else {
                dm = trunc(estm[4] / ex + estm[5]);
                m = (int) dm;
                d2 = dm - .5;
                d2 *= d2;
                d1 = dm + dm;
                for (i = 2; i <= m; ++i) {
                    dm -= 1.;
                    d1 -= 2.;
                    d2 -= d1;
                    ratio = (d3 + d2) / (twox + d1 - ratio);
                    blpha = (ratio + ratio * blpha) / dm;
                }
                bk1 = 1. / ((GEO_SQRT_2DPI + GEO_SQRT_2DPI * blpha) * sqrt(ex));
                if (*ize == 1) bk1 *= exp(-ex);
                wminf = estf[4] * (ex - fabs(ex - estf[6])) + estf[5];
            }
            bk2 = bk1 + bk1 * (nu + .5 - ratio) / ex;
        }

        *ncalc = *nb;
        bk[0] = bk1;
        if (iend == 0) return;
        j = 1 - k;
        if (j >= 0) bk[j] = bk2;
        if (iend == 1) return;

        m = GEO_MIN0((int) (wminf - nu), iend);
        for (i = 2; i <= m; ++i) {
            t1 = bk1;
            bk1 = bk2;
            twonu += 2.;
            if (ex < 1.) {
                if (bk1 >= DBL_MAX / twonu * ex) break;
            } else {
                if (bk1 / ex >= DBL_MAX / twonu) break;
            }
            bk2 = twonu / ex * bk1 + t1;
            ii = i;
            ++j;
            if (j >= 0) bk[j] = bk2;
        }
        m = ii;
        if (m == iend) return;
        ratio = bk2 / bk1;
        mplus1 = m + 1;
        *ncalc = -1;
        for (i = mplus1; i <= iend; ++i) {
            twonu += 2.;
            ratio = twonu / ex + 1. / ratio;
            ++j;
            if (j >= 1) {
                bk[j] = ratio;
            } else {
                if (bk2 >= DBL_MAX / ratio) return;
                bk2 *= ratio;
            }
        }
        *ncalc = GEO_MAX0(1, mplus1 - k);
        if (*ncalc == 1) bk[0] = bk2;
        if (*nb == 1) return;
L420:
        for (i = *ncalc; i < *nb; ++i) {
            if (bk[i - 1] >= DBL_MAX / bk[i]) return;
            bk[i] *= bk[i - 1];
            (*ncalc)++;
        }
    }
}

int GeoModels_bessel_k_scaled(double x, double nu, double *work, double *value)
{
    int nb, ncalc, ize = 2;
    double alpha;

    if (value == NULL || work == NULL || !isfinite(x) || !isfinite(nu) ||
        x <= 0.0 || nu < 0.0 || nu > (double) (INT_MAX - 2))
        return 0;

    nb = 1 + (int) floor(nu);
    if (nb < 1) return 0;
    alpha = nu - (double) (nb - 1);
    geo_K_bessel(&x, &alpha, &nb, &ize, work, &ncalc);
    if (ncalc != nb) return 0;
    if (!(work[nb - 1] > 0.0) || !isfinite(work[nb - 1])) return 0;

    *value = work[nb - 1];
    return 1;
}
