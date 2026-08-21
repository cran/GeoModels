#include "header.h"

/* Native bivariate densities implemented in Distributions.c. */
extern double log_biv_Norm(double corr, double zi, double zj,
                           double mi, double mj, double vari, double nugget);
extern double biv_skew(double corr, double z1, double z2,
                       double mi, double mj, double vari,
                       double skew, double nugget);
extern double biv_tukey_h(double corr, double data_i, double data_j,
                          double mean_i, double mean_j,
                          double tail, double sill);
extern double biv_T(double rho, double zi, double zj,
                    double nuu, double nugget);
extern double biv_sinh(double corr, double zi, double zj,
                       double mi, double mj, double skew,
                       double tail, double vari);
extern double d2lognorm(double x, double y, double sill, double nugget,
                        double mux, double muy, double rho);
extern double biv_gamma(double corr, double zi, double zj,
                        double mui, double muj, double shape);
extern double biv_Weibull(double corr, double zi, double zj,
                          double mui, double muj, double shape);

/* Copula densities implemented in Distributions.c. */
extern double biv_unif_CopulaGauss(double dat1, double dat2, double rho);
extern double biv_unif_CopulaSkewGauss(double dat1, double dat2,
                                       double rho, double alpha);
extern double biv_unif_CopulaClayton(double dat1, double dat2,
                                     double rho, double nu);

static R_INLINE double geo_density_clamp_corr(double rho)
{
    const double eps = 1e-10;
    if (!R_FINITE(rho)) return NA_REAL;
    if (rho >= 1.0) return 1.0 - eps;
    if (rho <= -1.0) return -1.0 + eps;
    return rho;
}

static R_INLINE double geo_density_clamp_prob(double u)
{
    const double eps = 1e-10;
    if (!R_FINITE(u)) return NA_REAL;
    if (u <= eps) return eps;
    if (u >= 1.0 - eps) return 1.0 - eps;
    return u;
}

static R_INLINE double geo_copula_density(int type_cop,
                                           double u1, double u2,
                                           double rho, double copula_par)
{
    double out = 0.0;

    if (type_cop == 1) {
        out = biv_unif_CopulaGauss(u1, u2, rho);
    } else if (type_cop == 2) {
        if (copula_par <= 0.0) return 0.0;
        const double logc = biv_unif_CopulaClayton(
            u1, u2, rho, copula_par
        );
        out = R_FINITE(logc) ? exp(logc) : 0.0;
    } else if (type_cop == 3) {
        if (fabs(copula_par) >= 1.0) return 0.0;
        const double alpha = copula_par / sqrt(1.0 - copula_par * copula_par);
        out = biv_unif_CopulaSkewGauss(u1, u2, rho, alpha);
    }

    return (R_FINITE(out) && out >= 0.0) ? out : 0.0;
}


/*
 * Direct vectorized evaluator of a fitted copula density. This routine is
 * used by GeoScatterplot() on the Uniform and Gaussian-score scales, where
 * the response marginal distribution must not enter the theoretical contour.
 */
void GeoCopulaDensity(int *type_cop,
                      double *copula_par,
                      double *u, double *v, int *n,
                      double *rho, double *res)
{
    /*
     * rho is evaluated from the fitted correlation model in R and passed
     * directly to this routine. This follows the contour construction used
     * in the skew-Gaussian copula paper and makes the value used by each
     * panel explicit and reproducible.
     */
    const double rho_use = geo_density_clamp_corr(*rho);

    if (!R_FINITE(rho_use)) {
        for (int i = 0; i < *n; ++i) res[i] = NA_REAL;
        return;
    }

    for (int i = 0; i < *n; ++i) {
        if (!R_FINITE(u[i]) || !R_FINITE(v[i])) {
            res[i] = NA_REAL;
            continue;
        }

        const double ui = geo_density_clamp_prob(u[i]);
        const double vi = geo_density_clamp_prob(v[i]);
        res[i] = geo_copula_density(
            *type_cop, ui, vi, rho_use, *copula_par
        );
    }
}

/*
 * Vectorized bivariate density evaluator for GeoScatterplot().
 *
 * type_cop:
 *   0 = native random-field bivariate density
 *   1 = Gaussian copula
 *   2 = Clayton copula
 *   3 = Skew-Gaussian copula
 *
 * nuis contains only the marginal nuisance parameters, in the order
 * returned by NuisParam2(model, FALSE, num_betas, copula = NULL), after
 * removing the mean parameters. The copula parameter is passed separately.
 */
void GeoBivDensity(int *model, int *cormod, int *type_cop,
                   double *copula_par,
                   double *x, double *y, int *n,
                   double *lag, double *lagt, double *par,
                   double *mean1, double *mean2,
                   double *nuis, int *NN1, int *NN2,
                   double *res)
{
    double rho0 = CorFct(cormod, *lag, *lagt, par, 0, 0);
    rho0 = geo_density_clamp_corr(rho0);

    if (!R_FINITE(rho0)) {
        for (int i = 0; i < *n; ++i) res[i] = NA_REAL;
        return;
    }

    for (int i = 0; i < *n; ++i) {
        const double xi = x[i];
        const double yi = y[i];
        double dens = 0.0;

        if (!R_FINITE(xi) || !R_FINITE(yi)) {
            res[i] = NA_REAL;
            continue;
        }

        if (*type_cop > 0) {
            double nugget = 0.0;
            double u1 = NA_REAL, u2 = NA_REAL;
            double g1 = 0.0, g2 = 0.0;

            switch (*model) {
            case 1: { /* Gaussian: nugget, sill */
                nugget = nuis[0];
                const double sill = nuis[1];
                if (sill > 0.0) {
                    const double sd = sqrt(sill);
                    const double z1 = (xi - *mean1) / sd;
                    const double z2 = (yi - *mean2) / sd;
                    u1 = pnorm(z1, 0.0, 1.0, 1, 0);
                    u2 = pnorm(z2, 0.0, 1.0, 1, 0);
                    g1 = dnorm(z1, 0.0, 1.0, 0) / sd;
                    g2 = dnorm(z2, 0.0, 1.0, 0) / sd;
                }
                break;
            }

            case 12: { /* StudentT: df, nugget, sill */
                const double df_par = nuis[0];
                nugget = nuis[1];
                const double sill = nuis[2];
                if (df_par > 0.0 && sill > 0.0) {
                    const double df = 1.0 / df_par;
                    const double sd = sqrt(sill);
                    const double z1 = (xi - *mean1) / sd;
                    const double z2 = (yi - *mean2) / sd;
                    u1 = pt(z1, df, 1, 0);
                    u2 = pt(z2, df, 1, 0);
                    g1 = dt(z1, df, 0) / sd;
                    g2 = dt(z2, df, 0) / sd;
                }
                break;
            }

            case 21: { /* Gamma: nugget, sill, shape */
                nugget = nuis[0];
                const double shape = nuis[2];
                if (xi > 0.0 && yi > 0.0 && shape > 0.0) {
                    const double a = shape / 2.0;
                    const double scale1 = 2.0 * exp(*mean1) / shape;
                    const double scale2 = 2.0 * exp(*mean2) / shape;
                    u1 = pgamma(xi, a, scale1, 1, 0);
                    u2 = pgamma(yi, a, scale2, 1, 0);
                    g1 = dgamma(xi, a, scale1, 0);
                    g2 = dgamma(yi, a, scale2, 0);
                }
                break;
            }

            case 22: { /* LogGaussian: nugget, sill */
                nugget = nuis[0];
                const double sill = nuis[1];
                if (xi > 0.0 && yi > 0.0 && sill > 0.0) {
                    const double sd = sqrt(sill);
                    const double lm1 = *mean1 - sill / 2.0;
                    const double lm2 = *mean2 - sill / 2.0;
                    u1 = plnorm(xi, lm1, sd, 1, 0);
                    u2 = plnorm(yi, lm2, sd, 1, 0);
                    g1 = dlnorm(xi, lm1, sd, 0);
                    g2 = dlnorm(yi, lm2, sd, 0);
                }
                break;
            }

            case 26: { /* Weibull: nugget, sill, shape */
                nugget = nuis[0];
                const double shape = nuis[2];
                if (xi > 0.0 && yi > 0.0 && shape > 0.0) {
                    const double den = gammafn(1.0 + 1.0 / shape);
                    const double scale1 = exp(*mean1) / den;
                    const double scale2 = exp(*mean2) / den;
                    u1 = pweibull(xi, shape, scale1, 1, 0);
                    u2 = pweibull(yi, shape, scale2, 1, 0);
                    g1 = dweibull(xi, shape, scale1, 0);
                    g2 = dweibull(yi, shape, scale2, 0);
                }
                break;
            }

            case 28: { /* Beta: nugget, sill, shape1, shape2, min, max */
                nugget = nuis[0];
                const double shape1 = nuis[2];
                const double shape2 = nuis[3];
                const double lower = nuis[4];
                const double upper = nuis[5];
                const double width = upper - lower;
                if (width > 0.0 && shape1 > 0.0 && shape2 > 0.0 &&
                    xi > lower && xi < upper && yi > lower && yi < upper) {
                    const double b1 = (xi - lower) / width;
                    const double b2 = (yi - lower) / width;
                    u1 = pbeta(b1, shape1, shape2, 1, 0);
                    u2 = pbeta(b2, shape1, shape2, 1, 0);
                    g1 = dbeta(b1, shape1, shape2, 0) / width;
                    g2 = dbeta(b2, shape1, shape2, 0) / width;
                }
                break;
            }

            case 50: { /* Beta2: nugget, sill, shape, min, max */
                nugget = nuis[0];
                const double shape = nuis[2];
                const double lower = nuis[3];
                const double upper = nuis[4];
                const double width = upper - lower;
                const double mu1 = 1.0 / (1.0 + exp(-(*mean1)));
                const double mu2 = 1.0 / (1.0 + exp(-(*mean2)));
                if (width > 0.0 && shape > 0.0 &&
                    xi > lower && xi < upper && yi > lower && yi < upper) {
                    const double b1 = (xi - lower) / width;
                    const double b2 = (yi - lower) / width;
                    const double a11 = shape * mu1;
                    const double a12 = shape * (1.0 - mu1);
                    const double a21 = shape * mu2;
                    const double a22 = shape * (1.0 - mu2);
                    u1 = pbeta(b1, a11, a12, 1, 0);
                    u2 = pbeta(b2, a21, a22, 1, 0);
                    g1 = dbeta(b1, a11, a12, 0) / width;
                    g2 = dbeta(b2, a21, a22, 0) / width;
                }
                break;
            }

            default:
                dens = NA_REAL;
            }

            if (R_FINITE(u1) && R_FINITE(u2) &&
                R_FINITE(g1) && R_FINITE(g2) && g1 >= 0.0 && g2 >= 0.0) {
                u1 = geo_density_clamp_prob(u1);
                u2 = geo_density_clamp_prob(u2);
                const double rho = geo_density_clamp_corr((1.0 - nugget) * rho0);
                if (R_FINITE(rho)) {
                    dens = geo_copula_density(
                        *type_cop, u1, u2, rho, *copula_par
                    ) * g1 * g2;
                }
            }

        } else {
            switch (*model) {
            case 1: { /* Gaussian: nugget, sill */
                const double nugget = nuis[0];
                const double sill = nuis[1];
                const double rho = geo_density_clamp_corr((1.0 - nugget) * rho0);
                if (sill > 0.0 && R_FINITE(rho)) {
                    dens = exp(log_biv_Norm(
                        rho, xi, yi, *mean1, *mean2, sill, 0.0
                    ));
                }
                break;
            }

            case 10: { /* SkewGaussian: nugget, sill, skew */
                const double nugget = nuis[0];
                const double sill = nuis[1];
                const double skew = nuis[2];
                if (sill > 0.0) {
                    dens = biv_skew(
                        rho0, xi, yi, *mean1, *mean2,
                        sill, skew, nugget
                    );
                }
                break;
            }

            case 12: { /* StudentT: df, nugget, sill */
                const double df = nuis[0];
                const double nugget = nuis[1];
                const double sill = nuis[2];
                if (df > 0.0 && sill > 0.0) {
                    const double sd = sqrt(sill);
                    dens = biv_T(
                        rho0,
                        (xi - *mean1) / sd,
                        (yi - *mean2) / sd,
                        df, nugget
                    ) / sill;
                }
                break;
            }

            case 20: { /* SinhAsinh: nugget, sill, skew, tail */
                const double nugget = nuis[0];
                const double sill = nuis[1];
                const double skew = nuis[2];
                const double tail = nuis[3];
                const double rho = geo_density_clamp_corr((1.0 - nugget) * rho0);
                if (sill > 0.0 && tail > 0.0 && R_FINITE(rho)) {
                    dens = biv_sinh(
                        rho, xi, yi, *mean1, *mean2,
                        skew, tail, sill
                    );
                }
                break;
            }

            case 21: { /* Gamma: nugget, sill, shape */
                const double nugget = nuis[0];
                const double shape = nuis[2];
                const double rho = geo_density_clamp_corr((1.0 - nugget) * rho0);
                if (xi > 0.0 && yi > 0.0 && shape > 0.0 && R_FINITE(rho)) {
                    const double logdens = biv_gamma(
                        rho, xi, yi, *mean1, *mean2, shape
                    );
                    dens = R_FINITE(logdens) ? exp(logdens) : 0.0;
                }
                break;
            }

            case 22: { /* LogGaussian: nugget, sill */
                const double nugget = nuis[0];
                const double sill = nuis[1];
                const double rho = geo_density_clamp_corr((1.0 - nugget) * rho0);
                if (xi > 0.0 && yi > 0.0 && sill > 0.0 && R_FINITE(rho)) {
                    dens = d2lognorm(
                        xi, yi, sill, nugget,
                        *mean1, *mean2, rho
                    );
                }
                break;
            }

            case 26: { /* Weibull: nugget, sill, shape */
                const double nugget = nuis[0];
                const double shape = nuis[2];
                const double rho = geo_density_clamp_corr((1.0 - nugget) * rho0);
                if (xi > 0.0 && yi > 0.0 && shape > 0.0 && R_FINITE(rho)) {
                    const double logdens = biv_Weibull(
                        rho, xi, yi, *mean1, *mean2, shape
                    );
                    dens = R_FINITE(logdens) ? exp(logdens) : 0.0;
                }
                break;
            }

            case 34: { /* Tukeyh: nugget, sill, tail */
                const double nugget = nuis[0];
                const double sill = nuis[1];
                const double tail = nuis[2];
                const double rho = geo_density_clamp_corr((1.0 - nugget) * rho0);
                if (sill > 0.0 && tail >= 0.0 && tail < 0.5 && R_FINITE(rho)) {
                    dens = biv_tukey_h(
                        rho, xi, yi, *mean1, *mean2, tail, sill
                    );
                }
                break;
            }

            default:
                dens = NA_REAL;
            }
        }

        if (!R_FINITE(dens) || dens < 0.0) dens = 0.0;
        res[i] = dens;
    }
}
