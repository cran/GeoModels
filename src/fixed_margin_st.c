/* Fixed-margin, nonseparable space-time correlations.
 *
 * 85: par = (scale_s, scale_t, power2_s, power2_t,
 *            smooth_s, smooth_t, sep).
 * 89: par = (scale_s, scale_t, smooth_s, smooth_t, sep).
 * sep is the normalized interaction rho in [-1,1], NOT eta or zeta.
 * The spatial dimension is set/restored at R call boundaries; it is read-only
 * during native evaluations.  Lag-only public calls default to dimension 2.
 */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "fixed_margin_st.h"

extern double CorFunWitMat(double lag, double scale, double nu);

static int geo_fixed_margin_st_dimension = 2;

int GeoFixedMarginSTDimension(void)
{
    return geo_fixed_margin_st_dimension;
}

SEXP GeoFixedMarginSTGetDimension(void)
{
    return ScalarInteger(geo_fixed_margin_st_dimension);
}

SEXP GeoFixedMarginSTSetDimension(SEXP dimension)
{
    if (XLENGTH(dimension) != 1)
        error("GeoModels: internal spatial dimension must be a scalar.");
    const double value = asReal(dimension);
    if (!R_FINITE(value) || value < 1.0 || value > 3.0 || value != floor(value))
        error("GeoModels: internal spatial dimension must be an integer in 1:3.");
    const int old = geo_fixed_margin_st_dimension;
    geo_fixed_margin_st_dimension = (int)value;
    return ScalarInteger(old);
}

int GeoFixedMarginSTValid(int code, const double *par, int dimension)
{
    const int n = code == 85 ? 7 : (code == 89 ? 5 : 0);
    if (!n || !par || dimension < 1 || dimension > 3) return 0;
    for (int i = 0; i < n; ++i) if (!R_FINITE(par[i])) return 0;
    if (par[0] <= 0.0 || par[1] <= 0.0 || fabs(par[n-1]) > 1.0) return 0;
    if (code == 89) return par[2] > 0.0 && par[3] > 0.0;
    /* Strong sufficient shape region of the sharp interaction theorem.
       The temporal dimension is always one.  Both inequalities are strict. */
    return par[4] > -0.5 && par[5] > -0.5 &&
           par[2] > 0.5 * (dimension + 1.0) + par[4] + 2.0 &&
           par[3] > par[5] + 3.0;
}

static double geo_st_interaction(double sep, double A, double B)
{
    /* A,B are positive.  In particular, the endpoints are valid. */
    if (sep == 0.0) return 0.0;
    return (2.0 * sep) / (A * (1.0 - sep) + B * (1.0 + sep));
}

/* Normalized GW evaluation for the fixed-margin family only.  The existing
   general-purpose 2F1 evaluator loses accuracy near integer nu=kappa+1/2 and
   when its continuation formula subtracts large terms.  Here the near-origin
   continuation is normalized analytically and the other branch is a positive
   series.  xd is x*dW/dx, allowing stable consecutive-shape differences. */
static int geo_st_series(long double a, long double b, long double c,
                         long double t, long double *v, long double *td)
{
    long double term = 1.0L, sum = 1.0L, deriv = 0.0L;
    for (int n = 1; n <= 20000; ++n) {
        const long double den = (c+n-1.0L)*n;
        if (den == 0.0L) return 0;
        term *= (a+n-1.0L)*(b+n-1.0L)*t/den;
        sum += term; deriv += n*term;
        if (!isfinite(sum) || !isfinite(deriv)) return 0;
        if (fabsl(term)*(n+1) <= 4.0L*LDBL_EPSILON*fmaxl(1.0L,fabsl(sum))) {
            *v=sum; *td=deriv; return 1;
        }
    }
    return 0;
}

/* Only shape-dependent constants are cached; x-dependent series remain exact.
 * Four slots cover both margins and their consecutive-shape evaluations. */
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
# define GEO_ST_TLS _Thread_local
#else
# define GEO_ST_TLS
#endif
typedef struct {
    int valid, have_normalizer, have_connection, have_integer;
    double mu, kappa;
    long double log_normalizer, log_connection, log_integer;
} GeoSTGWShape;
static GEO_ST_TLS GeoSTGWShape geo_st_gw_shapes[4];
static GEO_ST_TLS unsigned int geo_st_gw_next_shape;

static GeoSTGWShape *geo_st_gw_shape(double mu, double kappa)
{
    for (int i = 0; i < 4; ++i)
        if (geo_st_gw_shapes[i].valid && geo_st_gw_shapes[i].mu == mu &&
            geo_st_gw_shapes[i].kappa == kappa) return &geo_st_gw_shapes[i];
    GeoSTGWShape *ctx = &geo_st_gw_shapes[geo_st_gw_next_shape];
    geo_st_gw_next_shape = (geo_st_gw_next_shape + 1u) & 3u;
    ctx->valid = 1; ctx->mu = mu; ctx->kappa = kappa;
    ctx->have_normalizer = ctx->have_connection = ctx->have_integer = 0;
    return ctx;
}

static int geo_st_gw_eval(double x, double mu0, double k0,
                          double *value, double *xd)
{
    if (x == 0.0) { *value=1.0; *xd=0.0; return 1; }
    if (x >= 1.0) { *value=0.0; *xd=0.0; return 1; }
    const long double mu=mu0,k=k0,nu=k+0.5L;
    if (nu > (long double)INT_MAX-1) return 0;
    const long double t=(long double)x*x,z=1.0L-t,p=k+mu;
    const long double a=mu/2.0L,b=(mu+1.0L)/2.0L,c=k+mu+1.0L;
    GeoSTGWShape *ctx = geo_st_gw_shape(mu0, k0);
    long double sum, td;
    if (t < 0.25L && (mu+2.0L*k+1.0L)*x <= 3.0L) {
        const long double nearest=roundl(nu);
        if (nearest >= 1.0L && fabsl(nu-nearest) < 4.0L*DBL_EPSILON) {
            /* Logarithmic connection formula at integer nu. */
            const int m=(int)nearest;
            long double term=1.0L;sum=1.0L;td=0.0L;
            for(int n=1;n<m;++n) {
                term *= (a+n-1)*(b+n-1)*t/((long double)n*(n-m));
                sum += term; td += n*term;
            }
            long double f=0.0L,df=0.0L;term=1.0L;
            long double bracket=logl(t)-digamma(1.0)-digamma(m+1.0)+
                                 digamma((double)(a+m))+digamma((double)(b+m));
            for(int n=0;n<=20000;++n) {
                const long double add=term*bracket;
                f += add; df += n*add+term;
                if(n>2 && fabsl(term)*(fabsl(bracket)+n+1) <=
                   4.0L*LDBL_EPSILON*fmaxl(1.0L,fabsl(f))) break;
                if(n==20000) return 0;
                bracket += -1.0L/(n+1)-1.0L/(m+n+1)+1.0L/(a+m+n)+1.0L/(b+m+n);
                term *= (a+m+n)*(b+m+n)*t/((n+1.0L)*(m+n+1.0L));
            }
            if (!ctx->have_integer) {
                ctx->log_integer=lgammal(a+m)+lgammal(b+m)-lgammal(a)-lgammal(b)-
                                 lgammal(m)-lgammal(m+1);
                ctx->have_integer=1;
            }
            long double coef=expl(ctx->log_integer+m*logl(t));
            if(!(m&1)) coef=-coef;
            sum += coef*f;td += coef*(m*f+df);
        } else {
            long double f1,d1,f2,d2;
            if(!geo_st_series(a,b,1.0L-nu,t,&f1,&d1) ||
               !geo_st_series(c-a,c-b,1.0L+nu,t,&f2,&d2)) return 0;
            if (!ctx->have_connection) {
                ctx->log_connection=lgammal(c-a)+lgammal(c-b)+lgammal(-nu)-
                                    lgammal(nu)-lgammal(a)-lgammal(b);
                ctx->have_connection=1;
            }
            long double coef=expl(ctx->log_connection+nu*logl(t));
            /* Gamma(-nu) has alternating sign between its poles. */
            if (((long long)floorl(nu)) % 2 == 0) coef=-coef;
            sum=f1+coef*f2;td=d1+coef*(nu*f2+d2);
        }
        const long double pref=expl(p*log1pl(-t));
        *value=(double)(pref*sum);
        *xd=(double)(2.0L*pref*(td-p*t/z*sum));
    } else {
        /* a,b,c,z are positive: no subtraction of large analytic branches. */
        long double term=1.0L;sum=1.0L;td=0.0L;
        long double logscale=0.0L;
        int converged=0;
        for(int n=1;n<=2000000;++n) {
            term *= (a+n-1)*(b+n-1)*z/((c+n-1)*n);
            sum += term;td += n*term;
            if(sum > 1.0e150L) { sum*=1.0e-150L;td*=1.0e-150L;term*=1.0e-150L;
                                logscale+=150.0L*logl(10.0L); }
            if(term*(n+1) <= 4.0L*LDBL_EPSILON*sum) {converged=1;break;}
        }
        if(!converged || !isfinite(sum)) return 0;
        if (!ctx->have_normalizer) {
            ctx->log_normalizer=lgammal(c-a)+lgammal(c-b)-lgammal(nu)-lgammal(c);
            ctx->have_normalizer=1;
        }
        const long double w=expl(ctx->log_normalizer+p*logl(z)+logscale+logl(sum));
        *value=(double)w;
        *xd=(double)(-2.0L*t/z*w*(p+td/sum));
    }
    /* Correlations are nonnegative and decreasing; reject arithmetic failure.
       Only tiny positive roundoff at the origin is rounded to one. */
    if (!R_FINITE(*value) || !R_FINITE(*xd) || *value < -1e-12 ||
        *value > 1.0+1e-10 || *xd > 1e-10) return 0;
    if (*value < 0.0) *value=0.0;
    if (*value > 1.0) *value=1.0;
    if (*xd > 0.0) *xd=0.0;
    return 1;
}

static void geo_st_gw_parts(double h, double scale, double mu, double kappa,
                            double *w, double *difference)
{
    const double x = h / scale;
    if (h == 0.0) { *w = 1.0; *difference = 0.0; return; }
    if (x >= 1.0) { *w = 0.0; *difference = 0.0; return; }
    if (x == 0.0) { *w=1.0; *difference=0.0; return; }
    if (kappa == 0.0) {
        *w = exp(mu * log1p(-x));
        *difference = -x * (*w);
    } else if (kappa == 1.0) {
        const double p = exp((mu + 1.0) * log1p(-x));
        *w = p * (1.0 + (mu + 1.0) * x);
        *difference = -(mu + 2.0) * x * x * p;
    } else {
        double xd, p, pd;
        if (!geo_st_gw_eval(x,mu,kappa,w,&xd) ||
            !geo_st_gw_eval(x,mu+1.0,kappa,&p,&pd)) {
            *w=R_NaN;*difference=R_NaN;return;
        }
        /* Exact recurrence x P'(x)=(mu+2*kappa+1)(P-W), same radius. */
        *difference=pd/(mu+2.0*kappa+1.0);
    }
}

/* Separable evaluations need the margins but not the consecutive-shape
 * difference.  Keep the same robust GW evaluator for the actual values. */
static double geo_st_gw_value(double h, double scale, double mu, double kappa)
{
    const double x = h / scale;
    if (h == 0.0 || x == 0.0) return 1.0;
    if (x >= 1.0) return 0.0;
    if (kappa == 0.0) return exp(mu * log1p(-x));
    if (kappa == 1.0)
        return exp((mu+1.0)*log1p(-x))*(1.0+(mu+1.0)*x);
    double value, xd;
    return geo_st_gw_eval(x, mu, kappa, &value, &xd) ? value : R_NaN;
}

/* Q_nu(h) = h M_nu'(h).  A derivative identity avoids catastrophic
   cancellation in 2 nu (M_nu - M_{nu+1}) at short distances. */
static double geo_st_matern_q(double h, double scale, double nu)
{
    if (h == 0.0) return 0.0;
    const double x = h / scale;
    if (!R_FINITE(x)) return 0.0;
    const double log_x = x > 0.0 ? log(x) : log(h) - log(scale);
    if (nu > 1.0) {
        const double m = CorFunWitMat(h, scale, nu - 1.0);
        if (m == 0.0) return 0.0;
        if (!(m > 0.0) || !R_FINITE(m)) return R_NaN;
        return -exp(2.0 * log_x + log(m) - log(2.0 * (nu - 1.0)));
    }
    if (nu == 0.5) return x > 1.0e4 ? 0.0 : -exp(log_x - x);
    if (x <= 1.0e-10) {
        if (nu <= 0.9)
            return -exp((1.0-2.0*nu)*M_LN2 + lgammafn(1.0-nu) -
                        lgammafn(nu) + 2.0*nu*log_x);
        if (nu == 1.0)
            return -exp(2.0*log_x + log(M_LN2-log_x-0.57721566490153286061));
        const double delta=1.0-nu;
        const double aa=lgammafn(1.0+delta)-lgammafn(1.0-delta)+
                        2.0*delta*(M_LN2-log_x);
        const double logexpm1=aa>50.0 ? aa : log(expm1(aa));
        return -exp(2.0*log_x-log(2.0*delta)+logexpm1);
    }
    if (x == 0.0) {
        if (nu == 1.0) return 0.0;
        return -exp((1.0 - 2.0 * nu) * M_LN2 + lgammafn(1.0 - nu) -
                    lgammafn(nu) + 2.0 * nu * log_x);
    }
    const double bk = bessel_k(x, fabs(nu - 1.0), 2.0);
    if (!(bk > 0.0) || !R_FINITE(bk)) return R_NaN;
    return -exp((1.0 - nu) * M_LN2 - lgammafn(nu) +
                (nu + 1.0) * log_x + log(bk) - x);
}

double GeoFixedMarginSTCorr(int code, double h, double u,
                            const double *par, int dimension)
{
    if (!GeoFixedMarginSTValid(code, par, dimension) ||
        !R_FINITE(h) || !R_FINITE(u) || h < 0.0 || u < 0.0) return R_NaN;
    if (h == 0.0 && u == 0.0) return 1.0;
    if (code == 85) {
        if (h >= par[0] || u >= par[1]) return 0.0;
        if (par[6] == 0.0)
            return geo_st_gw_value(h, par[0], par[2], par[4]) *
                   geo_st_gw_value(u, par[1], par[3], par[5]);
        double ws, wt, ds, dt;
        geo_st_gw_parts(h, par[0], par[2], par[4], &ws, &ds);
        geo_st_gw_parts(u, par[1], par[3], par[5], &wt, &dt);
        if (h == 0.0) return wt;
        if (u == 0.0) return ws;
        const double Ls = dimension / (par[2] + 2.0*par[4] + dimension + 1.0);
        const double Lt = 1.0 / (par[3] + 2.0*par[5] + 2.0);
        const double Us = (2.0*par[4] + 1.0) / par[2];
        const double Ut = (2.0*par[5] + 1.0) / par[3];
        const double A = fmax(Ls * Lt, Us * Ut);
        const double B = fmax(Ls * Ut, Us * Lt);
        const double eta = geo_st_interaction(par[6], A, B);
        return ws * wt + eta * ds * dt;
    }
    const double ns = par[2], nt = par[3];
    const double ms = CorFunWitMat(h, par[0], ns);
    const double mt = CorFunWitMat(u, par[1], nt);
    if (h == 0.0) return mt;
    if (u == 0.0) return ms;
    if (par[4] == 0.0) return ms * mt;
    const double A = fmax((double)dimension, 4.0 * ns * nt);
    const double B = fmax(2.0 * dimension * nt, 2.0 * ns);
    const double zeta = geo_st_interaction(par[4], A, B);
    return ms * mt + zeta * geo_st_matern_q(h, par[0], ns) *
                              geo_st_matern_q(u, par[1], nt);
}

/* Registered internal vector evaluator for regression tests and diagnostics.
   It does not change the dimension context and does not apply a nugget. */
SEXP GeoFixedMarginSTEvaluate(SEXP model, SEXP hs, SEXP us, SEXP parameters,
                              SEXP dimension)
{
    if (XLENGTH(model) != 1 || XLENGTH(dimension) != 1)
        error("GeoModels: model and internal spatial dimension must be scalars.");
    const int code = asInteger(model);
    const double dd = asReal(dimension);
    const int n = code == 85 ? 7 : (code == 89 ? 5 : 0);
    if (!n || TYPEOF(hs) != REALSXP || TYPEOF(us) != REALSXP ||
        TYPEOF(parameters) != REALSXP || XLENGTH(hs) != XLENGTH(us) ||
        XLENGTH(parameters) != n || !R_FINITE(dd) || dd != floor(dd) ||
        dd < 1.0 || dd > 3.0)
        error("GeoModels: invalid fixed-margin space-time arguments.");
    if (!GeoFixedMarginSTValid(code, REAL(parameters), (int)dd))
        error("GeoModels: invalid fixed-margin space-time parameters.");
    const R_xlen_t count = XLENGTH(hs);
    SEXP out = PROTECT(allocVector(REALSXP, count));
    for (R_xlen_t i = 0; i < count; ++i) {
        if ((i & 8191) == 0) R_CheckUserInterrupt();
        REAL(out)[i] = GeoFixedMarginSTCorr(code, REAL(hs)[i], REAL(us)[i],
                                           REAL(parameters), (int)dd);
    }
    UNPROTECT(1);
    return out;
}
