#include "header.h"
/*************************************************/
//******************************************************************************
double polevl(double x, const double coef[], int N)

{
    double ans;
    int i;
    const double *p;
    p = coef;
    ans = *p++;
    i = N;
    do
        ans = ans * x + *p++;
    while (--i);   
    return (ans);
}
double p1evl(double x, const double coef[], int N)
{
    double ans;
    const double *p;
    int i;
    p = coef;
    ans = x + *p++;
    i = N - 1;
    do
        ans = ans * x + *p++;
    while (--i);
    return (ans);
}
/************************************************************************/
 int is_nonpos_int(double x)
{
    return x <= 0 && x == ceil(x) && fabs(x) < 1e13;
}
 
double poch(double a, double m)
{
    double r;
    r = 1.0;
    /* Recurse down */
    while (m >= 1.0) {
        if (a + m == 1) {
            break;
        }
        m -= 1.0;
        r *= (a + m);
        if (!R_FINITE(r) || r == 0) {
            break;
        }
    }

    /* Recurse up */
    while (m <= -1.0) {
        if (a + m == 0) {
            break;
        }
        r /= (a + m);
        m += 1.0;
        if (!R_FINITE(r) || r == 0) {
            break;
        }
    }

    if (m == 0) {
        /* Easy case */
        return r;
    }
    else if (a > 1e4 && fabs(m) <= 1) {
        /* Avoid loss of precision */
        return r * R_pow(a, m) * (
            1
            + m*(m-1)/(2*a)
            + m*(m-1)*(m-2)*(3*m-1)/(24*a*a)
            + m*m*(m-1)*(m-1)*(m-2)*(m-3)/(48*a*a*a)
            );
    }

    /* Check for infinity */
    if (is_nonpos_int(a + m) && !is_nonpos_int(a) && a + m != m) {
        return INFINITY;
    }

    /* Check for zero */
    if (!is_nonpos_int(a + m) && is_nonpos_int(a)) {
        return 0;
    }

    return(r * exp(lgammafn(a + m) - lgammafn(a)) * sign(gammafn(a + m)) * sign(gammafn(a)));
}
/************************************************************************/

// ===================================== START: Bivariate Normal  =====================================//



/*for bivariate t distributions*/
double A[] = {
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2
};

double B[] = {
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5
};

double C[] = {
    /* 1.00000000000000000000E0, */
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6
};


// Integrand function (derivatives of the Student-t cdf):
// bivariate of Gaussian bivariate rf
double log_biv2gauss(int *cormod, double dij,double *par, double data1, double data2, int first,int second)
{
double rhott,rhovv,rhotv,det,dens=0.0;
rhott=CorFct(cormod,0,0,par,first,first);
rhovv=CorFct(cormod,0,0,par,second,second);
rhotv=CorFct(cormod,dij,0,par,first,second);
det=rhott*rhovv-R_pow(rhotv,2);
dens=-0.5*(2*log(2*M_PI)+log(det)+(rhovv*R_pow(data1,2)+rhott*R_pow(data2,2)-2*(data1*data2)*rhotv)/det);
return dens;
}

// compute  bivariate normal standard pdf:
double d22norm(double x, double y,double v11,double v22,double v12)
{
  double res=0.0;
  double cc=sqrt(v11*v22);
  double rho=v12/cc;
  double omr=1-R_pow(rho,2);
  double aa= 2*M_PI*cc*sqrt(omr);
  double zz=R_pow(x,2)/v11  + R_pow(y,2)/v22-2*rho*x*y/cc;
  res=exp(-0.5*zz/omr)/aa;
  return(res);
}


double d2lognorm(double x, double y, double sill, double nugget, double mux, double muy, double rho)
{
    double sigma = sqrt(sill);
    double sigma2 = sill;
    double mu_i = exp(mux);double mu_j = exp(muy);
    double v_i = x / mu_i;double v_j = y / mu_j;
    double rho_eff = rho;  // Correlazione effettiva
    double one_minus_rho2 = 1.0 - rho_eff * rho_eff;
    double term_i = (log(v_i) + sigma2/2.0) / sigma;
    double term_j = (log(v_j) + sigma2/2.0) / sigma;
    double quad = (term_i*term_i + term_j*term_j - 2.0*rho_eff*term_i*term_j) / (2.0 * one_minus_rho2);
    double prefactor = 1.0 / (2.0 * M_PI * sigma2 * v_i * v_j * sqrt(one_minus_rho2));
    double density_V = prefactor * exp(-quad);
    double result = density_V / (mu_i * mu_j);
    return result;
}

static inline double logcosh(double b){
  double ab = fabs(b);
  return ab + log1p(exp(-2.0*ab)) - M_LN2; // log(cosh b)
}

double biv_log_sas(double corr, double zi, double zj,
                   double mi, double mj, double skew, double tail, double vari)
{
  const double sd = sqrt(vari);
  const double xi = (zi - mi) / sd;
  const double xj = (zj - mj) / sd;
  const double b1 = tail * asinh(xi) - skew;
  const double b2 = tail * asinh(xj) - skew;
  const double Z1 = sinh(b1);
  const double Z2 = sinh(b2);
  if (!isfinite(Z1) || !isfinite(Z2)) return -INFINITY;

  const double one_m_r2 = 1.0 - corr*corr;
  const double inv_one_m_r2 = 1.0 / one_m_r2;
  const double quad = Z1*Z1 + Z2*Z2 - 2.0*corr*Z1*Z2;
  const double exponent = -0.5 * quad * inv_one_m_r2;
  const double log_phi2 = -log(2.0*M_PI) - 0.5*log(one_m_r2) + exponent;
  const double log_jac =
      2.0*log(tail)
    + logcosh(b1) + logcosh(b2)
    - 0.5*log1p(xi*xi) - 0.5*log1p(xj*xj)
    - log(vari);

  return log_phi2 + log_jac;
}

double biv_sinh(double corr,double zi,double zj,double mi,double mj,double skew,double tail,double vari){
  double lp = biv_log_sas(corr, zi, zj, mi, mj, skew, tail, vari);
  return (lp == -INFINITY) ? 0.0 : exp(lp);
}



double biv_two_piece_bimodal(double rho,double zi,double zj,double sill,double nu,double delta,double eta,
             double p11,double mui,double muj)
{
double res=0.0;
double alpha=2*(delta+1)/nu;
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);

double nn=R_pow(2,1-alpha/2);

if(zi>=mui&&zj>=muj)
{res=          (R_pow(zistd*zjstd,alpha-1)*p11/R_pow(etamos,2*alpha))*biv_gamma_gen(rho,R_pow(zistd/etamos,alpha),R_pow(zjstd/etamos,alpha),0,0,nu,nn);  }
if(zi>=mui&&zj<muj)
{res=R_pow(-zistd*zjstd,alpha-1)*(1-eta-2*p11)/(R_pow(etamos,alpha)*R_pow(etamas,alpha))*biv_gamma_gen(rho,R_pow(zistd/etamos,alpha),R_pow(-zjstd/etamas,alpha),0,0,nu,nn);}
if(zi<mui&&zj>=muj)
{res=R_pow(-zistd*zjstd,alpha-1)*(1-eta-2*p11)/(R_pow(etamos,alpha)*R_pow(etamas,alpha))*biv_gamma_gen(rho,R_pow(-zistd/etamas,alpha),R_pow(zjstd/etamos,alpha),0,0,nu,nn);}
if(zi<mui&&zj<muj)
{res=     (R_pow(zistd*zjstd,alpha-1)*(p11+eta)/R_pow(etamas,2*alpha))*biv_gamma_gen(rho,R_pow(-zistd/etamas,alpha),R_pow(-zjstd/etamas,alpha),0,0,nu,nn);}
return(R_pow(alpha,2)*res/sill);
}



double biv_Weibull(double corr,double zi,double zj,double mui,double muj,double shape){

/*
 * Cache constants that depend only on the Weibull shape.  During one
 * composite-likelihood evaluation, shape is common to all pairs, so the
 * original code recomputed gammafn(), R_pow() and logarithms thousands of
 * times with identical arguments.  The cache is process-local, matching
 * the correlation-kernel context caches used elsewhere in the package.
 */
typedef struct {
    int valid;
    double shape;
    double half_shape;
    double shape_minus_one;
    double kpow;
    double log_prefactor;
} WeibullShapeContext;

static WeibullShapeContext ctx = {0, 0.0, 0.0, 0.0, 0.0, 0.0};

if (!ctx.valid || shape != ctx.shape) {
    const double k = 1.0 / gammafn(1.0 + 1.0 / shape);
    ctx.shape = shape;
    ctx.half_shape = 0.5 * shape;
    ctx.shape_minus_one = shape - 1.0;
    ctx.kpow = R_pow(k, -shape);
    ctx.log_prefactor = 2.0 * log(shape) - 2.0 * shape * log(k);
    ctx.valid = 1;
}

const double ci=exp(mui);
const double cj=exp(muj);
const double ui=zi/ci;
const double uj=zj/cj;
const double a=1.0-corr*corr;

double uipow, ujpow, uipow_half, ujpow_half;

/* Common integer shapes avoid generic power calls entirely. */
if (shape == 2.0) {
    uipow = ui * ui;
    ujpow = uj * uj;
    uipow_half = ui;
    ujpow_half = uj;
} else if (shape == 1.0) {
    uipow = ui;
    ujpow = uj;
    uipow_half = (ui >= 0.0) ? sqrt(ui) : R_pow(ui, ctx.half_shape);
    ujpow_half = (uj >= 0.0) ? sqrt(uj) : R_pow(uj, ctx.half_shape);
} else {
    uipow = R_pow(ui, shape);
    ujpow = R_pow(uj, shape);

    /* Weibull support is non-negative: sqrt(u^shape) = u^(shape/2). */
    uipow_half = (ui >= 0.0 && R_FINITE(uipow))
        ? sqrt(uipow) : R_pow(ui, ctx.half_shape);
    ujpow_half = (uj >= 0.0 && R_FINITE(ujpow))
        ? sqrt(ujpow) : R_pow(uj, ctx.half_shape);
}

const double z=2.0*fabs(corr)*uipow_half*ujpow_half*ctx.kpow/a;
const double A=ctx.log_prefactor+ctx.shape_minus_one*log(ui*uj)-log(a);
const double B=-ctx.kpow*(uipow+ujpow)/a;
const double res=A+B+log(bessel_i(z,0,2))+z-(mui+muj);
return res;
}

double psi(int i,int j,double p1,double p2,double p12,double**cache){
    if(i==0||j==0)return 0.0;
    if(i==1&&j==1)return(p1+p2-p1*p2)/(p2*p1*(p1+p2-p12));
    if(cache[i][j]!=-1.0)return cache[i][j];

    double aux=p1+p2-p12;
    double psi11=psi(i-1,j-1,p1,p2,p12,cache);
    double psi01=psi(i,j-1,p1,p2,p12,cache);
    double psi10=psi(i-1,j,p1,p2,p12,cache);
    double aux1=(psi11*p12+psi01*(p1-p12)+psi10*(p2-p12))/aux;
    double aux2=((i-p2/aux)/p2+(j-p1/aux)/p1)/aux;
    double aux3=(2.0-aux)/(aux*aux);
    cache[i][j]=aux1+aux2+aux3;
    return cache[i][j];
}

double** create_cache(int max_i,int max_j){
    const size_t nrow = (size_t) max_i + 1u;
    const size_t ncol = (size_t) max_j + 1u;
    double **cache = (double **) R_Calloc(nrow, double *);
    double *storage = (double *) R_Calloc(nrow * ncol, double);

    for (int i = 0; i <= max_i; i++) {
        cache[i] = storage + (size_t) i * ncol;
        for (int j = 0; j <= max_j; j++) {
            cache[i][j] = -1.0;
        }
    }
    return cache;
}

static void free_cache(double **cache){
    if (cache == NULL) return;
    R_Free(cache[0]);
    R_Free(cache);
}

double cov_binom_neg(int m,double p11,double p1, double p2) {
    double **cache;
    double a;

    /* Avoid allocating the recursion cache in the independence case. */
    if (fabs(p11 - p1 * p2) < 1e-32) return 0.0;

    cache = create_cache(m, m);
    a = psi(m, m, p1, p2, p11, cache) - ((double) m * (double) m) / (p1 * p2);
    free_cache(cache);
    return a;
}

/********************************************************************/
/********* poisson gamma correlation*********************************/
/********************************************************************/
/* Corrected PEM (Poisson-Erlang mixture) correlation.
 *
 * Theorem 5.1, with nu_i=a/mean_i and nu_j=a/mean_j:
 *
 * rho_M = rho^2 / sqrt((1+nu_i)(1+nu_j))
 *       + A * sum_{r,l>=0} C_{r,l} H_{r,l}(nu_i,nu_j),
 *
 * where H is the product of two regularized Gauss 2F1 functions.
 * hypergeo() returns the ordinary 2F1, so Gamma(r+2)^(-2) is included
 * explicitly in the logarithmic coefficient below.
 *
 * This core is used for the non-stationary case and as a numerical fallback.
 * The stationary case uses the closed Corollary 5.1 below.  Stationarity
 * means identical marginal distributions, not Lambda_i == Lambda_j.
 */

#define MIN_CORRELATION_THRESHOLD 1e-15
#define MIN_MEAN_THRESHOLD 1e-15
#define MIN_A_THRESHOLD 1e-15
#define MAX_NU_THRESHOLD 1e10
#ifndef RHO_EPS
#define RHO_EPS 1e-15
#endif

typedef struct {
    double sum;
    double compensation;
} KahanSum;

static R_INLINE void kahan_init(KahanSum *ks)
{
    ks->sum = 0.0;
    ks->compensation = 0.0;
}

static R_INLINE void kahan_add(KahanSum *ks, double value)
{
    const double y = value - ks->compensation;
    const double t = ks->sum + y;
    ks->compensation = (t - ks->sum) - y;
    ks->sum = t;
}

/* Poisson--Gamma numerical status.  These kernels are deliberately serial;
 * no OpenMP/thread-local state is used. */
typedef enum {
    PG_STATUS_OK = 0,
    PG_STATUS_INVALID_ARGUMENT = 1,
    PG_STATUS_SPECIAL_FUNCTION = 2,
    PG_STATUS_MAXITER = 3,
    PG_STATUS_INVALID_PROBABILITY = 4
} PGStatus;

static int pg_last_status = PG_STATUS_OK;

static R_INLINE void pg_set_status(PGStatus status)
{
    pg_last_status = (int)status;
}

static R_INLINE double pg_fail(PGStatus status)
{
    pg_set_status(status);
    return R_NaN;
}

void GeoPEMStatusGet(int *out)
{
    if (out) *out = pg_last_status;
}

#ifndef PG_PROB_TOL
#define PG_PROB_TOL 1e-10
#endif

#ifndef PG_POISSON_REL_TOL
#define PG_POISSON_REL_TOL 1e-8
#endif
#ifndef PG_MAX_OUTER
#define PG_MAX_OUTER 900
#endif
/*
 * The off-diagonal Poisson-Gamma series can require a longer outer tail
 * for large counts and high correlation.  Keep the global limit unchanged
 * and extend only the PGrt outer loop; ordinary cases still stop through the
 * existing relative convergence criterion.
 */
#ifndef PG_PGRT_MAX_OUTER
#define PG_PGRT_MAX_OUTER 1500
#endif
#ifndef PG_MAX_INNER
#define PG_MAX_INNER 900
#endif
#ifndef PG_SUM_TOL
#define PG_SUM_TOL 2e-14
#endif
#ifndef PG_MIN_QUIET
#define PG_MIN_QUIET 4
#endif

/* Safety guards for the non-stationary PEM correlation only.
 * The series itself is stopped by the relative convergence checks. */
#ifndef PEM_CORR_MAX_OUTER
#define PEM_CORR_MAX_OUTER 200000
#endif
#ifndef PEM_CORR_MAX_INNER
#define PEM_CORR_MAX_INNER 8192
#endif
#ifndef PEM_CORR_FAST_TOL
#define PEM_CORR_FAST_TOL 1e-10
#endif

static R_INLINE double pg_logspace_add_pos(double a, double b)
{
    if (a == R_NegInf) return b;
    if (b == R_NegInf) return a;
    if (!R_FINITE(a) || !R_FINITE(b)) return R_NaN;
    if (a < b) { const double tmp = a; a = b; b = tmp; }
    return a + log1p(exp(b - a));
}

static R_INLINE double pg_logspace_sub_pos(double a, double b)
{
    if (a == R_NegInf) return R_NaN;
    if (b == R_NegInf) return a;
    if (!R_FINITE(a) || !R_FINITE(b)) return R_NaN;
    if (b > a) {
        if (b - a <= 64.0 * DBL_EPSILON * fmax2(1.0, fabs(a)))
            return R_NegInf;
        return R_NaN;
    }
    if (a == b) return R_NegInf;
    {
        const double u = -expm1(b - a);
        if (!(u > 0.0) || !R_FINITE(u)) return R_NaN;
        return a + log(u);
    }
}

static R_INLINE double pg_exp_log_signed(double logabs, int sign)
{
    if (logabs == R_NegInf) return 0.0;
    if (!R_FINITE(logabs) || logabs > 700.0) return R_NaN;
    if (logabs < -745.0) return 0.0;
    return sign < 0 ? -exp(logabs) : exp(logabs);
}



/* q-sequence: I_w(p,q+l), B(p,q+l), l=0,1,... . */
typedef struct {
    double p, shape, w, logw, log1mw;
    double logI, logDelta, logBeta;
    int last_l1;
} PGBetaBSeq;

static R_INLINE int pg_beta_bseq_init_seeded(PGBetaBSeq *st,
                                              double p, double shape,
                                              double w, double seeded_logI,
                                              double seeded_logBeta)
{
    if (!(p > 0.0) || !(shape > 0.0) || !(w > 0.0 && w < 1.0) ||
        !R_FINITE(seeded_logBeta) || ISNAN(seeded_logI)) return 0;
    st->p = p;
    st->shape = shape;
    st->w = w;
    st->logw = log(w);
    st->log1mw = log1p(-w);
    st->logBeta = seeded_logBeta;
    st->logI = seeded_logI;
    st->logDelta = p * st->logw + shape * st->log1mw -
                   log(shape) - st->logBeta;
    st->last_l1 = 0;
    return !ISNAN(st->logDelta);
}


static R_INLINE double pg_beta_bseq_get(PGBetaBSeq *st, int l1)
{
    if (l1 < 0 || l1 < st->last_l1) return R_NaN;
    while (st->last_l1 < l1) {
        const double q = st->shape + (double)st->last_l1;
        st->logI = pg_logspace_add_pos(st->logI, st->logDelta);
        if (ISNAN(st->logI)) return R_NaN;
        st->logBeta += log(q) - log(st->p + q);
        st->logDelta += st->log1mw + log(st->p + q) - log(q + 1.0);
        ++st->last_l1;
    }
    return st->logI;
}

/* p-sequence: I_w(p+k,q), B(p+k,q), k=0,1,... . */
typedef struct {
    double p, q, w, logw, log1mw;
    double logI, logDelta, logBeta;
    int steps_since_anchor;
} PGBetaPSeq;

static R_INLINE int pg_beta_pseq_init(PGBetaPSeq *st, double p, double q,
                                      double w)
{
    if (!(p > 0.0) || !(q > 0.0) || !(w > 0.0 && w < 1.0)) return 0;
    st->p = p;
    st->q = q;
    st->w = w;
    st->logw = log(w);
    st->log1mw = log1p(-w);
    st->logBeta = lbeta(p, q);
    st->logI = pbeta(w, p, q, 1, 1);
    st->logDelta = p * st->logw + q * st->log1mw -
                   log(p) - st->logBeta;
    st->steps_since_anchor = 0;
    return !ISNAN(st->logI) && !ISNAN(st->logDelta) &&
           R_FINITE(st->logBeta);
}

static R_INLINE int pg_beta_pseq_advance(PGBetaPSeq *st)
{
    const double oldp = st->p;
    const double newp = oldp + 1.0;
    double next = pg_logspace_sub_pos(st->logI, st->logDelta);
    if (ISNAN(next) || st->logDelta > st->logI - 18.0 ||
        st->steps_since_anchor >= 31) {
        next = pbeta(st->w, newp, st->q, 1, 1);
        st->steps_since_anchor = 0;
    } else {
        ++st->steps_since_anchor;
    }
    if (ISNAN(next)) return 0;
    st->logDelta += st->logw + log(oldp + st->q) - log(newp);
    st->logBeta += log(oldp) - log(oldp + st->q);
    st->logI = next;
    st->p = newp;
    return !ISNAN(st->logDelta) && R_FINITE(st->logBeta);
}

static R_INLINE int pg_beta_bseq_from_pseq(PGBetaBSeq *st,
                                           const PGBetaPSeq *seed)
{
    return pg_beta_bseq_init_seeded(st, seed->p, seed->q, seed->w,
                                    seed->logI, seed->logBeta);
}

/* Fast positive-series evaluation of the non-stationary PEM correlation.
 *
 * In the stable representation of Theorem 5.1, for q=a+l and
 * w_i=mu_i/(a+mu_i),
 *
 *   I_{w_i}(r+1,q) = P(X_i > r),
 *
 * where X_i is negative binomial with size q and success probability
 * p_i=a/(a+mu_i).  Therefore the inner sum over r is
 *
 *   sum_{r>=0} I_{w_i}(r+1,q) I_{w_j}(r+1,q)
 *     = E[min(X_i,X_j)]
 *
 * for independent X_i and X_j.  Evaluating the two NB survival sequences
 * recursively avoids thousands of incomplete-beta evaluations per
 * correlation.  Occasional pbeta anchors prevent accumulated roundoff.
 *
 * This is algebraically the same positive-term representation used by
 * pem_corr_corrected_core; it only changes the order and numerical method
 * used to sum the terms.  R_NaN means "use the robust fallback".
 */
static double pem_nb_min_expect_fast(double q, double pi, double pj,
                                     double wi, double wj)
{
    if (!(q > 0.0) || !(pi > 0.0 && pi < 1.0) ||
        !(pj > 0.0 && pj < 1.0) ||
        !(wi > 0.0 && wi < 1.0) || !(wj > 0.0 && wj < 1.0))
        return R_NaN;

    const double logpmi0 = q * log(pi);
    const double logpmj0 = q * log(pj);

    /* If the mass at zero underflows, a forward PMF recurrence cannot
       recover reliably.  The general beta-based evaluator remains the
       safe fallback for such extreme parameter combinations. */
    if (logpmi0 < -700.0 || logpmj0 < -700.0)
        return R_NaN;

    double pmi = exp(logpmi0);
    double pmj = exp(logpmj0);
    double si = -expm1(logpmi0);  /* P(X_i > 0) */
    double sj = -expm1(logpmj0);  /* P(X_j > 0) */

    if (!R_FINITE(pmi) || !R_FINITE(pmj) ||
        !R_FINITE(si) || !R_FINITE(sj))
        return R_NaN;

    KahanSum total;
    kahan_init(&total);
    int quiet = 0;

    for (int r = 0; r < PEM_CORR_MAX_OUTER; ++r) {
        if (si < 0.0 && si > -1e-13) si = 0.0;
        if (sj < 0.0 && sj > -1e-13) sj = 0.0;
        if (!(si >= 0.0 && si <= 1.0 + 1e-12) ||
            !(sj >= 0.0 && sj <= 1.0 + 1e-12))
            return R_NaN;

        const double term = fmax2(0.0, si) * fmax2(0.0, sj);
        kahan_add(&total, term);

        if (r > 12) {
            const double scale = fmax2(1e-300, fabs(total.sum));
            if (term <= PEM_CORR_FAST_TOL * scale) {
                if (++quiet >= PG_MIN_QUIET)
                    return total.sum;
            } else {
                quiet = 0;
            }
        }

        /* Advance from r to r+1. */
        const double rr = (double)r;
        const double denom = rr + 1.0;
        pmi *= ((rr + q) / denom) * wi;
        pmj *= ((rr + q) / denom) * wj;

        if (!R_FINITE(pmi) || !R_FINITE(pmj) || pmi < 0.0 || pmj < 0.0)
            return R_NaN;

        si -= pmi;
        sj -= pmj;

        /* Re-anchor the tails periodically, and immediately if subtraction
           has produced a negative survival probability. */
        if (((r + 1) & 255) == 0 || si < 0.0 || sj < 0.0) {
            const double ppar = rr + 2.0;
            si = pbeta(wi, ppar, q, 1, 0);
            sj = pbeta(wj, ppar, q, 1, 0);
            if (!R_FINITE(si) || !R_FINITE(sj))
                return R_NaN;
        }
    }

    return R_NaN;
}

static double pem_corr_positive_fast(double corr, double mean_i,
                                     double mean_j, double a)
{
    const double rho2 = corr * corr;
    const double nu_i = a / mean_i;
    const double nu_j = a / mean_j;
    const double wi = 1.0 / (1.0 + nu_i);
    const double wj = 1.0 / (1.0 + nu_j);
    const double pi = nu_i / (1.0 + nu_i);
    const double pj = nu_j / (1.0 + nu_j);
    const double logD = log1p(nu_i) + log1p(nu_j);
    const double term_erlang = rho2 * exp(-0.5 * logD);

    if (!(rho2 > 0.0 && rho2 < 1.0) ||
        !(nu_i > 0.0) || !(nu_j > 0.0) ||
        !(wi > 0.0 && wi < 1.0) || !(wj > 0.0 && wj < 1.0))
        return R_NaN;

    const double log_pref = log(rho2) + (a + 1.0) * log1p(-rho2) +
                            log(nu_i) + log(nu_j) - lgammafn(a + 1.0) -
                            0.5 * logD;
    if (!R_FINITE(log_pref) || log_pref > 700.0)
        return R_NaN;
    if (log_pref < -745.0)
        return term_erlang;

    const double pref = exp(log_pref);
    KahanSum total;
    kahan_init(&total);
    int quiet_l = 0;
    double logcoef = lgammafn(a); /* log{rho2^l Gamma(a+l)/l!}, l=0 */

    for (int l = 0; l < PEM_CORR_MAX_INNER; ++l) {
        const double q = a + (double)l;
        const double inner = pem_nb_min_expect_fast(q, pi, pj, wi, wj);
        if (!R_FINITE(inner) || inner < 0.0)
            return R_NaN;

        double contribution = 0.0;
        if (inner > 0.0) {
            const double lc = logcoef + log(inner);
            if (!R_FINITE(lc) || lc > 700.0)
                return R_NaN;
            contribution = lc < -745.0 ? 0.0 : exp(lc);
        }
        kahan_add(&total, contribution);

        if (l > 10) {
            const double scale = fmax2(1e-300, fabs(total.sum));
            if (contribution <= PEM_CORR_FAST_TOL * scale) {
                if (++quiet_l >= PG_MIN_QUIET) {
                    const double ans = term_erlang + pref * total.sum;
                    if (R_FINITE(ans) && ans >= -PG_PROB_TOL &&
                        ans <= 1.0 + PG_PROB_TOL)
                        return fmax2(0.0, fmin2(1.0, ans));
                    return R_NaN;
                }
            } else {
                quiet_l = 0;
            }
        }

        if (l + 1 < PEM_CORR_MAX_INNER)
            logcoef += log(rho2) + log(a + (double)l) -
                       log((double)l + 1.0);
    }

    return R_NaN;
}

static double pem_corr_corrected_core(double corr, double mean_i,
                                      double mean_j, double a,
                                      int stationary)
{
    pg_set_status(PG_STATUS_OK);
    if (!R_FINITE(corr) || !R_FINITE(mean_i) || !R_FINITE(mean_j) ||
        !R_FINITE(a) || fabs(corr) > 1.0 + PG_PROB_TOL ||
        mean_i <= MIN_MEAN_THRESHOLD || mean_j <= MIN_MEAN_THRESHOLD ||
        a <= MIN_A_THRESHOLD)
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    const int equal_means =
        fabs(mean_i - mean_j) <=
        1e-12 * fmax2(1.0, fmax2(mean_i, mean_j));
    stationary = stationary || equal_means;

    if (fabs(corr) < MIN_CORRELATION_THRESHOLD) return 0.0;
    if (fabs(corr) >= 1.0 - RHO_EPS) {
        if (equal_means) return 1.0;
        corr = copysign(1.0 - RHO_EPS, corr);
    }

    /* The negative-binomial marginals converge to Poisson when mu/kappa
       is negligible.  This also avoids very large nu in the PEM series. */
    if (fmax2(mean_i, mean_j) / a <= PG_POISSON_REL_TOL) {
        const double out = corr_pois(corr, mean_i, mean_j);
        if (!R_FINITE(out)) return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
        return fmax2(0.0, fmin2(1.0, out));
    }

    const double rho2 = corr * corr;
    const double nu_i = a / mean_i;
    const double nu_j = a / mean_j;
    if (!R_FINITE(nu_i) || !R_FINITE(nu_j) ||
        !(nu_i > 0.0) || !(nu_j > 0.0))
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    /* The negative-binomial survival recurrence is substantially faster
       for non-stationary means, especially for small kappa.  Keep the
       beta-based positive series below as a robust fallback for extreme
       cases in which the fast recurrence cannot be initialized safely. */
    if (!stationary) {
        const double fast = pem_corr_positive_fast(corr, mean_i, mean_j, a);
        if (R_FINITE(fast)) return fast;
    }

    const double wi = 1.0 / (1.0 + nu_i);
    const double wj = 1.0 / (1.0 + nu_j);
    const double log_rho2 = log(rho2);
    const double log1m_rho2 = log1p(-rho2);
    const double logD = log1p(nu_i) + log1p(nu_j);

    /*
     * Stable positive-term form of Theorem 5.1.
     *
     * For p=r+1, q_l=a+l, and w=1/(1+nu), the regularized Gauss factor
     * appearing in the theorem satisfies
     *
     *  2F1~(1,1-q_l;p+1;-1/nu)
     *   = (1+nu)^(p+q_l-1) nu^(1-q_l)
     *     Gamma(q_l)/Gamma(p+q_l) I_w(p,q_l),
     *
     * where I_w is the regularized incomplete beta function.  Substitution
     * into the double series cancels the large Gamma and power factors
     * exactly and gives
     *
     *  rho_M = rho2/sqrt(D)
     *        + C sum_r sum_l rho2^l Gamma(a+l)/l!
     *              I_wi(r+1,a+l) I_wj(r+1,a+l),
     *
     *  C = rho2 (1-rho2)^(a+1) nu_i nu_j /
     *      {Gamma(a+1) sqrt(D)},
     *  D = (1+nu_i)(1+nu_j).
     *
     * Every summand is non-negative.  This form is algebraically identical
     * to the theorem but avoids the cancellation and slow hypergeometric
     * scaling that caused non-finite correlations for small kappa and large
     * non-stationary means.
     */
    const double term_erlang = exp(log_rho2 - 0.5 * logD);
    const double log_pref = log_rho2 + (a + 1.0) * log1m_rho2 +
                            log(nu_i) + log(nu_j) - lgammafn(a + 1.0) -
                            0.5 * logD;
    if (!R_FINITE(log_pref)) return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
    if (log_pref < -745.0) return term_erlang;
    if (log_pref > 700.0) return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
    const double pref = exp(log_pref);

    PGBetaPSeq pseq_i, pseq_j;
    if (!pg_beta_pseq_init(&pseq_i, 1.0, a, wi) ||
        (!stationary && !pg_beta_pseq_init(&pseq_j, 1.0, a, wj)))
        return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

    KahanSum total;
    kahan_init(&total);
    int quiet_r = 0;
    int outer_converged = 0;

    /* These are guards only.  The actual truncation is controlled by the
       relative contribution tests below.  The larger correlation-specific
       limits are needed because I_w(r+1,q) can decay slowly when nu is small
       (large local means), even though the series itself is well behaved. */
    for (int r = 0; r < PEM_CORR_MAX_OUTER; ++r) {
        KahanSum row;
        kahan_init(&row);
        int quiet_l = 0;
        int row_converged = 0;

        PGBetaBSeq bseq_i, bseq_j;
        if (!pg_beta_bseq_from_pseq(&bseq_i, &pseq_i) ||
            (!stationary && !pg_beta_bseq_from_pseq(&bseq_j, &pseq_j)))
            return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

        double logcoef_l = lgammafn(a);
        if (!R_FINITE(logcoef_l))
            return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

        for (int l = 0; l < PEM_CORR_MAX_INNER; ++l) {
            const double li = pg_beta_bseq_get(&bseq_i, l);
            const double lj = stationary ? li : pg_beta_bseq_get(&bseq_j, l);
            if (ISNAN(li) || ISNAN(lj))
                return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

            double term = 0.0;
            if (li != R_NegInf && lj != R_NegInf) {
                const double logterm = logcoef_l + li + lj;
                if (ISNAN(logterm) || logterm > 700.0)
                    return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
                term = (logterm < -745.0) ? 0.0 : exp(logterm);
            }

            if (!R_FINITE(term) || term < 0.0)
                return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
            kahan_add(&row, term);

            if (l > 10) {
                const double scale = fmax2(1e-300, fabs(row.sum));
                if (term <= PG_SUM_TOL * scale) {
                    if (++quiet_l >= PG_MIN_QUIET) {
                        row_converged = 1;
                        break;
                    }
                } else {
                    quiet_l = 0;
                }
            }

            if (l + 1 < PEM_CORR_MAX_INNER) {
                logcoef_l += log_rho2 + log(a + (double)l) -
                             log((double)l + 1.0);
            }
        }

        if (!row_converged)
            return pg_fail(PG_STATUS_MAXITER);

        kahan_add(&total, row.sum);

        if (r > 12) {
            const double scale = fmax2(1e-300, fabs(total.sum));
            if (fabs(row.sum) <= PG_SUM_TOL * scale) {
                if (++quiet_r >= PG_MIN_QUIET) {
                    outer_converged = 1;
                    break;
                }
            } else {
                quiet_r = 0;
            }
        }

        if (r + 1 < PEM_CORR_MAX_OUTER) {
            if (!pg_beta_pseq_advance(&pseq_i) ||
                (!stationary && !pg_beta_pseq_advance(&pseq_j)))
                return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
        }
    }

    if (!outer_converged)
        return pg_fail(PG_STATUS_MAXITER);
    if (total.sum < -PG_PROB_TOL || !R_FINITE(total.sum))
        return pg_fail(PG_STATUS_INVALID_PROBABILITY);

    const double result = term_erlang + pref * fmax2(0.0, total.sum);
    if (!R_FINITE(result) || result < -PG_PROB_TOL ||
        result > 1.0 + PG_PROB_TOL)
        return pg_fail(PG_STATUS_INVALID_PROBABILITY);

    return fmax2(0.0, fmin2(1.0, result));
}


/* Non-stationary PEM correlation from the simplified form of Theorem 5.1.
 *
 * Let q=rho^2, d=1-q, nu_i=kappa/mu_i, nu_j=kappa/mu_j, and define
 *
 *   A = d*nu_i*nu_j + nu_i + nu_j + 2,
 *   K = {A + sqrt[A^2-4(1+nu_i)(1+nu_j)]}/2,
 *   alpha = (1+nu_j)/K, beta = (1+nu_i)/K.
 *
 * Then
 *
 *   rho_M = q [ sqrt{(1+nu_i)(1+nu_j)}/K
 *             + d*nu_i*nu_j/(kappa*sqrt{(1+nu_i)(1+nu_j)})
 *               S_kappa(alpha,beta) ],
 *
 * where
 *
 *   S_kappa(alpha,beta)
 *     = sum_{h>=1} I_alpha(h,kappa) I_beta(h,kappa).
 *
 * This removes the outer Erlang-mixture sum used by the previous corrected
 * evaluator.  The beta tails are advanced by the existing stable p-sequence
 * recurrence, with periodic pbeta re-anchoring.  For kappa=1 and kappa=2 the
 * series is evaluated in closed form.  R_NaN means "use the legacy positive
 * double-series fallback".
 */
static double pem_corr_nonstationary_single_series(double corr,
                                                    double mean_i,
                                                    double mean_j,
                                                    double kappa)
{
    if (!R_FINITE(corr) || !R_FINITE(mean_i) || !R_FINITE(mean_j) ||
        !R_FINITE(kappa) || fabs(corr) > 1.0 + PG_PROB_TOL ||
        mean_i <= MIN_MEAN_THRESHOLD || mean_j <= MIN_MEAN_THRESHOLD ||
        kappa <= MIN_A_THRESHOLD)
        return R_NaN;

    if (fabs(corr) < MIN_CORRELATION_THRESHOLD) return 0.0;
    if (fabs(corr) >= 1.0 - RHO_EPS)
        corr = copysign(1.0 - RHO_EPS, corr);

    if (fmax2(mean_i, mean_j) / kappa <= PG_POISSON_REL_TOL) {
        const double out = corr_pois(corr, mean_i, mean_j);
        return R_FINITE(out) ? fmax2(0.0, fmin2(1.0, out)) : R_NaN;
    }

    const double q = corr * corr;
    const double d = 1.0 - q;
    const double nu_i = kappa / mean_i;
    const double nu_j = kappa / mean_j;
    if (!(q > 0.0 && q < 1.0) || !(d > 0.0) ||
        !R_FINITE(nu_i) || !R_FINITE(nu_j) ||
        !(nu_i > 0.0) || !(nu_j > 0.0))
        return R_NaN;

    const double B = 1.0 + nu_j;
    const double C = 1.0 + nu_i;
    const double uv = nu_i * nu_j;
    const double A = B + C + d * uv;

    /* Algebraically identical to A*A-4*B*C, but non-negative by
       construction and much safer when A is close to 2*sqrt(B*C). */
    const double BCdiff = B - C;
    const double disc = BCdiff * BCdiff +
                        d * uv * (2.0 * (B + C) + d * uv);
    if (!R_FINITE(A) || !R_FINITE(disc) || disc < 0.0)
        return R_NaN;

    const double K = 0.5 * (A + sqrt(disc));
    if (!R_FINITE(K) || !(K > 0.0)) return R_NaN;

    const double alpha = B / K;
    const double beta = C / K;
    if (!(alpha > 0.0 && alpha < 1.0) ||
        !(beta > 0.0 && beta < 1.0))
        return R_NaN;

    const double sqrtD = sqrt(B * C);
    if (!R_FINITE(sqrtD) || !(sqrtD > 0.0)) return R_NaN;

    double S = R_NaN;

    if (fabs(kappa - 1.0) <= 32.0 * DBL_EPSILON) {
        const double x = alpha * beta;
        const double den = 1.0 - x;
        if (!(x >= 0.0 && x < 1.0) || !(den > 0.0)) return R_NaN;
        S = x / den;
    } else if (fabs(kappa - 2.0) <= 64.0 * DBL_EPSILON) {
        const double x = alpha * beta;
        const double den = 1.0 - x;
        if (!(x >= 0.0 && x < 1.0) || !(den > 0.0)) return R_NaN;
        const double den2 = den * den;
        const double den3 = den2 * den;
        S = x / den +
            (2.0 - alpha - beta) * x / den2 +
            (1.0 - alpha) * (1.0 - beta) * x * (1.0 + x) / den3;
    } else {
        PGBetaPSeq seq_i, seq_j;
        if (!pg_beta_pseq_init(&seq_i, 1.0, kappa, alpha) ||
            !pg_beta_pseq_init(&seq_j, 1.0, kappa, beta))
            return R_NaN;

        KahanSum total;
        kahan_init(&total);
        const double asym_ratio = alpha * beta;
        double prev_term = R_NaN;
        int quiet = 0;
        int converged = 0;

        for (int h = 1; h <= PEM_CORR_MAX_OUTER; ++h) {
            const double li = seq_i.logI;
            const double lj = seq_j.logI;
            if (ISNAN(li) || ISNAN(lj)) return R_NaN;

            double term = 0.0;
            if (li != R_NegInf && lj != R_NegInf) {
                const double lt = li + lj;
                if (ISNAN(lt) || lt > 64.0 * DBL_EPSILON) return R_NaN;
                term = (lt < -745.0) ? 0.0 :
                       ((lt > 0.0) ? 1.0 : exp(lt));
            }
            if (!R_FINITE(term) || term < 0.0) return R_NaN;

            kahan_add(&total, term);

            if (term == 0.0) {
                converged = 1;
                break;
            }

            if (h > 8 && R_FINITE(prev_term) && prev_term > 0.0 &&
                term >= 0.0 && term <= prev_term) {
                const double ratio = term / prev_term;
                double rbound = fmax2(ratio, asym_ratio);
                if (rbound < 0.0) rbound = 0.0;
                if (rbound > 1.0 - 16.0 * DBL_EPSILON)
                    rbound = 1.0 - 16.0 * DBL_EPSILON;

                const double rem = (term == 0.0) ? 0.0 :
                                   term * rbound / (1.0 - rbound);
                const double scale = fmax2(1e-300, fabs(total.sum));
                if (R_FINITE(rem) && rem <= PG_SUM_TOL * scale) {
                    if (++quiet >= PG_MIN_QUIET) {
                        converged = 1;
                        break;
                    }
                } else {
                    quiet = 0;
                }
            } else {
                quiet = 0;
            }

            prev_term = term;
            if (h < PEM_CORR_MAX_OUTER) {
                if (!pg_beta_pseq_advance(&seq_i) ||
                    !pg_beta_pseq_advance(&seq_j))
                    return R_NaN;
            }
        }

        if (!converged || !R_FINITE(total.sum) || total.sum < 0.0)
            return R_NaN;
        S = total.sum;
    }

    if (!R_FINITE(S) || S < 0.0) return R_NaN;

    const double result = q *
        (sqrtD / K + d * nu_i * nu_j * S / (kappa * sqrtD));
    if (!R_FINITE(result) || result < -PG_PROB_TOL ||
        result > 1.0 + PG_PROB_TOL)
        return R_NaN;

    return fmax2(0.0, fmin2(1.0, result));
}

/* Stable special evaluation for kappa=3/2 in the stationary PEM
 * correlation.  In CorrPGs the hypergeometric factor is
 *
 *   F = 2F1(-1/2, 1/2; 2; z).
 *
 * At exactly kappa=3/2 the bundled generic 2F1 code sees
 * c-a-b=2 as an exact integer and, for z close to one, switches to its
 * integer-d digamma expansion.  That branch is unnecessarily delicate for
 * this very simple case and can make an optimizer behave discontinuously at
 * the exactly representable value 1.5.
 *
 * Use instead the exact elliptic-integral identity
 *
 *   F = 4[(z-1)K(z) + (z+1)E(z)] / (3*pi*z),
 *
 * evaluating K and E by the arithmetic-geometric mean.  For small z the
 * elliptic expression suffers cancellation, so use the defining convergent
 * series there.  No generic special-function implementation is modified.
 */
static double pem_corr_hyp_kappa15(double z)
{
    if (!R_FINITE(z) || z < 0.0 || z > 1.0) return R_NaN;
    if (z == 0.0) return 1.0;
    if (z >= 1.0 - 8.0 * DBL_EPSILON)
        return 8.0 / (3.0 * M_PI);

    if (z < 1.0e-2) {
        const double aa = -0.5;
        const double bb =  0.5;
        const double cc =  2.0;
        double sum = 1.0;
        double term = 1.0;
        for (int n = 0; n < 80; ++n) {
            term *= ((aa + (double)n) * (bb + (double)n) /
                     ((cc + (double)n) * ((double)n + 1.0))) * z;
            const double next = sum + term;
            if (!R_FINITE(next)) return R_NaN;
            sum = next;
            if (fabs(term) <= 4.0 * DBL_EPSILON * fmax2(1.0, fabs(sum)))
                break;
        }
        return sum;
    }

    double agm_a = 1.0;
    double agm_b = sqrt(1.0 - z);
    double corr_sum = 0.5 * z;  /* 2^(-1) c_0^2, c_0^2=z */
    double weight = 1.0;

    for (int it = 0; it < 32; ++it) {
        const double cnext = 0.5 * (agm_a - agm_b);
        const double anext = 0.5 * (agm_a + agm_b);
        const double bnext = sqrt(agm_a * agm_b);

        corr_sum += weight * cnext * cnext;
        weight *= 2.0;
        agm_a = anext;
        agm_b = bnext;

        if (fabs(cnext) <= 4.0 * DBL_EPSILON * agm_a)
            break;
    }

    const double K = M_PI / (2.0 * agm_a);
    const double E = K * (1.0 - corr_sum);
    const double num = (z - 1.0) * K + (z + 1.0) * E;
    const double F = 4.0 * num / (3.0 * M_PI * z);
    return R_FINITE(F) ? F : R_NaN;
}

#define PEM_CORR_INTEGER_MAX 4096

static int pem_corr_integer_degree(double a)
{
    const double aa = nearbyint(a);
    if (!(aa >= 1.0) || aa > (double)(PEM_CORR_INTEGER_MAX + 1) ||
        fabs(a - aa) > 1e-12 * fmax2(1.0, fabs(a)))
        return -1;
    return (int)aa - 1;
}

/* Stable contiguous recurrence for
 * F_m(z)=2F1(-m,1/2;2;z), m=0,1,... .  Direct Horner evaluation of the
 * terminating polynomial suffers severe coefficient cancellation near z=1
 * for large m; this recurrence remains positive and stable. */
static double pem_corr_hyp_integer(int degree, double z)
{
    if (degree < 0 || degree > PEM_CORR_INTEGER_MAX ||
        !R_FINITE(z) || z < 0.0 || z > 1.0)
        return R_NaN;
    if (degree == 0) return 1.0;

    double previous = 1.0;
    double current = 1.0 - 0.25 * z;
    if (degree == 1) return current;

    double apar = -1.0;
    for (int m = 1; m < degree; ++m) {
        const double den = 2.0 - apar;
        const double next =
            -(2.0 * apar - 2.0 - apar * z + 0.5 * z) * current / den -
            apar * (z - 1.0) * previous / den;
        if (!R_FINITE(next) || !(next > 0.0)) return R_NaN;
        previous = current;
        current = next;
        apar -= 1.0;
    }
    return current;
}

/* Stationary PEM correlation, Corollary 5.1.
 *
 * With q=rho^2 and nu=kappa/mu, the corrected stationary result is
 *
 *   rho_M = q [ 1 - nu*sqrt(1-q)/sqrt(nu^2(1-q)+4(1+nu))
 *                 * 2F1(1-kappa,1/2;2;z) ],
 *
 *   z = 4(1+nu)/(nu^2(1-q)+4(1+nu)).
 *
 * Gamma(2)=1, so the ordinary and regularized 2F1 coincide here.
 * This closed form replaces the double-series evaluation in the stationary
 * case.  If the bundled 2F1 evaluator fails, fall back to the corrected
 * non-stationary theorem specialized to equal means.
 */
static double pem_corr_stationary_eval(double corr, double mean, double a)
{
    pg_set_status(PG_STATUS_OK);
    if (!R_FINITE(corr) || !R_FINITE(mean) || !R_FINITE(a) ||
        fabs(corr) > 1.0 + PG_PROB_TOL ||
        mean <= MIN_MEAN_THRESHOLD || a <= MIN_A_THRESHOLD)
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    if (fabs(corr) < MIN_CORRELATION_THRESHOLD) return 0.0;
    if (fabs(corr) >= 1.0 - RHO_EPS) return 1.0;

    if (mean / a <= PG_POISSON_REL_TOL) {
        const double out = corr_pois_s(corr, mean);
        if (!R_FINITE(out)) return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
        return fmax2(0.0, fmin2(1.0, out));
    }

    const double q = corr * corr;
    const double d = 1.0 - q;
    const double nu = a / mean;
    if (!R_FINITE(nu) || !(nu > 0.0))
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    const double den = nu * nu * d + 4.0 * (1.0 + nu);
    if (!R_FINITE(den) || !(den > 0.0))
        return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

    double z = 4.0 * (1.0 + nu) / den;
    if (z < 0.0 && z > -64.0 * DBL_EPSILON) z = 0.0;
    if (z > 1.0 && z < 1.0 + 64.0 * DBL_EPSILON) z = 1.0;
    if (!(z >= 0.0 && z <= 1.0))
        return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

    const double scale = nu * sqrt(d / den);
    const int degree = pem_corr_integer_degree(a);
    double F;
    if (degree >= 0) {
        F = pem_corr_hyp_integer(degree, z);
    } else if (fabs(a - 1.5) <= 32.0 * DBL_EPSILON) {
        F = pem_corr_hyp_kappa15(z);
    } else {
        F = hypergeo(1.0 - a, 0.5, 2.0, z);
    }

    if (!R_FINITE(F) || !(F > 0.0))
        return pem_corr_corrected_core(corr, mean, mean, a, 1);

    const double result = q * (1.0 - scale * F);
    if (!R_FINITE(result) || result < -PG_PROB_TOL ||
        result > 1.0 + PG_PROB_TOL)
        return pem_corr_corrected_core(corr, mean, mean, a, 1);

    return fmax2(0.0, fmin2(1.0, result));
}

double CorrPGs(double corr, double mean, double a)
{
    return pem_corr_stationary_eval(corr, mean, a);
}

/* Non-stationary Theorem 5.1.  Keep the historical function name to avoid
 * changing callers.  The new single-series representation is the active
 * path; the previous positive double-series evaluator is retained only as a
 * robust fallback for extreme numerical cases. */
double CorrPGns_AdaptiveBlocks(double corr, double mean_i, double mean_j, double a)
{
    pg_set_status(PG_STATUS_OK);
    const double out =
        pem_corr_nonstationary_single_series(corr, mean_i, mean_j, a);
    if (R_FINITE(out)) return out;
    return pem_corr_corrected_core(corr, mean_i, mean_j, a, 0);
}

double corr_pois_gen(double corr, double mean_i, double mean_j, double a)
{
    pg_set_status(PG_STATUS_OK);
    if (!R_FINITE(corr) || !R_FINITE(mean_i) || !R_FINITE(mean_j) ||
        !R_FINITE(a) || fabs(corr) > 1.0 + PG_PROB_TOL ||
        mean_i <= MIN_MEAN_THRESHOLD || mean_j <= MIN_MEAN_THRESHOLD ||
        a <= MIN_A_THRESHOLD)
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    if (fabs(corr) < MIN_CORRELATION_THRESHOLD) return 0.0;
    const int equal_means =
        fabs(mean_i - mean_j) <=
        1e-12 * fmax2(1.0, fmax2(mean_i, mean_j));
    if (fabs(corr) >= 1.0 - RHO_EPS) {
        if (equal_means) return 1.0;
        corr = copysign(1.0 - RHO_EPS, corr);
    }
    if (equal_means)
        return CorrPGs(corr, 0.5 * (mean_i + mean_j), a);
    return CorrPGns_AdaptiveBlocks(corr, mean_i, mean_j, a);
}

void CorrPGsVec(double *corr, double *mean, double *shape, int *n, double *out)
{
    if (!corr || !mean || !shape || !n || !out || n[0] < 0) return;
    for (int i = 0; i < n[0]; ++i)
        out[i] = CorrPGs(corr[i], mean[0], shape[0]);
}

SEXP GeoCorrPGsVecCall(SEXP corrS, SEXP meanS, SEXP shapeS)
{
    SEXP corrR, meanR, shapeR, ans;
    int nprotect = 0;
    PROTECT(corrR = coerceVector(corrS, REALSXP)); ++nprotect;
    PROTECT(meanR = coerceVector(meanS, REALSXP)); ++nprotect;
    PROTECT(shapeR = coerceVector(shapeS, REALSXP)); ++nprotect;

    const R_xlen_t n = XLENGTH(corrR);
    const R_xlen_t nm = XLENGTH(meanR);
    const R_xlen_t ns = XLENGTH(shapeR);
    if (!((nm == 1 || nm == n) && (ns == 1 || ns == n))) {
        UNPROTECT(nprotect);
        error("GeoCorrPGsVecCall: mean and shape must have length 1 or length(corr)");
    }

    PROTECT(ans = allocVector(REALSXP, n)); ++nprotect;
    const double *cc = REAL(corrR);
    const double *mm = REAL(meanR);
    const double *aa = REAL(shapeR);
    double *out = REAL(ans);
    for (R_xlen_t i = 0; i < n; ++i)
        out[i] = CorrPGs(cc[i], mm[nm == 1 ? 0 : i], aa[ns == 1 ? 0 : i]);

    UNPROTECT(nprotect);
    return ans;
}

/******************************************/
/******** Poisson correlation *************/
/******************************************/
double corr_pois_s(double rho, double mi) {
    // Gestione casi limite
    if (fabs(rho) < 1e-15) return 0.0;  
    const double rho2 = rho * rho;      
    const double gg = 2.0 * mi / (1.0 - rho2);

    if (gg < 1e-8) return rho2 * gg * (0.25 * gg + 0.5);
    // Per valori molto grandi di gg, usa approssimazione asintotica
    if (gg > 700.0) { 
        const double sqrt_2pi_gg = sqrt(2.0 * M_PI * gg);
        return rho2 * (1.0 - 2.0 / sqrt_2pi_gg);
    }
    const double exp_neg_gg = exp(-gg);
    const double i0 = bessel_i(gg, 0, 1);
    const double i1 = bessel_i(gg, 1, 1);
    
    return rho2 * (1.0 - exp_neg_gg * (i0 + i1));
}
/******************************************/
double corr_pois_ns(double rho,double mi,double mj){
int r=0;
double res0=0.0,sum=0.0;
double rho2=rho*rho;
double ki=mi/(1.0-rho2);
double kj=mj/(1.0-rho2);
double K=rho2*(1.0-rho2)/sqrt(mi*mj);
while(r<8000){
    double aa=r+1;
    sum+=exp(log(igam(aa,ki))+log(igam(aa,kj)));
    if(fabs(sum-res0)<MACHEP) break;
    res0=sum;
    r++;
}
return sum*K;
}
/******************************************/
double corr_pois(double rho,double mi,double mj){
    if( (rho>(1-1e-6)) &&  rho<=1 &&  fabs(mi-mj)<1e-320){return(1.0);}
if(fabs(rho)<1e-10){return(0.0);}
    else{
    double res=0.0;
    if (fabs(mi-mj)<MACHEP)         res=corr_pois_s(rho,mi);
    else                            res=corr_pois_ns(rho,mi,mj);
    return(res);
}
}
/*****************************************/
/*****************************************/
double corr_tukeygh(double rho,double eta,double tail){

    if(fabs(rho)<1e-16) return 0.0;
    double rho2 = rho * rho;
    double tail2 = tail * tail;
    double eta2 = eta * eta;
    double u = 1.0 - tail;
    double a = 1.0 + rho;
    double mu, cova, vari, rho1 = 0.0;

    if(fabs(eta) > MACHEP){
        double A1 = exp(a * eta2 / (1.0 - tail * a));
        double A2 = 2.0 * exp(0.5 * eta2 * (1.0 - tail * (1.0 - rho2)) / (u * u - tail2 * rho2));
        double A3 = eta2 * sqrt(u * u - rho2 * tail2);
        mu = (exp(eta2 / (2.0 * u)) - 1.0) / (eta * sqrt(u));
        cova = ((A1 - A2 + 1.0) / A3 - mu * mu);
        vari = ((exp(2.0 * eta2 / (1.0 - 2.0 * tail)) - 2.0 * exp(eta2 / (2.0 * (1.0 - 2.0 * tail))) + 1.0) / 
               (eta2 * sqrt(1.0 - 2.0 * tail)) - mu * mu);
        rho1 = cova / vari;
    } else {
        double inv_tail = R_pow(1.0 - 2.0 * tail, -1.5);
        rho1 = (-rho / ((1.0 + tail * (rho - 1.0)) * (-1.0 + tail + tail * rho) * 
                        sqrt(1.0 + tail * (-2.0 + tail - tail * rho2)))) / inv_tail;
    }

    return rho1;
}

double corr_skewt(double corr,double df,double skew)
{
if(fabs(corr)<MACHEP){return(0.0);}
    else{
    double d,w,corr1,skew2,CorSkew,nu,l,y,KK;
    skew2=skew*skew; nu=df; l=df/2;  d=(df-1)/2; w=sqrt(1-skew2);
    y=corr;
    CorSkew=(2*skew2/(M_PI*w*w+skew2*(M_PI-2)))*(sqrt(1-y*y)+y*asin(y)-1)+w*w*y/(w*w+skew2*(1-2/M_PI)) ;
   
    if(df<99){
    KK=exp(log(M_PI)+log(nu-2)+2*lgammafn(d)
    - log(2)-log(M_PI*R_pow(gammafn(l),2)*(1+skew2) - skew2*(nu-2)*R_pow(gammafn(d),2)));
    corr1=KK * (hypergeo(0.5,0.5,l,y*y)*((1+skew2*(1-2/M_PI))*CorSkew+2*skew2/M_PI) - 2*skew2/M_PI);}
    else {corr1=CorSkew;}
return(corr1);
}
}

double biv_gamma(double corr, double zi, double zj,
                 double mui, double muj, double shape)
{
    // Input validation
    if (!R_FINITE(corr) || !R_FINITE(zi) || !R_FINITE(zj) ||
        !R_FINITE(mui)  || !R_FINITE(muj) || !R_FINITE(shape) ||
        zi <= 0.0 || zj <= 0.0 || shape <= 0.0)
        return R_NegInf;
    
    const double ac = fabs(corr);
    if (ac >= 1.0) return R_NegInf;
    
    // Lavora sempre in log: lci = log(ci), lcj = log(cj)
    const double lci = log(zi) - mui;
    const double lcj = log(zj) - muj;
    
    // Controllo overflow: se lci o lcj troppo grandi, densità ~ 0
    if (lci > 700.0 || lcj > 700.0) return R_NegInf;
    
    const double lgam = lgammafn(shape / 2.0);
    const double shape_half = shape / 2.0;
    
    if (ac > 1e-10) {  // Caso correlato
        const double a = (1.0 - ac) * (1.0 + ac);  // 1 - corr²
        const double loga = log1p(-ac * ac);
        
        // z = shape * |corr| * sqrt(ci*cj) / a
        // In log: log(z) = log(shape) + log(ac) + 0.5*(lci+lcj) - log(a)
        const double logz = log(shape) + log(ac) + 0.5*(lci + lcj) - log(a);
        if (logz > 700.0) return R_NegInf;  // z troppo grande
        
        const double z = exp(logz);
        
        // A = (shape/2 - 1)*log(ci*cj) - shape*(ci+cj)/(2a)
        // log(ci*cj) = lci + lcj (già stabile)
        // ci+cj = exp(lci) + exp(lcj): usa logsumexp ma poi devi exp()
        const double log_sum_c = (lci > lcj) 
            ? lci + log1p(exp(lcj - lci))
            : lcj + log1p(exp(lci - lcj));
        
        const double A = (shape_half - 1.0) * (lci + lcj) 
                       - shape * exp(log_sum_c) / (2.0 * a);
        
        // C = (1 - shape/2) * log(z/2) = (1 - shape/2) * (logz - log(2))
        const double C = (1.0 - shape_half) * (logz - M_LN2);
        
        // B = lgamma(shape/2) + (shape/2)*log(a) + shape*log(2) - shape*log(shape)
        const double B = lgam + shape_half * loga 
                       + shape * M_LN2 - shape * log(shape);
        
        // log(I_nu(z)) con Bessel scaled
        const double I_scaled = bessel_i(z, shape_half - 1.0, 2);
        if (!R_FINITE(I_scaled) || I_scaled <= 0.0) return R_NegInf;
        const double logI = log(I_scaled) + z;
        
        return A + C + logI - mui - muj - B;
        
    } else {  // Caso indipendente (corr ≈ 0)
        // Prodotto di due gamma indipendenti in log
        // log(f_i) = (shape/2)*log(shape/2) - (shape/2)*mui 
        //          + (shape/2 - 1)*log(zi) - shape*ci/2 - lgamma(shape/2)
        const double log_shape_half = log(shape_half);
        
        const double logB = shape_half * (log_shape_half - M_LN2 - mui)
                          + (shape_half - 1.0) * log(zi)
                          - shape * exp(lci) / 2.0
                          - lgam;
        
        const double logC = shape_half * (log_shape_half - M_LN2 - muj)
                          + (shape_half - 1.0) * log(zj)
                          - shape * exp(lcj) / 2.0
                          - lgam;
        
        return logB + logC;
    }
}

double biv_gamma_gen(double corr,double zi,double zj,double mui, double muj, double shape,double n)
{
    double a=0.0,A=0.0,D=0.0,z=0.0,res=0.0,B=0.0,C=0.0;
    double ci=zi/exp(mui);double cj=zj/exp(muj);
    double gam = gammafn(shape/2);
    a=1-R_pow(corr,2);  
  if (fabs(corr) > MACHEP) {
        z=n*fabs(corr)*sqrt(ci*cj)/a;
        A=(shape/2-1)*log(ci*cj) -n*(ci+cj)/(2*a); ///ok
        C=(1-shape/2)*log(z/2); //ok
        B=log(gam)+(shape/2)*log(a)+shape*log(2)-shape*log(n);  
        D=log(bessel_i(z,shape/2-1,2))+z; //ok
        res=(A+C+D)-(mui+muj+B);
        return(exp(res));
    }
    else
    {
        B=(R_pow((n/(2*exp(mui))),shape/2)*R_pow(zi,shape/2-1)*exp(-(n*ci/(2))))/gam;
        C=(R_pow((n/(2*exp(muj))),shape/2)*R_pow(zj,shape/2-1)*exp(-(n*cj/(2))))/gam;
        res=B*C;    
    }
return(res);
}

void biv_gamma_call(double *corr,double *zi,double *zj,double *mui, double *muj, double *shape, double *res)
{
    *res = biv_gamma(*corr,*zi,*zj,*mui,*muj,*shape);
}


    
double biv_skew(double corr, double z1, double z2, double mi, double mj, double vari, double skew, double nugget) {
    double zi = z1 - mi;
    double zj = z2 - mj;

    double skew2 = skew * skew;
    double skew3 = skew2 * skew;
    double skew4 = skew2 * skew2;
    double vari2 = vari;
    double vari4 = vari2 * vari2;
    double corr1 = (1 - nugget) * corr;
    double corr2 = corr * corr;
    double corr12 = corr1 * corr1;

    double s2v2 = skew2 * vari2;
    double pp1 = skew * vari2 * corr1;
    double pp = skew * vari2 * corr;
    double ff = pp1 + pp;
    double mm = pp1 - pp;
    double aux1 = vari2 + skew2;

    double aux2 = corr1 * vari2 + corr * skew2;
    double pdf1 = d22norm(zi, zj, aux1, aux1, aux2);

    double bb = vari4 * corr12 + 2 * s2v2 * corr * corr1 + skew4 * corr2 - vari4 - 2 * s2v2 - skew4;
    double lim1 = ((mm) * zj + (pp * corr1 + skew3 * corr2 - skew * vari2 - skew3) * zi) / bb;
    double lim2 = ((mm) * zi + (pp * corr1 + skew3 * corr2 - skew * vari2 - skew3) * zj) / bb;
    double a11 = (vari4 * corr12 + s2v2 * corr2 - vari4 - s2v2) / bb;
    double a12 = (vari4 * corr * corr12 + (s2v2 * corr2 - s2v2) * corr1 - vari4 * corr) / bb;
    double cdf1 = cdf_norm(lim1, lim2, a11, a12);

    // PDF 2
    aux2 = vari2 * corr1 - skew2 * corr;
    double pdf2 = d22norm(zi, zj, aux1, aux1, aux2);

    bb = vari4 * corr12 - 2 * s2v2 * corr * corr1 + skew4 * corr2 - vari4 - 2 * s2v2 - skew4;
    lim1 = ((ff) * zj + (-pp * corr1 + skew3 * corr2 - skew * vari2 - skew3) * zi) / bb;
    lim2 = -(-(ff) * zi + (pp * corr1 - skew3 * corr2 + skew * vari2 + skew3) * zj) / bb;
    a11 = (vari4 * corr12 + s2v2 * corr2 - vari4 - s2v2) / bb;
    a12 = -(vari4 * corr * corr12 + (s2v2 - s2v2 * corr2) * corr1 - vari4 * corr) / bb;
    double cdf2 = cdf_norm(lim1, lim2, a11, a12);

    double dens = 2 * (pdf1 * cdf1 + pdf2 * cdf2);
    return dens;
}

double biv_wrapped(double alfa, double u, double v, double mi, double mj, double nugget, double sill, double corr) {
    double x = u - 2 * atan(mi) - M_PI;
    double y = v - 2 * atan(mj) - M_PI;

    double s1 = sill;
    double s12 = sill * corr;
    double det = R_pow(s1, 2) - R_pow(s12, 2);
    double wrap_gauss = 0.0;

    double k1 = -alfa, k2 = -alfa;

    double common_factor = 1 / (2 * M_PI * sqrt(det)); // Precompute common factor

    while (k1 <= alfa) {
        while (k2 <= alfa) {
            double term1 = x + 2 * k2 * M_PI;
            double term2 = y + 2 * k1 * M_PI;
            double quadr = -0.5 * (1.0 / det) * (s1 * R_pow(term2, 2) + s1 * R_pow(term1, 2) - 2.0 * s12 * term1 * term2);

            wrap_gauss += common_factor * exp(quadr);
            k2++;
        }
        k1++;
        k2 = -alfa;
    }
    return wrap_gauss;
}

double biv_skew2(double corr, double zi, double zj, double vari1, double vari2, double rho, double skew1, double skew2) {
    // Calcoli preliminari
    double taul = R_pow(vari1, 0.5);
    double taur = R_pow(vari2, 0.5);
    double c1 = 1.0 / (1 - R_pow(corr, 2));
    
    // Variabili di calcolo del determinante
    double aux1 = R_pow(taul, 2) + R_pow(skew1, 2);
    double aux2 = taul * taur * corr + skew1 * skew2 * corr;
    double aux3 = R_pow(taur, 2) + R_pow(skew2, 2);
    double det = aux1 * aux3 - R_pow(aux2, 2);

    // Calcolo quadratico e PDF 1
    double quadr = (1.0 / det) * (aux1 * R_pow(zj, 2) + aux3 * R_pow(zi, 2) - 2.0 * aux2 * zi * zj);
    double pdf1 = (0.5 / M_PI) * (1 / sqrt(det)) * exp(-0.5 * quadr);

    // Calcoli intermedi per la CDF 1
    aux1 = (R_pow(skew1 / taul, 2) + 1) * c1;
    aux2 = (skew1 * skew2) / (taul * taur) + 1 * c1 * corr;
    aux3 = (R_pow(skew2 / taur, 2) + 1) * c1;
    det = aux1 * aux3 - R_pow(aux2, 2);
    double a11 = (1 / det) * aux3;
    double a12 = (1 / det) * aux2;
    double a22 = (1 / det) * aux1;
    double factor = c1 / det;
    double fact1 = factor * (aux3 * (skew1 / R_pow(taul, 2)) - aux2 * (skew2 / (taul * taur)) * corr);
    double fact2 = factor * (-aux3 * (skew1 / (taul * taur)) * corr + aux2 * (skew2 / R_pow(taur, 2)));
    double fact3 = factor * (aux2 * (skew1 / R_pow(taul, 2)) - aux1 * (skew2 / (taul * taur)) * corr);
    double fact4 = factor * (-aux2 * (skew1 / (taul * taur)) * corr + aux1 * (skew2 / R_pow(taur, 2)));

    double lim1 = fact1 * zi + fact2 * zj;
    double lim2 = fact3 * zi + fact4 * zj;
    double cdf1 = cdf_norm2(lim1, lim2, a11, a12, a22);

    // Calcoli per il secondo PDF
    aux1 = R_pow(taul, 2) + R_pow(skew1, 2);
    aux2 = taul * taur * corr - skew1 * skew2 * corr;
    aux3 = R_pow(taur, 2) + R_pow(skew2, 2);
    det = aux1 * aux3 - R_pow(aux2, 2);
    quadr = (1.0 / det) * (aux1 * R_pow(zj, 2) + aux3 * R_pow(zi, 2) - 2.0 * aux2 * zi * zj);
    double pdf2 = (0.5 / M_PI) * (1 / sqrt(det)) * exp(-0.5 * quadr);

    // Calcoli intermedi per la CDF 2
    aux1 = R_pow(skew1 / taul, 2) * c1 + c1;
    aux2 = (skew1 * skew2) / (taul * taur) * c1 * corr - c1 * corr;
    aux3 = R_pow(skew2 / taur, 2) * c1 + c1;
    det = aux1 * aux3 - R_pow(aux2, 2);
    a11 = (1 / det) * aux3;
    a12 = (1 / det) * aux2;
    a22 = (1 / det) * aux1;
    factor = c1 / det;
    fact1 = factor * (aux3 * (skew1 / R_pow(taul, 2)) - aux2 * (skew2 / (taul * taur)) * corr);
    fact2 = factor * (-aux3 * (skew1 / (taul * taur)) * corr + aux2 * (skew2 / R_pow(taur, 2)));
    fact3 = factor * (aux2 * (skew1 / R_pow(taul, 2)) - aux1 * (skew2 / (taul * taur)) * corr);
    fact4 = factor * (-aux2 * (skew1 / (taul * taur)) * corr + aux1 * (skew2 / R_pow(taur, 2)));

    lim1 = fact1 * zi + fact2 * zj;
    lim2 = fact3 * zi + fact4 * zj;
    double cdf2 = cdf_norm2(lim1, lim2, a11, a12, a22);

    // Calcolo finale della densità
    double dens = 2 * (pdf1 * cdf1 + pdf2 * cdf2);
    return dens;
}

double log_biv_Norm(double corr, double zi, double zj, double mi, double mj,
                    double vari, double nugget)
{
    const double r = corr * (1.0 - nugget);
    if (!isfinite(vari) || vari <= 0.0 || !isfinite(r) || fabs(r) >= 1.0 ||
        !isfinite(nugget) || nugget < 0.0 || nugget > 1.0)
        return NAN;

    const double sd = sqrt(vari);
    double u = zi - mi, v = zj - mj;
    /* Avoid overflowing an otherwise representable standardized residual. */
    u = isfinite(u) ? u / sd : zi / sd - mi / sd;
    v = isfinite(v) ? v / sd : zj / sd - mj / sd;
    /* Orthogonal coordinates give a sum of nonnegative terms, avoiding
       subtraction of nearly equal terms when |r| is close to one. */
    const double inv_sqrt2 = 0.70710678118654752440;
    const double zp = (u * inv_sqrt2 + v * inv_sqrt2) / sqrt(1.0 + r);
    const double zm = (u * inv_sqrt2 - v * inv_sqrt2) / sqrt(1.0 - r);
    const double quadratic = zp * zp + zm * zm;
    const double log_det = 2.0 * log(vari) + log1p(-r) + log1p(r);
    return -0.5 * (2.0 * log(2.0 * M_PI) + log_det + quadratic);
}



/*multivariate gaussian PDF*/
double dNnorm(int N, double **M, double *dat) {
    double dd, pdf = 0.0;
    int m, n, i;
    int *indx = (int *) R_Calloc(N, int);
    double *col = (double *) R_Calloc(N, double);
    double **inverse = (double **) R_Calloc(N, double *);

    for (i = 0; i < N; i++) {
        inverse[i] = (double *) R_Calloc(N, double);
    }
    ludcmp(M, N, indx, &dd);
    dd = 1.0;
    for (m = 0; m < N; m++) {
        dd *= M[m][m];
    }
    for (n = 0; n < N; n++) {
        for (m = 0; m < N; m++) {
            col[m] = 0.0;
        }
        col[n] = 1.0;
        lubksb(M, N, indx, col);
        for (m = 0; m < N; m++) {
            inverse[m][n] = col[m];
        }
    }
    pdf = R_pow(2 * M_PI, -N / 2) * R_pow(dd, -0.5) * exp(-0.5 * QFORM(inverse, dat, dat, N));
    R_Free(indx);
    R_Free(col);
    for (i = 0; i < N; i++) {
        R_Free(inverse[i]);
    }
    R_Free(inverse);
    return pdf;
}

static inline double inv_tukey_h(double z, double h) {
    if (h == 0.0) return z;
    const double arg = h * z * z;          
    double W = LambertW(arg);                  
    if (!isfinite(W) || W < 0.0) return NAN;  
    double u = copysign(sqrt(W / h), z);
    return u;
}

/* |du/dz| = exp(-0.5*W(h*z^2)) / (1 + W(h*z^2)),   h=0 -> 1 */
static inline double jac_tukey_h(double z, double h) {
    if (h == 0.0) return 1.0;
    const double arg = h * z * z;
    double W = LambertW(arg);
    if (!isfinite(W) || (1.0 + W) <= 0.0) return NAN;  // evita divisione per zero
    return exp(-0.5 * W) / (1.0 + W);
}

/* Bivariata Tukey-hh (Tukey-h a tratti: hr per z>=0, hl per z<0) */
double biv_tukey_hh(double rho,
                    double yi, double yj,
                    double mui, double muj,
                    double sill, double hl, double hr)
{
    /* Controlli base */
    if (!isfinite(rho)||!isfinite(yi)||!isfinite(yj)||
        !isfinite(mui)||!isfinite(muj)||!isfinite(sill)||
        !isfinite(hl) ||!isfinite(hr)) return NAN;

    if (sill <= 0.0 || hl < 0.0 || hr < 0.0 || hl >= 0.5 || hr >= 0.5) return NAN;

    /* Trasformazione y -> z = (y - mu)/sqrt(sill) */
    const double sqrt_sill = sqrt(sill);
    const double zi = (yi - mui) / sqrt_sill;
    const double zj = (yj - muj) / sqrt_sill;

    /* Scegli il parametro h a destra/sinistra */
    const double hi = (zi >= 0.0 ? hr : hl);
    const double hj = (zj >= 0.0 ? hr : hl);

    /* Inversa della trasformazione Tukey-h e jacobiani */
    const double ui = inv_tukey_h(zi, hi);
    const double uj = inv_tukey_h(zj, hj);
    if (!isfinite(ui) || !isfinite(uj)) return NAN;

    const double Ji = jac_tukey_h(zi, hi);
    const double Jj = jac_tukey_h(zj, hj);
    if (!isfinite(Ji) || !isfinite(Jj)) return NAN;

    /* Salvaguardia correlazione in (-1,1) */
    if (rho <= -1.0) rho = -1.0 + 1e-12;
    if (rho >=  1.0) rho =  1.0 - 1e-12;

    /* Densità normale bivariata standard in (ui,uj) con correlazione rho */
    double phi2 = dbnorm(ui, uj, 0.0, 0.0, 1.0, rho);
    if (!isfinite(phi2)) return NAN;

    /* Trasformazioni concatenate:
       f_Z(z_i,z_j) = phi2(ui,uj;r) * |du_i/dz_i| * |du_j/dz_j|
       f_Y(y_i,y_j) = f_Z(z_i,z_j) * |dz_i/dy_i| * |dz_j/dy_j| = ... / sill
    */
    double dens = (phi2 * Ji * Jj) / sill;

    return isfinite(dens) ? dens : NAN;
}



#define LOG_TOL   (-16.0)
#define EXP_CUTOFF (-600.0)
#define MAX_STREAK 4
#define MAX_TERMS 8192
#define KMAX_DEFAULT 50
#define KMAX_MAX 50
#define CONV_TOL 1e-6
#define CONV_STREAK_MIN 1

static inline double logsumexp_r_style(const double *log_terms, int n) {
    if (n <= 0) return R_NegInf;

    double m = R_NegInf;
    for (int i = 0; i < n; i++) {
        if (R_FINITE(log_terms[i]) && log_terms[i] > m) {
            m = log_terms[i];
        }
    }
    if (!R_FINITE(m)) return R_NegInf;

    double sum_exp = 0.0;
    for (int i = 0; i < n; i++) {
        if (R_FINITE(log_terms[i])) {
            double diff = log_terms[i] - m;
            if (diff > EXP_CUTOFF) {
                sum_exp += exp(diff);
            }
        }
    }
    // MIGLIORIA: Check su sum_exp con tolleranza numerica
    if (sum_exp <= 0.0) return R_NegInf;
    return m + log(sum_exp);
}
static inline double lchoose_fast_local(int n, int k, const double *lfact) {
    if (k < 0 || k > n) return R_NegInf;
    if (k == 0 || k == n) return 0.0;
    return lfact[n] - lfact[k] - lfact[n - k];
}

static double Ikell_log_r_style_local(int k, int ell, double z, double alpha,
                                      double log_alpha, const double *lfact) {
    int max_terms = (z <= 0.0) ? (ell + 1) : (k + 1);
    double *terms_log = (double*) R_alloc(max_terms, sizeof(double));
    
    for (int i = 0; i < max_terms; i++) terms_log[i] = R_NegInf;

    // MIGLIORIA: Gestione più robusta di z == 0
    double abs_z = fabs(z);
    const double Z_TOL = 1e-15;  // Tolleranza per considerare z ≈ 0

    if (z <= 0.0) {
        for (int m = 0; m <= ell; m++) {
            // Se z è vicino a 0 e ell-m > 0, il termine è trascurabile
            if (abs_z < Z_TOL && (ell - m) > 0) continue;

            double log_binom = lchoose_fast_local(ell, m, lfact);
            if (!R_FINITE(log_binom)) continue;

            double log_pow = 0.0;
            if (abs_z >= Z_TOL || (ell - m) == 0) {
                log_pow = (ell - m) * (abs_z < Z_TOL ? R_NegInf : log(abs_z));
            }
            
            double log_gamma = lfact[k + m];
            double log_alpha_term = -(k + m + 1.0) * log_alpha;

            terms_log[m] = log_binom + log_pow + log_gamma + log_alpha_term;
        }
        return logsumexp_r_style(terms_log, ell + 1);
    } else {
        for (int m = 0; m <= k; m++) {
            if (abs_z < Z_TOL && (k - m) > 0) continue;

            double log_binom = lchoose_fast_local(k, m, lfact);
            if (!R_FINITE(log_binom)) continue;

            double log_pow = 0.0;
            if (abs_z >= Z_TOL || (k - m) == 0) {
                log_pow = (k - m) * (abs_z < Z_TOL ? R_NegInf : log(z));
            }
            
            double log_gamma = lfact[ell + m];
            double log_alpha_term = -(ell + m + 1.0) * log_alpha;

            terms_log[m] = log_binom + log_pow + log_gamma + log_alpha_term;
        }
        double log_sum = logsumexp_r_style(terms_log, k + 1);
        if (!R_FINITE(log_sum)) return R_NegInf;
        return -alpha * z + log_sum;
    }
}

double biv_skewlaplace_kmax(double z1, double z2, double p, double rho, int Kmax) {
    if (p <= 0.0 || p >= 1.0 || fabs(rho) >= 1.0) return R_NegInf;
    if (!R_FINITE(z1) || !R_FINITE(z2)) return R_NegInf;
    if (Kmax < 0) Kmax = KMAX_DEFAULT;
    if (Kmax > KMAX_MAX) Kmax = KMAX_MAX;

    double Delta = 1.0 - rho * rho;
    // MIGLIORIA: Tolleranza più appropriata per Delta
    if (Delta <= 1e-14) return R_NegInf;

    double alpha = 1.0 / Delta;
    double log_p = log(p);
    double log_1mp = log1p(-p);
    double log_Delta = log(Delta);
    double log_abs_rho = log(fabs(rho));
    double log_alpha = log(alpha);

    double log_pref = 2.0 * log_p + 2.0 * log_1mp - 2.0 * log_Delta +
                      ((1.0 - p) * (z1 + z2)) / Delta;
    if (!R_FINITE(log_pref)) return R_NegInf;

    int max_fact_needed = 2 * Kmax + 2;
    double *lfact = (double*) R_alloc(max_fact_needed + 1, sizeof(double));
    for (int i = 0; i <= max_fact_needed; i++) {
        lfact[i] = lgammafn(i + 1.0);
    }

    int max_possible_terms = (Kmax + 1) * (Kmax + 1);
    double *log_terms_total = (double*) R_alloc(max_possible_terms, sizeof(double));
    int nterms = 0;

    for (int k = 0; k <= Kmax; k++) {
        for (int ell = 0; ell <= Kmax; ell++) {
            double log_coef = -2.0 * lfact[k] - 2.0 * lfact[ell] +
                              k * (2.0 * log_abs_rho + 2.0 * log_p - 2.0 * log_Delta) +
                              ell * (2.0 * log_abs_rho + 2.0 * log_1mp - 2.0 * log_Delta);

            // MIGLIORIA: Check di validità prima di calcolare
            if (!R_FINITE(log_coef)) continue;

            double log_I1 = Ikell_log_r_style_local(k, ell, z1, alpha, log_alpha, lfact);
            double log_I2 = Ikell_log_r_style_local(k, ell, z2, alpha, log_alpha, lfact);

            if (R_FINITE(log_I1) && R_FINITE(log_I2)) {
                double total_term = log_coef + log_I1 + log_I2;
                if (R_FINITE(total_term) && nterms < max_possible_terms) {
                    log_terms_total[nterms++] = total_term;
                }
            }
        }
    }
    
    if (nterms == 0) return R_NegInf;
    
    double log_density = log_pref + logsumexp_r_style(log_terms_total, nterms);
    return R_FINITE(log_density) ? log_density : R_NegInf;
}

double log_biv_skewlaplace(double rho, double z1, double z2, double m1, double m2, 
                           double sill, double p) {
    // MIGLIORIA: Validazione parametri all'inizio
    if (sill <= 0.0) return R_NegInf;
    if (p <= 0.0 || p >= 1.0) return R_NegInf;
    if (fabs(rho) >= 1.0) return R_NegInf;
    if (!R_FINITE(z1) || !R_FINITE(z2) || !R_FINITE(m1) || !R_FINITE(m2)) 
        return R_NegInf;

    // MIGLIORIA: Tolleranza più appropriata per rho ≈ 0
    const double RHO_TOL = 1e-10;
    
    if (fabs(rho) < RHO_TOL) {
        // Caso indipendente: somma di log-densità univariate
        return (one_log_SkewLaplace(z1, m1, sill, p) + 
                one_log_SkewLaplace(z2, m2, sill, p));
    }

    // Standardizzazione
    double sqrt_sill = sqrt(sill);
    z1 = (z1 - m1) / sqrt_sill;
    z2 = (z2 - m2) / sqrt_sill;

    double log_density_prev = R_NegInf;
    double log_density_curr = R_NegInf;
    int convergence_streak = 0;

    // MIGLIORIA: Inizio da Kmax più piccolo per casi semplici
    int Kmax_start = (fabs(z1) < 2.0 && fabs(z2) < 2.0) ? 5 : 10;
    int Kmax_step = 5;  // Incrementi più piccoli per convergenza più fine

    for (int Kmax = Kmax_start; Kmax <= KMAX_MAX; Kmax += Kmax_step) {
        log_density_curr = biv_skewlaplace_kmax(z1, z2, p, rho, Kmax);

        if (!R_FINITE(log_density_curr)) {
            log_density_prev = log_density_curr;
            continue;
        }

        if (Kmax > Kmax_start && R_FINITE(log_density_prev)) {
            double abs_diff = fabs(log_density_curr - log_density_prev);
            double rel_diff = abs_diff / (1.0 + fabs(log_density_curr));
            
            if (rel_diff < CONV_TOL) {
                convergence_streak++;
                if (convergence_streak >= CONV_STREAK_MIN) {
                    return (log_density_curr - log(sill));
                }
            } else {
                convergence_streak = 0;
            }
        }
        log_density_prev = log_density_curr;
    }

    // MIGLIORIA: Ritorna valore finale anche se non convergente
    if (R_FINITE(log_density_curr)) {
        return (log_density_curr - log(sill));
    }
    return R_NegInf;
}



static double logspace_add_binomneg(double logx, double logy)
{
    if (ISNAN(logx) || ISNAN(logy))
        return R_NaN;

    if (logx == R_NegInf)
        return logy;

    if (logy == R_NegInf)
        return logx;

    if (logx == R_PosInf || logy == R_PosInf)
        return R_PosInf;

    {
        double m = fmax(logx, logy);
        return m + log1p(exp(-fabs(logx - logy)));
    }
}


/*
 * Returns log(p^exponent), with the convention 0^0 = 1.
 * The exponents appearing in the bivariate pmf must be non-negative.
 */
static double log_prob_power_binomneg(double p, int exponent)
{
    if (exponent < 0 || !R_FINITE(p) || p < 0.0 || p > 1.0)
        return R_NaN;

    if (exponent == 0)
        return 0.0;

    if (p == 0.0)
        return R_NegInf;

    return ((double) exponent) * log(p);
}


/*
 * Add a log-term to a log-sum.  Terms equal to -Inf correspond to an
 * exact zero contribution and can be skipped.  NaN/+Inf signal failure.
 */
static int add_log_term_binomneg(double *logsum, double logterm)
{
    if (ISNAN(logterm) || logterm == R_PosInf)
        return 0;

    if (logterm != R_NegInf)
        *logsum = logspace_add_binomneg(*logsum, logterm);

    return 1;
}


/**********************************************************************/
/* Log-pmf for u < v: h1 in the paper                                  */
/**********************************************************************/

static double log_aux_biv_binomneg(
    int NN,
    int u,
    int v,
    double x,
    double y,
    double p11,
    const double *lgamma_cache,
    int cache_size)
{
    int a, i;
    double logdens1 = R_NegInf;
    double logdens2 = R_NegInf;

    double base_prob   = 1.0 + p11 - x - y;
    double diff_x      = x - p11;
    double diff_y      = y - p11;
    double one_minus_y = 1.0 - y;

    double lgamma_NN_u = lgammafn(NN + u);
    double lgamma_v_u  = lgammafn(v - u);

#define LG_BN(k) \
    (((k) >= 0 && (k) < cache_size) ? \
        lgamma_cache[(k)] : lgammafn((k) + 1.0))

    if (NN <= 0 || u < 0 || v <= u)
        return R_NaN;

    /*
     * The Frechet constraints imply that these four probabilities are
     * non-negative.  Allow only tiny round-off violations.
     */
    if (base_prob < -1e-12 ||
        diff_x < -1e-12 ||
        diff_y < -1e-12 ||
        one_minus_y < -1e-12)
        return R_NaN;

    base_prob   = fmax(0.0, base_prob);
    diff_x      = fmax(0.0, diff_x);
    diff_y      = fmax(0.0, diff_y);
    one_minus_y = fmax(0.0, one_minus_y);

    /* First sum in h1. */
    for (a = fmax_int(0, u - v + NN - 1);
         a <= NN - 2;
         a++) {

        int exp_y           = NN - a - 1;
        int exp_one_minus_y = v - u - NN + a + 1;

        double log_y =
            log_prob_power_binomneg(y, exp_y);

        double log_one_minus_y =
            log_prob_power_binomneg(one_minus_y, exp_one_minus_y);

        if (ISNAN(log_y) || ISNAN(log_one_minus_y))
            return R_NaN;

        for (i = fmax_int(0, a - u);
             i <= fmin_int(a, NN - 1);
             i++) {

            int exp_p11  = i + 1;
            int exp_base = u - a + i;
            int exp_dx   = NN - i - 1;
            int exp_dy   = a - i;

            double log_coeff1 =
                lgamma_NN_u
                - LG_BN(i)
                - LG_BN(NN - i - 1)
                - LG_BN(a - i)
                - lgammafn(u - a + i + 1.0);

            double log_coeff2 =
                lgamma_v_u
                - lgammafn(v + a - NN - u + 2.0)
                - LG_BN(NN - a - 2);

            double logterm =
                log_coeff1
                + log_coeff2
                + log_prob_power_binomneg(p11, exp_p11)
                + log_prob_power_binomneg(base_prob, exp_base)
                + log_prob_power_binomneg(diff_x, exp_dx)
                + log_prob_power_binomneg(diff_y, exp_dy)
                + log_one_minus_y
                + log_y;

            if (!add_log_term_binomneg(&logdens1, logterm))
                return R_NaN;
        }
    }

    /* Second sum in h1. */
    for (a = fmax_int(0, u - v + NN);
         a <= NN - 1;
         a++) {

        int exp_y           = NN - a;
        int exp_one_minus_y = v - u - NN + a;

        double log_y =
            log_prob_power_binomneg(y, exp_y);

        double log_one_minus_y =
            log_prob_power_binomneg(one_minus_y, exp_one_minus_y);

        if (ISNAN(log_y) || ISNAN(log_one_minus_y))
            return R_NaN;

        for (i = fmax_int(0, a - u);
             i <= fmin_int(a, NN - 1);
             i++) {

            int exp_p11  = i;
            int exp_base = u - a + i;
            int exp_dx   = NN - i;
            int exp_dy   = a - i;

            double log_coeff1 =
                lgamma_NN_u
                - LG_BN(i)
                - LG_BN(NN - i - 1)
                - LG_BN(a - i)
                - lgammafn(u - a + i + 1.0);

            double log_coeff2 =
                lgamma_v_u
                - lgammafn(v + a - NN - u + 1.0)
                - LG_BN(NN - a - 1);

            double logterm =
                log_coeff1
                + log_coeff2
                + log_prob_power_binomneg(p11, exp_p11)
                + log_prob_power_binomneg(base_prob, exp_base)
                + log_prob_power_binomneg(diff_x, exp_dx)
                + log_prob_power_binomneg(diff_y, exp_dy)
                + log_one_minus_y
                + log_y;

            if (!add_log_term_binomneg(&logdens2, logterm))
                return R_NaN;
        }
    }

#undef LG_BN

    return logspace_add_binomneg(logdens1, logdens2);
}


/**********************************************************************/
/* Log-pmf for u = v: h2 in the paper                                  */
/**********************************************************************/

static double log_aux_biv_binomneg_simple(
    int NN,
    int u,
    double p1,
    double p2,
    double p11,
    const double *lgamma_cache,
    int cache_size)
{
    int i;
    double logdens = R_NegInf;

    double p00 = 1.0 + p11 - p1 - p2;
    double p10 = p1 - p11;
    double p01 = p2 - p11;

    double lgamma_NN_u = lgammafn(NN + u);

#define LG_BN_SIMPLE(k) \
    (((k) >= 0 && (k) < cache_size) ? \
        lgamma_cache[(k)] : lgammafn((k) + 1.0))

    if (NN <= 0 || u < 0)
        return R_NaN;

    if (p00 < -1e-12 || p10 < -1e-12 || p01 < -1e-12)
        return R_NaN;

    p00 = fmax(0.0, p00);
    p10 = fmax(0.0, p10);
    p01 = fmax(0.0, p01);

    for (i = fmax_int(0, NN - u - 1);
         i <= NN - 1;
         i++) {

        int exp_p11 = i + 1;
        int exp_p00 = u - NN + 1 + i;
        int exp_p10 = NN - i - 1;
        int exp_p01 = NN - i - 1;

        double log_coeff =
            lgamma_NN_u
            - LG_BN_SIMPLE(i)
            - 2.0 * LG_BN_SIMPLE(NN - i - 1)
            - lgammafn(u - NN + i + 2.0);

        double logterm =
            log_coeff
            + log_prob_power_binomneg(p11, exp_p11)
            + log_prob_power_binomneg(p00, exp_p00)
            + log_prob_power_binomneg(p10, exp_p10)
            + log_prob_power_binomneg(p01, exp_p01);

        if (!add_log_term_binomneg(&logdens, logterm))
            return R_NaN;
    }

#undef LG_BN_SIMPLE

    return logdens;
}



/**********************************************************************/
/* Stable log-pmf dispatcher                                           */
/**********************************************************************/

double log_biv_binomneg(
    int NN,
    int u,
    int v,
    double p1,
    double p2,
    double p11,
    const double *lgamma_cache,
    int cache_size)
{
    double lower_p11;
    double upper_p11;

    if (NN <= 0 || u < 0 || v < 0)
        return R_NaN;

    if (!R_FINITE(p1) || !R_FINITE(p2) || !R_FINITE(p11))
        return R_NaN;

    if (p1 < 0.0 || p1 > 1.0 ||
        p2 < 0.0 || p2 > 1.0)
        return R_NaN;

    lower_p11 = fmax(0.0, p1 + p2 - 1.0);
    upper_p11 = fmin(p1, p2);

    /*
     * Reject genuine Frechet violations, but correct tiny numerical
     * deviations produced by the bivariate Gaussian cdf evaluation.
     */
    if (p11 < lower_p11 - 1e-12 ||
        p11 > upper_p11 + 1e-12)
        return R_NaN;

    p11 = fmin(upper_p11, fmax(lower_p11, p11));

    if (u < v) {
        return log_aux_biv_binomneg(
            NN, u, v, p1, p2, p11,
            lgamma_cache, cache_size
        );
    }

    if (u == v) {
        return log_aux_biv_binomneg_simple(
            NN, u, p1, p2, p11,
            lgamma_cache, cache_size
        );
    }

    return log_aux_biv_binomneg(
        NN, v, u, p2, p1, p11,
        lgamma_cache, cache_size
    );
}


/**********************************************************************/
/* Compatibility wrapper: returns the probability, as before           */
/**********************************************************************/

double biv_binomneg(
    int NN,
    int u,
    int v,
    double p1,
    double p2,
    double p11,
    double *lgamma_cache,
    int cache_size)
{
    double logp = log_biv_binomneg(
        NN, u, v, p1, p2, p11,
        lgamma_cache, cache_size
    );

    if (ISNAN(logp))
        return R_NaN;

    if (logp == R_NegInf || logp < log(DBL_MIN))
        return 0.0;

    if (logp == R_PosInf)
        return R_PosInf;

    return exp(logp);
}

/**********************************************************************/

void biv_binom_call(int *NN,int *u, int *v, double *p01, double *p10,double *p11,double *res)
{
    *res = biv_binom (*NN,*u,*v,*p01, *p10,*p11);
}


int check_biv_binom_params(int NN, int u, int v, double p1, double p2, double p11) {
    // Check basic bounds
    if (u < 0 || v < 0 || u > NN || v > NN) return 0;
    if (NN < 0) return 0;
    
    // Check probability bounds
    if (p1 < 0.0 || p1 > 1.0) return 0;
    if (p2 < 0.0 || p2 > 1.0) return 0;
    if (p11 < 0.0 || p11 > 1.0) return 0;
    
    // Check coherence constraints
    if (p1 < p11 || p2 < p11) return 0;  // p11 cannot exceed marginals
    if (1.0 + p11 - (p1 + p2) < 0.0) return 0;  // Double failure probability must be non-negative
    
    return 1;  // All checks passed
}

// Core bivariate binomial function for equal sample sizes
double biv_binom(int NN, int u, int v, double p1, double p2, double p11) {   
    int a;
    double kk = 0.0, dens = 0.0;
    
    // Parameter validation
    //if (!check_biv_binom_params(NN, u, v, p1, p2, p11)) {
    //    return 0.0;}
    
    // Edge case: if NN = 0
    if (NN == 0) {
        return (u == 0 && v == 0) ? 1.0 : 0.0;
    }
    
    double lgamma_NN_1 = lgammafn(NN + 1);
    
    for (a = fmax_int(0, u + v - NN); a <= fmin_int(u, v); a++) {
        if (u - a < 0 || v - a < 0 || NN - u - v + a < 0) continue;
        double lgamma_a1 = lgammafn(a + 1);
        double lgamma_u_a1 = lgammafn(u - a + 1);
        double lgamma_v_a1 = lgammafn(v - a + 1);
        double lgamma_NN_uv_a1 = lgammafn(NN - u - v + a + 1);

        kk = exp(lgamma_NN_1 - (lgamma_a1 + lgamma_u_a1 + lgamma_v_a1 + lgamma_NN_uv_a1));
        
        double term1 = (a == 0) ? 1.0 : R_pow(p11, a);
        double term2 = (u - a == 0) ? 1.0 : R_pow(p1 - p11, u - a);
        double term3 = (v - a == 0) ? 1.0 : R_pow(p2 - p11, v - a);
        double term4 = (NN - u - v + a == 0) ? 1.0 : R_pow(1.0 + p11 - (p1 + p2), NN - u - v + a);

        if (!R_FINITE(term1) || !R_FINITE(term2) || !R_FINITE(term3) || !R_FINITE(term4)) {
            continue;
        }
        dens += kk * term1 * term2 * term3 * term4;
    }
    return dens;
}

// Auxiliary function for different sample sizes
double aux_biv_binom(int n1, int n2, int u, int v, double p1, double p2, double p11) {
    int k;
    double dens = 0.0;
    int N = n1 - n2;
    
    // Basic parameter checks
    if (n1 <= 0 || n2 <= 0 || n1 <= n2) return 0.0;
    if (u < 0 || v < 0) return 0.0;
    if (p1 < 0.0 || p1 > 1.0) return 0.0;
    
    for (k = 0; k <= N; k++) {
        // Skip if u-k would be negative (invalid case)
        if (u - k < 0) continue;
        
        // Skip if u-k > n2 (impossible case)
        if (u - k > n2) continue;
        
        // Calculate binomial coefficient and probability for additional trials
        double log_binom_coeff = lgammafn(N + 1) - lgammafn(k + 1) - lgammafn(N - k + 1);
        double log_prob_k = k * log(p1) + (N - k) * log1p(-p1);
        
        // Calculate bivariate binomial for remaining trials
        double biv_prob = biv_binom(n2, u - k, v, p1, p2, p11);
        
        // Check for numerical issues
        if (!R_FINITE(biv_prob) || biv_prob <= 0.0) continue;
        
        double log_biv_prob = log(biv_prob);
        double log_term = log_binom_coeff + log_prob_k + log_biv_prob;
        
        // Check for numerical overflow/underflow
        if (R_FINITE(log_term)) {
            dens += exp(log_term);
        }
    }
    
    return dens;   
}

// Main bivariate binomial function handling all cases
double biv_binom222(int n1, int n2, int u, int v, double p1, double p2, double p11) {
    double res = 0.0;
    
    // Basic input validation
    if (n1 <= 0 || n2 <= 0) return 0.0;
    if (u < 0 || v < 0) return 0.0;
    if (u > n1 || v > n2) return 0.0;  // Impossible outcomes
    
    // Parameter validation for probabilities
    if (!check_biv_binom_params(fmax_int(n1, n2), u, v, p1, p2, p11)) {
        return 0.0;
    }
    
    if (n1 > n2) {
        res = aux_biv_binom(n1, n2, u, v, p1, p2, p11);
    } else if (n2 > n1) {
        res = aux_biv_binom(n2, n1, v, u, p2, p1, p11);
    } else {  // n1 == n2
        res = biv_binom(n1, u, v, p1, p2, p11);
    }
    
    // Final sanity check
    if (!R_FINITE(res) || res < 0.0) {
        return 0.0;
    }
    
    return res;
}



// bivariate pois-binomial 
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0, kk=0.0;
    a_u=fmin(u,v);
    dens=0;
    for(a=0;a<=a_u;a++){

    kk=exp(-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)));
        dens+=kk*(R_pow(NN*p11,a)*R_pow(NN*(p01-p11),u-a)*R_pow(NN*(p10-p11),v-a));
    }
    return(exp(-NN*(p01+p10-p11))*dens);
}

// bivariate pois-bineg 
double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0,pp=0.0,bb=0.0, kk=0.0;
    a_u=fmin_int(u,v);
    dens=0;
    pp=p11*p01*p10;
    bb=(1-p11)/p11;
    for(a=0;a<=a_u;a++){

       kk=exp(-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)));
       dens+=kk*(R_pow(NN*((p11*(p01+p10)-p01*p10*(p11+1))/pp),a)*R_pow(NN*((p10-p11)/(p10*p11)),u-a)*R_pow(NN*bb,v-a));
    }
    return(exp(-NN*bb)*dens);
}

// marginal binomial
//marginal geom
//marginal poisson
       
/* cdf inflated poisson*/
double ppoisinflated(double z, double m, double p)
{
 double res;
 double ans = ppois(z, m,1,0);
 if(z < 0) res=0;
 else     res=p + (1 - p) * ans;
 return(res);
}
/* cdf inflated binom neg*/
double pbneginflated(double z, int N,double m, double p)
{
 double res;   
 double ans = pnbinom(z,N,m,1,0);
 if(z < 0) res=0;
 else     res=p + (1 - p) * ans;
 return(res);
}



/**********************************************************/
double lgam_sgn(double x, int *sign)
{
    double p, q, u, w, z;
    int i;
    
    *sign = 1;
    
    if (!R_FINITE(x))
        return x;
    
    if (x < -34.0) {
        q = -x;
        w = lgam_sgn(q, sign);
        p = floor(q);
        if (p == q) {
        lgsing:
            //mtherr("lgam", SING);
           // printf("SIGN\n");
            return (NPY_INFINITY);
        }
        i = p;
        if ((i & 1) == 0)
            *sign = -1;
        else
            *sign = 1;
        z = q - p;
        if (z > 0.5) {
            p += 1.0;
            z = p - q;
        }
        z = q * sin(NPY_PI * z);
        if (z == 0.0)
            goto lgsing;

        z = LOGPI - log(z) - w;
        return (z);
    }
    
    if (x < 13.0) {
        z = 1.0;
        p = 0.0;
        u = x;
        while (u >= 3.0) {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while (u < 2.0) {
            if (u == 0.0)
                goto lgsing;
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if (z < 0.0) {
            *sign = -1;
            z = -z;
        }
        else
            *sign = 1;
        if (u == 2.0)
            return (log(z));
        p -= 2.0;
        x = x + p;
        p = x * polevl(x, B, 5) / p1evl(x, C, 6);
        return (log(z) + p);
    }
    
    if (x > MAXLGM) {
        return (*sign * NPY_INFINITY);
    }
    
    q = (x - 0.5) * log(x) - x + LS2PI;
    if (x > 1.0e8)
        return (q);
    
    p = 1.0 / (x * x);
    if (x >= 1000.0)
        q += ((7.9365079365079365079365e-4 * p
               - 2.7777777777777777777778e-3) * p
              + 0.0833333333333333333333) / x;
    else
        q += polevl(p, A, 4) / x;
    return (q);
}

/***********************************************************/
/*********** bivariate T distribution*********************/
double biv_T(double rho, double zi, double zj, double nuu, double nugget)
{
    int k = 0;
    double nu = 1 / nuu;
    double res0 = 0.0, RR = 0.0, pp1 = 0.0, pp2 = 0.0;
    double r2 = rho * rho;
    double bb1, bb2;
    double x = zi, y = zj;
    double cc = (nu + 1) / 2;
    double nu2 = nu / 2;
    double rho1 = rho * (1 - nugget);
    double rn = rho1 * rho1;
    double rho2 = pow1p(-r2, -nu2 - 1); // pow1p(x, y) computes (1 + x)^y
    double rho12 = pow1p(-rho1 * rho1, -nu - 0.5);
    double x1 = (x * x * (1 - r2) + nu * (1 - rn));
    double y1 = (y * y * (1 - r2) + nu * (1 - rn));
    double C, B;
    double b1 = exp(nu * log(nu) - cc * log(x1 * y1) + 2 * lgammafn(cc));
    double c1 = exp(log(M_PI) + 2 * lgammafn(nu / 2) + log(rho12) + log(rho2));
    double b2 = rho1 * x * y * R_pow(nu, nu + 2) * R_pow(x1 * y1, -nu2 - 1);
    double c2 = 2 * M_PI * pow1p(-rn, -nu - 0.5) * pow1p(-r2, -nu2 - 2);
    double a1 = 0, a2 = 0;
    double aux = R_pow(rho1 * x * y, 2) * pow1p(-r2, 2) / (x1 * y1);
    double aux1 = R_pow(rho * nu, 2) * pow1p(-rn, 2) / (x1 * y1);
    if (rho > DEPSILON)
    {
        while (k <= 3000)
        {
            // pp1 = hypergeo(cc + k, cc + k, 0.5, aux);
            pp1 = (0.5 - 2 * (cc + k)) * log1p(-aux) + log(hypergeo(0.5 - (cc + k), 0.5 - (cc + k), 0.5, aux)); // Euler
            bb1 = pp1 + k * log(aux1) + 2 * (lgammafn(cc + k) - lgammafn(cc)) - lgammafn(k + 1) - lgammafn(nu2 + k) + lgammafn(nu2);
            a1 += exp(bb1);
            // pp2 = hypergeo(nu2 + 1 + k, nu2 + 1 + k, 1.5, aux);
            pp2 = (1.5 - 2 * (nu2 + 1 + k)) * log1p(-aux) + log(hypergeo(1.5 - (nu2 + 1 + k), 1.5 - (nu2 + 1 + k), 1.5, aux)); // Euler
            bb2 = pp2 + k * log(aux1) + 2 * log((1 + k / nu2)) + lgammafn(nu2 + k) - lgammafn(k + 1) - lgammafn(nu2);
            a2 += exp(bb2);
            RR = (b1 / c1) * a1 + (b2 / c2) * a2;
            if (fabs(RR - res0) < 1e-10 || !R_FINITE(RR)) { break; }
            else { res0 = RR; }
            k++;
        }
        if (!R_finite(RR)) RR = 1e-320;
    }

    if (rho <= DEPSILON)
    {
        C = lgammafn(cc) + log(R_pow((1 + x * x / nu), -cc)) - log(sqrt(M_PI * nu)) - lgammafn(nu / 2);
        B = lgammafn(cc) + log(R_pow((1 + y * y / nu), -cc)) - log(sqrt(M_PI * nu)) - lgammafn(nu / 2);
        RR = exp(B) * exp(C);
    }

    return RR;
}

/*********** Appell F4 function ********/

void appellF4_call(double *a,double *b,double *c,double *d,double *x,double *y, double *res)
{
    *res = appellF4(*a,*b,*c,*d,*x,*y);
}

double appellF4(double a, double b, double c, double d, double x, double y)
{
    // Controlli di input
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(c) || !R_FINITE(d) || 
        !R_FINITE(x) || !R_FINITE(y)) {
        return R_NaN;
    }
    
    // Controlli per casi speciali
    if (fabs(y) < DBL_EPSILON) return hypergeo(a, b, c, x);
    if (fabs(x) >= 1.0 || fabs(y) >= 1.0) {
        // Fuori dal raggio di convergenza
        return R_NaN;
    }
    
    double sum = 0.0;
    double term = 0.0;
    double log_y = log(fabs(y));
    double log1p_neg_x = log1p(-x);
    
    // Pre-calcolo delle funzioni gamma invarianti
    double lgamma_a = lgammafn(a);
    double lgamma_b = lgammafn(b);
    double lgamma_d = lgammafn(d);
    
    // Parametri per il controllo della convergenza
    const int MAX_ITER = 1000;
    const double TOLERANCE = 1e-12;
    const double MIN_TERM = 1e-15;
    
    // Variabili per il controllo adattivo
    double prev_sum = 0.0;
    int stable_count = 0;
    const int STABLE_THRESHOLD = 3;
    
    for (int k = 0; k < MAX_ITER; k++) {
        // Calcolo logaritmico del termine per stabilità numerica
        double log_term = k * log_y + 
                         lgammafn(a + k) + lgammafn(b + k) + lgamma_d -
                         lgamma_a - lgamma_b - lgammafn(d + k) - lgammafn(k + 1) +
                         (c - a - k - b - k) * log1p_neg_x;
        
        // Controllo overflow/underflow
        if (log_term < -700.0) {  // exp(-700) ≈ 1e-304
            break;  // I termini successivi saranno trascurabili
        }
        if (log_term > 700.0) {   // exp(700) ≈ 1e+304
            return R_PosInf;  // Serie divergente
        }
        
        term = exp(log_term);
        
        // Calcolo della funzione ipergeometrica
        double hyp_val = hypergeo(c - a - k, c - b - k, c, x);
        if (!R_FINITE(hyp_val)) {
            if (k == 0) return R_NaN;
            break;  // Usa la somma parziale calcolata finora
        }
        
        term *= hyp_val;
        
        // Controllo convergenza
        if (fabs(term) < MIN_TERM) break;
        
        prev_sum = sum;
        sum += term;
        
        // Controllo stabilità della somma
        if (fabs(sum - prev_sum) < TOLERANCE * fabs(sum)) {
            stable_count++;
            if (stable_count >= STABLE_THRESHOLD) break;
        } else {
            stable_count = 0;
        }
        
        // Controllo per termine relativo
        if (k > 10 && fabs(term) < TOLERANCE * fabs(sum)) {
            break;
        }
    }
    
    if (!R_FINITE(sum)) {
        if (sum > 0) return R_PosInf;
        if (sum < 0) return R_NegInf;
        return R_NaN;
    }
    if (fabs(sum) < DBL_MIN) {
        return (sum >= 0) ? DBL_MIN : -DBL_MIN;
    }  
    return sum;
}



double appellF4_mod(double nu, double rho, double x, double y, double nugget)
{
    // Controlli di input
    if (!R_FINITE(nu) || !R_FINITE(rho) || !R_FINITE(x) || !R_FINITE(y) || !R_FINITE(nugget)) {
        return R_NaN;
    }
    
    if (nu <= 0) return R_NaN;
    if (fabs(rho) >= 1.0) return R_NaN;
    if (nugget < 0 || nugget >= 1.0) return R_NaN;
    
    // Casi speciali
    if (fabs(x) < DBL_EPSILON && fabs(y) < DBL_EPSILON) {
        // Limite per x,y -> 0
        double log_result = nu * log(nu) + 2 * lgammafn((nu + 1) / 2) -
                           log(M_PI) - 2 * lgammafn(nu / 2) - 
                           (nu + 1) * log(nu) - log(1 - rho * rho);
        return 4 * exp(log_result);
    }
    
    // Pre-calcolo di quantità riutilizzate
    const double rho_sq = rho * rho;
    const double one_minus_rho_sq = 1 - rho_sq;
    const double rho1 = rho * (1 - nugget);
    const double rho1_sq = rho1 * rho1;
    const double one_minus_rho1_sq = 1 - rho1_sq;
    
    // Controllo per evitare divisione per zero
    if (fabs(one_minus_rho_sq) < DBL_EPSILON || fabs(one_minus_rho1_sq) < DBL_EPSILON) {
        return R_NaN;
    }
    
    const double x_sq = x * x;
    const double y_sq = y * y;
    const double arg = (nu + 1) / 2;
    const double arg1 = nu / 2;
    
    // Calcolo di x2 e y2 con controllo underflow
    const double x2 = x_sq * one_minus_rho_sq + nu * one_minus_rho1_sq;
    const double y2 = y_sq * one_minus_rho_sq + nu * one_minus_rho1_sq;
    
    if (x2 <= 0 || y2 <= 0) return R_NaN;
    
    // Calcolo logaritmico per stabilità numerica
    double log_pp1, log_pp2;
    
    // log(pp1) = nu*log(nu) - arg*log(x2*y2) + 2*lgammafn(arg)
    log_pp1 = nu * log(nu) - arg * (log(x2) + log(y2)) + 2 * lgammafn(arg);
    
    // log(pp2) = log(π) + 2*lgammafn(arg1) + log(rho12) + log(rho22)
    // dove rho12 = (1-rho1²)^(-nu-0.5) e rho22 = (1-rho²)^(-nu/2-1)
    log_pp2 = log(M_PI) + 2 * lgammafn(arg1) - 
              (nu + 0.5) * log(one_minus_rho1_sq) - 
              (nu / 2 + 1) * log(one_minus_rho_sq);
    
    // Controllo overflow nei logaritmi
    if (!R_FINITE(log_pp1) || !R_FINITE(log_pp2)) {
        return R_NaN;
    }
    
    // Calcolo degli argomenti per appellF4
    const double x2_y2 = x2 * y2;
    const double denom_inv = 1.0 / x2_y2;
    
    // Primo argomento: (rho1 * x * y * (1-rho²))² / (x2 * y2)
    const double term1 = rho1 * x * y * one_minus_rho_sq;
    const double appell_arg1 = (term1 * term1) * denom_inv;
    
    // Secondo argomento: (rho * nu * (1-rho1²))² / (x2 * y2)
    const double term2 = rho * nu * one_minus_rho1_sq;
    const double appell_arg2 = (term2 * term2) * denom_inv;
    
    // Controllo che gli argomenti siano nel dominio di convergenza
    if (appell_arg1 >= 1.0 || appell_arg2 >= 1.0) {
        // Fuori dal raggio di convergenza
        return R_NaN;
    }
    
    // Chiamata alla funzione Appell F4 ottimizzata
    double app = appellF4(arg, arg, 0.5, arg1, appell_arg1, appell_arg2);
    
    if (!R_FINITE(app)) {
        return R_NaN;
    }
    
    // Calcolo finale usando aritmetica logaritmica quando possibile
    double log_result = log(4.0) + log_pp1 + log(fabs(app)) - log_pp2;
    
    // Controllo overflow nel risultato finale
    if (log_result > 700.0) return R_PosInf;
    if (log_result < -700.0) return 0.0;
    
    double result = exp(log_result);
    
    // Gestione del segno
    if (app < 0) result = -result;
    
    return result;
}
/*********** bivariate two piece-T distribution********************/ 

double biv_two_pieceT(double rho,double zi,double zj,double sill,double nuu,double eta,
             double p11,double mui,double muj,double nugget)
{
double res=0.0;
double nu=1/nuu;
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);

    if(rho>DBL_EPSILON){
    if(zistd>=0&&zjstd>=0)
{res=          (p11/R_pow(etamos,2))*appellF4_mod(nu,rho,zistd/etamos,zjstd/etamos,nugget);}
    if(zistd>=0&&zjstd<0)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho,zistd/etamos,-zjstd/etamas,nugget);}
      if(zistd<0&&zjstd>=0)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho,-zistd/etamas,zjstd/etamos,nugget);}
    if(zistd<0&&zjstd<0)
{res=    ((p11+eta)/R_pow(etamas,2))*appellF4_mod(nu,rho,-zistd/etamas,-zjstd/etamas,nugget);}

}
    if(rho<DBL_EPSILON){
    
    if(zi>=mui)
         {res=0.5*2*exp((nu/2)*log(nu)+lgammafn((nu+1)/2)-((nu+1)/2)*log(R_pow(zistd/etamos,2)+nu)-0.5*log(M_PI)-lgammafn(nu/2));}
         if(zj<muj)
         {res=0.5*2*exp((nu/2)*log(nu)+lgammafn((nu+1)/2)-((nu+1)/2)*log(R_pow(zjstd/etamas,2)+nu)-0.5*log(M_PI)-lgammafn(nu/2));}
      }
return(res/sill);
    
}

/***** bivariate half gaussian ****/     
double biv_half_Gauss(double rho,double zi,double zj)
{
double kk=0, dens=0,a=0,b=0,rho2=1-rho*rho;
  kk=(M_PI)*sqrt(rho2);
  a=exp(- (1/(2*(rho2)))*(R_pow(zi,2)+R_pow(zj,2)-2*rho*zi*zj));
  b=exp(- (1/(2*(rho2)))*(R_pow(zi,2)+R_pow(zj,2)+2*rho*zi*zj));
  dens=(a+b)/kk;
  return(dens);
}
/***** bivariate two piece gaussian *****/
double biv_two_pieceGaussian(double rho, double zi, double zj, double sill, double eta,
                              double p11, double mui, double muj)
{
    double res = 0.0;
    double etamas = 1 + eta;
    double etamos = 1 - eta;
    double zistd = (zi - mui) / sqrt(sill);
    double zjstd = (zj - muj) / sqrt(sill);
    
    double z1 = (zistd >= 0) ? zistd / etamos : -zistd / etamas;
    double z2 = (zjstd >= 0) ? zjstd / etamos : -zjstd / etamas;
    
    if (zistd >= 0 && zjstd >= 0) {
        res = (p11 / R_pow(etamos, 2)) * biv_half_Gauss(rho, z1, z2);
    } else if (zistd >= 0 && zjstd < 0) {
        res = ((1 - eta - 2 * p11) / (2 * (1 - eta * eta))) * biv_half_Gauss(rho, z1, -z2);
    } else if (zistd < 0 && zjstd >= 0) {
        res = ((1 - eta - 2 * p11) / (2 * (1 - eta * eta))) * biv_half_Gauss(rho, -z1, z2);
    } else if (zistd < 0 && zjstd < 0) {
        res = ((p11 + eta) / R_pow(etamas, 2)) * biv_half_Gauss(rho, -z1, -z2);
    }
    
    return res / sill;
}

double biv_beta(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2,double min,double max)
{
  double ki=0.0,kj=0.0,p1=0.0,p2=0.0,rho2=0.0,res=0.0,p3;
  double dd=max-min;
  zi=(zi-min)/dd;zj=(zj-min)/dd;
 ki=1-zi; kj=1-zj;
   double aa=0.5*(shape1+shape2);
if(rho) {
  rho2=rho*rho;

    p1=R_pow(zi*zj,shape1/2-1)*R_pow(ki*kj,shape2/2-1);
    p3=exp(2*lgammafn(aa)-(2*lgammafn(shape1/2)+2*lgammafn(shape2/2)-aa*log1p(-rho2)));
    p2= appellF4(aa,aa,shape1/2,shape2/2,rho2*zi*zj,rho2*ki*kj);
  res=p1*p2*p3;
} else  {p1=R_pow(zi,shape1/2-1)*R_pow(ki,shape2/2-1)*exp(lgammafn(aa)-lgammafn(shape1/2)-lgammafn(shape2/2));
         p2=R_pow(zj,shape1/2-1)*R_pow(kj,shape2/2-1)*exp(lgammafn(aa)-lgammafn(shape1/2)-lgammafn(shape2/2));
         res=p1*p2;}
return(res/R_pow(dd,2));
}


double biv_Kumara(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2,double min,double  max)
{
  double xx=0.0,yy=0.0,ki=0.0,kj=0.0,p1=0.0,p2=0.0,rho2=0.0,res=0.0;
    double dd=(max-min);
  zi=(zi-min)/dd;zj=(zj-min)/dd;

 ki=1-R_pow(zi,shape2); kj=1-R_pow(zj,shape2);
if(rho) {
  rho2=rho*rho;
  xx=rho2*R_pow(ki*kj,shape1);
  yy=rho2*(1-R_pow(ki,shape1))*(1-R_pow(kj,shape1));
  p1=R_pow(shape1*shape2,2)*R_pow(zi*zj,shape2-1)*R_pow(ki*kj,shape1-1)*R_pow(1-rho2,2);
  p2= appellF4(2,2,1,1,xx,yy);
  res=p1*p2;}
 else  {p1=shape1*shape2*R_pow(zi,shape2-1)* R_pow(ki,shape1-1);
         p2=shape1*shape2*R_pow(zj,shape2-1)* R_pow(kj,shape1-1);
         res=p1*p2;}
return(res/R_pow(dd,2));
}

double biv_Kumara2(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2,double min,double  max)
{
  double xx=0.0,yy=0.0,ki=0.0,kj=0.0,p1=0.0,p2=0.0,rho2=0.0,res=0.0;
  double mi=1/(1+exp(-ai)), mj=1/(1+exp(-aj));
  double dd=(max-min), shapei=log(0.5)/log1p(-R_pow(mi,shape2)),shapej=log(0.5)/log1p(R_pow(mj,shape2));
  zi=(zi-min)/dd;zj=(zj-min)/dd;

 ki=1-R_pow(zi,shape2); kj=1-R_pow(zj,shape2);
if(rho) {
  rho2=rho*rho;
  xx=rho2*R_pow(ki,shapei)*R_pow(kj,shapej);
  yy=rho2*(1-R_pow(ki,shapei))*(1-R_pow(kj,shapej));
  p1=R_pow(sqrt(shapei*shapej)*shape2,2)*R_pow(zi*zj,shape2-1)*R_pow(ki,shapei-1)*R_pow(kj,shapej-1)*R_pow(1-rho2,2);
  p2= appellF4(2,2,1,1,xx,yy);
  res=p1*p2;}
else  {p1=shapei*shape2*R_pow(zi,shape2-1)*R_pow(ki,shapei-1);
       p2=shapej*shape2*R_pow(zj,shape2-1)*R_pow(kj,shapej-1);
         res=p1*p2;}
return(res/R_pow(dd,2));
}

/***********************************************************************************/
/************ functions for binomial or negative binomial  two piece *********/
/***********************************************************************************/

/**********************************************************************/
double pblogi22(double lim1,double lim2,double corr)
{

double value=0.0,sum1=0.0,sum2=0.0,term=0.0,term2=0.0,kk=0,bb=0.0;
double corr21=1-R_pow(corr,2);
int m=0,n=0,p1,p2,p3;
double e1=exp(lim1);double e2=exp(lim2);
while(n<=400){
bb=exp(2*n*log(corr)+n*(lim1+lim2)- 2*log(n+1));
sum1=0.0;m=0;
p1=1+n;p3=p1+1;
while(m<=200)
{ 
      p2=2+n+m;
      term=exp(2*m*log(corr)+log(hypergeo(p1,p2,p3,-e1))+
                           log(hypergeo(p1,p2,p3,-e2))-2*lbeta(p1,m+1));
     sum1=sum1+ term;
     if(term<1e-7) {break;}
     m=m+1;      

}
term2=bb*sum1; 
sum2=sum2+term2;
if(term2<1e-7) {break;}
n=n+1;
}
kk=exp(2*log(corr21)+lim1+lim2);
value= sum2*kk;
 return(value);
}

// cdf bivariate half-normal distribution
//***********************************************//
// cdf univariate half-gaussian distribution
/***********************************************************************************/
/***********************************************************************************/
/********** bivariate log-logistic **********/
double biv_LogLogistic(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double c=gammafn(1+1/shape)*gammafn(1-1/shape);
    double A=0.0,res=0.0,B=0.0,C=0.0;
    double ci=exp(mui);double cj=exp(muj);
    double ki=R_pow(c*zi/ci,shape)+1;
    double kj=R_pow(c*zj/cj,shape)+1;
    double rho2=R_pow(corr,2);
    double kij=ki*kj;
   // if(corr)   {
        A=(R_pow(c*shape,2)/(ci*cj))*R_pow((c*c*zi*zj)/(ci*cj),shape-1)*R_pow(1-rho2,2);
        B=R_pow(kij,-2);
        C=appellF4(2,2,1,1,
        (rho2*R_pow(c*c*zi*zj,shape))/(R_pow(ci*cj,shape)*kij),
        rho2/(kij));
        res=A*B*C;
    //}else{    B=(c*shape/ci)*R_pow((c*zi/ci),shape-1)*R_pow(ki,-2);  C=(c*shape/cj)*R_pow((c*zj/cj),shape-1)*R_pow(kj,-2);
    //    res=B*C;}
    return(res);
}

/********** bivariate logistic **********/
double biv_Logistic(double corr,double zi,double zj,double mui, double muj, double sill)
{
    double a=0.0,A=0.0,res=0.0,B=0.0,C=0.0;
    double ci=mui;double cj=muj;
    double ki=exp((zi-ci)/sqrt(sill));
    double kj=exp((zj-cj)/sqrt(sill));
    double rho2=R_pow(corr,2);
    //if(corr)   {
        a=1-rho2;
        A=(ki*kj)/(R_pow(a,-2)*sill);
        B=R_pow((ki+1)*(kj+1),-2);
        C=appellF4(2,2,1,1,
        (rho2*ki*kj)/((ki+1)*(kj+1)),
        rho2/((ki+1)*(kj+1)));
        res=A*B*C;
    //} else{ B=ki*R_pow((ki+1),-2)/sqrt(sill);C=kj*R_pow((kj+1),-2)/sqrt(sill);res=B*C;}
    return(res);
}

//********** functions for bivariate tukey h *******************//

// compute lambert w function
double LambertW(double z) {
  const double eps = 4.0e-16;               // Precisione del calcolo
  const double em1 = 0.3678794411714423;    // e^(-1)
  double p, e, t, w;
  
  // Gestione dei casi speciali
  if (z < -em1 || isinf(z) || isnan(z)) {
    return R_NaN;  // Meglio restituire NAN piuttosto che 0.0 per input invalidi
  }
  
  if (z == 0.0) return 0.0;
  
  // Serie di approssimazione vicino a -e^(-1)
  if (z < -em1 + 1e-4) {
    double q = z + em1;
    double r = sqrt(q);
    double q2 = q * q;
    double q3 = q2 * q;
    return -1.0
           + 2.331643981597124203 * r
           - 1.812187885639363490 * q
           + 1.936631114492359755 * r * q
           - 2.353551201881614516 * q2
           + 3.066858901050631912 * r * q2
           - 4.175335600258177138 * q3
           + 5.858023729874774148 * r * q3
           - 8.401032217523977370 * q3 * q;
  }
  
  // Calcolo del valore iniziale di w per l'iterazione
  if (z < 1.0) {
    // Espansione in serie vicino a 0
    p = sqrt(2.0 * (M_E * z + 1.0)); // Uso di M_E da math.h
    w = -1.0 + p * (1.0 + p * (-1.0/3.0 + p * (1.0/72.0)));  // Migliorata precisione
  } else {
    // Approssimazione asintotica
    w = log(z);
    // Miglioramento per z > 3.0
    if (z > 3.0) {
      w -= log(w);
      // Ulteriore miglioramento per z molto grandi
      if (z > 30.0) {
        double l1 = log(w);
        w -= l1 * (1.0 - l1/(2.0 * w));
      }
    }
  }
  
  int i;
  for (i = 0; i < 6; i++) {  
    e = exp(w);
    t = w * e - z;
    p = w + 1.0;
    
    double f = e * p; 
    double f1 = (p + 1.0) * t / (2.0 * p);
    double d = t / (f - f1);
    
    w -= d;
    if (fabs(d) < eps * (1.0 + fabs(w))) {
      return w;
    }
  }
  return w;
}

// pdf bivariate gaussian distribution
double dbnorm(double x_i,double x_j,double mean_i,double mean_j,double sill,double corr)
{
    double  fraq = 1.0,dens = 0.0,aux1 = 1.0,z1 = 1.0,z2 = 1.0,z3 = 1.0,z = 1.0;
    fraq = 2*M_PI*sill*sqrt(1-corr*corr);
    z1   = (x_i - mean_i)*(x_i - mean_i);
    z2   = (x_j - mean_j)*(x_j - mean_j);
    z3   = 2*corr*(x_i - mean_i)*(x_j - mean_j);
    z    = (z1 + z2 - z3)/sill;
    aux1 = 2*(1-corr*corr); 
    dens = (1/fraq)*exp(-z/aux1);
    return(dens);
}

// pdf bivariate tukey h random field 
double biv_tukey_h(double corr,double data_i, double data_j, double mean_i, double mean_j, double tail, double sill)
{
    double z_i, z_j, x_i, x_j, jac_i, jac_j, phi2, dens;

    if (!R_FINITE(corr) || !R_FINITE(data_i) || !R_FINITE(data_j) ||
        !R_FINITE(mean_i) || !R_FINITE(mean_j) || !R_FINITE(tail) ||
        !R_FINITE(sill) || sill <= 0.0 || tail < 0.0) {
        return R_NaN;
    }

    z_i = (data_i - mean_i) / sqrt(sill);
    z_j = (data_j - mean_j) / sqrt(sill);

    x_i = inv_tukey_h(z_i, tail);
    x_j = inv_tukey_h(z_j, tail);
    jac_i = jac_tukey_h(z_i, tail);
    jac_j = jac_tukey_h(z_j, tail);

    if (!R_FINITE(x_i) || !R_FINITE(x_j) ||
        !R_FINITE(jac_i) || !R_FINITE(jac_j)) {
        return R_NaN;
    }

    if (corr <= -1.0) corr = -1.0 + 1e-12;
    if (corr >=  1.0) corr =  1.0 - 1e-12;

    phi2 = dbnorm(x_i, x_j, 0.0, 0.0, 1.0, corr);
    if (!R_FINITE(phi2)) return R_NaN;

    /* The two location-scale derivatives contribute the factor 1/sill. */
    dens = phi2 * jac_i * jac_j / sill;
    return R_FINITE(dens) ? dens : R_NaN;
}



/***** bivariate half tukey h ****/     
double biv_half_Tukeyh(double rho,double ti,double tj,double tail)
{
  double dens = 0.0;
  dens = biv_tukey_h(rho,ti,tj,0,0,tail,1) + biv_tukey_h(rho,-ti,-tj,0,0,tail,1) +
    biv_tukey_h(rho,-ti,tj,0,0,tail,1) +
    biv_tukey_h(rho,ti,-tj,0,0,tail,1);
  return(dens);
}
 
/***** bivariate two piece tukey h ****/ 
double biv_two_pieceTukeyh(double rho,double zi,double zj,double sill,double eta,double tail,
             double p11,double mui,double muj)
{
double res=0.0;  
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);
/*if(rho)   {*/
//if(zi>=mui&&zj>=muj)
    if(zistd>=0&&zjstd>=0)
{res=          (p11/R_pow(etamos,2))*biv_half_Tukeyh(rho,zistd/etamos,zjstd/etamos,tail);}
//if(zi>=mui&&zj<muj)
    if(zistd>=0&&zjstd<0)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Tukeyh(rho,zistd/etamos,-zjstd/etamas,tail);}
//if(zi<mui&&zj>=muj)
      if(zistd<0&&zjstd>=0)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Tukeyh(rho,-zistd/etamas,zjstd/etamos,tail);}
//if(zi<mui&&zj<muj)
    if(zistd<0&&zjstd<0)
{res=    ((p11+eta)/R_pow(etamas,2))*biv_half_Tukeyh(rho,-zistd/etamas,-zjstd/etamas,tail);}
/*}else{   if(zi>=mui)
         {res=dnorm(x_i,0,1,0)*x_i/((zistd/etamos)*(1+LambertW(tail*R_pow(zistd/etamos,2))));}
         if(zj<muj)
         {res=dnorm(x_j,0,1,0)*x_j/((zjstd/etamas)*(1+LambertW(tail*R_pow(zjstd/etamas,2))));}
      }*/
return(res/sill);
}
 


void biv_pois_call(double *corr,int *r, int *t, double *mean_i, double *mean_j,double *res)
{
    *res = biv_Poisson(*corr,*r,*t,*mean_i,*mean_j);
}


/*******************************************************************************/

/******************************************************************************
 * VERSIONE "LOG" (ritorna log(p)) delle funzioni che mi hai dato.
 *
 * NOTA:
 * - Tutte le funzioni Prt/Prr/Pr0/P00/biv_Poisson ora ritornano LOG-probabilità.
 * - In caso di probabilità 0 o input non valido, ritornano -INFINITY.
 *
 * Dipendenze attese (come nel tuo codice):
 * - R_pow, R_finite, lgammafn, igam, hyperg, aprox_reg_1F1
 ******************************************************************************/

/*******************************************************************************/
// log-sum-exp stabile (come tuo)
/*******************************************************************************/
static inline double log_sum_exp(double a, double b) {
    if (a == -INFINITY) return b;
    if (b == -INFINITY) return a;

    double max_val = (a > b) ? a : b;
    double min_val = (a > b) ? b : a;

    if (max_val - min_val > 700) return max_val;
    return max_val + log1p(exp(min_val - max_val));
}

/*******************************************************************************/
// log(1 - exp(x)) stabile per x <= 0
/*******************************************************************************/
static inline double log1mexp_local(double x) {
    // ritorna log(1 - exp(x)) per x<=0
    if (x > 0.0) return NAN;
    // soglia ~ log(0.5)
    if (x < -0.693147180559945309417232121458) {
        return log1p(-exp(x));
    } else {
        return log(-expm1(x));
    }
}

/*******************************************************************************/
// log(exp(a) - exp(b)) stabile, richiede a > b
/*******************************************************************************/
static inline double logdiffexp_local(double a, double b) {
    if (b >= a) return -INFINITY;
    return a + log1mexp_local(b - a);
}
/*******************************************************************************/
// Convergenza ultra-veloce (tuo)
/*******************************************************************************/
static inline bool check_ultra_fast_convergence(double current, double previous,
                                                double mean_scale, int iteration) {
    if (current == 0.0 && previous == 0.0) return true;

    int min_iter = (mean_scale > 20.0) ? 3 : (mean_scale > 5.0) ? 5 : 8;
    if (iteration < min_iter) return false;

    double abs_diff = fabs(current - previous);
    double scale = fmax(fabs(current), fabs(previous));

    double rel_tol = (mean_scale > 20.0) ? 1e-8 : (mean_scale > 5.0) ? 1e-10 : 1e-12;
    double abs_tol = (mean_scale > 20.0) ? 1e-10 : (mean_scale > 5.0) ? 1e-12 : 1e-15;

    return abs_diff < fmax(abs_tol, rel_tol * scale);
}
/*******************************************************************************/
// Cache lgamma (tuo): fast_lgamma(n) = lgamma(n+1)
/*******************************************************************************/
static double cached_lgamma[256];
static bool cache_initialized = false;

static inline double fast_lgamma(int n) {
    if (!cache_initialized) {
        for (int i = 0; i < 256; i++) cached_lgamma[i] = lgammafn(i + 1);
        cache_initialized = true;
    }
    if (n < 255) return cached_lgamma[n];
    return lgammafn(n + 1);
}
/*******************************************************************************/
// log(igam) safe (tuo)
/*******************************************************************************/
static inline double safe_log_igam_fast(double a, double x) {
    if (x <= 0.0 || a <= 0.0) return -INFINITY;

    if (x > a + 50.0 * sqrt(a) + 100.0) return -INFINITY;

    double igam_val = igam(a, x);
    if (igam_val <= 0.0 || !R_finite(igam_val)) return -INFINITY;
    if (igam_val < 1e-300) return -INFINITY;

    return log(igam_val);
}
/*******************************************************************************/
// Stima iterazioni (tuo)
/*******************************************************************************/
static inline int estimate_iterations(double mean_scale, int base_iter) {
    if (mean_scale > 50.0) return (int)fmin(base_iter / 4.0, 200.0);
    if (mean_scale > 20.0) return (int)fmin(base_iter / 3.0, 400.0);
    if (mean_scale > 10.0) return (int)fmin(base_iter / 2.0, 800.0);
    return base_iter;
}

/*******************************************************************************/
// helper: log(exp(a)+exp(b)+exp(c)) con log_sum_exp
/*******************************************************************************/
static inline double log_sum_exp3(double a, double b, double c) {
    return log_sum_exp(log_sum_exp(a, b), c);
}
/*******************************************************************************/
// helper: log( exp(a) + exp(b) - exp(c) )  (richiede risultato > 0)
/*******************************************************************************/
static inline double log_add_add_sub(double a, double b, double c) {
    // S = exp(a)+exp(b)-exp(c)
    // calcolo in scala con massimo per ridurre overflow
    double M = a;
    if (b > M) M = b;
    if (c > M) M = c;

    if (M == -INFINITY) return -INFINITY;

    double sa = (a == -INFINITY) ? 0.0 : exp(a - M);
    double sb = (b == -INFINITY) ? 0.0 : exp(b - M);
    double sc = (c == -INFINITY) ? 0.0 : exp(c - M);

    double s = sa + sb - sc;
    if (!(s > 0.0) || !isfinite(s)) return -INFINITY;

    return M + log(s);
}
/*******************************************************************************/
// helper: log( exp(a) + exp(b) + exp(c) - 1 )
/*******************************************************************************/
static inline double log_add3_minus1(double a, double b, double c) {
    // pos = exp(a)+exp(b)+exp(c)
    double log_pos = log_sum_exp3(a, b, c);
    // log(pos - 1) = logdiffexp(log_pos, 0) se pos>1
    return logdiffexp_local(log_pos, 0.0);
}
/*******************************************************************************/
/*
 * Off-diagonal Poisson probability: exact positive single-series form.
 *
 * For r > t >= 0, write rho2=corr^2, d=1-rho2,
 *   ai = mean_i/d, aj = mean_j/d,
 *   x  = rho2*ai, y = x*aj, n = r-t.
 * Then
 *
 * p_{r,t} = exp(-ai-aj) mean_i^r aj^t/t!
 *           * sum_{s>=0} y^s/s! * 1F1~(n+1; r+s+1; x),
 *
 * where 1F1~(a;b;x)=1F1(a;b;x)/Gamma(b).
 *
 * This is algebraically equivalent to the former Prt/Pr0 formulas, but it
 * contains only positive terms, removes the former double sum in Prt and the
 * subtraction of nearly equal positive quantities.  Hence it also avoids the
 * cancellation that could make logdiffexp_local() return -Inf in difficult
 * off-diagonal cases.
 */
static inline double log_reg_1F1_offdiag(int a, int b, double x) {
    const double log_gamma_b = (b < 255) ? fast_lgamma(b - 1) : lgammafn((double)b);

    /* Primary evaluation. */
    double h = hyperg((double)a, (double)b, x);
    if (R_finite(h) && h > 0.0) {
        return log(h) - log_gamma_b;
    }

    /* aprox_reg_1F1() already approximates the REGULARIZED 1F1. */
    double ar = aprox_reg_1F1(a, b, x);
    if (R_finite(ar) && ar > 0.0) {
        return log(ar);
    }

    /* Last fallback: Kummer transformation
     * 1F1(a;b;x)=exp(x) 1F1(b-a;b;-x). */
    const int ak = b - a;
    if (ak >= 0) {
        double hk = hyperg((double)ak, (double)b, -x);
        if (R_finite(hk) && hk > 0.0) {
            return x + log(hk) - log_gamma_b;
        }
    }

    return -INFINITY;
}

static double Prt_single(double corr, int r, int t,
                         double mean_i, double mean_j) {
    if (mean_i <= 0.0 || mean_j <= 0.0 || r < 0 || t < 0 || r <= t) {
        return -INFINITY;
    }

    const double rho2 = R_pow(corr, 2);
    const double one_minus_rho2 = 1.0 - rho2;
    if (one_minus_rho2 <= 1e-15 || rho2 <= 0.0) return -INFINITY;

    const double auxi = mean_i / one_minus_rho2;
    const double auxj = mean_j / one_minus_rho2;
    const double x = rho2 * auxi;
    const int n = r - t;

    /* y = rho2*auxi*auxj is used only through log(y), which avoids overflow. */
    const double log_y = log(rho2) + log(auxi) + log(auxj);

    /* A generous deterministic cap centred on the Bessel-like sqrt(y) scale.
     * The actual exit is controlled by the relative positive-term tail below. */
    double sqrt_y;
    if (0.5 * log_y > log(50000.0)) sqrt_y = 50000.0;
    else sqrt_y = exp(0.5 * log_y);

    int max_iter = (int)ceil(4.0 * sqrt_y + 4.0 * (r + t + 1) + 100.0);
    if (max_iter < 200) max_iter = 200;
    if (max_iter > 50000) max_iter = 50000;

    const double log_rel_tol = log(1e-13);
    double log_sum = -INFINITY;
    double prev_log_term = INFINITY;
    int tail_small = 0;
    bool converged = false;

    for (int s = 0; s < max_iter; ++s) {
        const int b = r + s + 1;
        const double log_reg = log_reg_1F1_offdiag(n + 1, b, x);

        if (log_reg == -INFINITY) {
            /* Do not terminate on one special-function failure: a later b can
             * still be evaluable.  Repeated failures trigger the legacy path. */
            ++tail_small;
            if (tail_small > 20 && s > 20) break;
            continue;
        }

        const double log_term = s * log_y - fast_lgamma(s) + log_reg;
        log_sum = (log_sum == -INFINITY) ? log_term : log_sum_exp(log_sum, log_term);

        /* All summands are positive.  Once terms are decreasing and several
         * consecutive terms are negligible relative to the accumulated sum,
         * truncation is safe at the requested scale. */
        if (s > 10 && log_term < prev_log_term && (log_term - log_sum) < log_rel_tol) {
            ++tail_small;
            if (tail_small >= 6) {
                converged = true;
                break;
            }
        } else {
            tail_small = 0;
        }
        prev_log_term = log_term;
    }

    if (log_sum == -INFINITY || !converged) return -INFINITY;

    const double log_base = -auxi - auxj
                          + r * log(mean_i)
                          + t * log(auxj)
                          - fast_lgamma(t);
    const double out = log_base + log_sum;
    return R_finite(out) ? out : -INFINITY;
}

/*******************************************************************************/
/*  Prt  (legacy exact representation; now used as numerical fallback) */
/*******************************************************************************/
double Prt(double corr, int r, int t, double mean_i, double mean_j) {
    if (mean_i <= 0.0 || mean_j <= 0.0 || r < 0 || t < 0 || r <= t) {
        return -INFINITY;
    }

    const double rho2 = R_pow(corr, 2);
    const double one_minus_rho2 = 1.0 - rho2;
    if (one_minus_rho2 <= 1e-15 || rho2 <= 0.0) return -INFINITY;

    const double mean_scale = mean_i + mean_j;

    const double log_rho2 = log(rho2);
    const double log_rho2_ratio = log_rho2 - log1p(-rho2);
    const double log_mean_i = log(mean_i);
    const double auxi = mean_i / one_minus_rho2;
    const double auxj = mean_j / one_minus_rho2;
    const double rho2auxi = rho2 * auxi;

    const int n = r - t;

    const int iter1 = estimate_iterations(mean_scale, 1000);
    const int iter2 = estimate_iterations(mean_scale, 500);

    const double log_gamma_t = fast_lgamma(t - 1);

    double log_sum  = -INFINITY;
    double log_sum1 = -INFINITY;
    double prev_log_sum  = -INFINITY;
    double prev_log_sum1 = -INFINITY;

    /*
     * Exact-control-flow cache for the expensive special functions used in
     * the first (m,k) sum.  Their arguments depend on s=m+k, so values are
     * identical along each diagonal.  We keep the original nested loops,
     * summation order and stopping rules unchanged and only reuse already
     * computed function values.  igam() and hyperg() themselves are untouched.
     *
     * iter1<=1000 and iter2<=500 by estimate_iterations(), hence the largest
     * offset used below is iter1+iter2+1 <= 1501.
     */
#define PRT_CACHE_MAX 1502
    double igam_cache[PRT_CACHE_MAX];
    double q1_cache[PRT_CACHE_MAX];
    unsigned char igam_ready[PRT_CACHE_MAX];
    unsigned char q1_ready[PRT_CACHE_MAX];
    memset(igam_ready, 0, sizeof(igam_ready));
    memset(q1_ready, 0, sizeof(q1_ready));

    for (int m = 0; m <= iter1; m++) {
        int tm = t + m;
        double log_gamma_tm = fast_lgamma(tm - 1);
        double log_gamma_m1 = fast_lgamma(m);

        // -------- prima somma interna su k --------
        double log_inner_sum = -INFINITY;
        double prev_log_inner = -INFINITY;
        int consecutive_small = 0;

        for (int k = 0; k <= iter2; k++) {
            int ss = m + k;
            int aa = t + ss + 1;
            int cc = r + ss + 1;
            int ig_idx = aa - t;  /* 1..iter1+iter2+1 */
            double log_gamma_cc = (cc < 255) ? fast_lgamma(cc - 1) : lgammafn(cc);

            double log_igam;
            if (ig_idx >= 0 && ig_idx < PRT_CACHE_MAX && igam_ready[ig_idx]) {
                log_igam = igam_cache[ig_idx];
            } else {
                log_igam = safe_log_igam_fast(aa, auxj);
                if (ig_idx >= 0 && ig_idx < PRT_CACHE_MAX) {
                    igam_cache[ig_idx] = log_igam;
                    igam_ready[ig_idx] = 1;
                }
            }

            if (log_igam == -INFINITY) {
                consecutive_small++;
                if (consecutive_small > 3) break;
                continue;
            }

            double log_q1;
            if (ss >= 0 && ss < PRT_CACHE_MAX && q1_ready[ss]) {
                log_q1 = q1_cache[ss];
            } else {
                double hyper = hyperg(n, cc, rho2auxi);

                if (R_finite(hyper) && hyper > 0.0) {
                    log_q1 = log(hyper) - log_gamma_cc;
                } else {
                    double approx = aprox_reg_1F1(n, cc, rho2auxi);
                    if (approx > 0.0 && R_finite(approx)) {
                        log_q1 = log(approx) - log_gamma_cc;
                    } else {
                        log_q1 = -INFINITY;
                    }
                }

                if (ss >= 0 && ss < PRT_CACHE_MAX) {
                    q1_cache[ss] = log_q1;
                    q1_ready[ss] = 1;
                }
            }

            if (log_q1 == -INFINITY) {
                consecutive_small++;
                if (consecutive_small > 3) break;
                continue;
            }

            double log_term = (k + m) * log_rho2_ratio
                            + log_gamma_tm - log_gamma_m1 - log_gamma_t
                            + (cc - 1) * log_mean_i
                            + log_igam + log_q1;

            if (log_term < -400) {
                consecutive_small++;
                if (consecutive_small > 5) break;
                continue;
            }

            consecutive_small = 0;
            log_inner_sum = (log_inner_sum == -INFINITY) ? log_term : log_sum_exp(log_inner_sum, log_term);

            if (k > 3 && check_ultra_fast_convergence(log_inner_sum, prev_log_inner, mean_scale, k)) break;
            prev_log_inner = log_inner_sum;
        }

        if (log_inner_sum != -INFINITY) {
            log_sum = (log_sum == -INFINITY) ? log_inner_sum : log_sum_exp(log_sum, log_inner_sum);
        }

        // -------- seconda somma (term1), control flow unchanged --------
        int cc2 = tm + n + 1;
        double hyper2 = hyperg(n + 1, cc2, rho2auxi);
        double log_q2;

        if (R_finite(hyper2) && hyper2 > 0.0) {
            log_q2 = log(hyper2) - ((cc2 < 255) ? fast_lgamma(cc2 - 1) : lgammafn(cc2));
        } else {
            double approx2 = aprox_reg_1F1(n + 1, cc2, rho2auxi);
            if (approx2 > 0.0 && R_finite(approx2)) {
                log_q2 = log(approx2) - ((cc2 < 255) ? fast_lgamma(cc2 - 1) : lgammafn(cc2));
            } else {
                continue;
            }
        }

        /* Reuse igam when the same (a, auxj) was already requested in the
         * first sum; this does not change any summation or stopping decision. */
        int ig2_idx = tm - t; /* m */
        double log_igam2;
        if (ig2_idx >= 0 && ig2_idx < PRT_CACHE_MAX && igam_ready[ig2_idx]) {
            log_igam2 = igam_cache[ig2_idx];
        } else {
            log_igam2 = safe_log_igam_fast(tm, auxj);
            if (ig2_idx >= 0 && ig2_idx < PRT_CACHE_MAX) {
                igam_cache[ig2_idx] = log_igam2;
                igam_ready[ig2_idx] = 1;
            }
        }
        if (log_igam2 == -INFINITY) continue;

        double log_term1 = m * log_rho2_ratio
                         + log_gamma_tm + (tm + n) * log_mean_i
                         - log_gamma_m1 - log_gamma_t
                         + log_q2 + log_igam2;

        if (log_term1 >= -400) {
            log_sum1 = (log_sum1 == -INFINITY) ? log_term1 : log_sum_exp(log_sum1, log_term1);
        }

        if (m > 3 &&
            check_ultra_fast_convergence(log_sum1, prev_log_sum1, mean_scale, m) &&
            check_ultra_fast_convergence(log_sum,  prev_log_sum,  mean_scale, m)) {
            break;
        }

        prev_log_sum  = log_sum;
        prev_log_sum1 = log_sum1;
    }

#undef PRT_CACHE_MAX

    // -------- finale in log: unchanged --------
    if (log_sum1 == -INFINITY) return -INFINITY;

    double a = -auxi + log_sum1;
    if (log_sum == -INFINITY) {
        return a;
    }

    double b = -auxi + log_sum;
    return logdiffexp_local(a, b);
}
/*******************************************************************************/
/*  Prr  (ora ritorna LOG(Prr)) */
/*******************************************************************************/
double Prr(double corr, int r, int t, double mean_i, double mean_j) {
    if (mean_i <= 0.0 || mean_j <= 0.0 || r < 0 || t < 0 || r != t) return -INFINITY;
    if (r == 0) return P00(corr, r, t, mean_i, mean_j);

    const double rho2 = R_pow(corr, 2);
    const double one_minus_rho2 = 1.0 - rho2;
    if (one_minus_rho2 <= 1e-15 || rho2 <= 0.0) return -INFINITY;

    const double auxi = mean_i / one_minus_rho2;
    const double auxj = mean_j / one_minus_rho2;
    const double rho2_auxi = rho2 * auxi;
    const double rho2_auxj = rho2 * auxj;

    const double log_rho2 = log(rho2);
    const double log_one_minus_rho2 = log1p(-rho2);
    const double mean_scale = mean_i + mean_j;

    const int iter1 = estimate_iterations(mean_scale, 600);
    const int iter2 = estimate_iterations(mean_scale, 300);
    const int iter_collapsed = iter1 + iter2 - 1;

    const double log_gamma_r = fast_lgamma(r - 1); /* log Gamma(r) */
    const double log_gamma_r1 = fast_lgamma(r);    /* log Gamma(r+1) */

    double log_sum  = -INFINITY; /* collapsed former double sum */
    double log_sum1 = -INFINITY; /* term1 */
    double log_sum2 = -INFINITY; /* term2+term3 */

    /*
     * Cache only repeated requests already present in opt1.
     * In the collapsed series, safe_log_igam_fast() is evaluated at
     *   a = r + ss + 1, x = auxi/auxj.
     * The following one-dimensional series asks again for
     *   a = r + k,      x = auxi/auxj.
     * Thus k>=1 reuses exactly the value computed at ss=k-1 whenever
     * that collapsed-series entry was reached.  igam() itself, all loop
     * bounds, summation order, and convergence tests are unchanged.
     */
#define PRR_AUX_CACHE_MAX 902
    double auxi_logigam_cache[PRR_AUX_CACHE_MAX];
    double auxj_logigam_cache[PRR_AUX_CACHE_MAX];
    int aux_cache_max = 0; /* largest valid index a-r; index 0 is uncached */

    {
        int stable_count = 0;
        int consecutive_small = 0;

        for (int ss = 0; ss < iter_collapsed; ++ss) {
            const int upper_m = (ss < iter2) ? ss : (iter2 - 1);
            const int lower_m = (ss >= iter1) ? (ss - iter1 + 1) : 0;
            if (lower_m > upper_m) continue;

            const double log_cum_upper =
                ((r + upper_m) < 255 ? fast_lgamma(r + upper_m) : lgammafn(r + upper_m + 1.0))
                - log_gamma_r1
                - ((upper_m < 255) ? fast_lgamma(upper_m) : lgammafn(upper_m + 1.0));

            double log_coeff = log_cum_upper;
            if (lower_m > 0) {
                const int lm1 = lower_m - 1;
                const double log_cum_lower =
                    ((r + lm1) < 255 ? fast_lgamma(r + lm1) : lgammafn(r + lm1 + 1.0))
                    - log_gamma_r1
                    - ((lm1 < 255) ? fast_lgamma(lm1) : lgammafn(lm1 + 1.0));
                log_coeff = logdiffexp_local(log_cum_upper, log_cum_lower);
                if (log_coeff == -INFINITY) continue;
            }

            const int aa = r + ss + 1;
            const double log_igam_i = safe_log_igam_fast(aa, auxi);
            const double log_igam_j = safe_log_igam_fast(aa, auxj);

            const int cache_idx = ss + 1;
            if (cache_idx < PRR_AUX_CACHE_MAX) {
                auxi_logigam_cache[cache_idx] = log_igam_i;
                auxj_logigam_cache[cache_idx] = log_igam_j;
                aux_cache_max = cache_idx;
            }

            /* For fixed x, once both tails underflow at increasing a they
             * cannot recover at larger a. */
            if (log_igam_i == -INFINITY || log_igam_j == -INFINITY) break;

            const double log_term = log_one_minus_rho2
                                  + ss * log_rho2
                                  + log_coeff
                                  + log_igam_i + log_igam_j;

            if (log_term < -500.0) {
                ++consecutive_small;
                if (consecutive_small > 8) break;
                continue;
            }
            consecutive_small = 0;

            log_sum = (log_sum == -INFINITY) ? log_term : log_sum_exp(log_sum, log_term);

            if (ss > 8 && (log_term - log_sum) < -29.933606208922594) {
                ++stable_count;
                if (stable_count >= 4) break;
            } else {
                stable_count = 0;
            }
        }
    }

    {
        double prev_log_sum1 = -INFINITY, prev_log_sum2 = -INFINITY;

        for (int k = 0; k < iter1; k++) {
            const double log_gamma_k1 = fast_lgamma(k);
            const double log_gamma_rk = (r + k < 255) ? fast_lgamma(r + k - 1) : lgammafn(r + k);

            double log_igam_auxi, log_igam_auxj;
            if (k >= 1 && k <= aux_cache_max) {
                log_igam_auxi = auxi_logigam_cache[k];
                log_igam_auxj = auxj_logigam_cache[k];
            } else {
                log_igam_auxi = safe_log_igam_fast(r + k, auxi);
                log_igam_auxj = safe_log_igam_fast(r + k, auxj);
            }
            if (log_igam_auxi == -INFINITY || log_igam_auxj == -INFINITY) break;

            const double log_term1 = k * log_rho2 + log_gamma_rk + log_igam_auxi + log_igam_auxj
                                   - log_gamma_k1 - log_gamma_r;

            if (log_term1 >= -500.0) {
                log_sum1 = (log_sum1 == -INFINITY) ? log_term1 : log_sum_exp(log_sum1, log_term1);
            }

            const double log_igam_rho_auxi = safe_log_igam_fast(r + k, rho2_auxi);
            const double log_igam_rho_auxj = safe_log_igam_fast(r + k, rho2_auxj);

            if (log_igam_rho_auxi != -INFINITY && log_igam_rho_auxj != -INFINITY) {
                const double log_term2 = -mean_i - r * log_rho2 + log_gamma_rk + log_igam_rho_auxi
                                       + log_igam_auxj - log_gamma_k1 - log_gamma_r;

                const double log_term3 = -mean_j - r * log_rho2 + log_gamma_rk + log_igam_auxi
                                       + log_igam_rho_auxj - log_gamma_k1 - log_gamma_r;

                const double log_terms23 = log_sum_exp(log_term2, log_term3);
                log_sum2 = (log_sum2 == -INFINITY) ? log_terms23 : log_sum_exp(log_sum2, log_terms23);
            }

            if (k > 3 &&
                check_ultra_fast_convergence(log_sum1, prev_log_sum1, mean_scale, k) &&
                check_ultra_fast_convergence(log_sum2, prev_log_sum2, mean_scale, k)) {
                break;
            }
            prev_log_sum1 = log_sum1;
            prev_log_sum2 = log_sum2;
        }
    }

#undef PRR_AUX_CACHE_MAX

    const double log_factor = r * log_one_minus_rho2;
    const double log_body = log_add_add_sub(log_sum2, log_sum, log_sum1);
    if (log_body == -INFINITY) return -INFINITY;

    return log_factor + log_body;
}

/*******************************************************************************/
/*  Pr0  (ora ritorna LOG(Pr0)) */
/*******************************************************************************/
double Pr0(double corr, int r, int t, double mean_i, double mean_j) {
    if (mean_i <= 0.0 || mean_j <= 0.0 || r < 0 || t < 0) return -INFINITY;

    const double rho2 = R_pow(corr, 2);
    const double one_minus_rho2 = 1.0 - rho2;
    if (one_minus_rho2 <= 1e-15 || rho2 <= 0.0) return -INFINITY;

    const double log_rho2 = log(rho2);
    const double log_rho_ratio = log_rho2 - log1p(-rho2);

    const double auxi = mean_i / one_minus_rho2;
    const double auxj = mean_j / one_minus_rho2;
    const double log_mean_i = log(mean_i);
    const double mean_scale = mean_i + mean_j;

    const int n = r - t;
    const int iter = estimate_iterations(mean_scale, 1500);

    const double rho2_auxi = rho2 * auxi;

    double log_sum = -INFINITY;
    double prev_log_sum = -INFINITY;
    int consecutive_fails = 0;

    for (int m = 0; m <= iter; m++) {
        double log_aux  = m * log_rho_ratio;
        double log_aux1 = (m + n) * log_mean_i;

        double hyper = hyperg(n, m + n + 1, rho2_auxi);
        if (!R_finite(hyper) || hyper <= 0.0) {
            consecutive_fails++;
            if (consecutive_fails > 10) break;
            continue;
        }

        double log_q2 = log(hyper) - ((m + n + 1 < 255) ? fast_lgamma(m + n) : lgammafn(m + n + 1));

        double log_igam = safe_log_igam_fast(m + 1, auxj);
        if (log_igam == -INFINITY) {
            consecutive_fails++;
            if (consecutive_fails > 10) break;
            continue;
        }

        double log_term = log_aux + log_aux1 + log_q2 + log_igam;

        if (log_term < -500) {
            consecutive_fails++;
            if (consecutive_fails > 15) break;
            continue;
        }

        consecutive_fails = 0;
        log_sum = (log_sum == -INFINITY) ? log_term : log_sum_exp(log_sum, log_term);

        if (m > 3 && check_ultra_fast_convergence(log_sum, prev_log_sum, mean_scale, m)) break;
        prev_log_sum = log_sum;
    }

    // pr0 = exp(log_part1) - exp(log_part2) (se log_sum finito e log_part1>log_part2)
    double log_part1 = -mean_i + n * log_mean_i - fast_lgamma(n);
    if (log_sum == -INFINITY) {
        // secondo termine ~ 0
        return log_part1;
    }

    double log_part2 = -auxi + log_sum;
    return logdiffexp_local(log_part1, log_part2);
}

/*******************************************************************************/
/*  P00  (ora ritorna LOG(P00)) */
/*******************************************************************************/
double P00(double corr, int r, int t, double mean_i, double mean_j) {
    (void)r; (void)t; // r=t=0 nel tuo dispatch, ma non lo uso qui

    if (mean_i <= 0.0 || mean_j <= 0.0) return -INFINITY;

    const double rho2 = R_pow(corr, 2);
    const double one_minus_rho2 = 1.0 - rho2;
    if (one_minus_rho2 <= 1e-15 || rho2 <= 0.0) return -INFINITY;

    const double log_rho2 = log(rho2);
    const double log_one_minus_rho2 = log1p(-rho2);
    const double auxi = mean_i / one_minus_rho2;
    const double auxj = mean_j / one_minus_rho2;
    const double mean_scale = mean_i + mean_j;

    const int max_iter = estimate_iterations(mean_scale, 1500);

    // sum = Σ_k rho2^k * igam(k+1,auxi) * igam(k+1,auxj)
    // in log: log_term = k log_rho2 + log_igam_i + log_igam_j
    double log_sum = -INFINITY;
    double prev_log_sum = -INFINITY;
    int consecutive_fails = 0;

    for (int k = 0; k < max_iter; ++k) {
        int kp1 = k + 1;

        double log_igam_i = safe_log_igam_fast((double)kp1, auxi);
        double log_igam_j = safe_log_igam_fast((double)kp1, auxj);

        if (log_igam_i == -INFINITY || log_igam_j == -INFINITY) {
            consecutive_fails++;
            if (consecutive_fails > 10) break;
            continue;
        }

        double log_term = k * log_rho2 + log_igam_i + log_igam_j;

        if (log_term < -500) {
            consecutive_fails++;
            if (consecutive_fails > 15) break;
            continue;
        }

        consecutive_fails = 0;
        log_sum = (log_sum == -INFINITY) ? log_term : log_sum_exp(log_sum, log_term);

        if (k > 3 && check_ultra_fast_convergence(log_sum, prev_log_sum, mean_scale, k)) break;
        prev_log_sum = log_sum;
    }

    // p00 = -1 + exp(-mi) + exp(-mj) + (1-rho2) * exp(log_sum)
    //     = ( exp(-mi) + exp(-mj) + exp(log1p(-rho2)+log_sum) ) - 1
    double a = -mean_i;
    double b = -mean_j;
    double c = (log_sum == -INFINITY) ? -INFINITY : (log_one_minus_rho2 + log_sum);

    // log(p00) = log( exp(a)+exp(b)+exp(c) - 1 )
    return log_add3_minus1(a, b, c);
}
/*******************************************************************************/
/*  Funzione principale: ora ritorna LOG densità */
/*******************************************************************************/
double biv_Poisson(double corr, int r, int t, double mean_i, double mean_j) {
    if (mean_i <= 0.0 || mean_j <= 0.0 || r < 0 || t < 0) return -INFINITY;
    if (fabs(corr) >= 1.0) return -INFINITY;

    double mean_scale = mean_i + mean_j;
    double corr_threshold = (mean_scale > 20.0) ? 1e-4 : (mean_scale > 5.0) ? 1e-5 : 1e-6;

    if (fabs(corr) <= corr_threshold) {
        // caso indipendente: log pmf prodotto
        double log_dens1 = -mean_i + r * log(mean_i) - fast_lgamma(r);
        double log_dens2 = -mean_j + t * log(mean_j) - fast_lgamma(t);
        double out = log_dens1 + log_dens2;
        return (R_finite(out)) ? out : -INFINITY;
    }

    double logdens = -INFINITY;

    if (r == t) {
        if (r == 0) {
            logdens = P00(corr, r, r, mean_i, mean_j);
        } else {
            logdens = Prr(corr, r, r, mean_i, mean_j);
        }
    } else {
        /* Exact positive single-series representation for every off-diagonal
         * case.  It unifies the former Pr0 and Prt branches.  The legacy
         * formulas are retained below as a numerical fallback only. */
        if (r > t) {
            logdens = Prt_single(corr, r, t, mean_i, mean_j);
            if (!R_finite(logdens)) {
                logdens = (t == 0) ? Pr0(corr, r, t, mean_i, mean_j)
                                   : Prt(corr, r, t, mean_i, mean_j);
            }
        } else {
            logdens = Prt_single(corr, t, r, mean_j, mean_i);
            if (!R_finite(logdens)) {
                logdens = (r == 0) ? Pr0(corr, t, r, mean_j, mean_i)
                                   : Prt(corr, t, r, mean_j, mean_i);
            }
        }
    }

    return (R_finite(logdens)) ? logdens : -INFINITY;
}
/******************************************************************************
 * (Opzionale) wrapper se ti serve ancora la densità “in scala normale”
 *
 * double biv_Poisson_dens(double corr, int r, int t, double mean_i, double mean_j) {
 *     double lp = biv_Poisson(corr, r, t, mean_i, mean_j);
 *     if (lp == -INFINITY) return 0.0;
 *     // clamp overflow
 *     if (lp > 700) return INFINITY;
 *     return exp(lp);
 * }
 ******************************************************************************/
double biv_PoissonZIP(double corr,int r, int t, double mean_i, double mean_j,double mup,double nugget1,double nugget2)
{
double dens=0.0,p,p00,p10,p01,p11;
p=pnorm(mup,0,1,1,0);
p00=pbnorm22(mup,mup,(1-nugget2)*corr);
p01=p-p00;
p10=p01;
p11=1-2*p+p00;

//Rprintf("%f %f \n",nugget1,nugget2);

if(r==0&&t==0)
     dens=p00  + p01*exp(-mean_i) + p10*exp(-mean_j)+p11*biv_Poisson((1-nugget1)*corr,0, 0, mean_i, mean_j);
if(r==0&&t>0)
      //dens=      p01*  exp(-mean_i+t*log(mean_i)-lgammafn(t+1))  + p11*biv_Poisson((1-nugget1)*corr,0, t, mean_i, mean_j);
      dens=      p01*  exp(-mean_j+t*log(mean_j)-lgammafn(t+1))  + p11*biv_Poisson((1-nugget1)*corr,0, t, mean_i, mean_j);
if(r>0&&t==0)
      //dens=      p10*  exp(-mean_j+r*log(mean_j)-lgammafn(r+1))  + p11*biv_Poisson((1-nugget1)*corr,r, 0, mean_i, mean_j);
      dens=      p10*  exp(-mean_i+r*log(mean_i)-lgammafn(r+1))  + p11*biv_Poisson((1-nugget1)*corr,r, 0, mean_i, mean_j);
if(r>0&&t>0)
      dens=      p11*biv_Poisson((1-nugget1)*corr,r, t, mean_i, mean_j);
return(dens);

}

double biv_Mis_PoissonZIP(double corr, double data_i, double data_j,
                          double mean_i, double mean_j, double mup,
                          double nugget1, double nugget2) {
    const int N = 2;
    double M[2][2];
    double dat[2];
    double dens = 0.0;

    double p = pnorm(mup, 0, 1, 1, 0);
    double one_minus_nug2_corr = (1.0 - nugget2) * corr;
    double one_minus_nug1_corr = (1.0 - nugget1) * corr;

    double p00 = pbnorm22(mup, mup, one_minus_nug2_corr);
    double p01 = p - p00;
    double p10 = p01;
    double p11 = 1.0 - 2.0 * p + p00;

    double corr1 = corr_pois(one_minus_nug1_corr, mean_i, mean_j);

    double sqrt_mean_i = sqrt(mean_i);
    double sqrt_mean_j = sqrt(mean_j);

    M[0][0] = mean_i;
    M[1][1] = mean_j;
    M[0][1] = M[1][0] = sqrt(mean_i * mean_j) * corr1;

    dat[0] = data_i - mean_i;
    dat[1] = data_j - mean_j;

    double pdf2 = dNnorm(N, (double **)M, dat); // Bivariate normal
    double pdf1i = dnorm(data_i, mean_i, sqrt_mean_i, 0); // Univariate marginals
    double pdf1j = dnorm(data_j, mean_j, sqrt_mean_j, 0);

    if (data_i > 0.0 && data_j > 0.0)
        dens = p11 * pdf2;
    else if (data_i > 0.0 && data_j == 0.0)
        dens = p11 * pdf2 + p10 * pdf1i;
    else if (data_i == 0.0 && data_j > 0.0)
        dens = p11 * pdf2 + p01 * pdf1j;
    else // data_i == 0.0 && data_j == 0.0
        dens = p11 * pdf2 + p01 * pdf1i + p10 * pdf1j + p00;

    return dens;
}

double biv_binomnegZINB(int N, double corr, int r, int t, double mean_i, double mean_j,
    double nugget1, double nugget2, double mup,
    double *lgamma_cache, int cache_size)
{
    double dens = 0.0, ap, ap00, ap01, ap10, ap11;
 
    ap    = pnorm(mup, 0, 1, 1, 0);
    ap00  = pbnorm22(mup, mup, (1 - nugget2) * corr);
    ap01  = ap - ap00;
    ap10  = ap01;
    ap11  = 1 - 2 * ap + ap00;
 
    double p11 = pbnorm22(mean_i, mean_j, (1 - nugget1) * corr);
    double p1  = pnorm(mean_i, 0, 1, 1, 0);
    double p2  = pnorm(mean_j, 0, 1, 1, 0);
 
    if (r == 0 && t == 0)
        dens = ap00
             + ap01 * R_pow(p1, N)
             + ap10 * R_pow(p2, N)
             + ap11 * biv_binomneg(N, 0, 0, p1, p2, p11, lgamma_cache, cache_size);
 
    if (r == 0 && t > 0)
        dens = ap01 * exp(lgammafn(N + t) - lgammafn(t + 1) - lgammafn(N)
                          + N * log(p2) + t * log1p(-p2))
             + ap11 * biv_binomneg(N, 0, t, p1, p2, p11, lgamma_cache, cache_size);
 
    if (r > 0 && t == 0)
        dens = ap10 * exp(lgammafn(N + r) - lgammafn(r + 1) - lgammafn(N)
                          + N * log(p1) + r * log1p(-p1))
             + ap11 * biv_binomneg(N, r, 0, p1, p2, p11, lgamma_cache, cache_size);
 
    if (r > 0 && t > 0)
        dens = ap11 * biv_binomneg(N, r, t, p1, p2, p11, lgamma_cache, cache_size);
 
    return dens;
}


/* binary misspecified negative binomial*/
double biv_binegbinary (int NN, int u, int v,double pu,double pv, double p11)
{
double res=0.0;double u_0,v_0,p00;

u_0=R_pow(pu,NN);
v_0=R_pow(pv,NN);
p00=R_pow(p11,NN);

if(u==0&&v==0)  res=p00;
if(u==0&&v>0)   res=u_0-p00;
if(v==0&&u>0)   res=v_0-p00;
if(u>0&&v>0)    res=1-(u_0-p00+v_0);
return(res);
}
/*
static double hyp2F1_trunc_recur(double a1, int m, double c, double z) {
    if (m < 0) return NA_REAL;
    // T0 = 1
    double term = 1.0;
    double sum  = 1.0;
    double comp = 0.0;

    for (int r = 1; r <= m; ++r) {
        // T_r = T_{r-1} * ((a+r-1)*(-m+r-1))/((c+r-1)*r) * z
        double num1 = a1 + (r - 1);
        double num2 = -m + (r - 1);   // alterna signo
        double den1 = c  + (r - 1);
        double den2 = (double)r;

        term *= (num1 * num2 / (den1 * den2)) * z;

        // Kahan summation
        double y = term - comp;
        double t = sum + y;
        comp = (t - sum) - y;
        sum = t;

        if (!R_finite(term)) break;
    }
    return sum;
}*/



// 2F1(a1, -m; c; z) con m >= 0 entero.
// Serie truncada con recurrencia:
// T0=1,  T_r = T_{r-1} * ((a1+r-1)*(-m+r-1))/((c+r-1)*r) * z
/*
static double hyp2F1_trunc_recur_general(double a1, int m, double c, double z) {
    if (m < 0) return NA_REAL;
    double term = 1.0;
    double sum  = 1.0;
    double comp = 0.0; // Kahan

    for (int r = 1; r <= m; ++r) {
        double num1 = a1 + (r - 1);
        double num2 = -m + (r - 1);
        double den1 = c  + (r - 1);
        double den2 = (double)r;
        term *= (num1 * num2 / (den1 * den2)) * z;

        // Kahan summation
        double y = term - comp;
        double t = sum + y;
        comp = (t - sum) - y;
        sum = t;

        if (!R_finite(term)) break;
    }
    return sum; 
}
*/

// ---------------------------------------------------------------------------



/********************************************************************/
/********* poisson gamma bivariate distribution *********************/
/********************************************************************/



#ifndef PG_PROB_MIN
#define PG_PROB_MIN 1e-300
#endif

#ifndef PG_RHO_EPS
#define PG_RHO_EPS 1e-15
#endif

#ifndef PG_TERM_TOL
#define PG_TERM_TOL 1e-14
#endif

#ifndef PG_CONV_TOL
#define PG_CONV_TOL 1e-13
#endif

#ifndef PG_INDEP_TOL
#define PG_INDEP_TOL 1e-8
#endif
/* 2F1(a1, -m; c; z), m >= 0 integer. Standard, not regularized. */
static double hyp2F1_trunc_recur_general(double a1, int m, double c, double z)
{
    if (m < 0 || !R_FINITE(a1) || !R_FINITE(c) || !R_FINITE(z)) return R_NaN;

    double term = 1.0;
    double sum  = 1.0;
    double comp = 0.0;

    for (int r = 1; r <= m; ++r) {
        const double num1 = a1 + (double)(r - 1);
        const double num2 = -((double)m) + (double)(r - 1);
        const double den1 = c + (double)(r - 1);
        const double den2 = (double)r;

        if (den1 == 0.0 || den2 == 0.0) return R_NaN;

        term *= (num1 * num2 / (den1 * den2)) * z;

        if (!R_FINITE(term)) return R_NaN;

        const double y = term - comp;
        const double tt = sum + y;
        comp = (tt - sum) - y;
        sum = tt;
    }

    return sum;
}

/* If you already have hyp2F1_trunc_recur, keep yours.
   This wrapper makes calls consistent. */
static R_INLINE double hyp2F1_trunc_recur_1(double a1, int m, double c, double z)
{
    return hyp2F1_trunc_recur_general(a1, m, c, z);
}

/* ------------------------------------------------------------------ */
/* Case r = t = 0                                                     */
/* ------------------------------------------------------------------ */

/* Corrected PEM bivariate probabilities, Theorem 5.2.
 *
 * The implementation below follows the corrected formulas directly.
 * The triple sums in the diagonal and off-diagonal cases are replaced by
 * their exact diagonal reductions:
 *
 *   sum_{l=0}^s (n)_l/l! = (n+1)_s/s!,
 *   sum_{l=0}^s (c)_l/l! = (c+1)_s/s!.
 *
 * No special-function implementation is modified.  hypergeo() is the
 * ordinary Gauss 2F1, hence Gamma(c)^(-1) factors from the regularized
 * functions are included explicitly in the logarithmic coefficients.
 */

/* Series limits are defined with the PEM correlation helpers above. */

static R_INLINE int pg_shape_is_integer(double a)
{
    if (!R_FINITE(a)) return 0;
    return fabs(a - nearbyint(a)) <= 1e-12 * fmax2(1.0, fabs(a));
}



static R_INLINE double pg_hyp_bparam(double a1, double bpar, double c,
                                     double z, int mpoly)
{
    if (mpoly >= 0)
        return hyp2F1_trunc_recur_1(a1, mpoly, c, z);
    return hypergeo(a1, bpar, c, z);
}

/*
 * Sequential evaluator for F(a1, 1-shape-l1; c; z).
 *
 * For moderately sized negative z, neighboring values in l1 are related by
 * the exact Gauss contiguous relation
 *
 *   (c-b) F(a,b-1;c;z)
 * + {2b-c+(a-b)z} F(a,b;c;z)
 * + b(z-1) F(a,b+1;c;z) = 0.
 *
 * We re-anchor every PG_HREC_BLOCK values with two calls to the unchanged
 * anchor path.  For integer shape these anchors are terminating 2F1
 * polynomials and are used throughout the audited negative-z range.  For
 * non-integer shape and sufficiently negative z, the first two values are
 * instead obtained from the positive Euler integral above, after which the
 * same exact contiguous recurrence is used through the row.  If a generated
 * value is non-finite, the current and following value are re-anchored by the
 * corresponding unchanged/direct or positive-integral path.
 */
#define PG_HREC_BLOCK 8
#define PG_HREC_ZMAX_NONINT 0.49
/*
 * Integer-shape rows use terminating 2F1 polynomials as exact periodic
 * anchors.  High-precision audits show that the same short 8-value
 * re-anchored contiguous blocks remain stable well beyond |z|=1.
 * Keep a conservative finite cap: outside it we preserve the previous
 * direct-polynomial path.  No special-function implementation is changed.
 */
#define PG_HREC_ZMAX_INT 50.0

typedef struct {
    double a1, shape, c, z;
    int a_int;
    int enabled;
    int use_integral_seeds;
    int force_direct_next;
    int last_l1;
    double fm2, fm1;
} PGHypBSeq;

/*
 * Robust positive-integral evaluator used only to seed the l1-contiguous
 * recurrence when the Gauss argument is too negative for the direct
 * hypergeo() call to be reliable for non-integer shape.  This does not
 * modify hypergeo() (or igam()) in any way.
 *
 * For A>0, C-A>0 and z<0, Euler's integral gives
 *
 *   2F1(A,B;C;z) = Gamma(C)/(Gamma(A)Gamma(C-A))
 *       integral_0^1 t^(A-1) (1-t)^(C-A-1) (1-z t)^(-B) dt.
 *
 * In PGrt, A and C-A are positive integers.  The integrand is therefore
 * positive and smooth on [0,1].  We integrate a mode-scaled version with
 * R's QUADPACK wrapper Rdqags, then restore the normalization in log scale.
 * Only the first two l1 values of a dangerous row are integrated; all
 * subsequent values are generated by the exact contiguous recurrence.
 */
#define PG_H1_QUAD_LIMIT 100

typedef struct {
    double A, N, B, z, shift;
} PGH1QuadCtx;

static R_INLINE double pg_h1_logkernel(double t, const PGH1QuadCtx *q)
{
    double out = 0.0;
    if (t <= 0.0) {
        if (q->A > 1.0) return R_NegInf;
    } else {
        out += (q->A - 1.0) * log(t);
    }
    if (t >= 1.0) {
        if (q->N > 1.0) return R_NegInf;
    } else {
        out += (q->N - 1.0) * log1p(-t);
    }
    out += (-q->B) * log1p(-q->z * t);
    return out;
}

static void pg_h1_quad_integrand(double *xx, int n, void *ex)
{
    PGH1QuadCtx *q = (PGH1QuadCtx *)ex;
    for (int i = 0; i < n; ++i) {
        const double lk = pg_h1_logkernel(xx[i], q);
        if (lk == R_NegInf) xx[i] = 0.0;
        else if (!R_FINITE(lk)) xx[i] = R_NaN;
        else {
            const double d = lk - q->shift;
            xx[i] = (d < -745.0) ? 0.0 : exp(d);
        }
    }
}

static double pg_hyp_integral_seed(double A, double B, double C, double z)
{
    const double N = C - A;
    if (!(A > 0.0) || !(N > 0.0) || !(z < 0.0) ||
        !R_FINITE(A) || !R_FINITE(B) || !R_FINITE(C) || !R_FINITE(z))
        return R_NaN;

    PGH1QuadCtx q;
    q.A = A; q.N = N; q.B = B; q.z = z; q.shift = R_NegInf;

    /* The derivative of the log-kernel has a quadratic numerator. */
    {
        const double aa = A - 1.0;
        const double nn = N - 1.0;
        const double mm = -B;
        const double qq = -z;
        const double a0 = aa;
        const double a1 = aa * (qq - 1.0) - nn + mm * qq;
        const double a2 = -qq * (aa + nn + mm);
        double cand[4];
        int nc = 0;
        if (aa == 0.0) cand[nc++] = 0.0;
        if (nn == 0.0) cand[nc++] = 1.0;
        if (fabs(a2) <= 1e-15 * fmax2(1.0, fabs(a1))) {
            if (a1 != 0.0) {
                const double r = -a0 / a1;
                if (r > 0.0 && r < 1.0) cand[nc++] = r;
            }
        } else {
            const double disc = a1*a1 - 4.0*a2*a0;
            if (disc >= 0.0 && R_FINITE(disc)) {
                const double sd = sqrt(disc);
                const double r1 = (-a1 + sd) / (2.0*a2);
                const double r2 = (-a1 - sd) / (2.0*a2);
                if (r1 > 0.0 && r1 < 1.0) cand[nc++] = r1;
                if (r2 > 0.0 && r2 < 1.0 && nc < 4) cand[nc++] = r2;
            }
        }
        /* Always include a safe interior point in case of roundoff in roots. */
        if (nc < 4) cand[nc++] = A / C;
        for (int i = 0; i < nc; ++i) {
            const double lk = pg_h1_logkernel(cand[i], &q);
            if (R_FINITE(lk) && (q.shift == R_NegInf || lk > q.shift)) q.shift = lk;
        }
    }
    if (!R_FINITE(q.shift)) return R_NaN;

    double lower = 0.0, upper = 1.0;
    double epsabs = 1e-13, epsrel = 2e-12;
    double result = 0.0, abserr = 0.0;
    int neval = 0, ier = 0, limit = PG_H1_QUAD_LIMIT, lenw = 4*PG_H1_QUAD_LIMIT, last = 0;
    int iwork[PG_H1_QUAD_LIMIT];
    double work[4*PG_H1_QUAD_LIMIT];
    Rdqags(pg_h1_quad_integrand, &q, &lower, &upper,
           &epsabs, &epsrel, &result, &abserr, &neval, &ier,
           &limit, &lenw, &last, iwork, work);
    if (ier > 1 || !(result > 0.0) || !R_FINITE(result)) return R_NaN;

    {
        const double logF = lgammafn(C) - lgammafn(A) - lgammafn(N) +
                            q.shift + log(result);
        if (!R_FINITE(logF) || logF > 700.0) return R_NaN;
        if (logF < -745.0) return 0.0;
        return exp(logF);
    }
}

static R_INLINE double pg_hyp_bseq_anchor(PGHypBSeq *st, int l1)
{
    const double bpar = 1.0 - (double)l1 - st->shape;
    const int mpoly = st->a_int ? (l1 + st->a_int - 1) : -1;
    if (st->use_integral_seeds)
        return pg_hyp_integral_seed(st->a1, bpar, st->c, st->z);
    return pg_hyp_bparam(st->a1, bpar, st->c, st->z, mpoly);
}

static R_INLINE void pg_hyp_bseq_init(PGHypBSeq *st, double a1, double shape,
                                      double c, double z, int a_int)
{
    st->a1 = a1;
    st->shape = shape;
    st->c = c;
    st->z = z;
    st->a_int = a_int;
    st->use_integral_seeds = (!a_int && R_FINITE(z) && z < -PG_HREC_ZMAX_NONINT);
    {
        const double zmax = a_int ? PG_HREC_ZMAX_INT : PG_HREC_ZMAX_NONINT;
        st->enabled = st->use_integral_seeds ||
                      (R_FINITE(z) && z <= 0.0 && fabs(z) <= zmax);
    }
    st->force_direct_next = 0;
    st->last_l1 = -1;
    st->fm2 = st->fm1 = R_NaN;
}

static R_INLINE double pg_hyp_bseq_get(PGHypBSeq *st, int l1)
{
    double val;

    /* The PGrt loops call this strictly in increasing l1 order. */
    if (!st->enabled || l1 != st->last_l1 + 1) {
        val = pg_hyp_bseq_anchor(st, l1);
        st->fm2 = st->fm1;
        st->fm1 = val;
        st->last_l1 = l1;
        return val;
    }

    if ((st->use_integral_seeds ? (l1 < 2) : ((l1 % PG_HREC_BLOCK) < 2)) ||
        st->force_direct_next || l1 < 2) {
        val = pg_hyp_bseq_anchor(st, l1);
        if (st->force_direct_next) st->force_direct_next = 0;
    } else {
        /* Current b is the parameter of F at l1-1. */
        const double bcur = 1.0 - (double)(l1 - 1) - st->shape;
        const double den = st->c - bcur;
        val = -((2.0*bcur - st->c + (st->a1 - bcur)*st->z) * st->fm1
                + bcur*(st->z - 1.0) * st->fm2) / den;
        if (!R_FINITE(val)) {
            val = pg_hyp_bseq_anchor(st, l1);
            st->force_direct_next = 1;
        }
    }

    st->fm2 = st->fm1;
    st->fm1 = val;
    st->last_l1 = l1;
    return val;
}



static R_INLINE double pg_finish_prob(double p)
{
    if (!R_FINITE(p)) return pg_fail(PG_STATUS_INVALID_PROBABILITY);
    if (p < -PG_PROB_TOL || p > 1.0 + PG_PROB_TOL)
        return pg_fail(PG_STATUS_INVALID_PROBABILITY);
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;
    return p;
}


static R_INLINE double pg_exp_logprob(double lp);

/* ================================================================== */
/* Positive PEM bivariate probabilities.                              */
/*                                                                    */
/* These kernels implement the cancellation-free representations used */
/* in Theorem 5.3 of the PEM manuscript:                              */
/*   - p00: positive double series in regularized incomplete beta;    */
/*   - pnn: positive double series, with the finite inner sum updated  */
/*          recursively in O(1) work per h;                           */
/*   - pnm: one positive double series for every n != m, including    */
/*          zero/nonzero pairs.                                       */
/*                                                                    */
/* The legacy kernels below are intentionally retained as an internal  */
/* numerical reference and are exposed only through a diagnostic       */
/* native entry point.                                                 */
/* ================================================================== */

#ifndef PG_POS_MAX_OUTER
#define PG_POS_MAX_OUTER 3000
#endif
#ifndef PG_POS_MAX_INNER
#define PG_POS_MAX_INNER 5000
#endif
#ifndef PG_POS_HYP_MAX
#define PG_POS_HYP_MAX 20000
#endif
#ifndef PG_POS_QUIET
#define PG_POS_QUIET 8
#endif

static R_INLINE int pg_pos_logterm_small(double logterm, double logsum)
{
    if (logsum == R_NegInf || ISNAN(logsum) || ISNAN(logterm)) return 0;
    if (logterm == R_NegInf) return 1;
    if (!R_FINITE(logsum) || !R_FINITE(logterm)) return 0;
    return logterm <= logsum + log(PG_SUM_TOL);
}

/* Direct positive-series evaluation of the regularized Gauss function
 *
 *   2F1~(aa,bb;cc;z)
 *     = sum_{j>=0} (aa)_j (bb)_j z^j / {Gamma(cc+j) j!}.
 *
 * The ordinary hypergeo() path is used when it is finite; otherwise the
 * positive series is accumulated in log scale.  This avoids overflow of the
 * ordinary 2F1 for large counts while preserving a fast path in routine
 * cases. */
static double pg_log_hyp2f1_reg_positive(double aa, double bb, double cc,
                                          double z, int *ok)
{
    if (ok) *ok = 0;
    if (!(aa > 0.0) || !(bb > 0.0) || !(cc > 0.0) ||
        !(z >= 0.0 && z < 1.0) ||
        !R_FINITE(aa) || !R_FINITE(bb) || !R_FINITE(cc) || !R_FINITE(z))
        return R_NaN;

    if (z == 0.0) {
        if (ok) *ok = 1;
        return -lgammafn(cc);
    }

    {
        const double direct = hypergeo(aa, bb, cc, z);
        if (direct > 0.0 && R_FINITE(direct)) {
            const double ans = log(direct) - lgammafn(cc);
            if (R_FINITE(ans)) {
                if (ok) *ok = 1;
                return ans;
            }
        }
    }

    const double logz = log(z);
    double logterm = -lgammafn(cc);
    double logsum = logterm;
    int quiet = 0;

    for (int j = 0; j < PG_POS_HYP_MAX; ++j) {
        logterm += log(aa + (double)j) + log(bb + (double)j) -
                   log(cc + (double)j) - log((double)j + 1.0) + logz;
        logsum = pg_logspace_add_pos(logsum, logterm);
        if (ISNAN(logsum)) return R_NaN;

        if (j > 10 && pg_pos_logterm_small(logterm, logsum)) {
            if (++quiet >= PG_POS_QUIET) {
                if (ok) *ok = 1;
                return logsum;
            }
        } else quiet = 0;
    }

    return R_NaN;
}

/* Positive p00 representation.  All summands are probabilities, so the
 * final double sum is accumulated in ordinary scale with compensated
 * summation.  The incomplete-beta sequence itself remains in log scale to
 * retain stability when z is small. */
static double PG00_positive(double corr, double mean_i, double mean_j,
                            double a)
{
    pg_set_status(PG_STATUS_OK);
    const double q = corr * corr;
    const double delta = 1.0 - q;
    const double nu_i = a / mean_i;
    const double nu_j = a / mean_j;
    const double zi = nu_i / (1.0 + nu_i);
    const double zj = nu_j / (1.0 + nu_j);
    if (!(q > 0.0 && q < 1.0) || !(nu_i > 0.0) || !(nu_j > 0.0) ||
        !(zi > 0.0 && zi < 1.0) || !(zj > 0.0 && zj < 1.0))
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    const double logq = log(q);
    const double logzi = log(zi);
    const double logzj = log(zj);
    const double logpref = (a + 1.0) * log(delta);
    const int stationary =
        fabs(zi - zj) <= 1e-14 * fmax2(1.0, fmax2(zi, zj));

    KahanSum total;
    kahan_init(&total);
    double logellw = 0.0; /* q^ell (a)_ell / ell! */
    int quiet_ell = 0;
    int outer_converged = 0;

    for (int ell = 0; ell < PG_POS_MAX_OUTER; ++ell) {
        const double cshape = a + (double)ell;
        PGBetaBSeq bi, bj;
        const double li0 = cshape * logzi;
        const double lj0 = stationary ? li0 : cshape * logzj;
        const double lbeta0 = -log(cshape); /* B(cshape,1) */

        if (!pg_beta_bseq_init_seeded(&bi, cshape, 1.0, zi, li0, lbeta0) ||
            (!stationary &&
             !pg_beta_bseq_init_seeded(&bj, cshape, 1.0, zj, lj0, lbeta0)))
            return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

        KahanSum row;
        kahan_init(&row);
        int quiet_k = 0;
        int row_converged = 0;
        const double rowbase = logpref + logellw;

        for (int k = 0; k < PG_POS_MAX_INNER; ++k) {
            const double li = pg_beta_bseq_get(&bi, k);
            const double lj = stationary ? li : pg_beta_bseq_get(&bj, k);
            if (ISNAN(li) || ISNAN(lj))
                return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

            const double logterm =
                (li == R_NegInf || lj == R_NegInf) ? R_NegInf :
                rowbase + (double)k * logq + li + lj;
            const double term = pg_exp_log_signed(logterm, +1);
            if (!R_FINITE(term))
                return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
            kahan_add(&row, term);

            if (k > 10) {
                const double scale = fmax2(1e-300, row.sum);
                if (term <= PG_SUM_TOL * scale) {
                    if (++quiet_k >= PG_POS_QUIET) {
                        row_converged = 1;
                        break;
                    }
                } else quiet_k = 0;
            }
        }
        if (!row_converged) return pg_fail(PG_STATUS_MAXITER);

        kahan_add(&total, row.sum);
        if (ell > 10) {
            const double scale = fmax2(1e-300, total.sum);
            if (row.sum <= PG_SUM_TOL * scale) {
                if (++quiet_ell >= PG_POS_QUIET) {
                    outer_converged = 1;
                    break;
                }
            } else quiet_ell = 0;
        }

        if (ell + 1 < PG_POS_MAX_OUTER) {
            logellw += logq + log(a + (double)ell) - log((double)ell + 1.0);
        }
    }

    if (!outer_converged) return pg_fail(PG_STATUS_MAXITER);
    return pg_finish_prob(total.sum);
}

/* Positive diagonal representation.  For fixed ell define
 *
 *   S_h = sum_{k=0}^h (n)_k/k! Delta_i(k,h,ell) Delta_j(k,h,ell).
 *
 * When h increases by one, every Delta receives the same new NB mass.
 * Maintaining S_h together with the two weighted one-dimensional sums and
 * W_h=sum_{k=0}^h (n)_k/k! eliminates the explicit finite k loop.
 *
 * The recurrence itself is kept in log scale because W_h can be very large.
 * The final positive h- and ell-sums are accumulated in ordinary scale with
 * compensated summation. */
static double PGrr_positive(double corr, int n, double mean_i, double mean_j,
                            double a)
{
    pg_set_status(PG_STATUS_OK);
    const double q = corr * corr;
    const double delta = 1.0 - q;
    const double nu_i = a / mean_i;
    const double nu_j = a / mean_j;
    const double zi = nu_i / (1.0 + nu_i);
    const double zj = nu_j / (1.0 + nu_j);
    if (n < 1 || !(q > 0.0 && q < 1.0) ||
        !(zi > 0.0 && zi < 1.0) || !(zj > 0.0 && zj < 1.0))
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    const double logq = log(q);
    const double log1mzi = log1p(-zi);
    const double log1mzj = log1p(-zj);
    const double logpref = ((double)n + a + 1.0) * log(delta);
    const int stationary =
        fabs(zi - zj) <= 1e-14 * fmax2(1.0, fmax2(zi, zj));

    KahanSum total;
    kahan_init(&total);
    double logellw = 0.0; /* q^ell (a)_ell / ell! */
    int quiet_ell = 0;
    int outer_converged = 0;

    for (int ell = 0; ell < PG_POS_MAX_OUTER; ++ell) {
        const double cshape = a + (double)ell;
        double logpi = dnbinom((double)n, cshape, zi, 1);
        double logpj = stationary ? logpi : dnbinom((double)n, cshape, zj, 1);
        if (ISNAN(logpi) || ISNAN(logpj))
            return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

        double logw = 0.0;  /* (n)_0/0! */
        double logW = 0.0;  /* sum_{k=0}^0 */
        double logAi = logpi;
        double logAj = stationary ? logAi : logpj;
        double logS = logpi + logpj;
        KahanSum row;
        kahan_init(&row);
        int quiet_h = 0;
        int row_converged = 0;
        const double rowbase = logpref + logellw;

        for (int h = 0; h < PG_POS_MAX_INNER; ++h) {
            const double logterm =
                logS == R_NegInf ? R_NegInf :
                rowbase + (double)h * logq + logS;
            const double term = pg_exp_log_signed(logterm, +1);
            if (!R_FINITE(term))
                return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
            kahan_add(&row, term);

            if (h > 10) {
                const double scale = fmax2(1e-300, row.sum);
                if (term <= PG_SUM_TOL * scale) {
                    if (++quiet_h >= PG_POS_QUIET) {
                        row_converged = 1;
                        break;
                    }
                } else quiet_h = 0;
            }

            if (h + 1 < PG_POS_MAX_INNER) {
                const double u = (double)n + (double)h;
                const double oldAi = logAi;
                const double oldAj = logAj;

                logpi += log(cshape + u) - log(u + 1.0) + log1mzi;
                if (stationary) logpj = logpi;
                else logpj += log(cshape + u) - log(u + 1.0) + log1mzj;

                logw += log((double)n + (double)h) - log((double)h + 1.0);
                const double logWnew = pg_logspace_add_pos(logW, logw);

                if (stationary) {
                    double inc = logpi + oldAi + M_LN2;
                    inc = pg_logspace_add_pos(inc,
                                              2.0 * logpi + logWnew);
                    logS = pg_logspace_add_pos(logS, inc);
                    logAi = pg_logspace_add_pos(oldAi, logpi + logWnew);
                    logAj = logAi;
                } else {
                    double inc = pg_logspace_add_pos(logpi + oldAj,
                                                      logpj + oldAi);
                    inc = pg_logspace_add_pos(inc, logpi + logpj + logWnew);
                    logS = pg_logspace_add_pos(logS, inc);
                    logAi = pg_logspace_add_pos(oldAi, logpi + logWnew);
                    logAj = pg_logspace_add_pos(oldAj, logpj + logWnew);
                }
                logW = logWnew;

                if (ISNAN(logS) || ISNAN(logAi) || ISNAN(logAj) || ISNAN(logW))
                    return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
            }
        }
        if (!row_converged) return pg_fail(PG_STATUS_MAXITER);

        kahan_add(&total, row.sum);
        if (ell > 10) {
            const double scale = fmax2(1e-300, total.sum);
            if (row.sum <= PG_SUM_TOL * scale) {
                if (++quiet_ell >= PG_POS_QUIET) {
                    outer_converged = 1;
                    break;
                }
            } else quiet_ell = 0;
        }

        if (ell + 1 < PG_POS_MAX_OUTER) {
            logellw += logq + log(a + (double)ell) - log((double)ell + 1.0);
        }
    }

    if (!outer_converged) return pg_fail(PG_STATUS_MAXITER);
    return pg_finish_prob(total.sum);
}

/* Positive off-diagonal representation, valid for n > m >= 0.
 *
 * The manuscript formula contains
 *
 *   2F1~(A, n+k+ell+a; n+k+1; z),  0<z<1.
 *
 * A Pfaff transformation gives
 *
 *   (1-z)^(-A) / Gamma(n+k+1)
 *   2F1(A, 1-ell-a; n+k+1; z/(z-1)).
 *
 * For fixed k, only the second hypergeometric parameter changes with ell.
 * We therefore make k the outer index and reuse the existing contiguous
 * PGHypBSeq recurrence.  This avoids a fresh hypergeo() call at every
 * (k,ell) pair and is substantially faster for the off-diagonal cases. */
static double PGrt_positive(double corr, int n, int m,
                            double mean_i, double mean_j, double a)
{
    pg_set_status(PG_STATUS_OK);
    const double q = corr * corr;
    const double delta = 1.0 - q;
    const double nu_i = a / mean_i;
    const double nu_j = a / mean_j;
    if (n <= m || m < 0 || !(q > 0.0 && q < 1.0) ||
        !(nu_i > 0.0) || !(nu_j > 0.0))
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    const double logq = log(q);
    const double logxy = log(nu_i) + log(nu_j);
    const double log1px = log1p(nu_i);
    const double log1py = log1p(nu_j);
    const double bx = 1.0 + nu_i - q;
    if (!(bx > 0.0)) return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    const double z = q / (1.0 + nu_i);
    const double z1 = -q / bx; /* z/(z-1) */
    const double A = (double)(n - m + 1);
    const int a_int = pg_shape_is_integer(a) ? (int)llround(a) : 0;

    /* Include the global prefactor and the Pfaff factor in every term.
       Because the series is positive and sums to a probability, individual
       terms are safely accumulated in ordinary scale. */
    const double logglobal =
        ((double)n + a) * log(delta) -
        lgammafn((double)m + 1.0) - lgammafn(a) -
        A * log1p(-z);

    double logcoef_k0 = logglobal +
        a * logxy +
        lgammafn((double)m + a) + lgammafn((double)n + a) -
        lgammafn(a) - lgammafn((double)n + 1.0) -
        ((double)n + a) * log1px -
        ((double)m + a) * log1py;
    if (!R_FINITE(logcoef_k0))
        return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

    KahanSum total;
    kahan_init(&total);
    int quiet_k = 0;
    int outer_converged = 0;

    for (int k = 0; k < PG_POS_MAX_OUTER; ++k) {
        const double C = (double)n + (double)k + 1.0;
        const int exact_ca = (m == 0 && k == 0);
        const int use_positive_seeds =
            (!a_int && fabs(z1) > PG_HREC_ZMAX_NONINT && !exact_ca);
        PGHypBSeq hseq;
        if (!use_positive_seeds && !exact_ca)
            pg_hyp_bseq_init(&hseq, A, a, C, z1, a_int);

        /* For large negative Pfaff arguments and half-integer shape, the
           legacy sequence uses two numerical quadratures per k-row.  Here
           the first two values are instead seeded from the original
           positive-z regularized series, then the same exact contiguous
           recurrence is used for all following ell. */
        double hm2 = R_NaN, hm1 = R_NaN;
        const double log_gamma_C = lgammafn(C);
        const double log_pfaff_back = A * log1p(-z);

        KahanSum row;
        kahan_init(&row);
        int quiet_ell = 0;
        int row_converged = 0;
        double logcoef = logcoef_k0;

        for (int ell = 0; ell < PG_POS_MAX_INNER; ++ell) {
            double logh;

            if (exact_ca) {
                const double B = 1.0 - (double)ell - a;
                logh = -B * log1p(-z1);
            } else if (use_positive_seeds) {
                double h1;
                if (ell < 2) {
                    int seed_ok = 0;
                    const double Borig = C + (double)ell + a - 1.0;
                    const double logreg = pg_log_hyp2f1_reg_positive(
                        A, Borig, C, z, &seed_ok);
                    if (!seed_ok || ISNAN(logreg))
                        return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
                    logh = logreg + log_gamma_C + log_pfaff_back;
                    if (!R_FINITE(logh) || logh > 700.0)
                        return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
                    h1 = logh < -745.0 ? 0.0 : exp(logh);
                } else {
                    const double bcur = 1.0 - (double)(ell - 1) - a;
                    const double den = C - bcur;
                    h1 = -((2.0*bcur - C + (A - bcur)*z1) * hm1 +
                           bcur*(z1 - 1.0) * hm2) / den;
                    if (!(h1 > 0.0) || !R_FINITE(h1)) {
                        int seed_ok = 0;
                        const double Borig = C + (double)ell + a - 1.0;
                        const double logreg = pg_log_hyp2f1_reg_positive(
                            A, Borig, C, z, &seed_ok);
                        if (!seed_ok || ISNAN(logreg))
                            return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
                        logh = logreg + log_gamma_C + log_pfaff_back;
                        if (!R_FINITE(logh) || logh > 700.0)
                            return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
                        h1 = logh < -745.0 ? 0.0 : exp(logh);
                    }
                }
                if (!(h1 > 0.0) || !R_FINITE(h1))
                    return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
                logh = log(h1);
                hm2 = hm1;
                hm1 = h1;
            } else {
                double h1 = pg_hyp_bseq_get(&hseq, ell);
                if (!(h1 > 0.0) || !R_FINITE(h1)) {
                    const double B = 1.0 - (double)ell - a;
                    h1 = hypergeo(A, B, C, z1);
                }
                if (!(h1 > 0.0) || !R_FINITE(h1))
                    return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
                logh = log(h1);
            }

            const double logterm = logcoef + logh;
            const double term = pg_exp_log_signed(logterm, +1);
            if (!R_FINITE(term))
                return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
            kahan_add(&row, term);

            if (ell > 10) {
                const double scale = fmax2(1e-300, row.sum);
                if (term <= PG_SUM_TOL * scale) {
                    if (++quiet_ell >= PG_POS_QUIET) {
                        row_converged = 1;
                        break;
                    }
                } else quiet_ell = 0;
            }

            if (ell + 1 < PG_POS_MAX_INNER) {
                const double ee = (double)ell;
                logcoef += logq + logxy +
                    log((double)m + (double)k + ee + a) +
                    log((double)n + (double)k + ee + a) -
                    log(ee + 1.0) - log(ee + a) -
                    log1px - log1py;
            }
        }
        if (!row_converged) return pg_fail(PG_STATUS_MAXITER);

        kahan_add(&total, row.sum);
        if (k > 10) {
            const double scale = fmax2(1e-300, total.sum);
            if (row.sum <= PG_SUM_TOL * scale) {
                if (++quiet_k >= PG_POS_QUIET) {
                    outer_converged = 1;
                    break;
                }
            } else quiet_k = 0;
        }

        if (k + 1 < PG_POS_MAX_OUTER) {
            const double kk = (double)k;
            logcoef_k0 += logq +
                log((double)m + kk + a) +
                log((double)n + kk + a) -
                log(kk + 1.0) - log((double)n + kk + 1.0) -
                log1px - log1py;
        }
    }

    if (!outer_converged) return pg_fail(PG_STATUS_MAXITER);
    return pg_finish_prob(total.sum);
}



static R_INLINE double pg_nb_logpmf(int z, double mean, double shape)
{
    if (z < 0 || !(mean > 0.0) || !(shape > 0.0) ||
        !R_FINITE(mean) || !R_FINITE(shape)) return R_NaN;
    const double den = shape + mean;
    return lgammafn(shape + (double)z) - lgammafn(shape) -
           lgammafn((double)z + 1.0) +
           shape * (log(shape) - log(den)) +
           (double)z * (log(mean) - log(den));
}

static R_INLINE double pg_exp_logprob(double lp)
{
    if (lp == R_NegInf) return 0.0;
    if (!R_FINITE(lp) || lp > PG_PROB_TOL)
        return pg_fail(PG_STATUS_INVALID_PROBABILITY);
    if (lp < -745.0) return 0.0;
    return pg_finish_prob(exp(lp));
}

double biv_PoissonGamma(double corr, int r, int t,
                        double mean_i, double mean_j, double a)
{
    pg_set_status(PG_STATUS_OK);
    if (!R_FINITE(corr) || !R_FINITE(mean_i) || !R_FINITE(mean_j) ||
        !R_FINITE(a) || fabs(corr) > 1.0 + PG_PROB_TOL ||
        r < 0 || t < 0 || !(mean_i > 0.0) || !(mean_j > 0.0) || !(a > 0.0))
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    const int equal_means =
        fabs(mean_i - mean_j) <=
        1e-12 * fmax2(1.0, fmax2(mean_i, mean_j));
    if (fabs(corr) >= 1.0 - PG_RHO_EPS) {
        if (equal_means) {
            if (r != t) return 0.0;
            return pg_exp_logprob(pg_nb_logpmf(r, 0.5 * (mean_i + mean_j), a));
        }
        corr = copysign(1.0 - PG_RHO_EPS, corr);
    }

    if (fabs(corr) <= PG_INDEP_TOL)
        return pg_exp_logprob(pg_nb_logpmf(r, mean_i, a) +
                              pg_nb_logpmf(t, mean_j, a));

    if (fmax2(mean_i, mean_j) / a <= PG_POISSON_REL_TOL) {
        const double lp = biv_Poisson(corr, r, t, mean_i, mean_j);
        if (lp == R_NegInf) return 0.0;
        if (!R_FINITE(lp)) return pg_fail(PG_STATUS_SPECIAL_FUNCTION);
        return pg_exp_logprob(lp);
    }

    if (r == t) {
        if (r == 0) return PG00_positive(corr, mean_i, mean_j, a);
        return PGrr_positive(corr, r, mean_i, mean_j, a);
    }
    return r > t ? PGrt_positive(corr, r, t, mean_i, mean_j, a)
                 : PGrt_positive(corr, t, r, mean_j, mean_i, a);
}



double biv_PoissonGammaZIP(double corr, int r, int t,
                           double mean_i, double mean_j, double mup,
                           double nugget1, double nugget2, double shape)
{
    pg_set_status(PG_STATUS_OK);
    if (!R_FINITE(corr) || !R_FINITE(mean_i) || !R_FINITE(mean_j) ||
        !R_FINITE(mup) || !R_FINITE(nugget1) || !R_FINITE(nugget2) ||
        !R_FINITE(shape) || fabs(corr) > 1.0 + PG_PROB_TOL ||
        r < 0 || t < 0 || !(mean_i > 0.0) || !(mean_j > 0.0) ||
        !(shape > 0.0) || nugget1 < 0.0 || nugget1 >= 1.0 ||
        nugget2 < 0.0 || nugget2 >= 1.0)
        return pg_fail(PG_STATUS_INVALID_ARGUMENT);

    double rho_z = (1.0 - nugget2) * corr;
    if (rho_z >= 1.0) rho_z = 1.0 - PG_RHO_EPS;
    if (rho_z <= -1.0) rho_z = -1.0 + PG_RHO_EPS;
    const double p = pnorm(mup, 0.0, 1.0, 1, 0);
    const double p00raw = pbnorm22(mup, mup, rho_z);
    if (!R_FINITE(p) || !R_FINITE(p00raw))
        return pg_fail(PG_STATUS_SPECIAL_FUNCTION);

    double p00 = p00raw;
    double p01 = p - p00;
    double p10 = p01;
    double p11 = 1.0 - 2.0 * p + p00;
    if (p00 < 0.0 && p00 > -PG_PROB_TOL) p00 = 0.0;
    if (p01 < 0.0 && p01 > -PG_PROB_TOL) p01 = 0.0;
    if (p10 < 0.0 && p10 > -PG_PROB_TOL) p10 = 0.0;
    if (p11 < 0.0 && p11 > -PG_PROB_TOL) p11 = 0.0;
    if (p00 < 0.0 || p01 < 0.0 || p10 < 0.0 || p11 < 0.0)
        return pg_fail(PG_STATUS_INVALID_PROBABILITY);

    double rho_count = (1.0 - nugget1) * corr;
    if (rho_count >= 1.0) rho_count = 1.0 - PG_RHO_EPS;
    if (rho_count <= -1.0) rho_count = -1.0 + PG_RHO_EPS;
    const double probi = shape / (shape + mean_i);
    const double probj = shape / (shape + mean_j);

    KahanSum dens;
    kahan_init(&dens);
    if (r == 0 && t == 0) {
        kahan_add(&dens, p00);
        kahan_add(&dens, p01 * exp(shape * log(probi)));
        kahan_add(&dens, p10 * exp(shape * log(probj)));
        const double joint = biv_PoissonGamma(rho_count, 0, 0,
                                               mean_i, mean_j, shape);
        if (!R_FINITE(joint)) return joint;
        kahan_add(&dens, p11 * joint);
    } else if (r == 0) {
        kahan_add(&dens, p01 * dnbinom((double)t, shape, probj, 0));
        const double joint = biv_PoissonGamma(rho_count, 0, t,
                                               mean_i, mean_j, shape);
        if (!R_FINITE(joint)) return joint;
        kahan_add(&dens, p11 * joint);
    } else if (t == 0) {
        kahan_add(&dens, p10 * dnbinom((double)r, shape, probi, 0));
        const double joint = biv_PoissonGamma(rho_count, r, 0,
                                               mean_i, mean_j, shape);
        if (!R_FINITE(joint)) return joint;
        kahan_add(&dens, p11 * joint);
    } else {
        const double joint = biv_PoissonGamma(rho_count, r, t,
                                               mean_i, mean_j, shape);
        if (!R_FINITE(joint)) return joint;
        kahan_add(&dens, p11 * joint);
    }
    return pg_finish_prob(dens.sum);
}

int GeoPGPreparePair(double data_i, double data_j,
                     double eta_i, double eta_j,
                     int *r, int *t, double *mean_i, double *mean_j)
{
    if (!r || !t || !mean_i || !mean_j) return -1;
    if (ISNAN(data_i) || ISNAN(data_j)) return 0;
    if (!R_FINITE(data_i) || !R_FINITE(data_j) ||
        !R_FINITE(eta_i) || !R_FINITE(eta_j) ||
        eta_i > 700.0 || eta_j > 700.0) return -1;
    const double ri = nearbyint(data_i), tj = nearbyint(data_j);
    if (ri < 0.0 || tj < 0.0 || ri > INT_MAX || tj > INT_MAX ||
        fabs(data_i - ri) > 1e-10 || fabs(data_j - tj) > 1e-10)
        return -1;
    *mean_i = exp(eta_i);
    *mean_j = exp(eta_j);
    if (!R_FINITE(*mean_i) || !R_FINITE(*mean_j) ||
        !(*mean_i > 0.0) || !(*mean_j > 0.0)) return -1;
    *r = (int)ri;
    *t = (int)tj;
    return 1;
}

int GeoPGPreparePairAt(int index, double data_i, double data_j,
                       double eta_i, double eta_j,
                       int *r, int *t, double *mean_i, double *mean_j)
{
    if (pg_count_n && pg_count1 && pg_count2 && index >= 0 &&
        index < pg_count_n[0]) {
        if (pg_count1[index] == INT_MIN || pg_count2[index] == INT_MIN)
            return 0;
        if (!R_FINITE(eta_i) || !R_FINITE(eta_j) ||
            eta_i > 700.0 || eta_j > 700.0) return -1;
        *mean_i = exp(eta_i);
        *mean_j = exp(eta_j);
        if (!R_FINITE(*mean_i) || !R_FINITE(*mean_j) ||
            !(*mean_i > 0.0) || !(*mean_j > 0.0)) return -1;
        *r = pg_count1[index];
        *t = pg_count2[index];
        return 1;
    }
    return GeoPGPreparePair(data_i, data_j, eta_i, eta_j,
                            r, t, mean_i, mean_j);
}

double GeoPGLogJoint(double corr, int r, int t,
                     double mean_i, double mean_j, double shape)
{
    const double p = biv_PoissonGamma(corr, r, t, mean_i, mean_j, shape);
    if (!R_FINITE(p) || p < 0.0 || p > 1.0 + PG_PROB_TOL) return R_NaN;
    return log(fmax2(p, PG_PROB_MIN));
}

double GeoPGLogMarginal(int z, double mean, double shape)
{
    return pg_nb_logpmf(z, mean, shape);
}

double GeoPGLogConditional(double corr, int r, int t,
                           double mean_i, double mean_j, double shape)
{
    const double lj = GeoPGLogJoint(corr, r, t, mean_i, mean_j, shape);
    const double lm = GeoPGLogMarginal(t, mean_j, shape);
    return (R_FINITE(lj) && R_FINITE(lm)) ? lj - lm : R_NaN;
}

double GeoPGZIPLogJoint(double corr, int r, int t,
                        double mean_i, double mean_j, double mup,
                        double nugget1, double nugget2, double shape)
{
    const double p = biv_PoissonGammaZIP(corr, r, t, mean_i, mean_j,
                                          mup, nugget1, nugget2, shape);
    if (!R_FINITE(p) || p < 0.0 || p > 1.0 + PG_PROB_TOL) return R_NaN;
    return log(fmax2(p, PG_PROB_MIN));
}

double GeoPGZIPLogMarginal(int z, double mean, double mup, double shape)
{
    if (z < 0 || !R_FINITE(mup)) return R_NaN;
    const double pzero = pnorm(mup, 0.0, 1.0, 1, 0);
    const double lnb = pg_nb_logpmf(z, mean, shape);
    if (!R_FINITE(pzero) || !R_FINITE(lnb)) return R_NaN;
    if (z == 0)
        return pg_logspace_add_pos(log(pzero), log1p(-pzero) + lnb);
    return log1p(-pzero) + lnb;
}

double GeoPGZIPLogConditional(double corr, int r, int t,
                              double mean_i, double mean_j, double mup,
                              double nugget1, double nugget2, double shape)
{
    const double lj = GeoPGZIPLogJoint(corr, r, t, mean_i, mean_j, mup,
                                       nugget1, nugget2, shape);
    const double lm = GeoPGZIPLogMarginal(t, mean_j, mup, shape);
    return (R_FINITE(lj) && R_FINITE(lm)) ? lj - lm : R_NaN;
}

double GeoPGGaussLogJoint(double corr, double data_i, double data_j,
                          double mean_i, double mean_j, double shape)
{
    if (!R_FINITE(data_i) || !R_FINITE(data_j) ||
        !(mean_i > 0.0) || !(mean_j > 0.0) || !(shape > 0.0))
        return R_NaN;
    const double vi = mean_i * (1.0 + mean_i / shape);
    const double vj = mean_j * (1.0 + mean_j / shape);
    double rr = corr_pois_gen(corr, mean_i, mean_j, shape);
    if (!R_FINITE(vi) || !R_FINITE(vj) || !(vi > 0.0) || !(vj > 0.0) ||
        !R_FINITE(rr)) return R_NaN;
    rr = fmax2(-1.0 + PG_RHO_EPS, fmin2(1.0 - PG_RHO_EPS, rr));
    const double zi = (data_i - mean_i) / sqrt(vi);
    const double zj = (data_j - mean_j) / sqrt(vj);
    const double omr2 = 1.0 - rr * rr;
    return -log(2.0 * M_PI) - 0.5 * (log(vi) + log(vj) + log(omr2)) -
           0.5 * (zi * zi - 2.0 * rr * zi * zj + zj * zj) / omr2;
}

double GeoPGGaussLogConditional(double corr, double data_i, double data_j,
                                double mean_i, double mean_j, double shape)
{
    const double lj = GeoPGGaussLogJoint(corr, data_i, data_j,
                                          mean_i, mean_j, shape);
    const double vj = mean_j * (1.0 + mean_j / shape);
    const double lm = dnorm(data_j, mean_j, sqrt(vj), 1);
    return (R_FINITE(lj) && R_FINITE(lm)) ? lj - lm : R_NaN;
}

/* Small diagnostics used by the external regression script. */
void GeoPEMCorrEval(double *corr, double *mean_i, double *mean_j,
                    double *shape, double *out)
{
    if (out) out[0] = corr_pois_gen(corr[0], mean_i[0], mean_j[0], shape[0]);
}

void GeoPEMPGrtBetaEval(double *corr, int *r, int *t,
                        double *mean_i, double *mean_j,
                        double *shape, double *out)
{
    if (out) out[0] = biv_PoissonGamma(corr[0], r[0], t[0],
                                       mean_i[0], mean_j[0], shape[0]);
}



void GeoPEMZIPEval(double *corr, int *r, int *t,
                   double *mean_i, double *mean_j, double *mup,
                   double *nugget1, double *nugget2,
                   double *shape, double *out)
{
    if (out) out[0] = biv_PoissonGammaZIP(corr[0], r[0], t[0],
        mean_i[0], mean_j[0], mup[0], nugget1[0], nugget2[0], shape[0]);
}

/*####*/
double int_kuma(double x,double eta, double gam,double k,double m)
{
    double res=0.0;

    res=R_pow( 1-R_pow(x,1/eta),1/gam)*R_pow(x,k-m)*R_pow(1-x,m);
    return (res);
}

/*####*/
void integr_kuma(double *x, int n, void *ex){
    int i;double eta,gam,k,m;
    eta =    ((double*)ex)[0];  
    gam =    ((double*)ex)[1];  
    k=       ((double*)ex)[2];  
    m=       ((double*)ex)[3];  

    for (i=0;i<n;i++) {
        x[i]=int_kuma(x[i],eta,gam,k,m);}
    return;
}

double kumaintegral(double *param) {
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;epsabs = R_pow(DBL_EPSILON, 0.25);epsrel = epsabs;    lenw = 4 * subdiv;           /* as instructed in WRE */
    iwork =   (int *) R_Calloc(subdiv, int);  /* idem */
    work = (double *) R_Calloc(lenw, double); /* idem */
    ex[0] = param[0]; //eta
    ex[1] = param[1]; //gam 
    ex[2] = param[2];  //k
    ex[3] = param[3];  //m
    lower=0;
    upper=1;
    Rdqags(integr_kuma, (void *) &ex,&lower, &upper, 
               &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);
    R_Free(iwork);R_Free(work);
    return(result);
}

/*######*/
  double mean_kuma(double eta,double gam){
    double out=eta*beta(1+(1/gam),eta);
    return(out);
  }
/*******/
   double var_kuma(double eta, double gam){
    double mm=mean_kuma(eta,gam);
    double out=eta*beta(1+2*(1/gam),eta)-mm*mm;
    return(out);
  }

/******/
/******/
double corr_kuma(double rho,double eta,double gam){

  double corr=0.0, tol=1e-6; int iter=0;
  double rho2=R_pow(rho,2);
  int k=0;int m=0;

  if(fabs(rho)< tol) return(0.0);

  /*###############*/
if (eta==1.0&&gam==1.0){
    corr=((2*(rho2*(3*rho2-1)-R_pow(rho2-1,2)*log1p(-rho2)))/R_pow(rho2,2))-3;
    return(corr);
  }
else
{
  if (eta==1.0&&gam!=1.0){
    double res_K=0.0,res_M=0.0,sum_M=0.0, bb,aa;
    iter=10000;
    tol=1e-14;
    while(k<=iter){
      res_M=0;m=0;
      bb=     2*(log1p(-rho2) + k*log(rho));
      while(m<=k){

        //A=exp(lbeta(1+k-m,1+(1/gam)+m));
        aa=    -2*lbeta(k-m+1,m+1);
        sum_M= exp(aa + bb + 2*lbeta(1+k-m,1+(1/gam)+m));
        res_M=res_M+  sum_M;
        if ((sum_M<tol)|(sum_M>1e300)){break;}
       m=m+1;
      }
      res_K=res_K+res_M;
      if (res_M<tol){ break;}
      k=k+1;
    }
    double mm=mean_kuma(eta,gam);
    double vv=var_kuma(eta,gam);
    corr=(res_K-R_pow(mm,2))/vv;
  }
/******/
if (eta!=1.0&&gam==1.0){
   double res_K=0.0,res_M=0.0,sum_M=0.0,bb,aa,c1,c2,c3;
    iter=10000;
    tol=1e-14;
    while(k<=iter){
      res_M=0;m=0;
        bb= 2*(log1p(-rho2) + k*log(rho));
      while(m<=k){
        aa=-2*lbeta(k-m+1,m+1);
        c1=lgammafn(1+m)+lgammafn(1+k-m)-lgammafn(2+k);
        c2=lgammafn(1+m)+lgammafn(1+(1/eta)+k-m)-lgammafn(2+(1/eta)+k);
        c3=aa + bb;
        sum_M= exp(2*c1+c3)+exp(2*c2+c3)-2*exp(c1+c2+c3);
        res_M=res_M+  sum_M;
        if ((sum_M<tol)|(sum_M>1e300)){ break;}
        m=m+1;
      }
      res_K=res_K+res_M;
      if (res_M<tol){break;}
      k=k+1;
    }
    double mm=mean_kuma(eta,gam);
    double vv=var_kuma(eta,gam);
    corr=(res_K-R_pow(mm,2))/vv;
  }
  
  /**+*********************/
  if ((eta!=1.0)&&(gam!=1.0)){
   
       double param[4] = {0.0, 0.0, 0.0, 0.0};
       param[0]=eta;param[1]=gam;  //mu,alpha //beta
       
 
    double res_K=0.0,res_M=0.0,sum_M=0.0,aa,bb,p2,p1;
    k=0;res_K=0;iter=10000;tol=1e-14;
     while (k<=iter){
      res_M=0.0;m=0;
      bb= 2*(log1p(-rho2) + k*log(rho));
       param[2]=k;
      while(m<=k){
        param[3]=m;
        p1=kumaintegral(param);
        aa=-2*lbeta(k-m+1,m+1);
        p2=exp(aa+bb);
        sum_M=p2*p1*p1;
        res_M=res_M+  sum_M;
        if ((sum_M<tol) | (sum_M>1e300)){break;}
        m=m+1;
      }
      res_K=res_K+res_M;
      if (res_M<tol){ break;}
      k=k+1;
    }
    double mm=mean_kuma(eta,gam);
    double vv=var_kuma(eta,gam);
    corr=(res_K-R_pow(mm,2))/vv;
  
  }
    return(corr);
 }
 
}
  /*********************/
void corr_kuma_vec(double *rho,double *eta,double *gam,double *res, int *n)
{
int i=0;
for(i=0;i<=*n;i++)  res[i]=corr_kuma(rho[i],eta[0],gam[0]);
}

 /*********************/

/*
 * Fast Gaussian-copula kernel.  The copula density has the closed form
 *
 *   log c(u,v;r) = -1/2 log(1-r^2)
 *                  - [r^2(x^2+y^2)-2rxy] / [2(1-r^2)],
 *
 * with x=Phi^{-1}(u), y=Phi^{-1}(v).  This avoids constructing three
 * Gaussian densities and then taking their ratio.  Pairwise likelihoods
 * revisit the same marginal uniforms many times, so the exact p -> qnorm(p)
 * map is also cached.  Cache keys are compared bit-for-bit; a hash collision
 * only causes a recomputation.
 */
#define GEO_GAUSSCOP_QCACHE_SIZE 65536u
typedef struct {
    uint64_t p_bits;
    double quantile;
} GeoGaussCopQCacheEntry;

static GeoGaussCopQCacheEntry geo_gausscop_qcache[GEO_GAUSSCOP_QCACHE_SIZE];

static R_INLINE double geo_gausscop_qnorm_cached(double p)
{
    uint64_t pb;
    memcpy(&pb, &p, sizeof(pb));
    uint64_t h = pb;
    h ^= h >> 30;
    h *= UINT64_C(0xbf58476d1ce4e5b9);
    h ^= h >> 27;
    h *= UINT64_C(0x94d049bb133111eb);
    h ^= h >> 31;

    GeoGaussCopQCacheEntry *e =
        &geo_gausscop_qcache[h & (GEO_GAUSSCOP_QCACHE_SIZE - 1u)];
    if (e->p_bits == pb) return e->quantile;

    const double q = qnorm(p, 0.0, 1.0, 1, 0);
    if (R_FINITE(q)) {
        e->p_bits = pb;
        e->quantile = q;
    }
    return q;
}

static R_INLINE double geo_gausscop_log_from_scores(double x, double y, double rho)
{
    rho = clamp_corr(rho);
    if (!R_FINITE(rho) || !R_FINITE(x) || !R_FINITE(y)) return R_NaN;
    const double r2 = rho * rho;
    const double den = 1.0 - r2;
    if (!(den > 0.0)) return R_NaN;
    return -0.5 * log1p(-r2)
           - 0.5 * (r2 * (x*x + y*y) - 2.0*rho*x*y) / den;
}

static R_INLINE double geo_gausscop_log(double dat1, double dat2, double rho)
{
    const double eps = 1e-12;
    dat1 = fmin(1.0-eps, fmax(eps, dat1));
    dat2 = fmin(1.0-eps, fmax(eps, dat2));
    const double a1 = geo_gausscop_qnorm_cached(dat1);
    const double a2 = geo_gausscop_qnorm_cached(dat2);
    return geo_gausscop_log_from_scores(a1, a2, rho);
}

void biv_unif_CopulaGauss_call(double *x,double *y,double *rho, double *res)
{
    *res = biv_unif_CopulaGauss(*x,*y,*rho);
}

// bivariate density Gaussian copula
double biv_unif_CopulaGauss(double dat1, double dat2, double rho)
{
    const double logc = geo_gausscop_log(dat1, dat2, rho);
    if (!R_FINITE(logc)) return R_NaN;
    return exp(logc);
}
/*
 * Pairwise skew-Gaussian copula fits revisit the same site-specific uniforms
 * many times within one objective evaluation.  The inverse skew-normal CDF
 * qsn() is substantially more expensive than the joint density itself, so
 * cache the exact (p, alpha) -> (quantile, density) map.  The cache is only
 * an execution optimization: keys are compared bit-for-bit and a collision
 * merely causes a recomputation.
 */
#define GEO_SKEWCOP_QCACHE_SIZE 65536u
typedef struct {
    uint64_t p_bits;
    uint64_t alpha_bits;
    double quantile;
    double density;
} GeoSkewCopQCacheEntry;

static GeoSkewCopQCacheEntry geo_skewcop_qcache[GEO_SKEWCOP_QCACHE_SIZE];

static inline uint64_t geo_double_bits(double x)
{
    uint64_t out;
    memcpy(&out, &x, sizeof(out));
    return out;
}

static inline uint64_t geo_skewcop_mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

static inline int geo_skewcop_quant_density(double p, double alpha,
                                            double *quantile, double *density)
{
    const double small = 1e-8;
    const double omega = sqrt(alpha * alpha + 1.0);

    /* eta = 0 is exactly the Gaussian copula. */
    if (alpha == 0.0) {
        const double q = qnorm(p, 0.0, 1.0, 1, 0);
        const double d = dnorm(q, 0.0, 1.0, 0);
        *quantile = q;
        *density = d;
        return R_FINITE(q) && R_FINITE(d) && d > 0.0;
    }

    const uint64_t pb = geo_double_bits(p);
    const uint64_t ab = geo_double_bits(alpha);
    const uint64_t h = geo_skewcop_mix64(pb ^ (geo_skewcop_mix64(ab) + UINT64_C(0x9e3779b97f4a7c15)));
    GeoSkewCopQCacheEntry *e = &geo_skewcop_qcache[h & (GEO_SKEWCOP_QCACHE_SIZE - 1u)];

    if (e->p_bits == pb && e->alpha_bits == ab) {
        *quantile = e->quantile;
        *density = e->density;
        return R_FINITE(*quantile) && R_FINITE(*density) && *density > 0.0;
    }

    const double q = qsn(p, omega, alpha, 0.0, small);
    const double d = dsn(q, omega, alpha, 0.0);
    if (R_FINITE(q) && R_FINITE(d) && d > 0.0) {
        e->p_bits = pb;
        e->alpha_bits = ab;
        e->quantile = q;
        e->density = d;
    }
    *quantile = q;
    *density = d;
    return R_FINITE(q) && R_FINITE(d) && d > 0.0;
}

// bivariate density skew-Gaussian copula
double biv_unif_CopulaSkewGauss(double dat1, double dat2, double rho, double alpha)
{
    const double eps = 1e-12;
    dat1 = fmin(1.0-eps, fmax(eps, dat1));
    dat2 = fmin(1.0-eps, fmax(eps, dat2));

    /* Preserve the exact Gaussian special case and avoid skew inversion work. */
    if (alpha == 0.0) return biv_unif_CopulaGauss(dat1, dat2, rho);

    double a1, a2, g1, g2;
    if (!geo_skewcop_quant_density(dat1, alpha, &a1, &g1) ||
        !geo_skewcop_quant_density(dat2, alpha, &a2, &g2))
        return R_NaN;

    const double joint = biv_skew(rho, a1, a2, 0, 0, 1, alpha, 0);
    if (!R_FINITE(joint) || joint <= 0.0) return R_NaN;
    return joint / (g1 * g2);
}


/*
 * Rectangle probability for the skew-Gaussian copula.
 *
 * Discrete copula likelihoods need
 *   P(u0 < U1 <= u1, v0 < U2 <= v1)
 * rather than the continuous copula density.  For the GeoModels
 * skew-Gaussian RF construction the copula density is available in closed
 * form, while its bivariate CDF is not.  We therefore integrate the copula
 * density directly on the bounded probability rectangle using tensor-product
 * Gauss-Legendre quadrature.  The marginal skew-normal quantiles are computed
 * only once per quadrature node (not once per tensor cell), which is critical
 * inside the pairwise likelihood optimizer.
 *
 * eta=0 is exactly the Gaussian copula and is evaluated by the existing
 * bivariate-normal rectangle formula.  rho=0 is exactly independence.
 */
static const double geo_gl10_x[10] = {
    -0.97390652851717172008, -0.86506336668898451073,
    -0.67940956829902440623, -0.43339539412924719080,
    -0.14887433898163121088,  0.14887433898163121088,
     0.43339539412924719080,  0.67940956829902440623,
     0.86506336668898451073,  0.97390652851717172008
};
static const double geo_gl10_w[10] = {
    0.06667134430868813759, 0.14945134915058059315,
    0.21908636251598204399, 0.26926671930999635509,
    0.29552422471475287017, 0.29552422471475287017,
    0.26926671930999635509, 0.21908636251598204399,
    0.14945134915058059315, 0.06667134430868813759
};
static const double geo_gl16_x[16] = {
    -0.98940093499164993260, -0.94457502307323257608,
    -0.86563120238783174388, -0.75540440835500303390,
    -0.61787624440264374845, -0.45801677765722738634,
    -0.28160355077925891323, -0.09501250983763744019,
     0.09501250983763744019,  0.28160355077925891323,
     0.45801677765722738634,  0.61787624440264374845,
     0.75540440835500303390,  0.86563120238783174388,
     0.94457502307323257608,  0.98940093499164993260
};
static const double geo_gl16_w[16] = {
    0.02715245941175409485, 0.06225352393864789286,
    0.09515851168249278481, 0.12462897125553387205,
    0.14959598881657673208, 0.16915651939500253819,
    0.18260341504492358887, 0.18945061045506849629,
    0.18945061045506849629, 0.18260341504492358887,
    0.16915651939500253819, 0.14959598881657673208,
    0.12462897125553387205, 0.09515851168249278481,
    0.06225352393864789286, 0.02715245941175409485
};

static double geo_skewcop_rect_rule(double ulo, double uhi,
                                    double vlo, double vhi,
                                    double rho, double alpha,
                                    const double *gx, const double *gw,
                                    int nq)
{
    const double hu = 0.5 * (uhi - ulo), cu = 0.5 * (uhi + ulo);
    const double hv = 0.5 * (vhi - vlo), cv = 0.5 * (vhi + vlo);
    double qu[16], du[16], qv[16], dv[16];

    if (!(hu > 0.0) || !(hv > 0.0)) return 0.0;
    for (int i = 0; i < nq; ++i) {
        const double u = cu + hu * gx[i];
        const double v = cv + hv * gx[i];
        if (!geo_skewcop_quant_density(u, alpha, &qu[i], &du[i]) ||
            !geo_skewcop_quant_density(v, alpha, &qv[i], &dv[i]))
            return R_NaN;
    }

    double acc = 0.0;
    for (int i = 0; i < nq; ++i) {
        for (int j = 0; j < nq; ++j) {
            const double joint = biv_skew(rho, qu[i], qv[j], 0, 0, 1, alpha, 0);
            if (!R_FINITE(joint) || joint < 0.0) return R_NaN;
            const double copdens = joint / (du[i] * dv[j]);
            if (!R_FINITE(copdens) || copdens < 0.0) return R_NaN;
            acc += gw[i] * gw[j] * copdens;
        }
    }
    return hu * hv * acc;
}

static double geo_skewcop_rect_prob(double ulo, double uhi,
                                    double vlo, double vhi,
                                    double rho, double alpha)
{
    if (!R_FINITE(ulo) || !R_FINITE(uhi) ||
        !R_FINITE(vlo) || !R_FINITE(vhi) ||
        !R_FINITE(rho) || !R_FINITE(alpha)) return R_NaN;

    ulo = fmin(1.0, fmax(0.0, ulo));
    uhi = fmin(1.0, fmax(0.0, uhi));
    vlo = fmin(1.0, fmax(0.0, vlo));
    vhi = fmin(1.0, fmax(0.0, vhi));
    if (uhi <= ulo || vhi <= vlo) return 0.0;

    rho = clamp_corr(rho);
    if (fabs(rho) < 1e-14)
        return (uhi - ulo) * (vhi - vlo);

    if (fabs(alpha) < 1e-14) {
        const double xh = qnorm(uhi, 0.0, 1.0, 1, 0);
        const double xl = qnorm(ulo, 0.0, 1.0, 1, 0);
        const double yh = qnorm(vhi, 0.0, 1.0, 1, 0);
        const double yl = qnorm(vlo, 0.0, 1.0, 1, 0);
        double p = safe_pbnorm22(xh, yh, rho) - safe_pbnorm22(xl, yh, rho)
                 - safe_pbnorm22(xh, yl, rho) + safe_pbnorm22(xl, yl, rho);
        return fmax(0.0, p);
    }

    double p = geo_skewcop_rect_rule(ulo, uhi, vlo, vhi, rho, alpha,
                                     geo_gl10_x, geo_gl10_w, 10);
    const double maxp = fmin(uhi - ulo, vhi - vlo);
    const double eta_abs = fabs(alpha) / sqrt(1.0 + alpha * alpha);
    const double edge = fmin(fmin(ulo, 1.0-uhi), fmin(vlo, 1.0-vhi));

    /* Strong dependence/skewness and extreme-tail rectangles are the cases
       where the density varies most within a cell.  Refine only there, or if
       the low-order result violates elementary probability bounds. */
    if (!R_FINITE(p) || p < 0.0 || p > maxp * (1.0 + 1e-7) ||
        fabs(rho) > 0.85 || eta_abs > 0.85 || edge < 1e-6) {
        p = geo_skewcop_rect_rule(ulo, uhi, vlo, vhi, rho, alpha,
                                  geo_gl16_x, geo_gl16_w, 16);
    }
    if (!R_FINITE(p)) return R_NaN;
    if (p < 0.0 && p > -1e-12) p = 0.0;
    if (p > maxp && p < maxp + 1e-10) p = maxp;
    return p;
}

void biv_unif_CopulaSkewGauss_rect_call(double *ulo, double *uhi,
                                        double *vlo, double *vhi,
                                        double *rho, double *nu, double *res)
{
    if (!R_FINITE(*nu) || fabs(*nu) >= 1.0) {
        *res = R_NaN;
        return;
    }
    const double alpha = (*nu) / sqrt(1.0 - (*nu) * (*nu));
    *res = geo_skewcop_rect_prob(*ulo, *uhi, *vlo, *vhi,
                                 clamp_corr(*rho), alpha);
}

 /*********************/
void biv_unif_CopulaClayton_call(double *x,double *y,double *rho, double *nu, double *res)
{
    *res = biv_unif_CopulaClayton(*x,*y,*rho,*nu);
}

/*
 * Fast log-Appell F4 evaluator for the constructive Clayton-like copula.
 *
 * In this caller the parameters always satisfy
 *   a = b = c + 1, d = 1,
 * with c = nu/2 > 0 and x,y >= 0.  Euler transformation of the inner
 * Gauss function gives
 *
 * F4(c+1,c+1;c,1;x,y)
 *   = sum_{k>=0} C_k P_{k+1}(x),
 * P_n(x) = 2F1(-n,-n;c;x),
 * C_0 = (1-x)^(-c-2),
 * C_{k+1}/C_k = y ((c+1+k)/(k+1))^2 / (1-x)^2.
 *
 * The polynomial sequence P_n is advanced with the three-term Jacobi
 * recurrence after the exact identity
 *   P_n(x)=(1-x)^n n!/(c)_n P_n^{(c-1,0)}((1+x)/(1-x)).
 * Thus the Clayton pair density no longer calls the generic hypergeo()
 * routine once per outer F4 term.  All F4 terms are positive in the
 * Clayton domain, so accumulation in log scale is cancellation-free.
 */
static R_INLINE double geo_clay_logaddexp(double a, double b)
{
    if (a == R_NegInf) return b;
    if (b == R_NegInf) return a;
    if (a < b) { const double t = a; a = b; b = t; }
    return a + log1p(exp(b - a));
}

static double geo_clayton_logf4_fast(double c, double x, double y, int *iters)
{
    const int MAX_ITER = 2000;
    const int MIN_K = 10;
    const int QUIET_NEED = 4;
    const double F4_LOG_TOL = -34.538776394910684; /* log(1e-15) */
    const double LOG_SCALE = 345.38776394910684; /* log(1e150) */

    if (iters) *iters = 0;
    if (!R_FINITE(c) || !R_FINITE(x) || !R_FINITE(y) ||
        c <= 0.0 || x < 0.0 || y < 0.0 || x >= 1.0 || y >= 1.0)
        return R_NaN;

    const double log1mx = log1p(-x);
    const double logy = (y > 0.0) ? log(y) : R_NegInf;

    /* P_0 and P_1, stored with a common dynamic scale. */
    double pnm1 = 1.0;
    double pn = 1.0 + x / c;
    double p_log_scale = 0.0;

    double logC = -(c + 2.0) * log1mx;
    double logsum = R_NegInf;
    int quiet = 0;

    for (int k = 0; k < MAX_ITER; ++k) {
        if (!(R_FINITE(pn) && pn > 0.0)) return R_NaN;

        const double logP = log(pn) + p_log_scale;
        const double logterm = logC + logP;
        const double newsum = geo_clay_logaddexp(logsum, logterm);

        if (k > MIN_K && logterm - newsum < F4_LOG_TOL) {
            if (++quiet >= QUIET_NEED) {
                if (iters) *iters = k + 1;
                return newsum;
            }
        } else {
            quiet = 0;
        }
        logsum = newsum;

        if (y == 0.0) {
            if (iters) *iters = k + 1;
            return logsum;
        }

        logC += logy +
                2.0 * (log(c + 1.0 + (double)k) - log((double)k + 1.0)) -
                2.0 * log1mx;

        /* Advance P_n -> P_{n+1}; here n=k+1 >= 1. */
        const double n = (double)k + 1.0;
        const double cn = c + n;
        const double den = cn * cn * (c + 2.0*n - 1.0);
        const double numer = c * (c - 1.0 + 2.0*n) + 2.0*n*n +
                             x * (2.0*n*cn + c - 1.0);
        const double A = (c + 2.0*n) * numer / den;
        const double B = n*n * (1.0-x)*(1.0-x) * (c + 2.0*n + 1.0) / den;
        double pnext = A * pn - B * pnm1;
        if (!(R_FINITE(pnext) && pnext > 0.0)) return R_NaN;

        const double mx = fmax(fabs(pn), fabs(pnext));
        if (mx > 1e150) {
            pn *= 1e-150;
            pnext *= 1e-150;
            p_log_scale += LOG_SCALE;
        } else if (mx > 0.0 && mx < 1e-150) {
            pn *= 1e150;
            pnext *= 1e150;
            p_log_scale -= LOG_SCALE;
        }
        pnm1 = pn;
        pn = pnext;
    }

    if (iters) *iters = MAX_ITER;
    return logsum;
}

// bivariate density clayton copula
double biv_unif_CopulaClayton(double dat1, double dat2, double rho, double nu)
{
    const double eps = 1e-12;
    dat1 = fmin(1.0-eps, fmax(eps, dat1));
    dat2 = fmin(1.0-eps, fmax(eps, dat2));
    if(!R_FINITE(nu) || nu <= 0.0 ||
       fabs(nu-nearbyint(nu)) > sqrt(DBL_EPSILON)) return R_NegInf;

    double nu2 = nu / 2;
    double rho2 = rho * rho;
    double a = nu2 + 1;

    // Calcola i termini a1 e a2
    double a1 = R_pow(dat1, 1 / nu2);
    double a2 = R_pow(dat2, 1 / nu2);

    // GeoModels' Clayton-like copula density from the RF construction.
    const double xx = rho2 * a1 * a2;
    const double yy = rho2 * (1.0 - a1) * (1.0 - a2);
    if (rho2 >= 1.0) return R_NegInf;

    double logf4 = geo_clayton_logf4_fast(nu2, xx, yy, NULL);
    if (!R_FINITE(logf4)) {
        /* Conservative fallback: preserve the historical generic evaluator
           outside the numerically validated fast recurrence path. */
        const double f4 = appellF4(a, a, nu2, 1, xx, yy);
        if (!R_FINITE(f4) || f4 <= 0.0) return R_NegInf;
        logf4 = log(f4);
    }
    return a * log1p(-rho2) + logf4;
}


double cdf_kuma(double y,double a, double b){
double res=1-R_pow(1-R_pow(y,a),b);
return(res);
}
double pdf_kuma(double y,double a, double b){
double res=(a*b)* R_pow(y,a-1)*R_pow(1-R_pow(y,a),b-1);
return(res);
}

/*
 * Exact site-level marginal cache for continuous copula likelihoods.
 *
 * In a pairwise composite likelihood the same site appears in several pairs.
 * Its marginal PIT and density therefore get recomputed repeatedly even within
 * a single objective evaluation.  Cache the exact endpoint map
 *
 *   (model, z, mu, marginal nuisance parameters) -> (CDF, PDF, score)
 *
 * for the continuous copula margins.  Nugget and copula parameters are
 * deliberately excluded from the key because they do not enter the marginal
 * distribution.  All key components are compared bit-for-bit; a hash
 * collision only causes a miss and a recomputation, never reuse of a wrong
 * value.  This is therefore an execution-only optimization.
 */
#define GEO_MARGIN_CACHE_SIZE 65536u
#define GEO_MARGIN_MAX_KEYPAR 4

typedef struct {
    int model;
    unsigned char npar;
    uint64_t key[2 + GEO_MARGIN_MAX_KEYPAR];
    double cdf;
    double pdf;
    double score;
} GeoMarginCacheEntry;

static GeoMarginCacheEntry geo_margin_cache[GEO_MARGIN_CACHE_SIZE];

static R_INLINE int geo_margin_keypars(int model, const double *nuis,
                                       uint64_t out[GEO_MARGIN_MAX_KEYPAR])
{
    int n = 0;
    switch(model) {
    case 1:   /* Gaussian: sill */
    case 22:  /* LogGaussian: sill */
    case 25:  /* Logistic: sill */
        out[n++] = geo_double_bits(nuis[1]);
        break;
    case 12:  /* StudentT: df, sill */
        out[n++] = geo_double_bits(nuis[0]);
        out[n++] = geo_double_bits(nuis[2]);
        break;
    case 21:  /* Gamma: shape */
    case 26:  /* Weibull: shape */
        out[n++] = geo_double_bits(nuis[2]);
        break;
    case 59:  /* SkewLaplace: sill, skew */
        out[n++] = geo_double_bits(nuis[1]);
        out[n++] = geo_double_bits(nuis[2]);
        break;
    case 28:  /* Beta: shape1, shape2, min, max */
    case 33:  /* Kumaraswamy: shape1, shape2, min, max */
        out[n++] = geo_double_bits(nuis[2]);
        out[n++] = geo_double_bits(nuis[3]);
        out[n++] = geo_double_bits(nuis[4]);
        out[n++] = geo_double_bits(nuis[5]);
        break;
    case 50:  /* Beta2 regression: shape, min, max */
    case 42:  /* Kumaraswamy2 regression: shape, min, max */
        out[n++] = geo_double_bits(nuis[2]);
        out[n++] = geo_double_bits(nuis[3]);
        out[n++] = geo_double_bits(nuis[4]);
        break;
    default:
        return 0;
    }
    return n;
}

static R_INLINE uint64_t geo_margin_hash(int model,
                                         const uint64_t *key, int nkey)
{
    uint64_t h = geo_skewcop_mix64((uint64_t)(unsigned int)model +
                                   UINT64_C(0x9e3779b97f4a7c15));
    for (int i = 0; i < nkey; ++i)
        h = geo_skewcop_mix64(h ^ (geo_skewcop_mix64(key[i]) +
                                   UINT64_C(0x9e3779b97f4a7c15)));
    return h;
}

static R_INLINE int geo_margin_make_key(int model, double z, double mu,
                                        const double *nuis,
                                        uint64_t key[2 + GEO_MARGIN_MAX_KEYPAR],
                                        int *nkey)
{
    key[0] = geo_double_bits(z);
    key[1] = geo_double_bits(mu);
    const int np = geo_margin_keypars(model, nuis, key + 2);
    if (np <= 0) return 0;
    *nkey = 2 + np;
    return np;
}

static R_INLINE int geo_margin_cache_lookup(int model, double z, double mu,
                                            const double *nuis,
                                            double *cdf, double *pdf,
                                            double *score)
{
    uint64_t key[2 + GEO_MARGIN_MAX_KEYPAR];
    int nkey = 0;
    const int np = geo_margin_make_key(model, z, mu, nuis, key, &nkey);
    if (np <= 0) return 0;

    const uint64_t h = geo_margin_hash(model, key, nkey);
    GeoMarginCacheEntry *e = &geo_margin_cache[h & (GEO_MARGIN_CACHE_SIZE - 1u)];
    if (e->model != model || e->npar != (unsigned char)np) return 0;
    for (int i = 0; i < nkey; ++i)
        if (e->key[i] != key[i]) return 0;

    *cdf = e->cdf;
    *pdf = e->pdf;
    *score = e->score;
    return 1;
}

static R_INLINE void geo_margin_cache_store(int model, double z, double mu,
                                            const double *nuis,
                                            double cdf, double pdf,
                                            double score)
{
    uint64_t key[2 + GEO_MARGIN_MAX_KEYPAR];
    int nkey = 0;
    const int np = geo_margin_make_key(model, z, mu, nuis, key, &nkey);
    if (np <= 0) return;
    if (!R_FINITE(cdf) || !R_FINITE(pdf) || !R_FINITE(score)) return;

    const uint64_t h = geo_margin_hash(model, key, nkey);
    GeoMarginCacheEntry *e = &geo_margin_cache[h & (GEO_MARGIN_CACHE_SIZE - 1u)];
    e->model = model;
    e->npar = (unsigned char)np;
    for (int i = 0; i < nkey; ++i) e->key[i] = key[i];
    e->cdf = cdf;
    e->pdf = pdf;
    e->score = score;
}

/*********************************************************************/
/*********************************************************************/
double biv_cop(double rho,int type_cop,int cond,
             double z1,double z2,double mu1,double mu2,double *nuis,int model, int NN1,int NN2)
             {
double dens=0.0,rho1=0.0,nu=0.0;
double g1=0.0,g2=0.0,a1=0.0,a2=0.0,b1=0.0,b2=0.0;
const double margin_mu1 = mu1;
const double margin_mu2 = mu2;

/*
 * The nuisance vector follows NuisParam2() after removal of the mean
 * coefficients.  StudentT is the only copula margin among the supported
 * models for which the nugget is not the first nuisance parameter:
 *   StudentT: df, nugget, sill, [nu]
 * All other supported copula margins start with nugget, sill, ... .
 */
double nugget_cop = (model == 12) ? nuis[1] : nuis[0];
if (!R_FINITE(rho) || !R_FINITE(nugget_cop) || nugget_cop < 0.0 || nugget_cop >= 1.0)
    return LOW;
rho1 = clamp_corr((1.0 - nugget_cop) * rho);

/* Gaussian/Gaussian has its own direct joint-density shortcut and must still
 * compute the standardized scores below.  Otherwise, if both endpoints are
 * already cached, the whole marginal switch can be skipped. */
const int gaussian_direct = (type_cop == 1 && model == 1);
int margin_hit1 = 0, margin_hit2 = 0;
if (!gaussian_direct) {
    margin_hit1 = geo_margin_cache_lookup(model, z1, margin_mu1, nuis,
                                          &a1, &g1, &b1);
    margin_hit2 = geo_margin_cache_lookup(model, z2, margin_mu2, nuis,
                                          &a2, &g2, &b2);
}

if (!(margin_hit1 && margin_hit2))
switch(model) // Correlation functions are in alphabetical order
    {

/******* models on the real line ****************/
    case 1: // gaussian
      if(!R_FINITE(nuis[1]) || nuis[1] <= 0.0) return LOW;
      b1=(z1-mu1)/sqrt(nuis[1]);b2=(z2-mu2)/sqrt(nuis[1]);
      /* Gaussian margin + Gaussian copula is exactly bivariate Gaussian.
       * Delay PIT/pdf work in that case and evaluate the joint log-density
       * directly below. */
      if(type_cop != 1) {
          a1=pnorm(b1,0,1,1,0);  a2=pnorm(b2,0,1,1,0);
          g1=dnorm(b1,0,1,0)/sqrt(nuis[1]);
          g2=dnorm(b2,0,1,0)/sqrt(nuis[1]);
      }
    break;
 case 34:  // Tukeyh
  
      b1=1;
      b2=1;
      a1=1;
      a2=1;
      g1=1;
      g2=1;
   break;
    case 25:  // Logistic
      if(!R_FINITE(nuis[1]) || nuis[1] <= 0.0) return LOW;
      b1=1;b2=1;
      a1=plogis(z1, mu1, sqrt(nuis[1]),1,0);
      a2=plogis(z2, mu2, sqrt(nuis[1]),1,0);
      g1=dlogis(z1, mu1, sqrt(nuis[1]),0);
      g2=dlogis(z2, mu2, sqrt(nuis[1]),0);
   break;
   case 12: // t student
      if(!R_FINITE(nuis[0]) || nuis[0] <= 0.0 ||
         !R_FINITE(nuis[2]) || nuis[2] <= 0.0) return LOW;
      b1=(z1-mu1)/sqrt(nuis[2]);b2=(z2-mu2)/sqrt(nuis[2]);
      a1=pt(b1,1/nuis[0],1,0);  a2=pt(b2,1/nuis[0],1,0);
      g1=dt(b1,1/nuis[0],0)/sqrt(nuis[2]); //marginal 1
      g2=dt(b2,1/nuis[0],0)/sqrt(nuis[2]); //marginal 2
   break;
/******* models on the  positive real line ****************/
     case 22: // lognormal
      if(!R_FINITE(nuis[1]) || nuis[1] <= 0.0 || z1 <= 0.0 || z2 <= 0.0) return LOW;
      b1=mu1-nuis[1]/2; b2=mu2-nuis[1]/2;
      a1=plnorm(z1,b1, sqrt(nuis[1]),1,0);
      a2=plnorm(z2,b2, sqrt(nuis[1]),1,0);
      g1=dlnorm(z1,b1, sqrt(nuis[1]),0);
      g2=dlnorm(z2,b2, sqrt(nuis[1]),0);
  break;
   case 21: // gamma
      if(!R_FINITE(nuis[2]) || nuis[2] <= 0.0 || z1 <= 0.0 || z2 <= 0.0) return LOW;
      b1=1;
      b2=1;
      a1=pgamma(z1,nuis[2]/2,R_pow(nuis[2]/(2*exp(mu1)),-1),1,0);  
      a2=pgamma(z2,nuis[2]/2,R_pow(nuis[2]/(2*exp(mu2)),-1),1,0);  
      g1=dgamma(z1,nuis[2]/2,R_pow(nuis[2]/(2*exp(mu1)),-1),0);  
      g2=dgamma(z2,nuis[2]/2,R_pow(nuis[2]/(2*exp(mu2)),-1),0); 
   break;
     case 26: // weibull
      if(!R_FINITE(nuis[2]) || nuis[2] <= 0.0 || z1 < 0.0 || z2 < 0.0) return LOW;
      b1=1;
      b2=1;
      a1=pweibull(z1,nuis[2],exp(mu1)/(gammafn(1+1/nuis[2])),1,0);
      a2=pweibull(z2,nuis[2],exp(mu2)/(gammafn(1+1/nuis[2])),1,0);
      g1=dweibull(z1,nuis[2],exp(mu1)/(gammafn(1+1/nuis[2])),0);
      g2=dweibull(z2,nuis[2],exp(mu2)/(gammafn(1+1/nuis[2])),0);
   break;
     case 59: // skewlaplace
      if(!R_FINITE(nuis[1]) || nuis[1] <= 0.0 ||
         !R_FINITE(nuis[2]) || nuis[2] <= 0.0 || nuis[2] >= 1.0) return LOW;
      b1=1;
      b2=1;
      a1=pSkewLaplace(z1,nuis[2],mu1,nuis[1]);
      a2=pSkewLaplace(z2,nuis[2],mu2,nuis[1]);
      g1=dSkewLaplace(z1,nuis[2],mu1,nuis[1]);
      g2=dSkewLaplace(z2,nuis[2],mu2,nuis[1]);
      //Rprintf("%f %f %f %f ",a1,a2,g1,g2);
   break;
/******* models on the a bounded support ****************/
   case 28: // Beta
      if(!R_FINITE(nuis[2]) || !R_FINITE(nuis[3]) || nuis[2] <= 0.0 || nuis[3] <= 0.0 ||
         !R_FINITE(nuis[4]) || !R_FINITE(nuis[5]) || nuis[4] >= nuis[5] ||
         z1 < nuis[4] || z1 > nuis[5] || z2 < nuis[4] || z2 > nuis[5]) return LOW;
      b1=(z1- nuis[4])/(nuis[5]-nuis[4]); b2=(z2- nuis[4])/(nuis[5]-nuis[4]);
      a1=pbeta(b1,nuis[2],nuis[3],1,0); a2=pbeta(b2,nuis[2],nuis[3],1,0);
      g1=dbeta(b1,nuis[2],nuis[3],0)/(nuis[5]-nuis[4]);//marginal 1
      g2=dbeta(b2,nuis[2],nuis[3],0)/(nuis[5]-nuis[4]);//marginal 2
   break;   
   case 50:  // Beta regression
      if(!R_FINITE(nuis[2]) || nuis[2] <= 0.0 ||
         !R_FINITE(nuis[3]) || !R_FINITE(nuis[4]) || nuis[3] >= nuis[4] ||
         z1 < nuis[3] || z1 > nuis[4] || z2 < nuis[3] || z2 > nuis[4]) return LOW;
      mu1=1/(1+exp(-mu1));
      mu2=1/(1+exp(-mu2));
    //  Rprintf("%f %f %f %f %f \n",mu1,rho,nuis[2],nuis[4],nuis[3]);
      b1=(z1- nuis[3])/(nuis[4]-nuis[3]); b2=(z2- nuis[3])/(nuis[4]-nuis[3]);
      a1=pbeta(b1, nuis[2]*mu1,(1-mu1)*nuis[2],1,0);
      a2=pbeta(b2, nuis[2]*mu2,(1-mu2)*nuis[2],1,0);
      g1=dbeta(b1, nuis[2]*mu1,(1-mu1)*nuis[2],0)/(nuis[4]-nuis[3]); //marginal 1
      g2=dbeta(b2, nuis[2]*mu2,(1-mu2)*nuis[2],0)/(nuis[4]-nuis[3]); //marginal 2
   break;
   case 33:  // Kumaraswamy
      if(!R_FINITE(nuis[2]) || !R_FINITE(nuis[3]) || nuis[2] <= 0.0 || nuis[3] <= 0.0 ||
         !R_FINITE(nuis[4]) || !R_FINITE(nuis[5]) || nuis[4] >= nuis[5] ||
         z1 < nuis[4] || z1 > nuis[5] || z2 < nuis[4] || z2 > nuis[5]) return LOW;
      b1=(z1- nuis[4])/(nuis[5]-nuis[4]);
      b2=(z2- nuis[4])/(nuis[5]-nuis[4]);
      a1= cdf_kuma(b1,nuis[2],nuis[3]);
      a2= cdf_kuma(b2,nuis[2],nuis[3]);
      g1= pdf_kuma(b1,nuis[2],nuis[3])/(nuis[5]-nuis[4]);
      g2= pdf_kuma(b2,nuis[2],nuis[3])/(nuis[5]-nuis[4]);
   break;
   case 42:  // Kumaraswamy regression
      if(!R_FINITE(nuis[2]) || nuis[2] <= 0.0 ||
         !R_FINITE(nuis[3]) || !R_FINITE(nuis[4]) || nuis[3] >= nuis[4] ||
         z1 < nuis[3] || z1 > nuis[4] || z2 < nuis[3] || z2 > nuis[4]) return LOW;
      mu1=1/(1+exp(-mu1));mu2=1/(1+exp(-mu2));
      b1=(z1- nuis[3])/(nuis[4]-nuis[3]);
      b2=(z2- nuis[3])/(nuis[4]-nuis[3]);
      mu1=log(1-R_pow(  mu1    ,nuis[2]))/log(0.5) ;
      mu2=log(1-R_pow(  mu2    ,nuis[2]))/log(0.5) ;
      a1= cdf_kuma(b1, nuis[2],1/mu1 );
      a2= cdf_kuma(b2, nuis[2],1/mu2 );
      g1= pdf_kuma(b1, nuis[2],1/mu1 )/(nuis[4]-nuis[3]);
      g2= pdf_kuma(b2, nuis[2],1/mu2 )/(nuis[4]-nuis[3]);
   break;
/************ discrete models *****************************/
case 30: // Poisson
  
      mu1=exp(mu1);mu2=exp(mu2);
      b2=dpois(z2,mu2,0);
      a1=qnorm(ppois(z1,  mu1,1,0),0,1,1,0);
      a2=qnorm(ppois(z1-1,mu1,1,0),0,1,1,0);
      g1=qnorm(ppois(z2,  mu2,1,0),0,1,1,0);      
      g2=qnorm(ppois(z2-1,mu2,1,0),0,1,1,0);  
      if(z1==0) {a2=-99;}
      if(z2==0) {g2=-99;}
   break;
      case 43: // Poisson inflated
 
      mu1=exp(mu1);mu2=exp(mu2);
      b2=dpois(z2,mu2,0);
      a1=qnorm(ppoisinflated(z1,  mu1,nuis[1]),0,1,1,0);
      a2=qnorm(ppoisinflated(z1-1,mu1,nuis[1]),0,1,1,0);
      g1=qnorm(ppoisinflated(z2,  mu2,nuis[1]),0,1,1,0);      
      g2=qnorm(ppoisinflated(z2-1,mu2,nuis[1]),0,1,1,0);  
      if(z1==0) {a2=-99;}
      if(z2==0) {g2=-99;}
   break;
     case 11: // Binomiale
   
    mu1 = pnorm(mu1, 0, 1, 1, 0);
    mu2 = pnorm(mu2, 0, 1, 1, 0);
    b2 = dbinom(z2, NN2, mu2, 0);
    a1 = (z1 >= NN1) ? R_PosInf : qnorm(pbinom(z1,   NN1, mu1, 1, 0), 0, 1, 1, 0);
    a2 = (z1 == 0)   ? R_NegInf : qnorm(pbinom(z1-1, NN1, mu1, 1, 0), 0, 1, 1, 0);
    g1 = (z2 >= NN2) ? R_PosInf : qnorm(pbinom(z2,   NN2, mu2, 1, 0), 0, 1, 1, 0);
    g2 = (z2 == 0)   ? R_NegInf : qnorm(pbinom(z2-1, NN2, mu2, 1, 0), 0, 1, 1, 0);  
   break;
       case 16: // binomial neg
  
    mu1 = pnorm(mu1, 0, 1, 1, 0);
    mu2 = pnorm(mu2, 0, 1, 1, 0);
    b2 = dnbinom(z2, NN2, mu2, 0);
    a1 = qnorm(pnbinom(z1,   NN1, mu1, 1, 0), 0, 1, 1, 0);
    a2 = (z1 == 0) ? R_NegInf : qnorm(pnbinom(z1-1, NN1, mu1, 1, 0), 0, 1, 1, 0);
    g1 = qnorm(pnbinom(z2,   NN2, mu2, 1, 0), 0, 1, 1, 0);
    g2 = (z2 == 0) ? R_NegInf : qnorm(pnbinom(z2-1, NN2, mu2, 1, 0), 0, 1, 1, 0);  
   break;
    case 45: // binomial neg  inflated   /// ok 
     
      mu1=pnorm(mu1,0,1,1,0);
      mu2=pnorm(mu2,0,1,1,0);
      b2=dnbinom(z2,NN2,mu2,0);
      a1=qnorm(pbneginflated(z1,  NN1,mu1,nuis[1]),0,1,1,0);
      a2=qnorm(pbneginflated(z1-1,NN1,mu1,nuis[1]),0,1,1,0);
      g1=qnorm(pbneginflated(z2,  NN2,mu2,nuis[1]),0,1,1,0);      
      g2=qnorm(pbneginflated(z2-1,NN2,mu2,nuis[1]),0,1,1,0);
      if(z1==0) {a2=-99;}
      if(z2==0) {g2=-99;}
   break;
       }

if (!(margin_hit1 && margin_hit2) && !gaussian_direct) {
    geo_margin_cache_store(model, z1, margin_mu1, nuis, a1, g1, b1);
    geo_margin_cache_store(model, z2, margin_mu2, nuis, a2, g2, b2);
}

/*********************** end cases ***************/
/******************copula gaussiana*************************/
if(type_cop==1)  {
   if(!(model==16 || model==11 || model==30)) {   // continuous
       if(model == 1) {
           /* Exact shortcut: Gaussian margins under a Gaussian copula are
            * a bivariate Gaussian law on the observed scale. */
           const double r2 = rho1 * rho1;
           const double den = 1.0 - r2;
           if (!(den > 0.0) || !R_FINITE(b1) || !R_FINITE(b2)) return LOW;
           dens = -log(2.0*M_PI) - log(nuis[1]) - 0.5*log1p(-r2)
                  - 0.5*(b1*b1 - 2.0*rho1*b1*b2 + b2*b2)/den;
       } else {
           if(!R_FINITE(g1) || !R_FINITE(g2) || g1 <= 0.0 || g2 <= 0.0)
               return LOW;
           const double logcop = geo_gausscop_log(a1, a2, rho1);
           if(!R_FINITE(logcop)) return LOW;
           dens = logcop + log(g1) + log(g2);
       }
   } else {                                      // discreti
       double prob = safe_pbnorm22(a1, g1, rho1) - safe_pbnorm22(a2, g1, rho1) - safe_pbnorm22(a1, g2, rho1) + safe_pbnorm22(a2, g2, rho1);
       prob = fmax(DBL_MIN, prob);
       dens = log(prob);
   }
}

/******************copula clayton and skew *************************/             
if(type_cop==2)
{
    if(!(model==16||model==11||model==30)) { // continuous models
        /* Copula parameter positions after removal of mean coefficients. */
        if(model==1 || model==22 || model==25) nu=nuis[2];
        if(model==12 || model==21 || model==26 || model==59) nu=nuis[3];
        if(model==50 || model==42) nu=nuis[5];
        if(model==28 || model==33) nu=nuis[6];

        /* GeoModels' Clayton-like RF construction uses positive integer nu. */
        if(!R_FINITE(nu) || nu <= 0.0 || fabs(nu-nearbyint(nu)) > sqrt(DBL_EPSILON))
            return LOW;
        if(!R_FINITE(g1) || !R_FINITE(g2) || g1 <= 0.0 || g2 <= 0.0)
            return LOW;

        double logcop = biv_unif_CopulaClayton(a1,a2,rho1,nu);
        if(!R_FINITE(logcop)) return LOW;
        dens = logcop + log(g1) + log(g2);
    } else {
        return LOW; /* discrete Clayton-like copula likelihood not implemented */
    }
}
/******************copula skew gaussian *************************/
if(type_cop==3)
{
    if(!(model==16||model==11||model==30)) { // continuous models
        /* nu is the bounded eta parameter in the skew-Gaussian paper. */
        if(model==1 || model==22 || model==25) nu=nuis[2];
        if(model==12 || model==21 || model==26 || model==59) nu=nuis[3];
        if(model==50 || model==42) nu=nuis[5];
        if(model==28 || model==33) nu=nuis[6];

        if(!R_FINITE(nu) || fabs(nu) >= 1.0)
            return LOW;
        if(!R_FINITE(g1) || !R_FINITE(g2) || g1 <= 0.0 || g2 <= 0.0)
            return LOW;

        const double copdens = biv_unif_CopulaSkewGauss(a1,a2,rho1,nu/sqrt(1-nu*nu));
        if(!R_FINITE(copdens) || copdens <= 0.0) return LOW;
        dens=log(copdens)+log(g1)+log(g2);
    } else {
        /* Poisson/Binomial/BinomialNeg.  Their switch entries above store
           Gaussian-score versions of the marginal CDF endpoints because the
           Gaussian copula uses them directly.  Transform back to the exact
           marginal uniforms and integrate the skew-Gaussian copula density
           over that probability rectangle. */
        nu=nuis[2];
        if(!R_FINITE(nu) || fabs(nu) >= 1.0) return LOW;
        const double uhi = pnorm(a1,0,1,1,0);
        const double ulo = pnorm(a2,0,1,1,0);
        const double vhi = pnorm(g1,0,1,1,0);
        const double vlo = pnorm(g2,0,1,1,0);
        const double alpha = nu/sqrt(1.0-nu*nu);
        double prob = geo_skewcop_rect_prob(ulo,uhi,vlo,vhi,rho1,alpha);
        if(!R_FINITE(prob) || prob <= 0.0) return LOW;
        prob = fmin(fmin(uhi-ulo, vhi-vlo), fmax(DBL_MIN, prob));
        dens=log(prob);
    }
}

if(cond) {
   if(!(model==16 || model==11 || model==30)) {   // continui
       if(type_cop == 1 && model == 1) {
           /* log marginal Gaussian density at site 2, evaluated without
            * materializing dnorm(). */
           const double logg2 = -0.5*log(2.0*M_PI*nuis[1]) - 0.5*b2*b2;
           dens -= logg2;
       } else {
           dens = dens - log(fmax(DBL_MIN, g2));
       }
   } else {                                      // discreti
       dens = dens - log(fmax(DBL_MIN, b2));
   }
}



if(!R_FINITE(dens)) return LOW;
return(dens);
}

/***********************/
/***********************/
/******* some marginals (log)pdf  *****************/
/***********************/
/***********************/

/* cdf skewlaplace */
double pSkewLaplace(double x, double sk, double mm, double vv) {
    double s = sqrt(vv);
    if (x >= mm) {
        return sk + (1.0 - sk) * (1.0 - exp(-sk * (x - mm) / s));
    } else {
        return sk * exp((1.0 - sk) * (x - mm) / s);
    }
}
/* pdf skewlaplace */
double dSkewLaplace(double x, double sk, double mm, double vv) {


    double s = sqrt(vv);
    double c = (sk * (1.0 - sk)) / s;

    if (x < mm) {
        return c * exp((1.0 - sk) * (x - mm) / s);
    } else {
        return c * exp(-sk * (x - mm) / s);
    }
}

/*******************************************************/

static const double GL_x[20] = {
    -0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188,
    -0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154195,
    -0.2277858511416451, -0.0765265211334973,  0.0765265211334973,  0.2277858511416451,
     0.3737060887154195,  0.5108670019508271,  0.6360536807265150,  0.7463319064601508,
     0.8391169718222188,  0.9122344282513259,  0.9639719272779138,  0.9931285991850949
};

static const double GL_w[20] = {
    0.0176140071391521, 0.0406014298003869, 0.0626720483341091, 0.0832767415767048,
    0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183820,
    0.1491729864726037, 0.1527533871307259, 0.1527533871307259, 0.1491729864726037,
    0.1420961093183820, 0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
    0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521
};

// Costante precomputata
static const double INV_2PI = 1.0 / (2.0 * M_PI);


/* Owen's T function ottimizzata */
double owens_t_optimized(double h, double a) {
    // Casi limite
    if (a == 0.0) return 0.0;
    if (fabs(a) < 1e-15) return 0.0;
    
    // Per |h| molto grande, Owen's T ≈ 0
    if (fabs(h) > 37.0) return 0.0;
    
    // Precomputa valori comuni
    const double h_squared = h * h;
    const double half_a = a * 0.5;
    const double neg_half_h2 = -0.5 * h_squared;
    
    double sum = 0.0;
    
    // Loop unrolling parziale per migliorare le performance
    for (int i = 0; i < 20; i += 2) {
        // Primo elemento
        double x1 = half_a * (GL_x[i] + 1.0);
        double x1_squared = x1 * x1;
        double denom1 = 1.0 + x1_squared;
        double fx1 = exp(neg_half_h2 * denom1) / denom1;
        sum += GL_w[i] * fx1;
        
        // Secondo elemento (se esiste)
        if (i + 1 < 20) {
            double x2 = half_a * (GL_x[i + 1] + 1.0);
            double x2_squared = x2 * x2;
            double denom2 = 1.0 + x2_squared;
            double fx2 = exp(neg_half_h2 * denom2) / denom2;
            sum += GL_w[i + 1] * fx2;
        }
    }
    
    return half_a * sum * INV_2PI;
}


/* Versione corretta senza simmetrie problematiche */
double psn(double x, double omega, double alpha, double tau) {
    if (omega <= 0.0) return NAN;
    
    const double z = x / omega;
    
    // Casi limite
    if (alpha == 0.0) return pnorm(z, 0.0, 1.0, 1, 0);
    if (fabs(z) > 8.0) return (z > 0) ? 1.0 : 0.0;
    
    // Calcola direttamente senza manipolare le simmetrie
    const double phi_z = pnorm(z, 0.0, 1.0, 1, 0);
    const double owens = owens_t_optimized(z, alpha);  // Usa z e alpha originali
    
    return phi_z - 2.0 * owens;
}


double dsn(double x, double omega, double alpha, double tau) {
    double z = x/ omega;
    double t = alpha * z + tau;
    return (2.0 / omega) * dnorm(z, 0.0, 1.0, 0) * pnorm(t, 0.0, 1.0, 1, 0);
}



double qsn(double p, double omega, double alpha, double tau, double tol) {
    if (omega <= 0) Rf_error("omega must be positive.");
    if (p <= 0.0) return -INFINITY;
    if (p >= 1.0) return INFINITY;

    int lower_tail = 1, log_p = 0;

    // Calcolo limiti iniziali basati sulla chi-quadro
    double max_q = sqrt(qchisq(p, 1, lower_tail, log_p)) + tau;
    double min_q = -sqrt(qchisq(1 - p, 1, lower_tail, log_p)) + tau;

    // Caso speciale tau == 0 e alpha infinito
    if (tau == 0 && isinf(alpha)) {
        return omega * (alpha > 0 ? max_q : min_q);
    }

    double abs_alpha = fabs(alpha);
    double p_adj = (alpha < 0) ? (1 - p) : p;

    // Inizializza estremi per la ricerca (bisezione + regula falsi alternata)
    double xa = qnorm(p_adj, 0.0, 1.0, 1, 0);
    double xb = sqrt(qchisq(p_adj, 1, 1, 0)) + fabs(tau);

    double fa = psn(xa, 1.0, abs_alpha, tau) - p_adj;
    double fb = psn(xb, 1.0, abs_alpha, tau) - p_adj;

    // Se i valori iniziali sono già abbastanza vicini, ritorna subito
    if (fabs(fa) < tol) {
        return omega * ((alpha < 0) ? -xa : xa);
    }
    if (fabs(fb) < tol) {
        return omega * ((alpha < 0) ? -xb : xb);
    }

    double xc = 0.0, fc = 0.0;
    int use_regula = 0;
    int iter = 0, max_iter = 100;

    while (iter < max_iter) {
        if (use_regula && (fb - fa) != 0.0) {
            // Regula falsi step
            xc = xb - fb * (xb - xa) / (fb - fa);
        } else {
            // Bisezione step
            xc = 0.5 * (xa + xb);
        }

        fc = psn(xc, 1.0, abs_alpha, tau) - p_adj;

        if (fabs(fc) < tol) {
            break;  // Convergenza raggiunta
        }

        // Aggiorna gli estremi dell'intervallo
        if (fc * fa < 0) {
            xb = xc;
            fb = fc;
        } else {
            xa = xc;
            fa = fc;
        }

        use_regula = !use_regula;
        iter++;
    }

    double x = (alpha < 0) ? -xc : xc;
    return omega * x;
}



double one_log_SkewGauss(double z, double m, double vari, double skew)
{
    const double LOG_2 = 0.6931471805599453; // log(2)
    double q = z - m;
    double skew2 = skew * skew;
    double omega = skew2 + vari;
    double sqrt_omega = sqrt(omega);
    double sqrt_var_omega = sqrt(vari * omega);

    double z1 = q / sqrt_omega;
    double z2 = skew * q / sqrt_var_omega;

    double res = LOG_2 - 0.5 * log(omega)
                 + dnorm(z1, 0.0, 1.0, 1)
                 + pnorm(z2, 0.0, 1.0, 1, 1);

    return res;
}

double one_log_SkewLaplace(double z, double mu, double sill, double skew) {
    double sigma = sqrt(sill);
    double z_std = (z - mu) / sigma; // Standardizzazione
    
    if (z_std >= 0) {
        return log(skew * (1 - skew)) - log(sigma) - (1 - skew) * z_std;
    } else {
        return log(skew * (1 - skew)) - log(sigma) + skew * z_std;
    }
}
double one_log_tukeyh(double z, double m, double sill, double tail)
{
    if (!isfinite(z) || !isfinite(m) || !isfinite(sill) || !isfinite(tail) ||
        sill <= 0.0 || tail < 0.0) return NAN;
    const double sqrt_sill = sqrt(sill);
    const double q = (z - m) / sqrt_sill;
    const double x = inv_tukey_h(q, tail);
    const double jac = jac_tukey_h(q, tail);
    if (!isfinite(x) || !isfinite(jac) || jac <= 0.0) return NAN;
    return dnorm(x, 0.0, 1.0, 1) + log(jac) - log(sqrt_sill);
}


double one_log_tukeyhh(double z,double m, double sill, double h1,double h2)
{
  double  res=0.0;
 if(z>=m){
    res=one_log_tukeyh(z,m,sill,h2);
          }
  if(z<m){
    res=one_log_tukeyh(z,m,sill,h1);
         }
  return(res);
}
double one_log_T(double z,double m, double sill, double df)
{
  double  res=0.0;
  double q=(z-m)/sqrt(sill);
    res=lgammafn(0.5*(df+1))-(0.5*(df+1))*log1p(q*q/df)-log(sqrt(M_PI*df))-lgammafn(df/2)-0.5*log(sill);
  return(res);
}

double one_log_sas(double z, double m, double skew, double tail, double vari)
{
  
  const double sigma = sqrt(vari);
  const double x = (z - m) / sigma;
  const double b = tail * asinh(x) - skew;
  const double S = sinh(b);
  const double logphi = dnorm(S, 0.0, 1.0, 1);
  const double log_num = log(tail) + 0.5 * log1p(S*S) + logphi;
  const double log_den = 0.5 * log1p(x*x) + log(sigma);
  return log_num - log_den;
}


double one_log_beta(double z, double shape1,double shape2,double min,double  max)
{
  double  res=0.0;
  double q=(z-min)/(max-min);
  res=(shape1/2-1)*log(q)+(shape2/2-1)*log1p(-q)+lgammafn(0.5*(shape1+shape2))-lgammafn(shape1/2)-lgammafn(shape2/2)-log(max-min);
  return(res);
}

double one_log_kumma2(double z,double m, double shape1,double shape2,double min,double  max)
{
  double  res=0.0,k;
  double q=(z-min)/(max-min);k=1-R_pow(q,shape2);
  double m1=1/(1+exp(-m));
  double shapei=log(0.5)/log1p(-R_pow(m1,shape2));
  res=log(shapei)+log(shape2)+(shape2-1)*log(q)+(shapei-1)*log(k)-log(max-min);
  return(res);
}

double one_log_kumma(double z,double m, double shape1,double shape2,double min,double  max)
{
  double  res=0.0,k;
  double q=(z-min)/(max-min);k=1-R_pow(q,shape2);
  res=log(shape1)+log(shape2)+(shape2-1)*log(q)+(shape1-1)*log(k)-log(max-min);
  return(res);
}

double one_log_loggaussian(double z,double m, double sill)
{
  double  res=0.0;
  double q=z*exp(sill/2);
  res=-0.5*R_pow((log(q)-m),2)/sill-log(q)-log(sqrt(sill))-0.5*log(2*M_PI)+sill/2;
  return(res);
}
double one_log_weibull(double z,double m, double shape)
{
  double  res=0.0;
  double scale1=exp(m)/(gammafn(1+1/shape));
  res=log(shape) -log(scale1)  +  (shape-1)*(log(z)-log(scale1))   -    R_pow(z/scale1,shape);
  return(res);
}
double one_log_gamma(double z,double m, double shape)
{
  double  res=0.0;
  res=(shape/2)*log(shape/(2*exp(m)))+(shape/2-1)*log(z)-(shape/(2*exp(m)))*z-log(gammafn(shape/2));
  return(res);
}
/************************************************************************/
double one_log_two_pieceTukey(double z,double m, double sill,double tail, double eta)
{
  double  res=0.0;
  double y=(z-m)/sqrt(sill);
 if(y>=0)res=one_log_tukeyh(y/(1-eta),0,1,tail);       
 if(y<0) res=one_log_tukeyh(y/(1+eta),0,1,tail); 
  return(res-log(sqrt(sill)));
}
double one_log_two_pieceT(double z,double m, double sill, double df, double eta)
{
  double  res=0.0;
  double y=(z-m)/sqrt(sill);
 if(y>=0) res=one_log_T(y/(1-eta),0,1,df);       
 if(y<0)  res=one_log_T(y/(1+eta),0,1,df); 
  return(res-log(sqrt(sill)));
}
double one_log_two_pieceGauss(double z,double m, double sill, double eta)
{
  double  res=0.0;
  double y=(z-m)/sqrt(sill);
 if(y>=0) res=dnorm(y/(1-eta),0,1,1);       
 if(y<0)  res=dnorm(y/(1+eta),0,1,1); 
  return(res-log(sqrt(sill)));
}

double one_log_wrapped(double alfa,double u,double mi,double sill)
{
    double x,wrap_gauss=0.0; 
    x=u-2*atan(mi)-M_PI; 
    double k1=-alfa;

     while(k1<=alfa){
      wrap_gauss +=  dnorm((x+2*k1*M_PI),0,sqrt(sill),0); 
     
     k1 = k1+1;
      }
return(log(wrap_gauss));
}

/************************************************************************/
double one_log_gammagem(double z,double shape,double n)
{
  double  res=0.0;
  res=(shape/2)*log(n/2)+(shape/2-1)*log(z)-0.5*n*z-lgammafn(shape/2);
  return(res);
}

double one_log_bomidal(double z,double m, double sill,double nu,double delta, double eta)
{
  double  res=0.0;
  double q=(z-m)/sqrt(sill);
  double alpha=2*(delta+1)/nu;
  double nn=R_pow(2,1-alpha/2);
 if(z>=m){
    res=log(alpha)+(alpha-1)*log(q)-(alpha-1)*log1p(-eta)-log(2)+one_log_gammagem(R_pow(z/(1-eta),alpha),nu,nn)-0.5*log(sill);    
          }
  if(z<m){
    res=log(alpha)+(alpha-1)*log(-q)-(alpha-1)*log1p(eta)-log(2)+one_log_gammagem(R_pow(-z/(1+eta),alpha),nu,nn)-0.5*log(sill);
         }
  return(res);
}
double one_log_BinomnegZIP(int z, double n, double mu, double mup)
{
    // Controlli di validità sui parametri
    if (z < 0 || n <= 0) {
        return -700.0;  // Parametri non validi
    }
    
    // CLIPPING DEI PARAMETRI LINEARI per evitare overflow in pnorm
    const double mu_safe = fmax(-12.0, fmin(12.0, mu));
    const double mup_safe = fmax(-12.0, fmin(12.0, mup));
    
    // Calcolo delle probabilità con clipping
    double p = pnorm(mup_safe, 0, 1, 1, 0);    // Probabilità di zero inflation
    double pp = pnorm(mu_safe, 0, 1, 1, 0);    // Parametro della binomiale negativa
    
    // CLIPPING DELLE PROBABILITÀ per stabilità numerica
    p = fmax(DBL_EPSILON, fmin(1- DBL_EPSILON , p));
    pp = fmax(DBL_EPSILON, fmin(1- DBL_EPSILON , pp));
    
    double res = 0.0;
    
    if (z == 0) {
        // Caso z = 0: P(Z=0) = p + (1-p) * P(NB=0)
        
        // Calcola P(NB=0) usando dnbinom per robustezza
        const double nb_zero = dnbinom(0, n, pp, 0);  // 0 = give_log=FALSE
        
        // Controllo per valori non finiti
        if (!R_FINITE(nb_zero) || nb_zero < 0) {
            return -700.0;
        }
        const double prob_zero = p + (1 - p) * nb_zero;
        if (prob_zero <= 1e-300) {
            return -700.0;
        }
        res = log(prob_zero);
        
    } else { // z > 0
        // Caso z > 0: P(Z=z) = (1-p) * P(NB=z)
        
        // Calcola log P(NB=z) direttamente
        const double log_nb_z = dnbinom(z, n, pp, 1);  // 1 = give_log=TRUE

        if (!R_FINITE(log_nb_z)) {
            return -700.0;
        }
        const double log_one_minus_p = log1p(-p);
        
        res = log_one_minus_p + log_nb_z;
    }
    
    // Controllo finale per valori non finiti
    if (!R_FINITE(res)) {
        return -700.0;
    }
    
    return res;
}
double one_log_PoisZIP(int z,double lambda, double mup)
{
  double  res=0.0;
  double  p=pnorm(mup,0,1,1,0);
 if(z==0){
    res=log(p+(1-p)*dpois(0,lambda,0));    
          }
  if(z>0){
    res=log1p(-p)+dpois(z,lambda,1);   
         }
  return(res);
}


double one_log_negbinom_marg(int u, int N, double p)
{
#define P_CLAMP 1e-10

    if (u < 0 || N <= 0) return R_NegInf;

    p = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, p));

    return dnbinom(u, N, p, 1);  /* give_log = TRUE */
}

double one_log_loglogistic(double z,double m, double shape)
{
  double  res;
  double c=gammafn(1+1/shape)*gammafn(1-1/shape);
  double k=R_pow(c*z/m,shape)+1;
  //res=log(c)+log(shape/m)+(shape-1)*log(c*z/m)-2*loh(k);
  res=log((c*shape/m)*R_pow((c*z/m),shape-1)*R_pow(k,-2));
  return(res);
}
double one_log_logistic(double z,double m, double sill)
{
  double  res;
  double k=exp((z-m)/sqrt(sill));
  res=log(k)-2*log(k+1)-0.5*log(sill);
  return(res);
}

// gentz
double bvu(double sh, double sk, double r) {
    const double two_pi = 6.28318530717958647693;
    double bvn = 0.0;
    double h = sh, k = sk;
    double hk = h * k;
    double hs, asr, sn, as, a, b, c, d, bs, rs, xs;
    int i, lg, ng;

    // Tabelle COMPLETAMENTE IDENTICHE - zero modifiche
    static const double W[10][3] = {
        {0.1713244923791705, 0.04717533638651177, 0.01761400713915212},
        {0.3607615730481384, 0.1069393259953183, 0.04060142980038694},
        {0.4679139345726904, 0.1600783285433464, 0.06267204833410906},
        {0.0,                0.2031674267230659, 0.08327674157670475},
        {0.0,                0.2334925365383547, 0.1019301198172404},
        {0.0,                0.2491470458134029, 0.1181945319615184},
        {0.0,                0.0,                0.1316886384491766},
        {0.0,                0.0,                0.1420961093183821},
        {0.0,                0.0,                0.1491729864726037},
        {0.0,                0.0,                0.1527533871307259}
    };
    static const double X[10][3] = {
        {-0.9324695142031522, -0.9815606342467191, -0.9931285991850949},
        {-0.6612093864662647, -0.9041172563704750, -0.9639719272779138},
        {-0.2386191860831970, -0.7699026741943050, -0.9122344282513259},
        {0.0,                 -0.5873179542866171, -0.8391169718222188},
        {0.0,                 -0.3678314989981802, -0.7463319064601508},
        {0.0,                 -0.1252334085114692, -0.6360536807265150},
        {0.0,                  0.0,                -0.5108670019508271},
        {0.0,                  0.0,                -0.3737060887154196},
        {0.0,                  0.0,                -0.2277858511416451},
        {0.0,                  0.0,                -0.07652652113349733}
    };

    // LOGICA IDENTICA - solo variabile ausiliaria per fabs(r)
    double abs_r = fabs(r);
    
    if (abs_r < 0.3) {
        ng = 0; lg = 3;
    } else if (abs_r < 0.75) {
        ng = 1; lg = 6;
    } else {
        ng = 2; lg = 10;
    }

    if (abs_r < 0.925) {
        hs = (h * h + k * k) / 2.0;
        asr = asin(r);
        
        // Loop IDENTICO - solo pre-computazione di 2.0 * two_pi
        double inv_4pi = 1.0 / (2.0 * two_pi);
        
        for (i = 0; i < lg; ++i) {
            sn = sin(asr * (X[i][ng] + 1.0) / 2.0);
            bvn += W[i][ng] * exp((sn * hk - hs) / (1.0 - sn * sn));
            
            sn = sin(asr * (-X[i][ng] + 1.0) / 2.0);
            bvn += W[i][ng] * exp((sn * hk - hs) / (1.0 - sn * sn));
        }
        
        // Formula IDENTICA - solo ottimizzazione della divisione finale
        bvn = bvn * asr * inv_4pi + pnorm(-h, 0.0, 1.0, 1, 0) * pnorm(-k, 0.0, 1.0, 1, 0);
        
    } else {
        if (r < 0.0) {
            k = -k;
            hk = -hk;
        }

        if (abs_r < 1.0) {
            as = (1.0 - r) * (1.0 + r);
            a = sqrt(as);
            bs = (h - k) * (h - k);
            c = (4.0 - hk) / 8.0;
            d = (12.0 - hk) / 16.0;
            
            // Pre-calcolo per evitare divisioni ripetute
            double bs_over_5 = bs / 5.0;
            double as_squared = as * as;
            
            bvn = a * exp(-(bs / as + hk) / 2.0) * 
                  (1.0 - c * (bs - as) * (1.0 - d * bs_over_5) / 3.0 + c * d * as_squared / 5.0);

            if (hk > -160.0) {
                b = sqrt(bs);
                bvn -= exp(-hk / 2.0) * sqrt(two_pi) * pnorm(-b / a, 0.0, 1.0, 1, 0) * b * 
                       (1.0 - c * bs * (1.0 - d * bs_over_5) / 3.0);
            }

            a = a / 2.0;
            for (i = 0; i < lg; ++i) {
                xs = (a * (X[i][ng] + 1.0));
                xs *= xs;
                rs = sqrt(1.0 - xs);
                bvn += a * W[i][ng] * (
                    exp(-bs / (2.0 * xs) - hk / (1.0 + rs)) / rs -
                    exp(-(bs / xs + hk) / 2.0) * (1.0 + c * xs * (1.0 + d * xs))
                );

                xs = as * (1.0 - X[i][ng]) * (1.0 - X[i][ng]) / 4.0;
                rs = sqrt(1.0 - xs);
                bvn += a * W[i][ng] * exp(-(bs / xs + hk) / 2.0) *
                    (exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs -
                     (1.0 + c * xs * (1.0 + d * xs)));
            }
            bvn = -bvn / two_pi;
        }

        if (r > 0.0)
            bvn += pnorm(-fmax(h, k), 0.0, 1.0, 1, 0);
        if (r < 0.0)
            bvn = -bvn + fmax(0.0, pnorm(-h, 0.0, 1.0, 1, 0) - pnorm(-k, 0.0, 1.0, 1, 0));
    }

    return bvn;
}
double bvn_hybrid(double h, double k, double rho) {
    double abs_rho = fabs(rho);
    // Gestione casi limite
    if (abs_rho >= 1.0) {
        if (rho >= 1.0) {
            return pnorm(fmin(h, k), 0.0, 1.0, 1, 0);
        } else {
            return fmax(0.0, pnorm(h, 0.0, 1.0, 1, 0) - pnorm(-k, 0.0, 1.0, 1, 0));
        }
    }
    // Usa Drezner-Wesolowsky solo per correlazioni molto basse
   // if (abs_rho < HYBRID_THRESHOLD) {
   //     return bvn_drezner_wesolowsky(h, k, rho);
   // } else {
        // versionec di gentz
        return bvu(h, k, rho);
    //}
}
// Wrapper IDENTICO al tuo bvnmvn
double bvnmvn(const double lower[2], const double upper[2], const int infin[2], double corr) {
    double result = 0.0;

    if (infin[0] == 2 && infin[1] == 2) {
        result = bvn_hybrid(lower[0], lower[1], corr)
               - bvn_hybrid(upper[0], lower[1], corr)
               - bvn_hybrid(lower[0], upper[1], corr)
               + bvn_hybrid(upper[0], upper[1], corr);
    } else if (infin[0] == 2 && infin[1] == 1) {
        result = bvn_hybrid(lower[0], lower[1], corr)
               - bvn_hybrid(upper[0], lower[1], corr);
    } else if (infin[0] == 1 && infin[1] == 2) {
        result = bvn_hybrid(lower[0], lower[1], corr)
               - bvn_hybrid(lower[0], upper[1], corr);
    } else if (infin[0] == 2 && infin[1] == 0) {
        result = bvn_hybrid(-upper[0], -upper[1], corr)
               - bvn_hybrid(-lower[0], -upper[1], corr);
    } else if (infin[0] == 0 && infin[1] == 2) {
        result = bvn_hybrid(-upper[0], -upper[1], corr)
               - bvn_hybrid(-upper[0], -lower[1], corr);
    } else if (infin[0] == 1 && infin[1] == 0) {
        result = bvn_hybrid(lower[0], -upper[1], -corr);
    } else if (infin[0] == 0 && infin[1] == 1) {
        result = bvn_hybrid(-upper[0], lower[1], -corr);
    } else if (infin[0] == 1 && infin[1] == 1) {
        result = bvn_hybrid(lower[0], lower[1], corr);
    } else if (infin[0] == 0 && infin[1] == 0) {
        result = bvn_hybrid(-upper[0], -upper[1], corr);
    }

    return result;
}

 double safe_pbnorm22(double x, double y, double rho)
{
    // casi esatti ai bordi — nessuna approssimazione
    if (!R_FINITE(x) && x < 0) return 0.0;          // Φ₂(-∞, y) = 0
    if (!R_FINITE(x) && x > 0) return pnorm(y,0,1,1,0); // Φ₂(+∞, y) = Φ(y)
    if (!R_FINITE(y) && y < 0) return 0.0;          // Φ₂(x, -∞) = 0
    if (!R_FINITE(y) && y > 0) return pnorm(x,0,1,1,0); // Φ₂(x, +∞) = Φ(x)
    return pbnorm22(x, y, rho);
}
double pbnorm22(double lim1,double lim2,double corr){
    double lower[2] = {-INFINITY, -INFINITY};
    double upper[2] = {lim1, lim2};
    int infin[2] = {0, 0};  // Entrambi (-inf, upper]
    double res=bvnmvn(lower, upper, infin, corr);
    //Rprintf("%f = %f %f %f  \n",res,lim1,lim2,corr);
    return res;
}

// CDF of a bivariate Gaussian distribution using bvnmvn
double cdf_norm(double lim1, double lim2, double a11, double a12) {
    double res = 0.0;
    double lower[2] = {0.0, 0.0};
    double upper[2] = {lim1 / sqrt(a11), lim2 / sqrt(a11)};
    int infin[2] = {0, 0};
    double corr = a12 / a11;
    double auxil = 1.0 - R_pow(corr, 2);
    double det = R_pow(a11, 2) - R_pow(a12, 2);
    res = a11 * sqrt(auxil / det) * bvnmvn(lower, upper, infin, corr);
    return res;
}
// CDF of a bivariate Gaussian distribution using general covariance
double cdf_norm2(double lim1, double lim2, double a11, double a12, double a22) {
    double res = 0.0;
    double lower[2] = {0.0, 0.0};
    double upper[2] = {lim1 / sqrt(a11), lim2 / sqrt(a22)};
    int infin[2] = {0, 0};
    double corr = a12 / sqrt(a11 * a22);
    double auxil = 1.0 - R_pow(corr, 2);
    double det = a11 * a22 - R_pow(a12, 2);
    res = sqrt(a11 * a22) * sqrt(auxil / det) * bvnmvn(lower, upper, infin, corr);

    return res;
}
