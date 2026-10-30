/*
 * Conditional simulation of the DIRECT asymmetric Laplace random field
 *
 *   Y_i = eta_i + sigma * [(Z1_i^2+Z2_i^2)/(2p)
 *                           - (Z3_i^2+Z4_i^2)/(2(1-p))]
 *
 * The four Z fields are independent Gaussian fields with common covariance R.
 * At an observed site and conditional on the other latent sites,
 * Z_i ~ N_4(mu_i, v_i I_4), constrained to Y_i = y_i.
 *
 * Parameterize the constraint using t >= 0, theta_plus, theta_minus:
 *   A = (Z1^2+Z2^2)/2 = p * (t + max(y,0)),
 *   B = (Z3^2+Z4^2)/2 = (1-p) * (t + max(-y,0)).
 *
 * With m_+ = ||mu_{1:2}||, m_- = ||mu_{3:4}||, the EXACT target
 * for the radial parameter, relative to dt, is proportional to
 *
 *   exp(-t/v) I_0(m_+ sqrt(2 A)/v) I_0(m_- sqrt(2 B)/v).
 *
 * Angles given t are independent von Mises with concentration
 * m_+ sqrt(2 A)/v and m_- sqrt(2 B)/v. This parameterization
 * AUTOMATICALLY includes the coarea/Jacobian factor for conditioning
 * on the quadratic constraint (no uncorrected surface-measure target).
 *
 * Radial updates use stepping-out/shrinkage slice sampling, and angular
 * updates use exact Best--Fisher rejection draws. Gibbs sweeps visit sites
 * in random order. The stationary target is the correct Gaussian latent
 * conditional distribution; finite burn-in/mixing diagnostics remain needed.
 *
 * The routine returns n x 4 x nrep latent Gaussian values, and R uses
 * standard Gaussian conditional simulation to obtain prediction values.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <limits.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define GEO_SL_MAX_STEPS 1000
#define GEO_SL_MAX_SHRINK 10000

static double geo_sl_u01(void) {
  double u = unif_rand();
  return (u > 0.0) ? u : 0x1p-53;
}

static double geo_sl_logi0(double x) {
  if (x == 0.0) return 0.0;
  if (x < 0.0 || !R_FINITE(x))
    Rf_error("Invalid SkewLaplace Bessel argument.");
  /* R C API: expo=2 returns exp(-x)*I0(x); expo=1 is unscaled. */
  const double scaled = bessel_i(x, 0.0, 2.0);
  if (!(scaled > 0.0) || !R_FINITE(scaled))
    Rf_error("SkewLaplace scaled Bessel evaluation failed.");
  return x + log(scaled);
}

typedef struct {
  double v, py, my, c_plus, c_minus, width;
} geo_sl_radial_info;

static double geo_sl_radial_logf(double t, const geo_sl_radial_info *r) {
  if (t < 0.0 || !R_FINITE(t)) return R_NegInf;
  const double kp = r->c_plus * sqrt(t + r->py);
  const double km = r->c_minus * sqrt(t + r->my);
  return -t / r->v + geo_sl_logi0(kp) + geo_sl_logi0(km);
}

static double geo_sl_radial_slice(double t, const geo_sl_radial_info *r) {
  const double old_log = geo_sl_radial_logf(t, r);
  if (!R_FINITE(old_log)) Rf_error("Non-finite SkewLaplace radial log-density.");
  const double level = old_log + log(geo_sl_u01());
  /* Fixed within each full conditional: state-dependent slice widths
     without an additional MH correction would destroy detailed balance. */
  double w = r->width;
  if (!R_FINITE(w) || w <= 0.0) w = r->v;
  double lo = t - w * geo_sl_u01();
  double hi = lo + w;
  /* Clip only the LEFT endpoint: moving hi after clipping lo changes the
     randomized bracket and is not the standard invariant slice kernel. */
  if (lo < 0.0) lo = 0.0;
  int it = 0;
  while (lo > 0.0 && geo_sl_radial_logf(lo, r) > level) {
    if (++it > GEO_SL_MAX_STEPS) Rf_error("SkewLaplace radial slice: left stepping-out failed.");
    lo = (lo > w) ? lo - w : 0.0;
  }
  it = 0;
  while (geo_sl_radial_logf(hi, r) > level) {
    if (++it > GEO_SL_MAX_STEPS) Rf_error("SkewLaplace radial slice: right stepping-out failed.");
    hi += w;
  }
  for (it = 0; it < GEO_SL_MAX_SHRINK; ++it) {
    const double proposed = lo + (hi - lo) * geo_sl_u01();
    if (geo_sl_radial_logf(proposed, r) >= level) return proposed;
    if (proposed < t) lo = proposed;
    else hi = proposed;
  }
  Rf_error("SkewLaplace radial slice: shrinkage failed to find a point.");
  return t;
}

/* Best and Fisher (1979) exact rejection algorithm for von Mises angles. */
static double geo_sl_vonmises(double kappa) {
  if (kappa < 1e-5) {
    /* Exact uniform-envelope rejection; inexpensive for tiny kappa. */
    for (;;) {
      const double angle = 2.0 * M_PI * geo_sl_u01() - M_PI;
      if (log(geo_sl_u01()) <= kappa * (cos(angle) - 1.0)) return angle;
    }
  }
  const double a = 1.0 + sqrt(1.0 + 4.0 * kappa * kappa);
  const double b = (a - sqrt(2.0 * a)) / (2.0 * kappa);
  const double r = (1.0 + b * b) / (2.0 * b);
  for (;;) {
    const double z = cos(M_PI * geo_sl_u01());
    double f = (1.0 + r * z) / (r + z);
    if (f < -1.0) f = -1.0;
    if (f > 1.0) f = 1.0;
    const double c = kappa * (r - f);
    const double u2 = geo_sl_u01();
    if (u2 < c * (2.0 - c) ||
        (c > 0.0 && log(c / u2) + 1.0 - c >= 0.0)) {
      return (geo_sl_u01() < 0.5 ? -1.0 : 1.0) * acos(f);
    }
  }
}

static void geo_sl_update_site(int i, int n, double p, const double *y,
                               const double *Q, double *z, double *Qz,
                               double *radial) {
  const double qii = Q[i + (size_t)n * i];
  const double v = 1.0 / qii;
  double mu[4], newz[4];
  for (int k = 0; k < 4; ++k) {
    mu[k] = z[i + (size_t)n * k] - Qz[i + (size_t)n * k] / qii;
  }
  const double mplus = hypot(mu[0], mu[1]);
  const double mminus = hypot(mu[2], mu[3]);
  const double py = y[i] > 0.0 ? y[i] : 0.0;
  const double my = y[i] < 0.0 ? -y[i] : 0.0;
  geo_sl_radial_info info;
  info.v = v;
  info.py = py;
  info.my = my;
  info.c_plus = mplus * sqrt(2.0 * p) / v;
  info.c_minus = mminus * sqrt(2.0 * (1.0 - p)) / v;
  info.width = 2.0 * sqrt(v * (v + 0.5*(mplus*mplus + mminus*mminus)));
  const double t = (mplus == 0.0 && mminus == 0.0)
                   ? -v * log(geo_sl_u01())
                   : geo_sl_radial_slice(radial[i], &info);
  radial[i] = t;
  const double rplus = sqrt(2.0 * p * (t + py));
  const double rminus = sqrt(2.0 * (1.0 - p) * (t + my));
  const double kp = rplus * mplus / v;
  const double km = rminus * mminus / v;
  const double theta_p = atan2(mu[1], mu[0]) + geo_sl_vonmises(kp);
  const double theta_m = atan2(mu[3], mu[2]) + geo_sl_vonmises(km);
  newz[0] = rplus * cos(theta_p);
  newz[1] = rplus * sin(theta_p);
  newz[2] = rminus * cos(theta_m);
  newz[3] = rminus * sin(theta_m);
  /* Maintain Qz incrementally. O(4n) arithmetic per updated site. */
  for (int k = 0; k < 4; ++k) {
    const size_t ik = i + (size_t)n * k;
    const double delta = newz[k] - z[ik];
    z[ik] = newz[k];
    if (delta != 0.0) {
      for (int j = 0; j < n; ++j) {
        Qz[j + (size_t)n * k] += Q[j + (size_t)n * i] * delta;
      }
    }
  }
}

/* Periodically rebuild Qz to prevent cumulative roundoff in long chains. */
static void geo_sl_recompute_Qz(int n, const double *Q, const double *z,
                                 double *Qz) {
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int j = 0; j < n; ++j)
        s += Q[i + (size_t)n*j] * z[j + (size_t)n*k];
      Qz[i + (size_t)n*k] = s;
    }
}

static void geo_sl_gibbs_sweep(int n, double p, const double *y, const double *Q,
                               double *z, double *Qz, double *radial, int *order) {
  for (int i = 0; i < n; ++i) order[i] = i;
  for (int i = n - 1; i > 0; --i) {
    const int j = (int) floor((i + 1.0) * geo_sl_u01());
    const int tmp = order[i]; order[i] = order[j]; order[j] = tmp;
  }
  for (int i = 0; i < n; ++i)
    geo_sl_update_site(order[i], n, p, y, Q, z, Qz, radial);
}

SEXP geo_skewlaplace_latent_gibbs(SEXP sy, SEXP sQ, SEXP sp,
                                  SEXP sburn, SEXP srep, SEXP sthin) {
  if (TYPEOF(sy) != REALSXP || TYPEOF(sQ) != REALSXP || TYPEOF(sp) != REALSXP ||
      TYPEOF(sburn) != INTSXP || TYPEOF(srep) != INTSXP || TYPEOF(sthin) != INTSXP)
    Rf_error("Invalid argument types for SkewLaplace latent Gibbs sampler.");
  const int n = (int) XLENGTH(sy);
  if (n < 1 || n > INT_MAX/4) Rf_error("Invalid observation size.");
  if (!Rf_isMatrix(sQ) || INTEGER(Rf_getAttrib(sQ, R_DimSymbol))[0] != n ||
      INTEGER(Rf_getAttrib(sQ, R_DimSymbol))[1] != n ||
      XLENGTH(sQ) != (R_xlen_t)n*n) Rf_error("Q must be an n by n precision matrix.");
  const double p = REAL(sp)[0];
  if (!R_FINITE(p) || p <= 0.0 || p >= 1.0)
    Rf_error("SkewLaplace skew must satisfy 0 < p < 1.");
  const int burn = INTEGER(sburn)[0], nr = INTEGER(srep)[0], thin = INTEGER(sthin)[0];
  if (burn < 0 || nr < 1 || thin < 1 ||
      ((double) burn + (double)(nr - 1) * thin) > INT_MAX)
    Rf_error("Invalid burn-in, thinning or replication count.");
  if ((double)n * 4.0 * nr > 1e8)
    Rf_error("Too many latent draws; decrease nrep or process in batches.");
  const double *Q = REAL(sQ), *y = REAL(sy);
  for (int i = 0; i < n; ++i) {
    if (!R_FINITE(y[i]) || !R_FINITE(Q[i + (size_t)n*i]) || Q[i + (size_t)n*i] <= 0.0)
      Rf_error("SkewLaplace observations or latent precision are invalid.");
  }
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)n*4*nr));
  SEXP dims = PROTECT(Rf_allocVector(INTSXP, 3));
  INTEGER(dims)[0] = n; INTEGER(dims)[1] = 4; INTEGER(dims)[2] = nr;
  Rf_setAttrib(out, R_DimSymbol, dims);
  double *z = (double *) R_alloc((size_t)n*4, sizeof(double));
  double *Qz = (double *) R_alloc((size_t)n*4, sizeof(double));
  double *radial = (double *) R_alloc((size_t)n, sizeof(double));
  int *order = (int *) R_alloc((size_t)n, sizeof(int));
  GetRNGstate();
  /* Independent marginal posterior initialization at each site (mu=0). */
  for (int i = 0; i < n; ++i) {
    const double t = -log(geo_sl_u01());
    radial[i] = t;
    const double rp = sqrt(2.0*p*(t + fmax(y[i], 0.0)));
    const double rm = sqrt(2.0*(1.0-p)*(t + fmax(-y[i], 0.0)));
    const double ap = 2.0*M_PI*geo_sl_u01();
    const double am = 2.0*M_PI*geo_sl_u01();
    z[i] = rp*cos(ap); z[i + (size_t)n] = rp*sin(ap);
    z[i + (size_t)n*2] = rm*cos(am); z[i + (size_t)n*3] = rm*sin(am);
  }
  geo_sl_recompute_Qz(n, Q, z, Qz);
  int sweep = 0;
  for (int b = 0; b < burn; ++b) {
    geo_sl_gibbs_sweep(n, p, y, Q, z, Qz, radial, order);
    ++sweep;
    if (sweep % 50 == 0) geo_sl_recompute_Qz(n, Q, z, Qz);
    if (sweep % 20 == 0) R_CheckUserInterrupt();
  }
  for (int r = 0; r < nr; ++r) {
    if (r > 0 || burn == 0) {
      /* Always move before storing the first draw when burn=0. */
      for (int t = 0; t < (r == 0 ? 1 : thin); ++t) {
        geo_sl_gibbs_sweep(n, p, y, Q, z, Qz, radial, order);
        ++sweep;
        if (sweep % 50 == 0) geo_sl_recompute_Qz(n, Q, z, Qz);
        if (sweep % 20 == 0) R_CheckUserInterrupt();
      }
    }
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < n; ++i)
        REAL(out)[i + (size_t)n*k + (size_t)n*4*r] = z[i + (size_t)n*k];
  }
  PutRNGstate();
  UNPROTECT(2);
  return out;
}
