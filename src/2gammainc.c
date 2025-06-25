


#include "header.h"


#define MAXITER1 2000
#define big 1.44115188075855872E+17
#define biginv 6.9388939039072284E-18

// Costanti ottimizzate
#define INV_SQRT_2PI 0.3989422804014327   // 1/sqrt(2*pi)
#define LOG_SQRT_2PI 0.9189385332046727   // log(sqrt(2*pi))
#define M_1_E 0.36787944117144232159       // 1/e

// Costanti per i regimi asintotici (valori tipici)
#ifndef SMALL
#define SMALL 20.0
#endif
#ifndef LARGE  
#define LARGE 200.0
#endif
#ifndef SMALLRATIO
#define SMALLRATIO 0.3
#endif
#ifndef LARGERATIO
#define LARGERATIO 4.5
#endif
#ifndef lanczos_g
#define lanczos_g 6.024680040776729583740234375
#endif
#ifndef NIC
#define NIC 8
#endif
#ifndef KIC  
#define KIC 25
#endif
#ifndef IGAM
#define IGAM 1
#endif
#ifndef IGAMC
#define IGAMC 2
#endif

// Cache per valori calcolati frequentemente
static struct {
    double a_cached;
    double lgamma_cached;
    int valid;
} lgamma_cache = {0, 0, 0};





//************************************** incomplete gamma*****************************************//

/**
 */
static inline double lgamma_cached(double a) {
    if (lgamma_cache.valid && lgamma_cache.a_cached == a) {
        return lgamma_cache.lgamma_cached;
    }
    
    lgamma_cache.a_cached = a;
    lgamma_cache.lgamma_cached = lgammafn(a);
    lgamma_cache.valid = 1;
    
    return lgamma_cache.lgamma_cached;
}

/**
 *  
 */
static inline double handle_edge_cases(double a, double x, int lower_part) {
    // Caso più comune prima per branch prediction
    if (__builtin_expect(x > 0 && a > 0 && R_FINITE(a) && R_FINITE(x), 1)) {
        return -1.0; // Non è un caso limite
    }
    
    if (x < 0 || a < 0) {
        return R_NaN;
    } else if (a == 0) {
        return (x > 0) ? (lower_part ? 1.0 : 0.0) : R_NaN;
    } else if (x == 0) {
        return lower_part ? 0.0 : 1.0;
    } else if (!R_FINITE(a)) {
        if (!R_FINITE(x)) return R_NaN;
        return lower_part ? 0.0 : 1.0;
    } else if (!R_FINITE(x)) {
        return lower_part ? 1.0 : 0.0;
    }
    
    return -1.0;
}

/**
 * 
 */
double ratevl(double x, const double num[], int M,
              const double denom[], int N) {
    int i, dir;
    double y, num_ans, denom_ans;
    double absx = fabs(x);
    const double *p;
    
    if (absx > 1) {
        dir = -1;
        p = num + M;
        y = 1 / x;
    } else {
        dir = 1;
        p = num;
        y = x;
    }
    
    num_ans = *p;
    p += dir;
    for (i = 1; i <= M; i++) {
        num_ans = num_ans * y + *p;
        p += dir;
    }
    
    if (absx > 1) {
        p = denom + N;
    } else {
        p = denom;
    }
    
    denom_ans = *p;
    p += dir;
    for (i = 1; i <= N; i++) {
        denom_ans = denom_ans * y + *p;
        p += dir;
    }
    
    if (absx > 1) {
        i = N - M;
        return pow(x, i) * num_ans / denom_ans;
    } else {
        return num_ans / denom_ans;
    }
}

/**
 * 
 */
double lanczos_sum_expg_scaled(double x) {
    return ratevl(x, lanczos_sum_expg_scaled_num,
                  sizeof(lanczos_sum_expg_scaled_num) / sizeof(lanczos_sum_expg_scaled_num[0]) - 1,
                  lanczos_sum_expg_scaled_denom,
                  sizeof(lanczos_sum_expg_scaled_denom) / sizeof(lanczos_sum_expg_scaled_denom[0]) - 1);
}

/**
 *  35: Versione ottimizzata di igam_fac
 */
double igam_fac(double a, double x) {
    double ax, fac, res, num;
    const double threshold = 0.4;

    if (__builtin_expect(fabs(a - x) <= threshold * fabs(a), 0)) {
        fac = a + lanczos_g - 0.5;
        res = sqrt(fac * M_1_E) / lanczos_sum_expg_scaled(a);

        if ((a < 200) && (x < 200)) {
            res *= exp(a - x) * pow(x / fac, a);
        } else {
            num = x - a - lanczos_g + 0.5;
            res *= exp(a * log1pmx(num / fac) + x * (0.5 - lanczos_g) / fac);
        }
        return res;
    }
    
    ax = a * log(x) - x - lgamma_cached(a);
    return (ax < -MAXLOG) ? 0.0 : exp(ax);
}

/**
 *  
 */
double igamc_continued_fraction(double a, double x) {
    int i;
    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;

    ax = igam_fac(a, x);
    if (__builtin_expect(ax == 0.0, 0)) {
        return 0.0;
    }

    y = 1.0 - a;
    z = x + y + 1.0;
    c = 0.0;
    
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = x + 1.0;
    qkm1 = z * x;
    ans = pkm1 / qkm1;

    for (i = 0; i < MAXITER1 - 1; i += 2) {
        // Prima iterazione
        c += 1.0;
        y += 1.0;
        z += 2.0;
        yc = y * c;
        pk = pkm1 * z - pkm2 * yc;
        qk = qkm1 * z - qkm2 * yc;
        
        if (__builtin_expect(qk != 0, 1)) {
            r = pk / qk;
            t = fabs((ans - r) / r);
            ans = r;
            if (__builtin_expect(t <= DBL_EPSILON, 0)) break;
        } else {
            t = 1.0;
        }
        
        if (__builtin_expect(fabs(pk) > big, 0)) {
            pkm2 *= biginv; pkm1 *= biginv;
            qkm2 *= biginv; qkm1 *= biginv;
        }
        
        pkm2 = pkm1; pkm1 = pk;
        qkm2 = qkm1; qkm1 = qk;
        
        // Seconda iterazione
        c += 1.0;
        y += 1.0;
        z += 2.0;
        yc = y * c;
        pk = pkm1 * z - pkm2 * yc;
        qk = qkm1 * z - qkm2 * yc;
        
        if (__builtin_expect(qk != 0, 1)) {
            r = pk / qk;
            t = fabs((ans - r) / r);
            ans = r;
            if (__builtin_expect(t <= DBL_EPSILON, 0)) break;
        } else {
            t = 1.0;
        }
        
        if (__builtin_expect(fabs(pk) > big, 0)) {
            pkm2 *= biginv; pkm1 *= biginv;
            qkm2 *= biginv; qkm1 *= biginv;
        }
        
        pkm2 = pkm1; pkm1 = pk;
        qkm2 = qkm1; qkm1 = qk;
    }

    return ans * ax;
}

/**
 *  
 */
double igam_series(double a, double x) {
    int i;
    double ans, ax, c, r;
    double compensation = 0.0, y, t;

    ax = igam_fac(a, x);
    if (__builtin_expect(ax == 0.0, 0)) {
        return 0.0;
    }

    r = a;
    c = 1.0;
    ans = 1.0;

    for (i = 0; i < MAXITER1; i++) {
        r += 1.0;
        c *= x / r;
        
        // Kahan summation
        y = c - compensation;
        t = ans + y;
        compensation = (t - ans) - y;
        ans = t;
        
        if (__builtin_expect(c <= DBL_EPSILON * ans, 0)) {
            break;
        }
    }

    return (ans * ax) / a;
}

/**
 *  40: Serie igamc ottimizzata
 */
double igamc_series(double a, double x) {
    int n;
    double fac = 1.0;
    double sum = 0.0;
    double term, logx;
    double compensation = 0.0, y, t;
    
    logx = log(x);
    const double log_term_base = a * logx - lgamma_cached(a);

    for (n = 1; n < MAXITER1; n++) {
        fac *= -x / n;
        term = fac / (a + n);
        
        // Kahan summation
        y = term - compensation;
        t = sum + y;
        compensation = (t - sum) - y;
        sum = t;
        
        if (__builtin_expect(fabs(term) <= DBL_EPSILON * fabs(sum), 0)) {
            break;
        }
    }

    term = -expm1(a * logx - lgamma1p(a));
    return term - exp(log_term_base) * sum;
}

/**
 *  
 */
double asymptotic_series(double a, double x, int func) {
    int k, n, sgn;
    int maxpow = 0;
    double lambda = x / a;
    double sigma = (x - a) / a;
    double eta, res, ck, ckterm, term, absterm;
    double absoldterm = R_PosInf;
    double etapow[NIC] = {1};
    double sum = 0;
    double afac = 1;
    double compensation = 0.0, y, t;
    
    sgn = (func == IGAM) ? -1 : 1;
    
    const double log_sigma = log1pmx(sigma);
    eta = (lambda > 1) ? sqrt(-2 * log_sigma) : 
          (lambda < 1) ? -sqrt(-2 * log_sigma) : 0;
    
    const double sqrt_a_half = sqrt(0.5 * a);
    res = 0.5 * erfc(sgn * eta * sqrt_a_half);
    
    for (k = 0; k < KIC; k++) {
        ck = d[k][0];
        
        for (n = 1; n < NIC; n++) {
            if (n > maxpow) {
                etapow[n] = eta * etapow[n-1];
                maxpow += 1;
            }
            
            ckterm = d[k][n] * etapow[n];
            ck += ckterm;
            
            if (__builtin_expect(fabs(ckterm) < DBL_EPSILON * fabs(ck), 0)) {
                break;
            }
        }
        
        term = ck * afac;
        absterm = fabs(term);
        
        if (__builtin_expect(absterm > absoldterm, 0)) {
            break;
        }
        
        // Kahan summation
        y = term - compensation;
        t = sum + y;
        compensation = (t - sum) - y;
        sum = t;
        
        if (__builtin_expect(absterm < DBL_EPSILON * fabs(sum), 0)) {
            break;
        }
        
        absoldterm = absterm;
        afac /= a;
    }
    
    const double exp_factor = exp(-0.5 * a * eta * eta);
    const double sqrt_factor = INV_SQRT_2PI / sqrt(a);
    
    res += sgn * exp_factor * sum * sqrt_factor;
    
    return res;
}

/**
 *  
 */
double igam(double a, double x) {
    double absxma_a;
    double result;
    
    result = handle_edge_cases(a, x, 1);
    if (__builtin_expect(result >= 0.0, 0)) {
        return result;
    }

    absxma_a = fabs(x - a) / a;
    
    if (__builtin_expect((x <= 1.0) || (x <= a), 1)) {
        return igam_series(a, x);
    }
    
    if (((a > SMALL) && (a < LARGE) && (absxma_a < SMALLRATIO)) ||
        ((a > LARGE) && (absxma_a < LARGERATIO / sqrt(a)))) {
        return asymptotic_series(a, x, IGAM);
    }

    return 1.0 - igamc(a, x);
}

/**
 *  
 */
double igamc(double a, double x) {
    double absxma_a;
    double result;
    
    result = handle_edge_cases(a, x, 0);
    if (__builtin_expect(result >= 0.0, 0)) {
        return result;
    }

    absxma_a = fabs(x - a) / a;
    
    if (((a > SMALL) && (a < LARGE) && (absxma_a < SMALLRATIO)) ||
        ((a > LARGE) && (absxma_a < LARGERATIO / sqrt(a)))) {
        return asymptotic_series(a, x, IGAMC);
    }

    if (__builtin_expect(x > 1.1, 1)) {
        return (x < a) ? 1.0 - igam_series(a, x) : 
                        igamc_continued_fraction(a, x);
    } else if (x <= 0.5) {
        return (-0.4 / log(x) < a) ? 1.0 - igam_series(a, x) : 
                                    igamc_series(a, x);
    } else {
        return (x * 1.1 < a) ? 1.0 - igam_series(a, x) : 
                              igamc_series(a, x);
    }
}

/**
 *  
 */
void igam_call(double *a, double *x, double *res) {
    __builtin_prefetch(a, 0, 3);
    __builtin_prefetch(x, 0, 3);
    __builtin_prefetch(res, 1, 3);
    
    *res = igam(*a, *x);
}

/**
 *  
 */
void igam_batch(double *a, double *x, double *res, int n) {
    int i;
    
    lgamma_cache.valid = 0;
    
    for (i = 0; i < n; i++) {
        res[i] = igam(a[i], x[i]);
    }
}

//************************************* END igam.c*****************************************
