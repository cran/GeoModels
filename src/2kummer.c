//VERSIONE OTTIMIZZATA COMPLETA - SOLO VELOCITÀ E CACHE

#include "header.h"

// Cache ottimizzata
#define CACHE_SIZE 4096
#define CACHE_MASK 4095

typedef struct {
    double a, b, c, result;
    uint64_t hash;
    uint32_t timestamp;
    uint16_t valid;
    uint16_t hits;
} KummerCacheEntry;

static KummerCacheEntry cache[CACHE_SIZE];
static int cache_initialized = 0;
static uint32_t global_timestamp = 0;
static uint32_t cache_hits = 0, cache_misses = 0;

// Hash function veloce
static inline uint64_t fast_hash_v2(double a, double b, double c) {
    uint64_t ia = *(uint64_t*)&a;
    uint64_t ib = *(uint64_t*)&b; 
    uint64_t ic = *(uint64_t*)&c;
    
    ia &= 0x7FFFFFFFFFFFFFFFULL;
    ic &= 0x7FFFFFFFFFFFFFFFULL;
    
    uint64_t hash = ia * 0x9E3779B97F4A7C15ULL;
    hash ^= ib + 0x9E3779B97F4A7C15ULL + (hash << 6) + (hash >> 2);
    hash ^= ic + 0x9E3779B97F4A7C15ULL + (hash << 6) + (hash >> 2);
    
    return hash;
}

// Inizializzazione cache veloce
static inline void init_cache_optimized() {
    if (!cache_initialized) {
        memset(cache, 0, sizeof(cache));
        cache_initialized = 1;
        global_timestamp = 1;
        cache_hits = cache_misses = 0;
    }
}

// Cache lookup ottimizzato
static inline int cache_lookup_fast(double a, double b, double c, double *result) {
    init_cache_optimized();
    uint64_t hash = fast_hash_v2(a, b, c);
    int idx = hash & CACHE_MASK;
    
    KummerCacheEntry *entry = &cache[idx];
    
    if (entry->valid && entry->hash == hash) {
        if (entry->a == a && entry->b == b && entry->c == c) {
            *result = entry->result;
            entry->timestamp = ++global_timestamp;
            entry->hits++;
            cache_hits++;
            return 1;
        }
    }
    
    cache_misses++;
    return 0;
}

// Cache store con eviction intelligente
static inline void cache_store_smart(double a, double b, double c, double result) {
    uint64_t hash = fast_hash_v2(a, b, c);
    int idx = hash & CACHE_MASK;
    
    KummerCacheEntry *entry = &cache[idx];
    
    if (entry->valid) {
        uint32_t age = global_timestamp - entry->timestamp;
        if (entry->hits > 3 && age < 1000) {
            for (int i = 1; i < 8; i++) {
                int alt_idx = (idx + i) & CACHE_MASK;
                KummerCacheEntry *alt = &cache[alt_idx];
                if (!alt->valid || (global_timestamp - alt->timestamp) > age) {
                    entry = alt;
                    break;
                }
            }
        }
    }
    
    entry->a = a;
    entry->b = b;
    entry->c = c;
    entry->result = result;
    entry->hash = hash;
    entry->timestamp = ++global_timestamp;
    entry->valid = 1;
    entry->hits = 0;
}

// Parametri ottimizzati
typedef struct {
    double a, b, x, c;
    double a_minus_1, b_minus_a_minus_1, neg_x;
    double inv_gamma_a;
} OptimizedParams;

static inline void setup_params(OptimizedParams *p, double a, double b, double x, double c) {
    p->a = a;
    p->b = b;
    p->x = x;
    p->c = c;
    p->a_minus_1 = a - 1.0;
    p->b_minus_a_minus_1 = b - a - 1.0;
    p->neg_x = -x;
    p->inv_gamma_a = exp(-lgammafn(a));
}

// Convergenza veloce
static inline int fast_converged(double current, double previous, double threshold) {
    double diff = fabs(current - previous);
    if (diff <= 1e-300) return 1;
    if (current == 0.0) return 0;
    return (diff <= threshold * fabs(current));
}

// Potenza ottimizzata
static inline double optimized_pow(double base, double exp) {
    if (exp == 0.0) return 1.0;
    if (exp == 1.0) return base;
    if (exp == 2.0) return base * base;
    if (exp == 0.5) return sqrt(base);
    if (exp == -1.0) return 1.0 / base;
    if (exp == -0.5) return 1.0 / sqrt(base);
    return R_pow(base, exp);
}

// Integrand 1 ottimizzato
void integrand1_opt(double *t, int n, void *params) {
    OptimizedParams *p = (OptimizedParams*)params;
    const double a_m1 = p->a_minus_1;
    const double b_ma_m1 = p->b_minus_a_minus_1;
    const double neg_x = p->neg_x;
    
    for(int i = 0; i < n; i++) {
        const double ti = t[i];
        const double exp_term = exp(neg_x * ti);
        const double pow_term1 = optimized_pow(ti, a_m1);
        const double pow_term2 = optimized_pow(1.0 + ti, b_ma_m1);
        t[i] = exp_term * pow_term1 * pow_term2;
    }
}

// Integrand 2 ottimizzato
void integrand2_opt(double *u, int n, void *params) {
    OptimizedParams *p = (OptimizedParams*)params;
    const double a_m1 = p->a_minus_1;
    const double b_ma_m1 = p->b_minus_a_minus_1;
    const double neg_x = p->neg_x;
    const double c = p->c;
    
    for(int i = 0; i < n; i++) {
        const double ui = u[i];
        const double one_minus_u = 1.0 - ui;
        const double t = c / one_minus_u;
        const double jacobian = c / (one_minus_u * one_minus_u);
        const double exp_term = exp(neg_x * t);
        const double pow_term1 = optimized_pow(t, a_m1);
        const double pow_term2 = optimized_pow(1.0 + t, b_ma_m1);
        u[i] = jacobian * exp_term * pow_term1 * pow_term2;
    }
}

// Serie ottimizzata con unrolling
static double optimized_series_v2(double a, double b, double x) {
    double result = 1.0;
    double term = 1.0;
    const double tol = 1e-15;
    
    const double x_val = x;
    double a_k = a;
    double b_k = b;
    double k_double = 1.0;
    
    // Unroll prime iterazioni
    term *= (a_k / (k_double * b_k)) * x_val;
    result += term;
    if (fabs(term) <= tol * fabs(result)) return result;
    
    a_k += 1.0; b_k += 1.0; k_double = 2.0;
    
    term *= (a_k / (k_double * b_k)) * x_val;
    result += term;
    if (fabs(term) <= tol * fabs(result)) return result;
    
    a_k += 1.0; b_k += 1.0; k_double = 3.0;
    
    for (int k = 3; k <= 500; k++) {
        term *= (a_k / (k_double * b_k)) * x_val;
        result += term;
        
        if (fast_converged(result, result - term, tol)) break;
        if (fabs(term) > 1e100) return R_NaN;
        
        a_k += 1.0;
        b_k += 1.0; 
        k_double += 1.0;
    }
    
    return result;
}

// CHGUIT ottimizzato
void chguit_opt(double a, double b, double x, double *hu, int *id) {
    *id = 9;
    *hu = 0.0;
    
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(x) || x <= 0.0) {
        *id = -1;
        return;
    }
    
    double c = fmax(15.0 / sqrt(x), 8.0);
    OptimizedParams params;
    setup_params(&params, a, b, x, c);
    
    if (!R_FINITE(params.inv_gamma_a)) {
        *id = -1;
        return;
    }
    
    double zero = 0.0, one = 1.0;
    double epsabs = 1e-12;
    double epsrel = 1e-12;
    int last = 0;
    int limit = 150;
    int lenw = 600;
    
    int *iwork = (int*)calloc(limit, sizeof(int));
    double *work = (double*)calloc(lenw, sizeof(double));
    
    if (!iwork || !work) {
        free(iwork);
        free(work);
        *id = -1;
        return;
    }
    
    double result1, abserr1;
    int neval1, ier1;
    
    Rdqags(integrand1_opt, &params, &zero, &c, &epsabs, &epsrel, &result1, &abserr1, 
           &neval1, &ier1, &limit, &lenw, &last, iwork, work);
    
    if (ier1 > 1) {
        free(iwork);
        free(work);
        *id = -1;
        return;
    }
    
    result1 *= params.inv_gamma_a;
    
    double result2, abserr2;
    int neval2, ier2;
    
    Rdqags(integrand2_opt, &params, &zero, &one, &epsabs, &epsrel, &result2, &abserr2,
           &neval2, &ier2, &limit, &lenw, &last, iwork, work);
    
    if (ier2 > 1) {
        free(iwork);
        free(work);
        *id = -1;
        return;
    }
    
    result2 *= params.inv_gamma_a;
    *hu = result1 + result2;
    
    double total_error = abserr1 + abserr2;
    if (*hu != 0.0 && total_error > 0.0) {
        *id = (int)fmax(1.0, -log10(total_error / fabs(*hu)));
    }
    
    free(iwork);
    free(work);
}

// CHGM ottimizzato
void chgm_opt(double a, double b, double x, double *hg) {
    const double pi_val = M_PI;
    double a0 = a, a1 = a, x0 = x;
    *hg = 0.0;
    
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(x)) {
        *hg = R_NaN;
        return;
    }
    
    if (x < 0.0) {
        a = b - a;
        a0 = a;
        x = fabs(x);
    }
    
    int nl = 0, la = 0;
    if (a >= 2.0) {
        nl = 1;
        la = (int)a;
        a = a - la - 1.0;
    }
    
    double y0 = 0.0, y1 = 0.0;
    
    for (int n = 0; n <= nl; n++) {
        if (a0 >= 2.0) a = a + 1.0;
        
        if (x <= 35.0 + fabs(b) || a < 0.0) {
            *hg = optimized_series_v2(a, b, x);
            if (x0 < 0.0) *hg *= exp(x0);
        } else {
            double sum1 = 1.0, sum2 = 1.0;
            double r1 = 1.0, r2 = 1.0;
            const double inv_x = 1.0 / x;
            
            for (int i = 1; i <= 10; i++) {
                r1 *= -(a + i - 1.0) * (a - b + i) * inv_x / i;
                r2 *= -(b - a + i - 1.0) * (a - i) * inv_x / i;
                sum1 += r1;
                sum2 += r2;
                
                if (fabs(r1) + fabs(r2) < 1e-15 * (fabs(sum1) + fabs(sum2))) break;
            }
            
            double hg1, hg2;
            if (x0 >= 0.0) {
                hg1 = exp(lgammafn(b) - lgammafn(b-a)) * optimized_pow(x, -a) * cos(pi_val * a) * sum1;
                hg2 = exp(lgammafn(b) - lgammafn(a) + x) * optimized_pow(x, a - b) * sum2;
            } else {
                hg1 = exp(lgammafn(b) - lgammafn(b-a) + x0) * optimized_pow(x, -a) * cos(pi_val * a) * sum1;
                hg2 = exp(lgammafn(b) - lgammafn(a)) * optimized_pow(x, a - b) * sum2;
            }
            *hg = hg1 + hg2;
        }
        
        if (n == 0) y0 = *hg;
        if (n == 1) y1 = *hg;
    }
    
    if (a0 >= 2.0) {
        for (int i = 1; i <= la - 1; i++) {
            double temp = ((2.0 * a - b + x) * y1 + (b - a) * y0) / a;
            y0 = y1;
            y1 = temp;
            a += 1.0;
        }
        *hg = y1;
    }
    
    a = a1;
    x = x0;
}

// CHGUS ottimizzato
void chgus_opt(double a, double b, double x, double *hu, int *id) {
    *id = -100;
    *hu = 0.0;
    
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(x) || x <= 0.0) {
        *id = -1;
        return;
    }
    
    const double pi_val = M_PI;
    
    double ga = gammafn(a);
    double gb = gammafn(b);
    double gab = gammafn(1.0 + a - b);
    double gb2 = gammafn(2.0 - b);
    
    if (!R_FINITE(ga) || !R_FINITE(gb) || !R_FINITE(gab) || !R_FINITE(gb2)) {
        *id = -1;
        return;
    }
    
    double sin_pi_b = sin(pi_val * b);
    if (fabs(sin_pi_b) < 1e-15) {
        *id = -1;
        return;
    }
    
    double hu0 = pi_val / sin_pi_b;
    double r1 = hu0 / (gab * gb);
    double r2 = hu0 * optimized_pow(x, 1.0 - b) / (ga * gb2);
    *hu = r1 - r2;
    
    double hmax = fabs(*hu), hmin = fabs(*hu);
    
    double a_k = a;
    double b_k = b;
    double ab_k = a - b + 1.0;
    double b1_k = 2.0 - b;
    
    for (int j = 1; j <= 150; j++) {
        r1 *= (a_k / (j * b_k)) * x;
        r2 *= (ab_k / (j * b1_k)) * x;
        
        double prev_hu = *hu;
        *hu += r1 - r2;
        
        double hua = fabs(*hu);
        if (hua > hmax) hmax = hua;
        if (hua < hmin && hua > 0.0) hmin = hua;
        
        if (fast_converged(*hu, prev_hu, 1e-15)) break;
        
        a_k += 1.0;
        b_k += 1.0;
        ab_k += 1.0;
        b1_k += 1.0;
    }
    
    if (hmin > 0.0 && hmax > 0.0) {
        double d1 = log10(hmax);
        double d2 = log10(hmin);
        *id = (int)(15 - fabs(d1 - d2));
    } else {
        *id = 10;
    }
}

// CHGUL ottimizzato
void chgul_opt(double a, double b, double x, double *hu, int *id) {
    *id = -100;
    *hu = 0.0;
    
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(x) || x <= 0.0) {
        *id = -1;
        return;
    }
    
    double aa = a - b + 1.0;
    int il1 = (a == floor(a) && a <= 0.0);
    int il2 = (aa == floor(aa) && aa <= 0.0);
    
    int nm = 0;
    if (il1) nm = (int)fabs(a);
    if (il2) nm = (int)fabs(aa);
    
    if (il1 || il2) {
        *hu = 1.0;
        double r = 1.0;
        const double inv_x = 1.0 / x;
        double a_k = a;
        double ab_k = a - b + 1.0;
        
        for (int k = 1; k <= nm; k++) {
            r *= -a_k * ab_k * inv_x / k;
            *hu += r;
            a_k += 1.0;
            ab_k += 1.0;
        }
        *hu = optimized_pow(x, -a) * (*hu);
        *id = 10;
    } else {
        *hu = 1.0;
        double r = 1.0, r0 = 0.0;
        const double inv_x = 1.0 / x;
        double a_k = a;
        double ab_k = a - b + 1.0;
        
        for (int k = 1; k <= 40; k++) {
            r *= -a_k * ab_k * inv_x / k;
            double ra = fabs(r);
            
            if (k > 5 && (ra >= r0 || ra < 1e-15 * fabs(*hu))) break;
            
            r0 = ra;
            *hu += r;
            a_k += 1.0;
            ab_k += 1.0;
        }
        *id = (r0 > 0.0) ? (int)fmax(1.0, -log10(r0)) : 15;
        *hu = optimized_pow(x, -a) * (*hu);
    }
}

// CHGUBI ottimizzato
void chgubi_opt(double a, double b, double x, double *hu, int *id) {
    *id = -100;
    *hu = 0.0;
    
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(x) || x <= 0.0) {
        *id = -1;
        return;
    }
    
    const double EL = 0.5772156649015329;
    
    int n = (int)fabs(b - 1);
    double rn1 = 1.0, rn = 1.0;
    for (int j = 1; j <= n; j++) {
        rn *= j;
        if (j == n - 1) rn1 = rn;
    }
    
    double ps = digamma(a);
    double ga = gammafn(a);
    
    if (!R_FINITE(ps) || !R_FINITE(ga) || ga == 0.0) {
        *id = -1;
        return;
    }
    
    double a0, a1, a2, ua, ub, ga1;
    
    if (b > 0.0) {
        a0 = a;
        a1 = a - n;
        a2 = a1;
        ga1 = gammafn(a1);
        if (!R_FINITE(ga1) || ga1 == 0.0) {
            *id = -1;
            return;
        }
        ua = optimized_pow(-1, n - 1) / (rn * ga1);
        ub = rn1 / ga * optimized_pow(x, -n);
    } else {
        a0 = a + n;
        a1 = a0;
        a2 = a;
        ga1 = gammafn(a1);
        if (!R_FINITE(ga1) || ga1 == 0.0) {
            *id = -1;
            return;
        }
        ua = optimized_pow(-1, n - 1) / (rn * ga) * optimized_pow(x, n);
        ub = rn1 / ga1;
    }
    
    double hm1 = 1.0;
    double r = 1.0;
    double hmax = 1.0, hmin = 1.0;
    
    for (int k = 1; k <= 150; k++) {
        r *= (a0 + k - 1.0) * x / ((n + k) * k);
        double prev_hm1 = hm1;
        hm1 += r;
        
        double hu1 = fabs(hm1);
        if (hu1 > hmax) hmax = hu1;
        if (hu1 < hmin && hu1 > 0.0) hmin = hu1;
        
        if (fast_converged(hm1, prev_hm1, 1e-15)) break;
    }
    
    *id = (hmin > 0.0 && hmax > 0.0) ? (int)(15 - fabs(log10(hmax) - log10(hmin))) : 10;
    
    hm1 *= log(x);
    
    double s0 = 0.0;
    for (int m = 1; m <= n; m++) {
        if (b >= 0.0) s0 -= 1.0 / m;
        if (b < 0.0) s0 += (1.0 - a) / (m * (a + m - 1.0));
    }
    
    double hm2 = ps + 2.0 * EL + s0;
    r = 1.0;
    
    for (int k = 1; k <= 150; k++) {
        double s1 = 0.0, s2 = 0.0;
        
        if (b > 0.0) {
            for (int m = 1; m <= k; m++) {
                s1 -= (m + 2.0 * a - 2.0) / (m * (m + a - 1.0));
            }
            for (int m = 1; m <= n; m++) {
                s2 += 1.0 / (k + m);
            }
        } else {
            for (int m = 1; m <= k + n; m++) {
                s1 += (1.0 - a) / (m * (m + a - 1.0));
            }
            for (int m = 1; m <= k; m++) {
                s2 += 1.0 / m;
            }
        }
        
        double hw = 2.0 * EL + ps + s1 - s2;
        r *= (a0 + k - 1.0) * x / ((n + k) * k);
        double prev_hm2 = hm2;
        hm2 += r * hw;
        
        if (fast_converged(hm2, prev_hm2, 1e-15)) break;
    }
    
    double hm3 = (n == 0) ? 0.0 : 1.0;
    r = 1.0;
    for (int k = 1; k <= n - 1; k++) {
        r *= (a2 + k - 1.0) * x / ((k - n) * k);
        hm3 += r;
    }
    
    double sa = ua * (hm1 + hm2);
    double sb = ub * hm3;
    *hu = sa + sb;
    
    if (sa * sb < 0.0 && sa != 0.0 && *hu != 0.0) {
        int id1 = (int)log10(fabs(sa));
        int id2 = (int)log10(fabs(*hu));
        *id = fmax(1, *id - abs(id1 - id2));
    }
}

// CHGU ottimizzato
void chgu_opt(double a, double b, double x, double *hu, int *md, int *isfer) {
    *isfer = 0;
    *hu = 0.0;
    *md = 0;
    
    double cached_result;
    if (cache_lookup_fast(a, b, x, &cached_result)) {
        *hu = cached_result;
        *md = 0;
        return;
    }
    
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(x) || x <= 0.0) {
        *isfer = -1;
        *hu = R_NaN;
        return;
    }
    
    double aa = a - b + 1.0;
    
    int il1 = (a == floor(a) && a <= 0.0);
    int il2 = (aa == floor(aa) && aa <= 0.0);
    int il3 = (fabs(a * aa) / x <= 2.0);
    int bn = (b == floor(b) && b != 0.0);
    
    double hu_best = 0.0;
    int id_best = -100;
    int md_best = 0;
    
    if (b != floor(b) && x <= 50.0) {
        int id;
        double hu_temp;
        chgus_opt(a, b, x, &hu_temp, &id);
        if (id > id_best) {
            id_best = id;
            hu_best = hu_temp;
            md_best = 1;
        }
    }
    
    if (il1 || il2 || il3) {
        int id;
        double hu_temp;
        chgul_opt(a, b, x, &hu_temp, &id);
        if (id > id_best) {
            id_best = id;
            hu_best = hu_temp;
            md_best = 2;
        }
    }
    
    if (a >= 1.0 && bn && x <= 30.0) {
        int id;
        double hu_temp;
        chgubi_opt(a, b, x, &hu_temp, &id);
        if (id > id_best) {
            id_best = id;
            hu_best = hu_temp;
            md_best = 3;
        }
    }
    
    if (a >= 1.0 && id_best < 5) {
        int id;
        double hu_temp;
        chguit_opt(a, b, x, &hu_temp, &id);
        if (id > id_best) {
            id_best = id;
            hu_best = hu_temp;
            md_best = 4;
        }
    }
    
    if (a < 1.0 && id_best < 5) {
        if (b <= a) {
            double a00 = a, b00 = b;
            a = a - b + 1.0;
            b = 2.0 - b;
            int id;
            double hu_temp;
            chguit_opt(a, b, x, &hu_temp, &id);
            hu_temp = optimized_pow(x, 1.0 - b00) * hu_temp;
            if (id > id_best) {
                id_best = id;
                hu_best = hu_temp;
                md_best = 4;
            }
            a = a00;
            b = b00;
        } else if (bn && !il1) {
            int id;
            double hu_temp;
            chgubi_opt(a, b, x, &hu_temp, &id);
            if (id > id_best) {
                id_best = id;
                hu_best = hu_temp;
                md_best = 3;
            }
        }
    }
    
    *hu = hu_best;
    *md = md_best;
    
    if (id_best <= 0 || !R_FINITE(hu_best)) {
        *isfer = 1;
        if (!R_FINITE(hu_best)) *hu = R_NaN;
    } else {
        cache_store_smart(a, b, x, *hu);
    }
}

// Funzione principale kummer ottimizzata
double kummer(double a, double b, double c) {
    // Casi immediati
    if (a == 0.0) return 1.0;
    if (c == 0.0) return R_NaN;
    if (a == 1.0 && b == 1.0) return exp(c);
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(c) || c <= 0.0) return R_NaN;
    
    // Check cache
    double cached_result;
    if (cache_lookup_fast(a, b, c, &cached_result)) {
        return cached_result;
    }
    
    // Casi molto grandi
    if (c > 1e50) {
        double result = optimized_pow(c, -a);
        cache_store_smart(a, b, c, result);
        return result;
    }
    
    // Computation principale
    double result;
    int method, status;
    chgu_opt(a, b, c, &result, &method, &status);
    
    if (status == -1 || !R_FINITE(result)) {
        return R_NaN;
    }
    
    // Store in cache
    cache_store_smart(a, b, c, result);
    return result;
}
