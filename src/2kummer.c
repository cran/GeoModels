#include "header.h"
// Optimized Kummer function implementation with improved performance and robustness


// Structure for integration parameters
typedef struct {
    double a, b, x, c;
} IntegrandParams;


// Optimized convergence check with relative and absolute tolerances
static inline int converged(double current, double previous, double tol_rel, double tol_abs) {
    if (current == 0.0) return (fabs(previous) <= tol_abs);
    double rel_err = fabs((current - previous) / current);
    return (rel_err <= tol_rel) || (fabs(current - previous) <= tol_abs);
}



// Integrand function for first integral (0 to C) - optimized
void integrand1_opt(double *t, int n, void *params) {
    IntegrandParams *p = (IntegrandParams*)params;
    const double a_minus_1 = p->a - 1.0;
    const double b_minus_a_minus_1 = p->b - p->a - 1.0;
    const double neg_x = -p->x;
    
    for(int i = 0; i < n; i++) {
        const double ti = t[i];
        t[i] = exp(neg_x * ti) * R_pow(ti, a_minus_1) * R_pow(1.0 + ti, b_minus_a_minus_1);
    }
}

// Integrand function for second integral (C to infinity) - optimized
void integrand2_opt(double *u, int n, void *params) {
    IntegrandParams *p = (IntegrandParams*)params;
    const double a_minus_1 = p->a - 1.0;
    const double b_minus_a_minus_1 = p->b - p->a - 1.0;
    const double neg_x = -p->x;
    const double c = p->c;
    
    for(int i = 0; i < n; i++) {
        const double ui = u[i];
        const double one_minus_u = 1.0 - ui;
        const double t = c / one_minus_u;
        const double jacobian = c / (one_minus_u * one_minus_u);
        u[i] = jacobian * exp(neg_x * t) * R_pow(t, a_minus_1) * R_pow(1.0 + t, b_minus_a_minus_1);
    }
}

// Optimized CHGUIT with better memory management and error handling
void chguit_opt(double a, double b, double x, double *hu, int *id) {
    *id = 9;
    *hu = 0.0;
    
    // Input validation
    if (!isfinite(a) || !isfinite(b) || !isfinite(x) || x <= 0.0) {
        *id = -1;
        return;
    }
    
    // Optimize c selection based on x
    double c = fmax(12.0 / x, 5.0);  // Ensure minimum c for stability
    IntegrandParams params = {a, b, x, c};
    
    // Pre-compute gamma factor
    double log_gamma_a = lgammafn(a);
    if (!isfinite(log_gamma_a)) {
        *id = -1;
        return;
    }
    double inv_gamma_a = exp(-log_gamma_a);
    
    // Fixed variables for address-taking
     double zero = 0.0, one = 1.0;
    double epsabs = 1e-12;  // Increased precision
     double epsrel = 1e-12;
    int last = 0;
    
    // Memory allocation with error checking
    int limit = 200;  // Increased limit for better accuracy
    int lenw = 800;
    int *iwork = (int*)calloc(limit, sizeof(int));
    double *work = (double*)calloc(lenw, sizeof(double));
    
    if (!iwork || !work) {
        free(iwork);
        free(work);
        *id = -1;
        return;
    }
    
    // First integral: 0 to C
    double result1, abserr1;
    int neval1, ier1;
    
    Rdqags(integrand1_opt, &params, &zero, &c, &epsabs, &epsrel, &result1, &abserr1, 
           &neval1, &ier1, &limit, &lenw, &last, iwork, work);
    
    if (ier1 > 1) {  // Integration failed
        free(iwork);
        free(work);
        *id = -1;
        return;
    }
    
    result1 *= inv_gamma_a;
    
    // Second integral: C to infinity (transformed to 0 to 1)
    double result2, abserr2;
    int neval2, ier2;
    
    Rdqags(integrand2_opt, &params, &zero, &one, &epsabs, &epsrel, &result2, &abserr2,
           &neval2, &ier2, &limit, &lenw, &last, iwork, work);
    
    if (ier2 > 1) {  // Integration failed
        free(iwork);
        free(work);
        *id = -1;
        return;
    }
    
    result2 *= inv_gamma_a;
    
    *hu = result1 + result2;
    
    // Estimate precision based on integration errors
    double total_error = abserr1 + abserr2;
    if (*hu != 0.0 && total_error > 0.0) {
        *id = (int)fmax(1.0, -log10(total_error / fabs(*hu)));
    }
    
    free(iwork);
    free(work);
}

// Optimized CHGM with better convergence detection
void chgm_opt(double a, double b, double x, double *hg) {
    const double pi_val = M_PI;
    double a0 = a, a1 = a, x0 = x;
    *hg = 0.0;
    
    // Input validation
    if (!isfinite(a) || !isfinite(b) || !isfinite(x)) {
        *hg = R_NaN;
        return;
    }
    
    // DLMF 13.2.39
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
    const double tol_rel = 1e-15;
    const double tol_abs = 1e-300;
    
    for (int n = 0; n <= nl; n++) {
        if (a0 >= 2.0) a = a + 1.0;
        
        if (x <= 30.0 + fabs(b) || a < 0.0) {
            // Series expansion
            *hg = 1.0;
            double rg = 1.0;
            double prev_hg = 0.0;
            
            // Pre-compute constants to avoid repeated calculations
            const double x_factor = x;
            double a_term = a;
            double b_term = b;
            
            for (int j = 1; j <= 1000; j++) {  // Increased max iterations
                rg *= (a_term / (j * b_term)) * x_factor;
                prev_hg = *hg;
                *hg += rg;
                
                // Optimized convergence check
                if (converged(*hg, prev_hg, tol_rel, tol_abs) || fabs(rg) < tol_abs) {
                    break;
                }
                
                a_term += 1.0;
                b_term += 1.0;
            }
            
            if (x0 < 0.0) *hg *= exp(x0);
        } else {
            // Asymptotic expansion for large x
            double sum1 = 1.0, sum2 = 1.0;
            double r1 = 1.0, r2 = 1.0;
            const double inv_x = 1.0 / x;
            
            // Use more terms for better accuracy
            for (int i = 1; i <= 12; i++) {
                r1 *= -(a + i - 1.0) * (a - b + i) * inv_x / i;
                r2 *= -(b - a + i - 1.0) * (a - i) * inv_x / i;
                sum1 += r1;
                sum2 += r2;
                
                // Early termination if terms become negligible
                if (fabs(r1) + fabs(r2) < tol_abs * (fabs(sum1) + fabs(sum2))) break;
            }
            
            double hg1, hg2;
            if (x0 >= 0.0) {
                hg1 = exp(lgammafn(b) - lgammafn(b-a)) * R_pow(x, -a) * cos(pi_val * a) * sum1;
                hg2 = exp(lgammafn(b) - lgammafn(a) + x) * R_pow(x, a - b) * sum2;
            } else {
                hg1 = exp(lgammafn(b) - lgammafn(b-a) + x0) * R_pow(x, -a) * cos(pi_val * a) * sum1;
                hg2 = exp(lgammafn(b) - lgammafn(a)) * R_pow(x, a - b) * sum2;
            }
            *hg = hg1 + hg2;
        }
        
        if (n == 0) y0 = *hg;
        if (n == 1) y1 = *hg;
    }
    
    // Recurrence relation for a0 >= 2
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

// Optimized CHGUS with better numerical stability
void chgus_opt(double a, double b, double x, double *hu, int *id) {
    *id = -100;
    *hu = 0.0;
    
    // Input validation
    if (!isfinite(a) || !isfinite(b) || !isfinite(x) || x <= 0.0) {
        *id = -1;
        return;
    }
    
    const double pi_val = M_PI;
    const double tol_rel = 1e-15;
    const double tol_abs = 1e-300;
    
    // Pre-compute gamma functions with error checking
    double ga = gammafn(a);
    double gb = gammafn(b);
    double gab = gammafn(1.0 + a - b);
    double gb2 = gammafn(2.0 - b);
    
    if (!isfinite(ga) || !isfinite(gb) || !isfinite(gab) || !isfinite(gb2)) {
        *id = -1;
        return;
    }
    
    double sin_pi_b = sin(pi_val * b);
    if (fabs(sin_pi_b) < 1e-15) {  // Near singularity
        *id = -1;
        return;
    }
    
    double hu0 = pi_val / sin_pi_b;
    double r1 = hu0 / (gab * gb);
    double r2 = hu0 * R_pow(x, 1.0 - b) / (ga * gb2);
    *hu = r1 - r2;
    
    double hmax = fabs(*hu), hmin = fabs(*hu);
    double prev_hu = 0.0;
    
    // Pre-compute x factor
    const double x_factor = x;
    double a_term = a;
    double b_term = b;
    double ab_term = a - b + 1.0;
    double b1_term = 1.0 - b + 1.0;
    
    for (int j = 1; j <= 200; j++) {  // Increased iterations
        r1 *= (a_term / (j * b_term)) * x_factor;
        r2 *= (ab_term / (j * b1_term)) * x_factor;
        
        prev_hu = *hu;
        *hu += r1 - r2;
        
        double hua = fabs(*hu);
        if (hua > hmax) hmax = hua;
        if (hua < hmin && hua > 0.0) hmin = hua;
        
        // Optimized convergence check
        if (converged(*hu, prev_hu, tol_rel, tol_abs)) break;
        
        a_term += 1.0;
        b_term += 1.0;
        ab_term += 1.0;
        b1_term += 1.0;
    }
    
    // Precision estimation
    if (hmin > 0.0 && hmax > 0.0) {
        double d1 = log10(hmax);
        double d2 = log10(hmin);
        *id = (int)(15 - fabs(d1 - d2));
    } else {
        *id = 10;
    }
}

// Optimized CHGUL with early termination
void chgul_opt(double a, double b, double x, double *hu, int *id) {
    *id = -100;
    *hu = 0.0;
    
    // Input validation
    if (!isfinite(a) || !isfinite(b) || !isfinite(x) || x <= 0.0) {
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
        double a_term = a;
        double ab_term = a - b + 1.0;
        
        for (int k = 1; k <= nm; k++) {
            r *= -a_term * ab_term * inv_x / k;
            *hu += r;
            a_term += 1.0;
            ab_term += 1.0;
        }
        *hu = R_pow(x, -a) * (*hu);
        *id = 10;
    } else {
        *hu = 1.0;
        double r = 1.0, r0 = 0.0;
        const double inv_x = 1.0 / x;
        const double tol = 1e-15;
        double a_term = a;
        double ab_term = a - b + 1.0;
        
        for (int k = 1; k <= 50; k++) {  // Increased iterations
            r *= -a_term * ab_term * inv_x / k;
            double ra = fabs(r);
            
            // Improved convergence criteria
            if (k > 5 && (ra >= r0 || ra < tol * fabs(*hu))) break;
            
            r0 = ra;
            *hu += r;
            a_term += 1.0;
            ab_term += 1.0;
        }
        *id = (r0 > 0.0) ? (int)fmax(1.0, -log10(r0)) : 15;
        *hu = R_pow(x, -a) * (*hu);
    }
}

// Optimized CHGUBI with better memory usage
void chgubi_opt(double a, double b, double x, double *hu, int *id) {
    *id = -100;
    *hu = 0.0;
    
    // Input validation
    if (!isfinite(a) || !isfinite(b) || !isfinite(x) || x <= 0.0) {
        *id = -1;
        return;
    }
    
    const double EL = 0.5772156649015329; // Euler's constant
    const double tol_rel = 1e-15;
    const double tol_abs = 1e-300;
    
    int n = (int)fabs(b - 1);
    double rn1 = 1.0, rn = 1.0;
    for (int j = 1; j <= n; j++) {
        rn *= j;
        if (j == n - 1) rn1 = rn;
    }
    
    double ps = digamma(a);
    double ga = gammafn(a);
    
    if (!isfinite(ps) || !isfinite(ga) || ga == 0.0) {
        *id = -1;
        return;
    }
    
    double a0, a1, a2, ua, ub, ga1;
    
    if (b > 0.0) {
        a0 = a;
        a1 = a - n;
        a2 = a1;
        ga1 = gammafn(a1);
        if (!isfinite(ga1) || ga1 == 0.0) {
            *id = -1;
            return;
        }
        ua = R_pow(-1, n - 1) / (rn * ga1);
        ub = rn1 / ga * R_pow(x, -n);
    } else {
        a0 = a + n;
        a1 = a0;
        a2 = a;
        ga1 = gammafn(a1);
        if (!isfinite(ga1) || ga1 == 0.0) {
            *id = -1;
            return;
        }
        ua = R_pow(-1, n - 1) / (rn * ga) * R_pow(x, n);
        ub = rn1 / ga1;
    }
    
    // First sum - optimized
    double hm1 = 1.0, r = 1.0;
    double hmax = 1.0, hmin = 1.0, prev_hm1 = 0.0;
    const double x_factor = x;
    const double inv_denom = 1.0 / (n + 1);
    
    for (int k = 1; k <= 200; k++) {
        r *= (a0 + k - 1.0) * x_factor * inv_denom / k;
        prev_hm1 = hm1;
        hm1 += r;
        
        double hu1 = fabs(hm1);
        if (hu1 > hmax) hmax = hu1;
        if (hu1 < hmin && hu1 > 0.0) hmin = hu1;
        
        if (converged(hm1, prev_hm1, tol_rel, tol_abs)) break;
    }
    
    // Precision estimation
    double da1 = (hmax > 0.0) ? log10(hmax) : 0.0;
    double da2 = (hmin > 0.0) ? log10(hmin) : 0.0;
    *id = (int)(15 - fabs(da1 - da2));
    
    hm1 *= log(x);
    
    // Second sum with precomputed constants
    double s0 = 0.0;
    for (int m = 1; m <= n; m++) {
        if (b >= 0.0) s0 -= 1.0 / m;
        if (b < 0.0) s0 += (1.0 - a) / (m * (a + m - 1.0));
    }
    
    double hm2 = ps + 2.0 * EL + s0;
    r = 1.0;
    double prev_hm2 = hm2;
    
    for (int k = 1; k <= 200; k++) {
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
        r *= (a0 + k - 1.0) * x_factor / ((n + k) * k);
        prev_hm2 = hm2;
        hm2 += r * hw;
        
        if (converged(hm2, prev_hm2, tol_rel, tol_abs)) break;
    }
    
    // Third sum
    double hm3 = (n == 0) ? 0.0 : 1.0;
    r = 1.0;
    for (int k = 1; k <= n - 1; k++) {
        r *= (a2 + k - 1.0) * x_factor / ((k - n) * k);
        hm3 += r;
    }
    
    double sa = ua * (hm1 + hm2);
    double sb = ub * hm3;
    *hu = sa + sb;
    
    // Adjust precision if terms cancel
    if (sa * sb < 0.0 && sa != 0.0 && *hu != 0.0) {
        int id1 = (int)log10(fabs(sa));
        int id2 = (int)log10(fabs(*hu));
        *id = fmax(1, *id - abs(id1 - id2));
    }
}

// Optimized main CHGU function with better algorithm selection
void chgu_opt(double a, double b, double x, double *hu, int *md, int *isfer) {
    *isfer = 0;
    *hu = 0.0;
    *md = 0;
    
    // Comprehensive input validation
    if (!isfinite(a) || !isfinite(b) || !isfinite(x) || x <= 0.0) {
        *isfer = -1;
        *hu = R_NaN;
        return;
    }
    
    double aa = a - b + 1.0;
    
    // Pre-compute conditions
    int il1 = (a == floor(a) && a <= 0.0);
    int il2 = (aa == floor(aa) && aa <= 0.0);
    int il3 = (fabs(a * (a - b + 1.0)) / x <= 2.0);
    int bn = (b == floor(b) && b != 0.0);
    
    // Improved algorithm selection criteria
    int bl1 = (x <= 8.0 || (x <= 15.0 && a <= 3.0));  // Extended range
    int bl2 = ((x > 8.0 && x <= 20.0) && (a >= 1.0 && b >= a + 4.0));
    int bl3 = (x > 20.0 && a >= 5.0 && b >= a + 5.0);
    
    int id_best = -100;
    double hu_best = 0.0;
    int md_best = 0;
    
    // Try multiple methods and select the best result
    
    // Method 1: CHGUS for non-integer b
    if (b != floor(b)) {
        int id;
        double hu_temp;
        chgus_opt(a, b, x, &hu_temp, &id);
        if (id > id_best) {
            id_best = id;
            hu_best = hu_temp;
            md_best = 1;
        }
    }
    // Method 2: CHGUL for special cases or large x
    if (il1 || il2 || il3 || x > 25.0) {
        int id;
        double hu_temp;
        chgul_opt(a, b, x, &hu_temp, &id);
        if (id > id_best) {
            id_best = id;
            hu_best = hu_temp;
            md_best = 2;
        }
    }
    // Method 3: CHGUBI for integer b
    if (a >= 1.0 && bn && (bl1 || bl2 || bl3)) {
        int id;
        double hu_temp;
        chgubi_opt(a, b, x, &hu_temp, &id);
        if (id > id_best) {
            id_best = id;
            hu_best = hu_temp;
            md_best = 3;
        }
    }
    
    // Method 4: CHGUIT for numerical integration
    if (a >= 1.0 && (!bn || !(bl1 || bl2 || bl3))) {
        int id;
        double hu_temp;
        chguit_opt(a, b, x, &hu_temp, &id);
        if (id > id_best) {
            id_best = id;
            hu_best = hu_temp;
            md_best = 4;
        }
    }
    // Handle a < 1 cases with transformation
    if (a < 1.0) {
        if (b <= a) {
            double a00 = a, b00 = b;
            a = a - b + 1.0;
            b = 2.0 - b;
            int id;
            double hu_temp;
            chguit_opt(a, b, x, &hu_temp, &id);
            hu_temp = R_pow(x, 1.0 - b00) * hu_temp;
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
    
    // Set final results
    *hu = hu_best;
    *md = md_best;
    
    // Check if we got a valid result
    if (id_best <= 0 || !isfinite(hu_best)) {
        *isfer = 1;
        if (!isfinite(hu_best)) *hu = R_NaN;
    }
}

// Optimized wrapper function with enhanced error handling
double kummer(double a, double b, double c) {
    // Comprehensive input validation
  //  if (!isfinite(a) || !isfinite(b) || !isfinite(c)) {
   //     return R_NaN;
   // }
    
    if (c <= 0.0) {
        return R_NaN;
    }
    // Handle special cases for improved performance
    if (a == 0.0) return 1.0;
    if (fabs(c) < 1e-300) return R_NaN;
    // Handle very large arguments
    if (c > 1e100) {
        // Use asymptotic expansion
        return R_pow(c, -a);
    }
    double result;
    int method, status;
    // Call optimized computation
    chgu_opt(a, b, c, &result, &method, &status);
    // Enhanced error reporting
    if (status == -1) {
        Rprintf("Error: Invalid input parameters for Kummer function\n");
        return R_NaN;
    } else if (status == 1) {
        Rprintf("Warning: Kummer function computation may have reduced accuracy\n");
    }
    
    // Final validation
    if (!isfinite(result)) {
        return R_NaN;
    }
    
    return result;
}
