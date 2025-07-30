#include "header.h"




/**
 * Approssimazione della forma regolarizzata di 1F1
 * Ottimizzazioni:
 * - Limite di iterazioni esplicito
 * - Condizione di uscita più precisa
 * - Pre-calcolo di valori comuni
 */
double aprox_reg_1F1(int n, int m, double z) {
    // Pre-calcolo del termine costante
    double p1 = exp(z + (n-m) * log(z) - lgammafn(n));
    
    double s1 = 0.0;
    double prev_s1 = -1.0; 
    const int MAX_ITER = 1000;
    const double TOLERANCE = 1e-10;
    
    for (int k = 0; k < MAX_ITER; k++) {
        double term = poch(1-n, k) * poch(m-n, k) * R_pow(z, -k) / gammafn(k+1);
        s1 += term;
        
        // Ottimizzazione: usa sia il valore assoluto del termine che il cambio relativo
        if (fabs(term) < TOLERANCE || fabs(s1 - prev_s1) < TOLERANCE * fabs(s1)) {
            break;
        }
        
        prev_s1 = s1;
    }
    
    return p1 * s1;
}



/* ----------------------------------------------------------
   FUNCTION 1F1– optmized version
 ---------------------------------------------------------- */

/* ---------- cache locale per lgamma -------------------------------- */
static double lgamma_cache_key = 0.0;
static double lgamma_cache_val = 0.0;
static double hy1f1p(double a, double b, double x, double *err);
static double hy1f1a(double a, double b, double x, double *err);
static double hyp2f0(double a, double b, double x, int type, double *err);

static inline double cached_lgamma(double x) {
    if (x == lgamma_cache_key) return lgamma_cache_val;
    lgamma_cache_key = x;
    lgamma_cache_val = lgammafn(x);
    return lgamma_cache_val;
}

/* ---------- funzione principale 1F1 ------------------------------- */
double hyperg(double a, double b, double x) {
    /* casi speciali – risposte immediate */
    if (fabs(x) < 1e-8) return 1.0 + a * x / b;
    if (fabs(a) < 1e-15) return 1.0;
    if (fabs(a - b) < 1e-15) return exp(x);

    int use_asymptotic = fabs(x) > (10.0 + fabs(a) + fabs(b));
    double err;
    double res = use_asymptotic ?
                 hy1f1a(a, b, x, &err) :
                 hy1f1p(a, b, x, &err);

    if (err > 1e-12 && fabs(b - a) > 0.001 * fabs(a)) {
        double alt = exp(x) * hyperg(b - a, b, -x);
        return (err < 1e-6) ? res : alt;
    }
    return res;
}


/* ---------- serie di potenze 1F1 ---------------------------------- */
static double hy1f1p(double a, double b, double x, double *err) {
    if (b == 0.0)           { *err = 1.0; return R_PosInf; }
    if (a == 0.0)           { *err = 0.0; return 1.0; }

    double sum  = 1.0;
    double term = 1.0;
    double n    = 1.0;
    double maxn = 200.0 + 2.0 * fabs(a) + 2.0 * fabs(b);

    while (n <= maxn && fabs(term) > 1e-16) {
        term *= x * (a + n - 1.0) / ((b + n - 1.0) * n);
        sum  += term;
        n    += 1.0;
    }
    *err = fabs(term / (sum != 0.0 ? sum : 1.0));
    return sum;
}

/* ---------- espansione asintotica 1F1 ----------------------------- */
static double hy1f1a(double a, double b, double x, double *err) {
    if (x == 0.0) { *err = 1.0; return R_PosInf; }

    double lgb  = cached_lgamma(b);
    double lgba = cached_lgamma(b - a);
    double lga  = (a <= 0.0) ? 0.0 : lgammafn(a);

    double t = x + log(fabs(x)) * (a - b) + lgb;
    double u = -log(fabs(x)) * a + lgb;

    double err1, err2;
    double h1 = hyp2f0(a, a - b + 1.0, -1.0 / x, 1, &err1);
    h1 *= exp(u - lgba);
    err1 *= exp(u - lgba);

    double h2 = hyp2f0(b - a, 1.0 - a, 1.0 / x, 2, &err2);
    h2 *= exp(t - lga);
    err2 *= exp(t - lga);

    double res = (x < 0.0) ? h1 : h2;
    *err = fabs(err1) + fabs(err2);
    return res;
}
/* ---------- 2F0 (serie asintotica ausiliaria) ---------------------- */
static double hyp2f0(double a, double b, double x, int type, double *err) {
    double an = a;
    double bn = b;
    double a0 = 1.0;
    double sum = 0.0;
    double n = 1.0;
    const int MAX_ITER = 200;
    *err = 1.0;

    do {
        if (an == 0.0 || bn == 0.0) {
            sum += a0;
            *err = 0.0;
            return sum;
        }
        a0 *= an * bn * x / n;
        if (fabs(a0) > 1e200) break;     /* evita overflow */
        sum += a0;
        an += 1.0;
        bn += 1.0;
        n  += 1.0;
    } while (n <= MAX_ITER && fabs(a0) > 1e-16);

    *err = fabs(a0);
    return sum;
}




void  hyperg_call(double *a,double *b,double *x,double *res)
{
    *res = hyperg(*a,*b,*x);
}


