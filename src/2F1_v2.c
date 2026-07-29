#include "header.h"
#define CONV_TOL 1.0e-15
static double hyp2f1_series_fast(double a, double b, double c, double x);
static double hyp2f1_near_one(double a, double b, double c, double x);
static double hyp2f1_negative_x(double a, double b, double c, double x);
//static double hyp2f1_integer_params(double a, double b, double c, double x);
static double hyp2f1_half_integer_smooth(double a, double b, double c, double x);
double hypergeom_2f1(double a, double b, double c, double x);



static int is_half_integer(double x) {
    double twice = 2.0 * x;
    double rounded = round(twice);
    if (fabs(twice - rounded) < 1e-8) {
        int int_val = (int)rounded;
        return (int_val % 2 != 0);
    }
    return 0;
}

/*
static double hyp2f1_integer_params(double a, double b, double c, double x) {
    int n_terms;

    if (a <= 0.0 && a == floor(a)) {
        n_terms = (int)(-a);
    } else if (b <= 0.0 && b == floor(b)) {
        n_terms = (int)(-b);
        // Usa simmetria: scambia a e b
        double temp = a;
        a = b;
        b = temp;
    } else {
        return NA_REAL; // non dovrebbe succedere
    }

    double sum = 1.0;
    double term = 1.0;
    double a_n = a;
    double b_n = b;
    double c_n = c;

    for (int n = 1; n <= n_terms; n++) {
        term *= (a_n * b_n * x) / (c_n * n);
        sum += term;
        a_n += 1.0;
        b_n += 1.0;
        c_n += 1.0;
    }

    return sum;
}*/


// Serie ottimizzata con somma compensata
static double hyp2f1_series_fast(double a, double b, double c, double x) {
    double sum = 1.0;
    double term = 1.0;
    double c_err = 0.0; // compensazione Kahan

    double a_n = a;
    double b_n = b;
    double c_n = c;

    for (int n = 1; n < MAX_ITERATIONS; n++) {
        term *= (a_n * b_n * x) / (c_n * n);

        // Somma compensata (Kahan)
        double y = term - c_err;
        double t = sum + y;
        c_err = (t - sum) - y;
        sum = t;

        // Controllo convergenza
        if (fabs(term) < fabs(sum) * CONV_TOL) break;

        // Stima del resto geometrico
        double r_est = fabs((a_n * b_n * x) / ((c_n + 1.0) * (n + 1.0)));
        if (r_est < 1.0) {
            double remainder_bound = fabs(term) * r_est / (1.0 - r_est);
            if (remainder_bound < 1e-15 * fabs(sum)) break;
        }

        a_n += 1.0;
        b_n += 1.0;
        c_n += 1.0;
    }

    return sum;
}

// La tua versione (classica, con gammafn e R_pow)
static double hyp2f1_near_one(double a, double b, double c, double x) {
    double delta = c - a - b;
    
    if (fabs(delta) < MACHEP * 10) {
        return R_pow(1.0 - x, -a);
    }

    if (delta > ETHRESH) {
        double y = 1.0 - x;

        double gamma_c   = gammafn(c);
        double gamma_d   = gammafn(delta);
        double gamma_c_a = gammafn(c - a);
        double gamma_c_b = gammafn(c - b);
        if (!R_FINITE(gamma_d) || !R_FINITE(gamma_c_a) || !R_FINITE(gamma_c_b)) {
            return NA_REAL;
        }
        double coeff1 = gamma_c * gamma_d / (gamma_c_a * gamma_c_b);
        double term1  = coeff1 * hyp2f1_series_fast(a, b, a + b - c + 1.0, y);

        double gamma_neg_d = gammafn(-delta);
        double gamma_a     = gammafn(a);
        double gamma_b     = gammafn(b);
        if (!R_FINITE(gamma_neg_d) || !R_FINITE(gamma_a) || !R_FINITE(gamma_b)) {
            return term1;
        }
        double coeff2 = gamma_c * gamma_neg_d / (gamma_a * gamma_b);
        double term2  = coeff2 * R_pow(y, delta) *
                        hyp2f1_series_fast(c - a, c - b, c - a - b + 1.0, y);
        
        return term1 + term2;
    } else if (delta < -ETHRESH) {
        return R_pow(1.0 - x, -a) *
               hyp2f1_series_fast(a, c - b, c, x / (x - 1.0));
    } else {
        return R_pow(1.0 - x, -a);
    }
}

// x < 0 (Pfaff)
static double hyp2f1_negative_x(double a, double b, double c, double x) {
    double y = x / (x - 1.0);
    double factor = R_pow(1.0 - x, -a);
    return factor * hyp2f1_series_fast(a, c - b, c, y);
}

// parametri interi negativi: lasciata invariata
// ...

// half integer smooth: basta cambiare le chiamate
static double hyp2f1_half_integer_smooth(double a, double b, double c, double x) {
    double delta = c - a - b;
    if (delta > 0.5 && delta < 10.5 && fabs(delta - round(delta)) < 1e-12) {
        if (x > 0.9) {
            double z = x / (x - 1.0);
            double factor = R_pow(1.0 - x, -a);
            double series_result = hyp2f1_series_fast(a, c - b, c, z);
            if (R_FINITE(series_result)) {
                return factor * series_result;
            }
            return hyp2f1_series_fast(a, b, c, x);
        }
        return hyp2f1_series_fast(a, b, c, x);
    }
    return hyp2f1_series_fast(a, b, c, x);
}

// Dispatcher principale
double hypergeom_2f1(double a, double b, double c, double x) {
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(c) || !R_FINITE(x)) return NA_REAL;
    if (fabs(x) > 1.0) return NA_REAL;

    if (x == 0.0) return 1.0;
    if (a == 0.0 || b == 0.0) return 1.0;
    if (c <= 0.0 && c == floor(c)) return R_PosInf;

    if (fabs(x - 1.0) < MACHEP) {
        double delta = c - a - b;
        if (delta > 0.0) {
            return gammafn(c) * gammafn(delta) / (gammafn(c - a) * gammafn(c - b));
        } else if (delta == 0.0) {
            return R_PosInf;
        } else {
            return NA_REAL;
        }
    }

   /* if ((a <= 0.0 && a == floor(a)) || (b <= 0.0 && b == floor(b))) {
        return hyp2f1_integer_params(a, b, c, x);
    }*/

    double delta = c - a - b;
    if (fabs(delta - floor(delta + 0.5)) < 1e-10 && delta > 0.5) {
        return hyp2f1_half_integer_smooth(a, b, c, x);
    }

    if (fabs(x) <= 0.7) {
        return hyp2f1_series_fast(a, b, c, x);
    } else if (x < 0.0) {
        return hyp2f1_negative_x(a, b, c, x);
    } else if (x < 0.95) {
        if (c - a - b > 0.5) {
            return hyp2f1_series_fast(a, b, c, x);
        } else {
            return hyp2f1_near_one(a, b, c, x);
        }
    } else {
        return hyp2f1_near_one(a, b, c, x);
    }
}

double hypergeo2(double a, double b, double c, double x) {
    if (is_half_integer(c)) {return hypergeo(a, b, c, x);}
    else return hypergeom_2f1(a, b, c, x);
}
