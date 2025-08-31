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






double hy1f1p(double a, double b, double x, double *err) {
    // Check precoci invariati
    if (b == 0) {
        *err = 1.0;
        return NPY_INFINITY;
    }
    
    if (a == 0) {
        *err = 0;
        return 1.0;
    }
    
    double sum = 1.0;
    double c = 0.0;  // Correzione Kahan
    double term = 1.0;
    double max_term = 1.0;
    
    // Pre-calcolo per evitare divisioni ripetute nel loop
    const double x_val = x;
    const double abs_a = fabs(a);
    const double abs_b = fabs(b);
    const int maxn = (int)(200.0 + 2 * abs_a + 2 * abs_b);
    
    // Parametri del loop ottimizzati
    double an = a;
    double bn = b;
    
    for (int n = 1; n <= maxn; n++) {
        // Calcolo ottimizzato del fattore moltiplicativo
        const double factor = (an * x_val) / (bn * n);
        term *= factor;
        
        // Check overflow anticipato
        const double abs_term = fabs(term);
        if (abs_term > DBL_MAX * 0.1) {
            *err = 1.0;
            return sum;
        }
        
        // Aggiornamento max_term
        if (abs_term > max_term) {
            max_term = abs_term;
        }
        
        // Kahan summation
        const double y = term - c;
        const double temp_sum = sum + y;
        c = (temp_sum - sum) - y;
        sum = temp_sum;
        
        // Test di convergenza ottimizzato
        if (abs_term <= MACHEP) {
            break;
        }
        
        // Incremento dei parametri
        an += 1.0;
        bn += 1.0;
    }
    
    // Stima errore
    *err = (sum != 0.0) ? fabs(c / sum) : fabs(c);
    if (*err != *err) {
        *err = 1.0;
    }
    
    return sum;
}

double hyp2f0(double a, double b, double x, int type, double *err) {
    double an = a;
    double bn = b;
    double sum = 1.0;  // Inizializza con il primo termine
    double term = 1.0;
    double prev_abs_term = 1e9;
    double max_term = 1.0;
    
    const int MAX_ITER = 200;
    const double overflow_limit = DBL_MAX * 0.1;
    
    for (int n = 1; n <= MAX_ITER; n++) {
        // Check terminazione anticipata
        if (an == 0 || bn == 0) {
            goto pdone;
        }
        
        // Calcolo del fattore ottimizzato
        const double factor = (an * bn * x) / n;
        term *= factor;
        
        const double abs_term = fabs(term);
        
        // Check overflow
        if (abs_term > overflow_limit) {
            goto error;
        }
        
        // Test divergenza
        if (abs_term > prev_abs_term) {
            goto ndone;
        }
        
        prev_abs_term = abs_term;
        sum += term;
        
        if (abs_term > max_term) {
            max_term = abs_term;
        }
        
        // Test convergenza
        if (abs_term <= MACHEP) {
            break;
        }
        
        an += 1.0;
        bn += 1.0;
    }
    
pdone:  // Serie convergente
    *err = fabs(MACHEP * max_term);
    return sum;
    
ndone:  // Serie non convergente
    {
        const double inv_x = 1.0 / x;
        const double n_val = an - a;  // Numero di iterazioni effettuate
        
        // Fattori di convergenza ottimizzati
        switch (type) {
            case 1:
                term *= (0.5 + (0.125 + 0.25 * b - 0.5 * a + 0.25 * x - 0.25 * n_val) * inv_x);
                break;
            case 2:
                term *= (2.0/3.0 - b + 2.0 * a + x - n_val);
                break;
        }
        
        sum += term;
        *err = MACHEP * max_term + fabs(term);
        return sum;
    }
    
error:  // Overflow
    *err = NPY_INFINITY;
    return sum;
}

double hy1f1a(double a, double b, double x, double *err) {
    // Caso speciale
    if (x == 0) {
        *err = 1.0;
        return NPY_INFINITY;
    }
    
    // Pre-calcoli ottimizzati
    const double abs_x = fabs(x);
    const double log_abs_x = log(abs_x);
    const double inv_x = 1.0 / x;
    
    // Primo termine: exp(u) / Γ(b-a) * 2F0(a, a-b+1, -1/x)
    double u = -log_abs_x * a;
    if (b > 0) {
        u += lgammafn(b);
    }
    u += x + log_abs_x * (a - b);
    
    double err1, err2;
    double h1 = hyp2f0(a, a - b + 1, -inv_x, 1, &err1);
    
    // Gestione dei casi di underflow/overflow
    double temp1;
    if (b - a > 0) {
        temp1 = exp(u - lgammafn(b - a));
    } else {
        temp1 = exp(u) / gammafn(b - a);
    }
    
    h1 *= temp1;
    err1 *= fabs(temp1);
    
    // Secondo termine: exp(t) / Γ(a) * 2F0(b-a, 1-a, 1/x)
    double t = x + log_abs_x * (a - b);
    if (b > 0) {
        t += lgammafn(b);
    }
    
    double h2 = hyp2f0(b - a, 1.0 - a, inv_x, 2, &err2);
    
    double temp2;
    if (a > 0) {
        temp2 = exp(t - lgammafn(a));
    } else {
        temp2 = exp(t) / gammafn(a);
    }
    
    h2 *= temp2;
    err2 *= fabs(temp2);
    
    // Selezione del termine dominante
    double asum = (x < 0.0) ? h1 : h2;
    double acanc = fabs(err1) + fabs(err2);
    
    // Aggiustamento per b < 0
    if (b < 0) {
        const double gamma_b = gammafn(b);
        asum *= gamma_b;
        acanc *= fabs(gamma_b);
    }
    
    // Normalizzazione errore
    if (asum != 0.0) {
        acanc /= fabs(asum);
    }
    if (acanc != acanc) {
        acanc = 1.0;  // NaN
    }
    
    if (asum == NPY_INFINITY || asum == -NPY_INFINITY) {
        acanc = 0;
    }
    
    *err = acanc * 30.0;
    return asum;
}

double hyperg(double a, double b, double x) {
    // Ottimizzazione per parametri molto grandi
    const double alpha = (b - a - 1) / x;
    
    if (alpha > 0 && b > 1000000 && x > 1000000) {
        // Calcoli ottimizzati per grandi parametri
        const double log_b = lgammafn(b);
        const double log_diff = lgammafn(b - a);
        
        if (b < x + a + 1) {
            // Approssimazione asintotica ottimizzata
            const double log_alpha = log(alpha);
            const double one_minus_alpha = 1.0 - alpha;
            
            const double aux1 = log_b + x + (a - 1) * log1p(-alpha) - lgammafn(a) - log_diff;
            const double aux2 = alpha * x * (log_alpha - 1) + 0.5 * (log(2 * alpha * M_PI) - log(x));
            const double aux3 = ((2 - a * alpha) * (1 - a)) / (2 * one_minus_alpha * one_minus_alpha) + 1.0 / (12 * alpha);
            
            return exp(aux1 + aux2) * (1 + aux3 / x);
        } else {
            const double denom = b - a - x - 1;
            const double inv_denom_pow_a = exp(-a * log(denom));
            const double aux1 = exp(log_b - log_diff) * inv_denom_pow_a;
            const double aux2 = 1.0 - (a * (a + 1) * (a + 1 - b)) / (2 * denom * denom);
            
            return aux1 * aux2;
        }
    }

    // Trasformazione di Kummer ottimizzata
    const double temp = b - a;
    if (fabs(temp) < 0.001 * fabs(a)) {
        return exp(x) * hyperg(temp, b, -x);
    }
    
    // Strategia di calcolo ottimizzata
    const double abs_x = fabs(x);
    const double abs_a = fabs(a);
    const double abs_b = fabs(b);
    const bool use_power_series = abs_x < (10 + abs_a + abs_b);
    
    double psum, asum, pcanc, acanc;
    
    if (use_power_series) {
        psum = hy1f1p(a, b, x, &pcanc);
        if (pcanc < 1.0e-15) {
            return psum;
        }
        asum = hy1f1a(a, b, x, &acanc);
    } else {
        asum = hy1f1a(a, b, x, &acanc);
        if (acanc < 1.0e-15) {
            return asum;
        }
        psum = hy1f1p(a, b, x, &pcanc);
    }
    
    return (acanc < pcanc) ? asum : psum;
}

void  hyperg_call(double *a,double *b,double *x,double *res)
{
    *res = hyperg(*a,*b,*x);
}


