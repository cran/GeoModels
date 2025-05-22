#include "header.h"


/**
 * Calcola il logaritmo della forma regolarizzata di 1F1
 * Ottimizzazioni: 
 * - Evitare calcoli ripetuti
 * - Pre-calcolare valori comuni
 * - Migliorare la leggibilità
 */
double log_hyp1F1_reg(int n, int m, double z) {
    if (n < m) {
        // Pre-calcolo dei valori comuni
        double poch_2_m = poch(2-m, n-1);
        double z_1_m = R_pow(z, 1-m);
        double exp_z = exp(z);
        
        // Calcolo dei termini
        double term1 = 0.0;
        for (int k = 0; k <= (n-1); k++) {
            term1 += R_pow(-z, k) * poch(1-n, k) / (gammafn(k+1) * poch(2-m, k));
        }
        
        double term2 = 0.0;
        for (int k = 0; k <= (m-n-1); k++) {
            term2 += R_pow(z, k) * poch(1-m+n, k) / (gammafn(k+1) * poch(2-m, k));
        }
        
        // Risultato finale
        return log(poch_2_m) + log(z_1_m) + log(exp_z * term1 - term2) - lgammafn(n);
    } else {
        double term3 = 0.0;
        for (int k = 0; k <= (n-m); k++) {
            term3 += exp(log(poch(m-n, k)) + log(R_pow(-z, k)) - 
                        (lgammafn(k+1) + log(poch(m, k))));
        }
        return z + log(term3) - lgammafn(m);
    }
}

/**
 * Calcola il logaritmo della forma regolarizzata di 1F1 per casi speciali
 * Ottimizzazione: casi speciali pre-calcolati per n = 1,2,3,4,5
 */
double log_regularized1F1(int n, int m, double z) {
    // Ottimizzazione: gestione separata dei casi particolari
    switch (n) {
        case 1:
            return z + (1-m) * log(z) + log(igam(-1 + m, z));
            
        case 2: {
            double res = exp(-lgammafn(m-1)) + exp(z) * R_pow(z, 1-m) * (2 - m + z) * igam(-1 + m, z);
            return log(res);
        }
            
        case 3: {
            double term1 = (4-m+z) / gammafn(-1+m);
            double term2 = exp(z) * R_pow(z,1-m) * (6-5*m + m*m + 6*z-2*m*z+z*z) * igam(-1 + m, z);
            return log(0.5 * (term1 + term2));
        }
            
        case 4: {
            double inv_gammafn_m1 = 1.0 / gammafn(-1 + m);
            double term1 = (18 - 8*m + m*m + 10*z - 2*m*z + z*z) * inv_gammafn_m1;
            
            // Pre-calcolo per evitare ripetizioni
            double exp_z_pow = exp(z) * R_pow(z, 1-m);
            double m2 = m*m;
            double m3 = m2*m;
            double z2 = z*z;
            double z3 = z2*z;
            
            double term2 = exp_z_pow * (24 - 26*m + 9*m2 - m3 + 36*z - 21*m*z + 
                          3*m2*z + 12*z2 - 3*m*z2 + z3) * igam(-1 + m, z);
            
            return log((term1 + term2) / 6.0);
        }
            
        case 5: {
            double inv_gammafn_m1 = 1.0 / gammafn(-1 + m);
            
            // Pre-calcolo delle potenze
            double m2 = m*m;
            double m3 = m2*m;
            double m4 = m3*m;
            double z2 = z*z;
            double z3 = z2*z;
            double z4 = z3*z;
            
            double term1 = (96 - 58*m + 13*m2 - m3 + 86*z - 31*m*z + 3*m2*z + 
                          18*z2 - 3*m*z2 + z3) * inv_gammafn_m1;
            
            double exp_z_pow = exp(z) * R_pow(z, 1-m);
            double term2 = exp_z_pow * (120 - 154*m + 71*m2 - 14*m3 + m4 + 
                          240*z - 188*m*z + 48*m2*z - 4*m3*z + 
                          120*z2 - 54*m*z2 + 6*(z*m)*(z*m) + 
                          20*z3 - 4*m*z3 + z4) * igam(-1 + m, z);
            
            return log((term1 + term2) / 24.0);
        }
        
        default:
            // Caso generico
            return log(hyperg(n, m, z)) - lgammafn(m);
    }
}

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

/**
 * Wrapper per la funzione log_regularized1F1
 */
void reghyperg_call(int *a, int *b, double *x, double *res) {
    *res = log_regularized1F1(*a, *b, *x);
}

/**
 * Calcola la funzione ipergeometrica confluente 1F1
 * Ottimizzazioni:
 * - Migliorata la struttura del codice
 * - Migliorato il branch prediction
 * - Evitate alcune divisioni e calcoli ripetuti
 */
double hyperg(double a, double b, double x) {
    double alpha = (b - a - 1) / x;
    
    if (alpha > 0 && b > 1000000 && x > 1000000) {
        double res;
        
        if (b < x + a + 1) {
            // Approssimazione per b < x + a + 1
            double aux1 = exp(lgammafn(b) + x + (a-1) * log1p(-alpha) - lgammafn(a) - lgammafn(b-a));
            double aux2 = exp((alpha * x) * (log(alpha) - 1) + 0.5 * (log(2 * alpha * M_PI) - log(x)));
            double aux3 = ((2 - a * alpha) * (1 - a)) / ((2 * pow(1 - alpha, 2))) + (1 / (12 * alpha));
            
            res = aux1 * aux2 * (1 + aux3 / x);
        } else {
            double inv_denom = 1.0 / pow(b - a - x - 1, a);
            double aux1 = exp(lgammafn(b) - lgammafn(b-a)) * inv_denom;
            double aux2 = 1.0 - (a * (a + 1) * (a + 1 - b)) / (2 * pow(b - a - x - 1, 2));
            
            res = aux1 * aux2;
        }
        
        return res;
    }

    // Trasformazione di Kummer (se applicabile)
    double temp = b - a;
    if (fabs(temp) < 0.001 * fabs(a)) {
        return exp(x) * hyperg(temp, b, -x);
    }
    double psum, asum, pcanc, acanc;
    bool use_power_series = fabs(x) < (10 + fabs(a) + fabs(b));
    
    if (use_power_series) {
        psum = hy1f1p(a, b, x, &pcanc);
        
        // Se la precisione è buona, usala subito
        if (pcanc < 1.0e-15) {
            return psum;
        }
        
        asum = hy1f1a(a, b, x, &acanc);
    } else {
        psum = hy1f1a(a, b, x, &pcanc);
        if (pcanc < 1.0e-15) {
            return psum;
        }
        
        asum = hy1f1p(a, b, x, &acanc);
    }
    return (acanc < pcanc) ? asum : psum;
}

/**
 * Calcola 1F1 usando serie di potenze
 * Ottimizzazioni:
 * - Rilevamento più rapido della convergenza
 * - Migliore stima dell'errore
 * - Gestione ottimizzata dei casi speciali
 */
double hy1f1p(double a, double b, double x, double *err) {
    // Check precoce per i valori singolari
    if (b == 0) {
        *err = 1.0;
        return NPY_INFINITY;
    }
    
    if (a == 0) {
        *err = 0;  // Nessun errore
        return 1.0;  // 1F1(0,b,x) = 1
    }
    
    double an = a;
    double bn = b;
    double sum = 1.0;
    double c = 0.0;  // Correzione per Kahan summation
    double n = 1.0;
    double t = 1.0;
    double a0 = 1.0;
    double maxt = 0.0;
    
    // Limite massimo di iterazioni basato sui parametri
    double maxn = 200.0 + 2 * fabs(a) + 2 * fabs(b);
    
    // Criterio di convergenza migliorato
    while (t > MACHEP && n <= maxn) {
        // Calcolo del termine attuale
        double u = x * (an / (bn * n));
        
        // Controllo per overflow
        double temp = fabs(u);
        if ((temp > 1.0) && (maxt > (DBL_MAX / temp))) {
            *err = 1.0;  // Blowup: stima errore 100%
            return sum;
        }
        
        // Aggiornamento del termine e della somma
        a0 *= u;
        
        // Kahan summation per ridurre errori di arrotondamento
        double y = a0 - c;
        double sumc = sum + y;
        c = (sumc - sum) - y;
        sum = sumc;
        
        t = fabs(a0);
        if (t > maxt) {
            maxt = t;
        }
        
        // Aggiornamento dei parametri per il prossimo termine
        an += 1.0;
        bn += 1.0;
        n += 1.0;
    }
    
    // Stima dell'errore
    *err = (sum != 0.0) ? fabs(c / sum) : fabs(c);
    
    // Controllo NaN
    if (*err != *err) {
        *err = 1.0;
    }
    
    return sum;
}

/**
 * Calcola 1F1 usando serie asintotiche
 * Ottimizzazioni:
 * - Migliore gestione dei casi speciali
 * - Riutilizzo di calcoli intermedi
 */
double hy1f1a(double a, double b, double x, double *err) {
    // Caso speciale x = 0
    if (x == 0) {
        *err = 1.0;
        return NPY_INFINITY;
    }
    
    // Pre-calcoli per i termini comuni
    double temp = log(fabs(x));
    double t = x + temp * (a - b);
    double u = -temp * a;
    
    if (b > 0) {
        temp = lgammafn(b);
        t += temp;
        u += temp;
    }

    double err1, err2;
    double h1 = hyp2f0(a, a - b + 1, -1.0 / x, 1, &err1);
    
    temp = exp(u) / gammafn(b - a);
    h1 *= temp;
    err1 *= temp;
    
    double h2 = hyp2f0(b - a, 1.0 - a, 1.0 / x, 2, &err2);
    if (a < 0) {
        temp = exp(t) / gammafn(a);
    } else {
        temp = exp(t - lgammafn(a));
    }
    
    h2 *= temp;
    err2 *= temp;
    
    // Selezione del termine in base al segno di x
    double asum = (x < 0.0) ? h1 : h2;
    double acanc = fabs(err1) + fabs(err2);
    
    // Aggiustamento per b < 0
    if (b < 0) {
        temp = gammafn(b);
        asum *= temp;
        acanc *= fabs(temp);
    }
    if (asum != 0.0) {
        acanc /= fabs(asum);
    }
    if (acanc != acanc) {
        acanc = 1.0;  // NaN
    }
    
    if (asum == NPY_INFINITY || asum == -NPY_INFINITY) {
        acanc = 0;  // Infinito
    }
    acanc *= 30.0;
    
    *err = acanc;
    return asum;
}

/**
 * Calcola la funzione ipergeometrica 2F0
 * Ottimizzazioni:
 * - Condizione di uscita migliorata
 * - Pre-calcolo di valori comuni
 * - Riduzione del numero di divisioni
 */
double hyp2f0(double a, double b, double x, int type, double *err) {
    double an = a;
    double bn = b;
    double a0 = 1.0;
    double alast = 1.0;
    double sum = 0.0;
    double n = 1.0;
    double t = 1.0;
    double tlast = 1.0e9;
    double maxt = 0.0;
    
    const int MAX_ITER = 200;
    
    do {
        if (an == 0 || bn == 0) {
            goto pdone;
        }
        
        // Calcolo del termine attuale
        double u = an * (bn * x / n);
        
        // Controllo overflow
        double temp = fabs(u);
        if ((temp > 1.0) && (maxt > (DBL_MAX / temp))) {
            goto error;
        }
        
        a0 *= u;
        t = fabs(a0);
        if (t > tlast) {
            goto ndone;
        }
        
        tlast = t;
        sum += alast;  
        alast = a0;
        
        if (n > MAX_ITER) {
            goto ndone;
        }
        an += 1.0;
        bn += 1.0;
        n += 1.0;
        if (t > maxt) {
            maxt = t;
        }
    } while (t > MACHEP);
    
pdone:  // Serie convergente
    // Stima dell'errore dovuto all'arrotondamento
    *err = fabs(MACHEP * (n + maxt));
    alast = a0;
    goto done;
    
ndone:  // Serie non convergente
    n -= 1.0;
    x = 1.0 / x;
    
    // "Fattori di convergenza" per migliorare l'accuratezza
    switch (type) {
        case 1:
            alast *= (0.5 + (0.125 + 0.25 * b - 0.5 * a + 0.25 * x - 0.25 * n) / x);
            break;
            
        case 2:
            alast *= 2.0 / 3.0 - b + 2.0 * a + x - n;
            break;
            
        default:
            // Nessuna modifica per tipi non riconosciuti
            break;
    }
    // Stima dell'errore
    *err = MACHEP * (n + maxt) + fabs(a0);
done:
    sum += alast;
    return sum;
    
error:  // Gestione dell'errore di overflow
    *err = NPY_INFINITY;
    return sum;
}

void  hyperg_call(double *a,double *b,double *x,double *res)
{
    *res = hyperg(*a,*b,*x);
}


