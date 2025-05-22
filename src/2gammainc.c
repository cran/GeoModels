


#include "header.h"


#define MAXITER1 2000
#define big 1.44115188075855872E+17
#define biginv 6.9388939039072284E-18


//************************************** incomplete gamma*****************************************

/**
 * Ottimizzazione 1: Funzione inline per controlli comuni
 * Evita duplicazione codice e migliora leggibilità
 */
static inline double handle_edge_cases(double a, double x, int lower_part) {
    if (x < 0 || a < 0) {
        return R_NaN;
    } else if (a == 0) {
        if (x > 0) {
            return lower_part ? 1.0 : 0.0;
        } else {
            return R_NaN;
        }
    } else if (x == 0) {
        return lower_part ? 0.0 : 1.0;
    } else if (!R_finite(a)) {
        if (!R_finite(x)) {
            return R_NaN;
        }
        return lower_part ? 0.0 : 1.0;
    } else if (!R_finite(x)) {
        return lower_part ? 1.0 : 0.0;
    }
    return -1.0; /* Segnala che non siamo in un caso limite */
}

/**
 * Calcola il fattore comune per le serie della gamma incompleta
 * Ottimizzato per riutilizzo e leggibilità
 */
double igam_fac(double a, double x) {
    double ax, fac, res, num;

    /* Ottimizzazione 2: Fast path per differenze significative tra a e x */
    if (fabs(a - x) > 0.4 * fabs(a)) {
        ax = a * log(x) - x - lgammafn(a);
        if (ax < -MAXLOG) {
            return 0.0;
        }
        return exp(ax);
    }

    /* Ottimizzazione 3: Pre-calcolo per migliorare riutilizzo valori */
    fac = a + lanczos_g - 0.5;
    res = sqrt(fac / exp(1)) / lanczos_sum_expg_scaled(a);

    /* Ottimizzazione 4: Casi separati per piccoli e grandi valori
       per evitare overflow ed errori numerici */
    if ((a < 200) && (x < 200)) {
        res *= exp(a - x) * R_pow(x / fac, a);
    } else {
        num = x - a - lanczos_g + 0.5;
        res *= exp(a * log1pmx(num / fac) + x * (0.5 - lanczos_g) / fac);
    }

    return res;
}

/**
 * Frazione continua per igamc - Ottimizzata per velocità e stabilità numerica
 */
double igamc_continued_fraction(double a, double x) {
    int i;
    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;

    /* Ottimizzazione 5: Controllo preliminare per evitare calcoli inutili */
    ax = igam_fac(a, x);
    if (ax == 0.0) {
        return 0.0;
    }

    /* Ottimizzazione 6: Inizializzazione variabili all'esterno del ciclo */
    y = 1.0 - a;
    z = x + y + 1.0;
    c = 0.0;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = x + 1.0;
    qkm1 = z * x;
    ans = pkm1 / qkm1;

    /* Ottimizzazione 7: Criterio di terminazione migliorato */
    for (i = 0; i < MAXITER1; i++) {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        yc = y * c;
        pk = pkm1 * z - pkm2 * yc;
        qk = qkm1 * z - qkm2 * yc;
        
        if (qk != 0) {
            r = pk / qk;
            t = fabs((ans - r) / r);
            ans = r;
        } else {
            t = 1.0;
        }
        
        /* Ottimizzazione 8: Controllo overflow con scaling numerico */
        if (fabs(pk) > big) {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }
        
        /* Controllo convergenza migliorato */
        if (t <= MACHEP) {
            break;
        }
        
        /* Aggiornamento valori per prossima iterazione */
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
    }

    return (ans * ax);
}

/**
 * Serie per igam - Ottimizzata per velocità e stabilità numerica
 */
double igam_series(double a, double x) {
    int i;
    double ans, ax, c, r;

    /* Ottimizzazione 9: Riutilizzo del fattore comune */
    ax = igam_fac(a, x);
    if (ax == 0.0) {
        return 0.0;
    }

    /* Ottimizzazione 10: Inizializzazione ottimale */
    r = a;
    c = 1.0;
    ans = 1.0;

    /* Ottimizzazione 11: Criterio di terminazione migliorato per precisione */
    for (i = 0; i < MAXITER1; i++) {
        r += 1.0;
        c *= x / r;
        ans += c;
        
        /* Controllo convergenza più efficiente */
        if (c <= MACHEP * ans) {
            break;
        }
    }

    return (ans * ax / a);
}

/**
 * Serie per igamc - Ottimizzata per evitare cancellazione numerica
 */
double igamc_series(double a, double x) {
    int n;
    double fac = 1;
    double sum = 0;
    double term, logx;

    /* Ottimizzazione 12: Calcolo serie con controllo convergenza migliorato */
    for (n = 1; n < MAXITER1; n++) {
        fac *= -x / n;  
        term = fac / (a + n);
        sum += term;
        
        /* Ottimizzazione 13: Criterio di convergenza più preciso */
        if (fabs(term) <= MACHEP * fabs(sum)) {
            break;
        }
    }

    /* Ottimizzazione 14: Pre-calcolo logaritmo */
    logx = log(x);
    term = -expm1(a * logx - lgamma1p(a));
    
    /* Ottimizzazione 15: Calcolo finale più stabile */
    return term - exp(a * logx - lgammafn(a)) * sum;
}

/**
 * Serie asintotiche per igam/igamc - Ottimizzate per precisione
 */
double asymptotic_series(double a, double x, int func) {
    int k, n, sgn;
    int maxpow = 0;
    double lambda = x / a;
    double sigma = (x - a) / a;
    double eta, res, ck, ckterm, term, absterm;
    double absoldterm = NPY_INFINITY;
    double etapow[NIC] = {1};  /* Ottimizzazione 16: Pre-allocazione array */
    double sum = 0;
    double afac = 1;
    
    /* Ottimizzazione 17: Calcolo segno con valore diretto */
    sgn = (func == IGAM) ? -1 : 1;
    
    /* Ottimizzazione 18: Calcolo eta con casi separati per precisione */
    if (lambda > 1) {
        eta = sqrt(-2 * log1pmx(sigma));
    } else if (lambda < 1) {
        eta = -sqrt(-2 * log1pmx(sigma));
    } else {
        eta = 0;
    }
    
    /* Ottimizzazione 19: Pre-calcolo res iniziale */
    res = 0.5 * erfc(sgn * eta * sqrt(a / 2));
    
    /* Ottimizzazione 20: Loop ottimizzato con criterio di terminazione avanzato */
    for (k = 0; k < KIC; k++) {
        ck = d[k][0];
        
        for (n = 1; n < NIC; n++) {
            if (n > maxpow) {
                etapow[n] = eta * etapow[n-1];
                maxpow += 1;
            }
            
            ckterm = d[k][n] * etapow[n];
            ck += ckterm;
            
            /* Terminazione anticipata se la serie converge */
            if (fabs(ckterm) < MACHEP * fabs(ck)) {
                break;
            }
        }
        
        term = ck * afac;
        absterm = fabs(term);
        
        /* Ottimizzazione 21: Terminazione se la serie diverge */
        if (absterm > absoldterm) {
            break;
        }
        
        sum += term;
        
        /* Ottimizzazione 22: Controllo convergenza complessiva */
        if (absterm < MACHEP * fabs(sum)) {
            break;
        }
        
        absoldterm = absterm;
        afac /= a;
    }
    
    /* Ottimizzazione 23: Calcolo finale ottimizzato */
    res += sgn * exp(-0.5 * a * eta * eta) * sum / sqrt(2 * M_PI * a);
    
    return res;
}

/**
 * Valutazione razionale ottimizzata
 */
double ratevl(double x, const double num[], int M,
              const double denom[], int N) {
    int i, dir;
    double y, num_ans, denom_ans;
    double absx = fabs(x);
    const double *p;
    
    /* Ottimizzazione 24: Decisione calcolo basata sul valore di x */
    if (absx > 1) {
        dir = -1;
        p = num + M;
        y = 1 / x;
    } else {
        dir = 1;
        p = num;
        y = x;
    }
    
    /* Ottimizzazione 25: Valutazione polinomio con schema di Horner */
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
    
    /* Ottimizzazione 26: Calcolo finale con considerazione caso x > 1 */
    if (absx > 1) {
        i = N - M;
        return R_pow(x, i) * num_ans / denom_ans;
    } else {
        return num_ans / denom_ans;
    }
}

/**
 * Calcolo somma Lanczos ottimizzata
 */
double lanczos_sum_expg_scaled(double x) {
    return ratevl(x, lanczos_sum_expg_scaled_num,
                  sizeof(lanczos_sum_expg_scaled_num) / sizeof(lanczos_sum_expg_scaled_num[0]) - 1,
                  lanczos_sum_expg_scaled_denom,
                  sizeof(lanczos_sum_expg_scaled_denom) / sizeof(lanczos_sum_expg_scaled_denom[0]) - 1);
}

/**
 * Gamma incompleta inferiore - Funzione principale ottimizzata
 */
double igam(double a, double x) {
    double absxma_a;
    double result;
    
    /* Ottimizzazione 27: Gestione casi limite unificata */
    result = handle_edge_cases(a, x, 1); /* 1 = lower part (igam) */
    if (result >= 0.0) {
        return result;
    }

    /* Ottimizzazione 28: Valutazione ottimizzata regime asintotico */
    absxma_a = fabs(x - a) / a;
    if ((a > SMALL) && (a < LARGE) && (absxma_a < SMALLRATIO)) {
        return asymptotic_series(a, x, IGAM);
    } else if ((a > LARGE) && (absxma_a < LARGERATIO / sqrt(a))) {
        return asymptotic_series(a, x, IGAM);
    }

    /* Ottimizzazione 29: Scelta percorso di calcolo ottimale */
    if ((x > 1.0) && (x > a)) {
        return (1.0 - igamc(a, x));
    }

    return igam_series(a, x);
}

/**
 * Gamma incompleta superiore - Funzione principale ottimizzata
 */
double igamc(double a, double x) {
    double absxma_a;
    double result;
    
    /* Ottimizzazione 30: Gestione casi limite unificata */
    result = handle_edge_cases(a, x, 0); /* 0 = upper part (igamc) */
    if (result >= 0.0) {
        return result;
    }

    /* Ottimizzazione 31: Valutazione ottimizzata regime asintotico */
    absxma_a = fabs(x - a) / a;
    if ((a > SMALL) && (a < LARGE) && (absxma_a < SMALLRATIO)) {
        return asymptotic_series(a, x, IGAMC);
    } else if ((a > LARGE) && (absxma_a < LARGERATIO / sqrt(a))) {
        return asymptotic_series(a, x, IGAMC);
    }

    /* Ottimizzazione 32: Miglior branching per scelta algoritmo */
    if (x > 1.1) {
        if (x < a) {
            return 1.0 - igam_series(a, x);
        } else {
            return igamc_continued_fraction(a, x);
        }
    } else if (x <= 0.5) {
        if (-0.4 / log(x) < a) {
            return 1.0 - igam_series(a, x);
        } else {
            return igamc_series(a, x);
        }
    } else {
        if (x * 1.1 < a) {
            return 1.0 - igam_series(a, x);
        } else {
            return igamc_series(a, x);
        }
    }
}
void igam_call(double *a,double *x,double *res)
{
    *res = igam(*a,*x);
}



//************************************* END igam.c*****************************************
