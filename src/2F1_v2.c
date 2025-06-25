#include "header.h"

#define CONV_TOL 1.0e-15
static double hyp2f1_series(double a, double b, double c, double x) ;
static double hyp2f1_near_one(double a, double b, double c, double x) ;
static double hyp2f1_negative_x(double a, double b, double c, double x);
static double hyp2f1_integer_params(double a, double b, double c, double x);
static double hyp2f1_half_integer_smooth(double a, double b, double c, double x);
static double hyp2f1_series_extended(double a, double b, double c, double x);
double hypergeom_2f1(double a, double b, double c, double x);

static double hyp2f1_series(double a, double b, double c, double x) {
    double sum = 1.0;
    double term = 1.0;
    double a_n = a;
    double b_n = b;
    double c_n = c;
    double ratio_old = 1.0;
    
    for (int n = 1; n < MAX_ITERATIONS; n++) {
        // Calcolo del termine usando le proprietà dei simboli di Pochhammer
        term *= (a_n * b_n * x) / (c_n * n);
        sum += term;
        
        // Test di convergenza migliorato
        double ratio = fabs(term / sum);
        if (ratio < CONV_TOL) break;
        
        // Test di stagnazione (se il rapporto non migliora)
        if (n > 50 && ratio >= ratio_old * 0.99) break;
        ratio_old = ratio;
        
        a_n += 1.0;
        b_n += 1.0;
        c_n += 1.0;
    }
    
    return sum;
}

// Trasformazione per x vicino a 1 usando funzioni R
static double hyp2f1_near_one(double a, double b, double c, double x) {
    double delta = c - a - b;
    
    // Caso speciale: c = a + b (evita problemi numerici)
    if (fabs(delta) < MACHEP * 10) {
        return R_pow(1.0 - x, -a);
    }
    
    // Usa la trasformazione standard con funzioni gamma R
    if (delta > ETHRESH) {
        // 2F1(a,b;c;x) = Γ(c)Γ(c-a-b)/[Γ(c-a)Γ(c-b)] * 2F1(a,b;a+b-c+1;1-x)
        //                + Γ(c)Γ(a+b-c)/[Γ(a)Γ(b)] * (1-x)^(c-a-b) * 2F1(c-a,c-b;c-a-b+1;1-x)
        
        double y = 1.0 - x;
        
        // Primo termine
        double gamma_c = gammafn(c);
        double gamma_delta = gammafn(delta);
        double gamma_c_a = gammafn(c - a);
        double gamma_c_b = gammafn(c - b);
        
        if (!R_FINITE(gamma_delta) || !R_FINITE(gamma_c_a) || !R_FINITE(gamma_c_b)) {
            return NA_REAL;
        }
        
        double coeff1 = gamma_c * gamma_delta / (gamma_c_a * gamma_c_b);
        double term1 = coeff1 * hyp2f1_series(a, b, a + b - c + 1.0, y);
        
        // Secondo termine
        double gamma_neg_delta = gammafn(-delta);
        double gamma_a = gammafn(a);
        double gamma_b = gammafn(b);
        
        if (!R_FINITE(gamma_neg_delta) || !R_FINITE(gamma_a) || !R_FINITE(gamma_b)) {
            return term1; // Restituisce solo il primo termine
        }
        
        double coeff2 = gamma_c * gamma_neg_delta / (gamma_a * gamma_b);
        double term2 = coeff2 * R_pow(y, delta) * 
                      hyp2f1_series(c - a, c - b, c - a - b + 1.0, y);
        
        return term1 + term2;
        
    } else if (delta < -ETHRESH) {
        // Caso delta negativo - usa ricorrenza
        return R_pow(1.0 - x, -a) * hyp2f1_series(a, c - b, c, x / (x - 1.0));
    } else {
        // Caso limite delta ≈ 0
        return R_pow(1.0 - x, -a);
    }
}

// Trasformazione per x < 0 usando identità di Pfaff
static double hyp2f1_negative_x(double a, double b, double c, double x) {
    // 2F1(a,b;c;x) = (1-x)^(-a) * 2F1(a,c-b;c;x/(x-1))
    double y = x / (x - 1.0);
    double factor = R_pow(1.0 - x, -a);
    return factor * hyp2f1_series(a, c - b, c, y);
}

// Gestione di parametri interi negativi (serie termina)
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
        return NA_REAL; // Non dovrebbe succedere
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
}

// Caso speciale per smooth = n + 0.5 (mezzi interi) - versione migliorata
static double hyp2f1_half_integer_smooth(double a, double b, double c, double x) {
    double delta = c - a - b;
    
    // Verifica che delta sia effettivamente un piccolo intero positivo
    if (delta > 0.5 && delta < 10.5 && fabs(delta - round(delta)) < 1e-12) {

        
        // Per x molto vicino a 1, usa la trasformazione di Euler
        if (x > 0.9) {
            // 2F1(a,b;c;x) = (1-x)^(c-a-b) * 2F1(c-a,c-b;c;x)
            // Ma questo può essere instabile, quindi usiamo un approccio alternativo
            
            // Usa la relazione: 2F1(a,b;c;x) = 2F1(a,c-b;c;x/(x-1)) * (1-x)^(-a)
            double z = x / (x - 1.0);
            
            if (fabs(z) < 0.8) { // Assicurati che z sia nel dominio di convergenza
                double factor = R_pow(1.0 - x, -a);
                if (R_FINITE(factor) && factor > 0) {
                    double series_result = hyp2f1_series(a, c - b, c, z);
                    if (R_FINITE(series_result)) {
                        return factor * series_result;
                    }
                }
            }
            
            // Fallback: usa serie con precisione estesa
            return hyp2f1_series_extended(a, b, c, x);
        }
        
        // Per x non troppo vicino a 1, usa la serie standard
        return hyp2f1_series(a, b, c, x);
    }
    
    return hyp2f1_series(a, b, c, x);
}

// Serie con precisione estesa per casi difficili
static double hyp2f1_series_extended(double a, double b, double c, double x) {
    double sum = 1.0;
    double term = 1.0;
    double a_n = a;
    double b_n = b;
    double c_n = c;
    
    // Usa tolleranza più stringente e più iterazioni
    const double ext_tol = 1.0e-16;
    const int max_iter = 5000;
    
    for (int n = 1; n < max_iter; n++) {
        // Calcola il termine con controllo di overflow
        double new_term = term * (a_n * b_n * x) / (c_n * n);
        
        // Controllo di overflow/underflow
        if (!R_FINITE(new_term)) {
            break;
        }
        
        term = new_term;
        sum += term;
        
        // Test di convergenza più stringente
        double ratio = fabs(term);
        double sum_abs = fabs(sum);
        
        if (sum_abs > 0 && ratio/sum_abs < ext_tol) {
            break;
        }
        
        // Aggiorna i parametri
        a_n += 1.0;
        b_n += 1.0;
        c_n += 1.0;
        
        // Test aggiuntivo per convergenza lenta
        if (n > 100 && n % 50 == 0) {
            if (ratio/sum_abs > 1e-10) { // Se converge troppo lentamente
                // Prova un approccio di accelerazione (trasformazione di Euler)
                break;
            }
        }
    }
    
    return sum;
}

// Funzione principale usando Rmath
double hypergeom_2f1(double a, double b, double c, double x) {
    
    // Validazione input usando funzioni R
    if (!R_FINITE(a) || !R_FINITE(b) || !R_FINITE(c) || !R_FINITE(x)) {
        return NA_REAL;
    }
    
    if (fabs(x) > 1.0) {
        return NA_REAL; // Fuori dal dominio supportato
    }
    
    // Casi banali
    if (x == 0.0) return 1.0;
    if (a == 0.0 || b == 0.0) return 1.0;
    
    // Gestione c = 0, -1, -2, ... (poli)
    if (c <= 0.0 && c == floor(c)) {
        return R_PosInf;
    }
    
    // Caso x = 1 usando la formula di Gauss
    if (fabs(x - 1.0) < MACHEP) {
        double delta = c - a - b;
        if (delta > 0.0) {
            return gammafn(c) * gammafn(delta) / (gammafn(c - a) * gammafn(c - b));
        } else if (delta == 0.0) {
            return R_PosInf;
        } else {
            return NA_REAL; // Serie divergente
        }
    }
    
    // Gestione parametri interi negativi (serie polinomiale)
    if ((a <= 0.0 && a == floor(a)) || (b <= 0.0 && b == floor(b))) {
        return hyp2f1_integer_params(a, b, c, x);
    }
    
    // NUOVO: Controllo speciale per il caso smooth = n + 0.5
    double delta = c - a - b;
    if (fabs(delta - floor(delta + 0.5)) < 1e-10 && delta > 0.5) {
        return hyp2f1_half_integer_smooth(a, b, c, x);
    }
    
    // Selezione algoritmo basata su x e parametri (ORIGINALE)
    if (fabs(x) <= 0.7) {
        // Serie standard per |x| piccolo
        return hyp2f1_series(a, b, c, x);
    } else if (x < 0.0) {
        // Trasformazione per x negativo
        return hyp2f1_negative_x(a, b, c, x);
    } else if (x < 0.95) {
        // Zona intermedia - scelta basata sui parametri
        if (c - a - b > 0.5) {
            return hyp2f1_series(a, b, c, x);
        } else {
            return hyp2f1_near_one(a, b, c, x);
        }
    } else {
        // x vicino a 1
        return hyp2f1_near_one(a, b, c, x);
    }
}

double hypergeo2(double a, double b, double c, double x) {    
    return hypergeom_2f1(a, b, c, x);
}
