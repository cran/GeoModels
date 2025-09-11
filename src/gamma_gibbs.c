#include "header.h"
#define IDX(i,j,nrow) ((i) + (int)(j) * (int)(nrow))





/* ====================== Helpers internos (static) ====================== */

static inline void sample_sphere(double Rlen, int v, double *out) {
  double sum_sq = 0.0;
  for (int i = 0; i < v; ++i) {
    double z = norm_rand();
    out[i] = z;
    sum_sq += z * z;
  }
  double inv_len = 1.0 / sqrt(sum_sq);
  for (int i = 0; i < v; ++i) out[i] *= inv_len * Rlen;  // z/||z|| * R
}

static inline double mvn_iso_pdf(const double *x, const double *u,
                                 double s, int v) {
  // f(x) = (2π s)^(-v/2) * exp( -||x-u||^2 / (2s) )
  double norm2 = 0.0;
  for (int i = 0; i < v; ++i) {
    double d = x[i] - u[i];
    norm2 += d * d;
  }
  double coeff = pow(2.0 * M_PI * s, -0.5 * (double)v);
  return coeff * exp(-norm2 / (2.0 * s));
}

static inline void simple_krige_from_inv_core(const double *inv,
                                              int n, int iii, double *w_out, double *var_out) {
  if (iii < 0 || iii >= n) error("iii fuera de rango.");
  double diag_val = inv[IDX(iii, iii, n)];
  if (!(diag_val > 0.0)) error("Inv[ii,ii] no positivo.");
  int pos = 0;
  for (int c = 0; c < n; ++c) {
    if (c == iii) continue;
    w_out[pos++] = -inv[IDX(iii, c, n)] / diag_val;
  }
  *var_out = 1.0 / diag_val;
}

/* ====================== Funciones EXPORTADAS (C64) ====================== */

/* ---- 1) Kriging simple a partir de la inversa ----
 inv:    matriz n×n (col-major, como en R)
 n:      dimensión (n)
 iii:    índice 1-based (como en R); internamente se pasa a 0-based
 w_out:  salida de largo n-1
 var_out: salida escalar
 */
void C64_simple_krige_from_inv(double *inv,
                               long long *n,
                               long long *iii,
                               double *w_out,
                               double *var_out) {
  int nloc = (int)(*n);
  int iloc = (int)((*iii) - 1);  // R->C
  if (nloc <= 1) error("n debe ser >= 2.");
  simple_krige_from_inv_core(inv, nloc, iloc, w_out, var_out);
}

/* ---- 2) Densidad normal multivariante isótropa (s*I) ----
 x, u:  vectores de largo v
 s:     varianza escalar (>0)
 v:     dimensión
 out:   escalar con la densidad
 */
void C64_mult_dist_u(double *x,
                     double *u,
                     double *s,
                     long long *v,
                     double *out) {
  if (*s <= 0.0) error("s debe ser positivo.");
  if (*v <= 0)   error("v debe ser positivo.");
  *out = mvn_iso_pdf(x, u, *s, (int)(*v));
}

/* ---- 3) Bucle condicional (Gibbs + MH) IN-PLACE sobre U ----
 SigmaInv: n×n inversa (col-major)
 y:        vector de largo n
 n:        número de sitios
 v:        dimensión del vector U_i (columnas de U)
 U:        matriz n×v (col-major), se ACTUALIZA IN-PLACE
 nIte:     número de sweeps de Gibbs
 nRep:     iteraciones de Metropolis-Hastings por sitio
 */
void gamma_gibbs_sampler(double *SigmaInv,
                         double *y,
                         int *n,
                         int *v,
                         double *U,
                         int *nIte,
                         int *nRep) {
  const int nloc   = *n;
  const int vloc   = *v;
  const int nIteL  = *nIte;
  const int nRepL  = *nRep;
  
  if (nloc <= 0 || vloc <= 0) error("n y v deben ser positivos.");
  if (nIteL < 0 || nRepL < 0) error("nIte y nRep no pueden ser negativos.");
  
  // Buffers - simplified without int/long long switching
  int *perm = (int*) R_alloc(nloc, sizeof(int));
  for (int i = 0; i < nloc; ++i) perm[i] = i;
  
  double *u_mean = (double*) R_alloc(vloc, sizeof(double));
  double *us     = (double*) R_alloc(vloc, sizeof(double));
  double *u_sph  = (double*) R_alloc(vloc, sizeof(double));
  double *wbuf   = (double*) R_alloc((nloc>1? nloc-1:1), sizeof(double)); // w de largo n-1
  
  GetRNGstate();
  
  for (int ite = 0; ite < nIteL; ++ite) {
    
    // 1) Barajar índices (Fisher–Yates)
    for (int ii = nloc - 1; ii > 0; --ii) {
      // j ~ U{0..ii}
      int j = (int) floor(unif_rand() * (double)(ii + 1));
      int tmp     = perm[ii];
      perm[ii]    = perm[j];
      perm[j]     = tmp;
    }
    
    // 2) Barrido de Gibbs
    for (int idx = 0; idx < nloc; ++idx) {
      int x_alpha = perm[idx];
      
      // 2a) Kriging (w y var desde SigmaInv)
      double var_alpha;
      simple_krige_from_inv_core(SigmaInv, nloc, x_alpha, wbuf, &var_alpha);
      
      // 2b) u_mean = U[-x_alpha,]^T %*% w
      for (int j = 0; j < vloc; ++j) {
        double acc = 0.0;
        int p = 0;
        // Recorre filas k != x_alpha
        for (int k = 0; k < nloc; ++k) {
          if (k == x_alpha) continue;
          acc += U[IDX(k, j, nloc)] * wbuf[p++];
        }
        u_mean[j] = acc;
      }
      
      // 2c) Inicializar us en la esfera ||u|| = sqrt(2*y[x_alpha])
      double Rlen = sqrt(2.0 * y[x_alpha]);
      us[0] = Rlen;
      for (int j = 1; j < vloc; ++j) us[j] = 0.0;
      
      // 2d) Metropolis–Hastings (nRepL iteraciones)
      for (int rep = 0; rep < nRepL; ++rep) {
        sample_sphere(Rlen, vloc, u_sph);
        
        double fu1 = mvn_iso_pdf(us,   u_mean, var_alpha, vloc);
        double fu2 = mvn_iso_pdf(u_sph,u_mean, var_alpha, vloc);
        
        if (unif_rand() * fu1 < fu2) {
          // aceptar
          memcpy(us, u_sph, (size_t)vloc * sizeof(double));
        }
      }
      
      // 2e) Actualizar U (fila x_alpha en cada columna j)
      for (int j = 0; j < vloc; ++j) {
        U[IDX(x_alpha, j, nloc)] = us[j];
      }
    } // fin barrido
  } // fin iteraciones
  
  PutRNGstate();
}
