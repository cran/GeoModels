#include "header.h"


void krige_from_inv(double *Inv_mat, int *n, int *iii, double *mean, double *var){
  int j, count = 0;
  double diag_val = Inv_mat[(*iii-1)*(*n) + (*iii-1)];
  
  for (j = 0; j < *n; j++) {
    if (j != (*iii-1)) {
      mean[count] = -Inv_mat[(*iii-1)*(*n) + j] / diag_val;
      count++;
    }
  }
  *var = 1.0 / diag_val;
}

// Function to generate multivariate normal with diagonal covariance
// For 2D case with diagonal covariance matrix
void mvrnorm_diag_2d(double *mu, double *sigma_diag, double *result) {
  // Generate two independent normal variates
  result[0] = mu[0] + sigma_diag[0] * norm_rand();
  result[1] = mu[1] + sigma_diag[0] * norm_rand();
}

// Main projection function
void rnorm_constraint_simple(double *A, double *b, double *mu, 
                             double *sigma, double *result) {
  double a1 = A[0], a2 = A[1];
  double mu1 = mu[0], mu2 = mu[1];
  double sig = *sigma;
  
  GetRNGstate();
  // Calculate the unit normal vector n.vec = A / ||A||
  double A_norm = sqrt(a1*a1 + a2*a2);
  
  if (A_norm == 0.0) {
    // Handle the case where A is zero vector
    result[0] = NAN;
    result[1] = NAN;
  PutRNGstate();
    return;
  }
  
  double n1 = a1 / A_norm;
  double n2 = a2 / A_norm;
  
  // Generate unconstrained sample from N(mu, sigma²I)
  double x1_sim = mu1 + sig * norm_rand();
  double x2_sim = mu2 + sig * norm_rand();
  
  // Calculate the projection: result = x.sim - n.vec %*% t(n.vec) %*% (x.sim - c(b,0))
  // Note: The R code uses c(b,0) but this seems incorrect for the constraint a1*x1 + a2*x2 = b
  // The correct translation point should be on the constraint line
  
  // Find the closest point on the constraint line to the origin that satisfies a1*x1 + a2*x2 = b
  double t = (*b) / (A_norm * A_norm);
  double proj_x1 = a1 * t;
  double proj_x2 = a2 * t;
  
  // Calculate the vector from projection point to sample
  double dx1 = x1_sim - proj_x1;
  double dx2 = x2_sim - proj_x2;
  
  // Calculate the dot product with normal vector
  double dot_product = dx1 * n1 + dx2 * n2;
  
  // Project orthogonally onto the constraint line
  result[0] = x1_sim - dot_product * n1;
  result[1] = x2_sim - dot_product * n2;
  PutRNGstate();  
}



void skew_gaussian_gibbs_sampler(double *data_obs, int *n, 
                                 double *Sigma_mat_inv, double *eta, 
                                 int *n_iter, double *data_x, double *data_y)
 {
  int i, k, ii, jj, accepted, trials;
  double c_i, mu_x, mu_y, sigma_val, w_plus;
  double *mean = (double *)malloc((*n-1)*sizeof(double));
  double var;
  double A_mat[2], mu_vec[2], result[2];
  
  // Controllo allocazione memoria
  if (!mean) {
    error("Errore di allocazione memoria per mean");
    return;
  }
  
  GetRNGstate();
  int *random_ind = (int *)malloc(*n * sizeof(int));
  
  if (!random_ind) {
    free(mean);
    error("Errore di allocazione memoria per random_ind");
    return;
  }
  
  for (jj = 0; jj < *n_iter; jj++) {
    // Initialize indices
    for (i = 0; i < *n; i++) random_ind[i] = i;
    
    // Fisher-Yates shuffle
    for (i = *n - 1; i > 0; i--) {
      int swap_idx = (int)(unif_rand() * (i + 1));
      int temp = random_ind[i];
      random_ind[i] = random_ind[swap_idx];
      random_ind[swap_idx] = temp;
    }
    
    for (k = 0; k < *n; k++) {
      ii = random_ind[k];
      c_i = data_obs[ii];
      
      // CORREZIONE: Indice Fortran corretto
      int iii = ii + 1; // Conversione a 1-based
      krige_from_inv(Sigma_mat_inv, n, &iii, mean, &var);
      
      // CORREZIONE: Calcolo corretto delle medie condizionali
      // mean[] contiene già i pesi/coefficienti dalla funzione krige_from_inv
      mu_x = 0.0;
      mu_y = 0.0;
      int count = 0;
      
      for (i = 0; i < *n; i++) {
        if (i != ii) {
          mu_x += mean[count] * data_x[i];
          mu_y += mean[count] * data_y[i];
          count++;
        }
      }
      
      // Controllo validità varianza
      if (var <= 0) {
        var = 1e-6; // Valore minimo per evitare problemi numerici
      }
      sigma_val = sqrt(var);
      
      // Calculate proposal weights
      w_plus = pnorm(0, mu_y, sigma_val, 0, 0);
      
      // Controllo che w_plus sia valido
      if (!R_finite(w_plus) || w_plus < 0 || w_plus > 1) {
        w_plus = 0.5; // Valore di default
      }
      
      // Prepare mu vector
      mu_vec[0] = mu_x;
      mu_vec[1] = mu_y;
      
      // Rejection sampling loop
      accepted = 0;
      trials = 0;
      const int MAX_TRIALS = 100; // Aumentato per maggiore robustezza
      
      while(!accepted && trials < MAX_TRIALS){
        trials++;
        
        if (unif_rand() < w_plus) {
          // Y_i >= 0 case: η₁X + η₂|Y| = c, with Y >= 0
          A_mat[0] = eta[0];
          A_mat[1] = eta[1];  
          
          rnorm_constraint_simple(A_mat, &c_i, mu_vec, &sigma_val, result);
          
          // Check if result is valid and Y >= 0
          if(R_finite(result[0]) && R_finite(result[1]) && result[1] >= 0){
            accepted = 1;
          }
        } else {
          // Y_i <= 0 case: η₁X + η₂|Y| = c, with Y <= 0
          A_mat[0] = eta[0];
          A_mat[1] = -eta[1]; // Negative because |Y| = -Y when Y <= 0
          
          rnorm_constraint_simple(A_mat, &c_i, mu_vec, &sigma_val, result);
          
          // Check if result is valid and Y <= 0
          if(R_finite(result[0]) && R_finite(result[1]) && result[1] <= 0){
            accepted = 1;
          }
        }
      }
      
      // Update if accepted, otherwise keep current values
      if (accepted) {
        data_x[ii] = result[0];
        data_y[ii] = result[1];
      }
      // Se non accettato, mantieni i valori correnti (implicito)
    }
  }
  
  free(random_ind);
  free(mean);
  PutRNGstate();
}

