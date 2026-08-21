#include "header.h"

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


static double rtruncnorm_positive_exact(double mu, double sd)
{
    double a, log_surv_a, u, log_surv_z, z;

    a = -mu / sd;

    /*
     * log_surv_a = log P{N(0,1) >= a}.
     */
    log_surv_a = pnorm5(a, 0.0, 1.0, 0, 1);

    u = unif_rand();
    if (u <= 0.0)
        u = DBL_MIN;

    /*
     * If Z is standard normal truncated to [a,+Inf), then
     *
     * P(Z >= z) = U P(Z >= a).
     */
    log_surv_z = log_surv_a + log(u);
    z = qnorm5(log_surv_z, 0.0, 1.0, 0, 1);

    return mu + sd * z;
}


/*
 * Draw from N(mu,sd^2) truncated according to the selected branch:
 *
 * branch = +1: X >= 0;
 * branch = -1: X <= 0.
 */
static double rtruncnorm_branch_exact(double mu, double sd, int branch)
{
    if (branch > 0)
        return rtruncnorm_positive_exact(mu, sd);

    /*
     * X <= 0 is equivalent to -X >= 0.
     */
    return -rtruncnorm_positive_exact(-mu, sd);
}


/*
 * Exact single-site Gibbs sampler for
 *
 *     z_i = gamma |G1_i| + G2_i,
 *
 * where
 *
 *     G1 ~ N_n(0,Sigma),
 *     G2 ~ N_n(0,Sigma),
 *
 * independently, and Sigma_mat_inv = Sigma^{-1}.
 *
 * Historical .C interface:
 *
 * eta[0] = gamma;
 * eta[1] is ignored and is retained only to preserve compatibility
 * with the existing R-to-C interface.
 *
 * data_x contains G1;
 * data_y contains G2.
 */
void skew_gaussian_gibbs_sampler(double *data_obs,
                                 int *n,
                                 double *Sigma_mat_inv,
                                 double *eta,
                                 int *n_iter,
                                 double *data_x,
                                 double *data_y)
{
    int N, niter;
    int iter, k, i, j;
    int swap_index, tmp, branch;
    int *random_ind;

    double gamma, one_plus_gamma2;
    double qii, var, sd_linear, sd_trunc;
    double sum_x, sum_y, mu_x, mu_y;
    double m_plus, m_minus;
    double mu_trunc_plus, mu_trunc_minus;
    double log_q_plus, log_q_minus;
    double prob_plus;
    double x_new, y_new;
    double log_ratio;

    N = *n;
    niter = *n_iter;
    gamma = eta[0];

    if (N <= 0)
        error("n must be positive");

    if (niter <= 0)
        error("n_iter must be positive");

    if (!R_FINITE(gamma))
        error("The latent skew coefficient gamma must be finite");

    /*
     * Validate quantities before opening the RNG state.
     */
    for (i = 0; i < N; i++) {
        if (!R_FINITE(data_obs[i]))
            error("Non-finite latent observation");

        qii = Sigma_mat_inv[i + N * i];

        if (!R_FINITE(qii) || qii <= 0.0)
            error("Non-positive diagonal element in the precision matrix");
    }

    random_ind = (int *) R_alloc((size_t) N, sizeof(int));

    one_plus_gamma2 = 1.0 + gamma * gamma;

    GetRNGstate();

    for (iter = 0; iter < niter; iter++) {

        /*
         * Random ordering within each Gibbs sweep.
         */
        for (i = 0; i < N; i++)
            random_ind[i] = i;

        for (i = N - 1; i > 0; i--) {
            swap_index = (int) floor(unif_rand() * (i + 1));

            if (swap_index > i)
                swap_index = i;

            tmp = random_ind[i];
            random_ind[i] = random_ind[swap_index];
            random_ind[swap_index] = tmp;
        }

        for (k = 0; k < N; k++) {

            i = random_ind[k];
            qii = Sigma_mat_inv[i + N * i];

            /*
             * Conditional distributions:
             *
             * G1_i | G1_-i ~ N(mu_x,var),
             * G2_i | G2_-i ~ N(mu_y,var).
             */
            sum_x = 0.0;
            sum_y = 0.0;

            for (j = 0; j < N; j++) {
                if (j == i)
                    continue;

                /*
                 * R matrices are stored column-wise:
                 * Sigma_mat_inv[i + N*j] = Q_{ij}.
                 */
                sum_x += Sigma_mat_inv[i + N * j] * data_x[j];
                sum_y += Sigma_mat_inv[i + N * j] * data_y[j];
            }

            var = 1.0 / qii;
            mu_x = -sum_x / qii;
            mu_y = -sum_y / qii;

            /*
             * For branch s in {-1,+1},
             *
             * L_s = G2_i + s gamma G1_i.
             *
             * L_s | rest ~ N(m_s,var(1+gamma^2)).
             */
            sd_linear = sqrt(var * one_plus_gamma2);
            sd_trunc = sqrt(var / one_plus_gamma2);

            m_plus = mu_y + gamma * mu_x;
            m_minus = mu_y - gamma * mu_x;

            /*
             * Conditional means of G1_i given L_s = z_i.
             */
            mu_trunc_plus =
                mu_x
                + gamma / one_plus_gamma2
                * (data_obs[i] - m_plus);

            mu_trunc_minus =
                mu_x
                - gamma / one_plus_gamma2
                * (data_obs[i] - m_minus);

            /*
             * Exact unnormalized branch probabilities:
             *
             * q_s =
             * f_{L_s}(z_i)
             * P(s G1_i >= 0 | L_s = z_i).
             */
            log_q_plus =
                dnorm4(data_obs[i], m_plus, sd_linear, 1)
                + pnorm5(mu_trunc_plus / sd_trunc,
                         0.0, 1.0, 1, 1);

            log_q_minus =
                dnorm4(data_obs[i], m_minus, sd_linear, 1)
                + pnorm5(-mu_trunc_minus / sd_trunc,
                         0.0, 1.0, 1, 1);

            /*
             * Stable evaluation of
             *
             * q_plus / (q_plus + q_minus).
             */
            if (!R_FINITE(log_q_plus) &&
                !R_FINITE(log_q_minus)) {

                prob_plus = 0.5;

            } else if (!R_FINITE(log_q_minus)) {

                prob_plus = 1.0;

            } else if (!R_FINITE(log_q_plus)) {

                prob_plus = 0.0;

            } else {

                log_ratio = log_q_minus - log_q_plus;

                if (log_ratio > 40.0)
                    prob_plus = 0.0;
                else if (log_ratio < -40.0)
                    prob_plus = 1.0;
                else
                    prob_plus = 1.0 / (1.0 + exp(log_ratio));
            }

            branch = (unif_rand() < prob_plus) ? 1 : -1;

            /*
             * Rejection-free sampling from the appropriate truncated
             * Gaussian distribution.
             */
            if (branch > 0) {
                x_new =
                    rtruncnorm_branch_exact(mu_trunc_plus,
                                            sd_trunc,
                                            1);
            } else {
                x_new =
                    rtruncnorm_branch_exact(mu_trunc_minus,
                                            sd_trunc,
                                            -1);
            }

            /*
             * Enforce the nonlinear restriction exactly:
             *
             * z_i = gamma |G1_i| + G2_i
             *     = G2_i + branch gamma G1_i.
             */
            y_new =
                data_obs[i]
                - ((double) branch) * gamma * x_new;

            data_x[i] = x_new;
            data_y[i] = y_new;
        }
    }

    PutRNGstate();

}
