#include "header.h"
#include <limits.h>

/*
 * Count-constrained latent Gaussian Gibbs sampler used by GeoSimcond for
 * direct Binomial and Negative-Binomial random fields.
 *
 * The observed count at site i constrains the signs of independent latent
 * Gaussian copies Z_l(s_i).  A site block is sampled exactly by first drawing
 * the sign configuration from its conditional-Bernoulli law and then drawing
 * independent univariate truncated Gaussian magnitudes.  Prediction is done
 * inside the same native call, avoiding storage of L x m x nrep latent arrays.
 */

static double cc_logspace_add(double a, double b)
{
    if (!R_FINITE(a)) return b;
    if (!R_FINITE(b)) return a;
    if (a < b) {
        double tmp = a; a = b; b = tmp;
    }
    return a + log1p(exp(b - a));
}

/* Standard normal truncated below at alpha. */
static double rnorm_lower_trunc(double alpha)
{
    if (alpha <= 0.0) {
        double z;
        do { z = norm_rand(); } while (z <= alpha);
        return z;
    }

    /* Robert's exponential rejection sampler for a positive lower bound. */
    {
        double rate = 0.5 * (alpha + sqrt(alpha * alpha + 4.0));
        for (;;) {
            double z = alpha + exp_rand() / rate;
            double logacc = -0.5 * (z - rate) * (z - rate);
            if (log(unif_rand()) <= logacc) return z;
        }
    }
}

static double rnorm_sign_trunc(double mean, double sd, int positive)
{
    double alpha = -mean / sd;
    double z;
    if (positive) {
        z = rnorm_lower_trunc(alpha);
    } else {
        z = -rnorm_lower_trunc(-alpha);
    }
    return mean + sd * z;
}

/*
 * Draw q Bernoulli signs conditional on exactly c successes.
 * logp/logq contain log P(S=1) and log P(S=0).  To reduce work, when c>q/2
 * the dynamic program samples the q-c failures and then complements the draw.
 */
static void conditional_bernoulli_draw(int q, int c,
                                       const double *logp,
                                       const double *logq,
                                       int *signs,
                                       double *logF,
                                       int max_k)
{
    int sample_failures = (c > q - c);
    int target = sample_failures ? (q - c) : c;
    int stride = max_k + 1;
    int j, k;

    if (target == 0) {
        int success_value = sample_failures ? 1 : 0;
        for (j = 0; j < q; ++j) signs[j] = success_value;
        return;
    }

    for (j = 0; j <= q; ++j) {
        for (k = 0; k <= target; ++k) {
            logF[j * stride + k] = R_NegInf;
        }
    }
    logF[0] = 0.0;

    for (j = 1; j <= q; ++j) {
        double la = sample_failures ? logq[j - 1] : logp[j - 1];
        double lb = sample_failures ? logp[j - 1] : logq[j - 1];
        int kmax = (j < target) ? j : target;
        for (k = 0; k <= kmax; ++k) {
            double t0 = logF[(j - 1) * stride + k] + lb;
            double t1 = R_NegInf;
            if (k > 0) t1 = logF[(j - 1) * stride + (k - 1)] + la;
            logF[j * stride + k] = cc_logspace_add(t0, t1);
        }
    }

    k = target;
    for (j = q; j >= 1; --j) {
        int chosen;
        if (k == 0) {
            chosen = 0;
        } else if (k == j) {
            chosen = 1;
        } else {
            double la = sample_failures ? logq[j - 1] : logp[j - 1];
            double lp = la + logF[(j - 1) * stride + (k - 1)]
                        - logF[j * stride + k];
            double pr = exp(lp);
            if (!R_FINITE(pr)) pr = 0.0;
            if (pr < 0.0) pr = 0.0;
            if (pr > 1.0) pr = 1.0;
            chosen = (unif_rand() < pr);
        }
        if (chosen) --k;
        signs[j - 1] = sample_failures ? !chosen : chosen;
    }
}

static void latent_sweep(double *Z, double *H,
                         const double *Q, const double *mu,
                         const int *free_len, const int *success_count,
                         const int *forced_success,
                         int L, int m,
                         double *cond_mean, double *delta,
                         double *logp, double *logq, int *signs,
                         double *logF, int max_k)
{
    int i, l;
    int inc = 1;
    double one = 1.0;

    for (i = 0; i < m; ++i) {
        int q = free_len[i];
        int c = success_count[i];
        int forced = forced_success[i] - 1; /* -1 means none */
        double qii = Q[i + (size_t)m * i];
        double sd = sqrt(1.0 / qii);
        double *Zi = Z + (size_t)L * i;
        double *Hi = H + (size_t)L * i;

        for (l = 0; l < L; ++l) {
            cond_mean[l] = Zi[l] - Hi[l] / qii;
            delta[l] = 0.0;
        }

        if (q > 0) {
            for (l = 0; l < q; ++l) {
                double zeta = cond_mean[l] / sd;
                logp[l] = pnorm5(zeta, 0.0, 1.0, 1, 1);
                logq[l] = pnorm5(zeta, 0.0, 1.0, 0, 1);
            }
            conditional_bernoulli_draw(q, c, logp, logq, signs, logF, max_k);
        }

        for (l = 0; l < L; ++l) {
            double oldz = Zi[l];
            double newz;
            if (l < q) {
                newz = rnorm_sign_trunc(cond_mean[l], sd, signs[l]);
            } else if (l == forced) {
                newz = rnorm_sign_trunc(cond_mean[l], sd, 1);
            } else {
                newz = cond_mean[l] + sd * norm_rand();
            }
            Zi[l] = newz;
            delta[l] = newz - oldz;
        }

        /* H = (Z-mu) Q.  Updating site i changes H by delta %*% Q[i,].
         * Q is symmetric, so its i-th row equals the contiguous i-th column. */
        F77_CALL(dger)(&L, &m, &one, delta, &inc,
                       (double *)(Q + (size_t)m * i), &inc,
                       H, &L);

        if ((i & 31) == 0) R_CheckUserInterrupt();
    }
}

static void mvn_add_upper_chol(double *x, const double *chol, int p)
{
    char uplo = 'U', trans = 'T', diag = 'N';
    int inc = 1;
    F77_CALL(dtrmv)(&uplo, &trans, &diag, &p,
                    (double *)chol, &p, x, &inc FCONE FCONE FCONE);
}

/*
 * .Call interface.
 * model_code: 0 = Binomial, 1 = Negative Binomial.
 */
typedef struct {
    SEXP args[15];
    int rng_active;
} GeoCountSimContext;

static void geo_count_rng_cleanup(void *data)
{
    GeoCountSimContext *ctx = (GeoCountSimContext *) data;
    if (ctx->rng_active) {
        ctx->rng_active = 0;
        PutRNGstate();
    }
}

static SEXP geo_count_cond_sim_body(void *data)
{
    GeoCountSimContext *ctx = (GeoCountSimContext *) data;
    SEXP QSEXP = ctx->args[0];
    SEXP muObsSEXP = ctx->args[1];
    SEXP freeLenSEXP = ctx->args[2];
    SEXP successCountSEXP = ctx->args[3];
    SEXP forcedSuccessSEXP = ctx->args[4];
    SEXP WSEXP = ctx->args[5];
    SEXP muLocSEXP = ctx->args[6];
    SEXP cholCondSEXP = ctx->args[7];
    SEXP cholLocSEXP = ctx->args[8];
    SEXP sizePredSEXP = ctx->args[9];
    SEXP modelCodeSEXP = ctx->args[10];
    SEXP burnSEXP = ctx->args[11];
    SEXP thinSEXP = ctx->args[12];
    SEXP nrepSEXP = ctx->args[13];
    SEXP LSEXP = ctx->args[14];

    SEXP qdim, wdim, ans;
    int m, p, L, model_code, burn, thin, nrep;
    int i, j, l, rep, sweep;
    int max_q = 0, max_k = 0;
    int max_pred_size = 0;
    int r_nb = 0;
    double *Q, *mu, *W, *muLoc, *cholCond, *cholLoc, *out;
    int *free_len, *success_count, *forced_success, *size_pred;

    if (!isReal(QSEXP) || !isMatrix(QSEXP)) error("Q must be a numeric matrix.");
    if (!isReal(WSEXP) || !isMatrix(WSEXP)) error("W must be a numeric matrix.");

    qdim = getAttrib(QSEXP, R_DimSymbol);
    wdim = getAttrib(WSEXP, R_DimSymbol);
    m = INTEGER(qdim)[0];
    if (INTEGER(qdim)[1] != m) error("Q must be square.");
    if (INTEGER(wdim)[0] != m) error("W has incompatible dimensions.");
    p = INTEGER(wdim)[1];

    L = asInteger(LSEXP);
    model_code = asInteger(modelCodeSEXP);
    burn = asInteger(burnSEXP);
    thin = asInteger(thinSEXP);
    nrep = asInteger(nrepSEXP);
    if (L < 1 || p < 1 || m < 1 || burn < 0 || thin < 1 || nrep < 1)
        error("Invalid dimensions or MCMC controls in GeoCountCondSim.");
    if (!(model_code == 0 || model_code == 1)) error("Unknown count model code.");

    if (XLENGTH(muObsSEXP) != m || XLENGTH(freeLenSEXP) != m ||
        XLENGTH(successCountSEXP) != m || XLENGTH(forcedSuccessSEXP) != m)
        error("Observation-side vectors have incompatible lengths.");
    if (XLENGTH(muLocSEXP) != p || XLENGTH(sizePredSEXP) != p)
        error("Prediction-side vectors have incompatible lengths.");
    if (XLENGTH(cholCondSEXP) != (R_xlen_t)p * p ||
        XLENGTH(cholLocSEXP) != (R_xlen_t)p * p)
        error("Prediction Cholesky factors have incompatible dimensions.");

    Q = REAL(QSEXP);
    mu = REAL(muObsSEXP);
    W = REAL(WSEXP);
    muLoc = REAL(muLocSEXP);
    cholCond = REAL(cholCondSEXP);
    cholLoc = REAL(cholLocSEXP);
    free_len = INTEGER(freeLenSEXP);
    success_count = INTEGER(successCountSEXP);
    forced_success = INTEGER(forcedSuccessSEXP);
    size_pred = INTEGER(sizePredSEXP);

    for (i = 0; i < m; ++i) {
        int q = free_len[i], c = success_count[i];
        if (q < 0 || q > L || c < 0 || c > q)
            error("Invalid count constraint supplied to GeoCountCondSim.");
        if (forced_success[i] < 0 || forced_success[i] > L)
            error("Invalid forced-success index supplied to GeoCountCondSim.");
        if (forced_success[i] > 0 && forced_success[i] != q + 1)
            error("A forced success must immediately follow the free count-constrained trials.");
        if (q > max_q) max_q = q;
        {
            int kk = (c < q - c) ? c : q - c;
            if (kk > max_k) max_k = kk;
        }
        if (!(Q[i + (size_t)m * i] > 0.0) || !R_FINITE(Q[i + (size_t)m * i]))
            error("The latent precision matrix has an invalid diagonal.");
    }

    for (j = 0; j < p; ++j) {
        if (size_pred[j] < 1) error("Prediction sizes must be positive integers.");
        if (size_pred[j] > max_pred_size) max_pred_size = size_pred[j];
    }
    if (model_code == 1) {
        r_nb = size_pred[0];
        for (j = 1; j < p; ++j)
            if (size_pred[j] != r_nb)
                error("Negative-Binomial prediction requires one common r.");
    }

    {
        double *Z = (double *) R_alloc((size_t)L * m, sizeof(double));
        double *H = (double *) R_alloc((size_t)L * m, sizeof(double));
        double *centered = (double *) R_alloc((size_t)L * m, sizeof(double));
        double *cond_mean = (double *) R_alloc(L, sizeof(double));
        double *delta = (double *) R_alloc(L, sizeof(double));
        double *logp = (double *) R_alloc((max_q > 0 ? max_q : 1), sizeof(double));
        double *logq = (double *) R_alloc((max_q > 0 ? max_q : 1), sizeof(double));
        int *signs = (int *) R_alloc((max_q > 0 ? max_q : 1), sizeof(int));
        double *logF = (double *) R_alloc((size_t)(max_q + 1) * (max_k + 1), sizeof(double));
        int *perm = (int *) R_alloc((max_q > 0 ? max_q : 1), sizeof(int));
        double *pred_mean = (double *) R_alloc((size_t)L * p, sizeof(double));
        double *base_pred = (double *) R_alloc(p, sizeof(double));
        double *work_p = (double *) R_alloc(p, sizeof(double));
        int *succ = (int *) R_alloc(p, sizeof(int));
        int *done = (int *) R_alloc(p, sizeof(int));

        PROTECT(ans = allocMatrix(REALSXP, p, nrep));
        out = REAL(ans);

        for (j = 0; j < p; ++j) {
            base_pred[j] = muLoc[j];
            for (i = 0; i < m; ++i)
                base_pred[j] -= mu[i] * W[i + (size_t)m * j];
        }

        GetRNGstate();
        ctx->rng_active = 1;

        /* Random valid initialization. */
        for (i = 0; i < m; ++i) {
            int q = free_len[i];
            int c = success_count[i];
            int forced = forced_success[i] - 1;
            double *Zi = Z + (size_t)L * i;
            for (l = 0; l < L; ++l) Zi[l] = mu[i] + norm_rand();
            for (l = 0; l < q; ++l) perm[l] = l;
            for (l = q - 1; l > 0; --l) {
                int k = (int) floor(unif_rand() * (l + 1));
                int tmp = perm[l]; perm[l] = perm[k]; perm[k] = tmp;
            }
            for (l = 0; l < q; ++l) {
                int pos = 0;
                int k;
                for (k = 0; k < c; ++k) {
                    if (perm[k] == l) { pos = 1; break; }
                }
                Zi[l] = pos ? (fabs(Zi[l]) + 0.1) : -(fabs(Zi[l]) + 0.1);
            }
            if (forced >= 0) Zi[forced] = fabs(Zi[forced]) + 0.1;
        }

        /* H = (Z-mu) %*% Q. */
        for (i = 0; i < m; ++i) {
            for (l = 0; l < L; ++l) {
                centered[l + (size_t)L * i] = Z[l + (size_t)L * i] - mu[i];
            }
        }
        {
            char N = 'N';
            double one = 1.0, zero = 0.0;
            F77_CALL(dgemm)(&N, &N, &L, &m, &m,
                            &one, centered, &L, Q, &m,
                            &zero, H, &L FCONE FCONE);
        }

        for (sweep = 0; sweep < burn; ++sweep) {
            latent_sweep(Z, H, Q, mu, free_len, success_count, forced_success,
                         L, m, cond_mean, delta, logp, logq, signs, logF, max_k);
        }

        for (rep = 0; rep < nrep; ++rep) {
            int nsweeps = (rep == 0 && burn == 0) ? 1 : thin;
            for (sweep = 0; sweep < nsweeps; ++sweep) {
                latent_sweep(Z, H, Q, mu, free_len, success_count, forced_success,
                             L, m, cond_mean, delta, logp, logq, signs, logF, max_k);
            }

            /* pred_mean = Z %*% W + muLoc - muObs %*% W. */
            {
                char N = 'N';
                double one = 1.0, zero = 0.0;
                F77_CALL(dgemm)(&N, &N, &L, &p, &m,
                                &one, Z, &L, W, &m,
                                &zero, pred_mean, &L FCONE FCONE);
            }
            for (j = 0; j < p; ++j) {
                for (l = 0; l < L; ++l)
                    pred_mean[l + (size_t)L * j] += base_pred[j];
            }

            for (j = 0; j < p; ++j) {
                out[j + (size_t)p * rep] = 0.0;
                succ[j] = 0;
                done[j] = 0;
            }

            for (l = 0; l < L; ++l) {
                for (j = 0; j < p; ++j) work_p[j] = norm_rand();
                mvn_add_upper_chol(work_p, cholCond, p);

                if (model_code == 0) {
                    for (j = 0; j < p; ++j) {
                        if (l < size_pred[j] && pred_mean[l + (size_t)L * j] + work_p[j] > 0.0)
                            out[j + (size_t)p * rep] += 1.0;
                    }
                } else {
                    for (j = 0; j < p; ++j) if (!done[j]) {
                        if (pred_mean[l + (size_t)L * j] + work_p[j] > 0.0) {
                            succ[j] += 1;
                            if (succ[j] >= r_nb) {
                                out[j + (size_t)p * rep] = (double)((l + 1) - r_nb);
                                done[j] = 1;
                            }
                        }
                    }
                }
            }

            if (model_code == 0 && max_pred_size > L) {
                int trial;
                for (trial = L; trial < max_pred_size; ++trial) {
                    if ((trial & 1023) == 0) R_CheckUserInterrupt();
                    for (j = 0; j < p; ++j) work_p[j] = norm_rand();
                    mvn_add_upper_chol(work_p, cholLoc, p);
                    for (j = 0; j < p; ++j) {
                        if (trial < size_pred[j] && muLoc[j] + work_p[j] > 0.0)
                            out[j + (size_t)p * rep] += 1.0;
                    }
                }
            }

            if (model_code == 1) {
                int n_done = 0;
                int trial = L;
                for (j = 0; j < p; ++j) n_done += done[j];
                while (n_done < p) {
                    if ((trial & 1023) == 0) R_CheckUserInterrupt();
                    if ((double)trial > (double)L + 10000000.0 || trial == INT_MAX)
                        error("Negative-Binomial conditional simulation exceeded the safety trial limit.");
                    for (j = 0; j < p; ++j) work_p[j] = norm_rand();
                    mvn_add_upper_chol(work_p, cholLoc, p);
                    ++trial;
                    for (j = 0; j < p; ++j) if (!done[j]) {
                        if (muLoc[j] + work_p[j] > 0.0) {
                            succ[j] += 1;
                            if (succ[j] >= r_nb) {
                                out[j + (size_t)p * rep] = (double)(trial - r_nb);
                                done[j] = 1;
                                ++n_done;
                            }
                        }
                    }
                }
            }

            R_CheckUserInterrupt();
        }

        UNPROTECT(1);
    }

    return ans;
}

SEXP GeoCountCondSim(SEXP QSEXP, SEXP muObsSEXP,
                     SEXP freeLenSEXP, SEXP successCountSEXP,
                     SEXP forcedSuccessSEXP,
                     SEXP WSEXP, SEXP muLocSEXP,
                     SEXP cholCondSEXP, SEXP cholLocSEXP,
                     SEXP sizePredSEXP, SEXP modelCodeSEXP,
                     SEXP burnSEXP, SEXP thinSEXP, SEXP nrepSEXP,
                     SEXP LSEXP)
{
    GeoCountSimContext ctx = {{QSEXP, muObsSEXP, freeLenSEXP, successCountSEXP, forcedSuccessSEXP, WSEXP, muLocSEXP, cholCondSEXP, cholLocSEXP, sizePredSEXP, modelCodeSEXP, burnSEXP, thinSEXP, nrepSEXP, LSEXP}, 0};
    return R_ExecWithCleanup(geo_count_cond_sim_body, &ctx,
                             geo_count_rng_cleanup, &ctx);
}
