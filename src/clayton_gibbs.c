#include "header.h"

#define IDX(i,j,nrow) ((i) + (j) * (nrow))

/*
 * Latent Gibbs sampler for the constructive Clayton-like random field
 *
 *   U_nu(s) = { A(s) / [A(s) + B(s)] }^(nu/2),
 *   A(s) = sum_{k=1}^nu Z_k(s)^2,
 *   B(s) = W_1(s)^2 + W_2(s)^2,
 *
 * where the nu + 2 latent fields are independent standard Gaussian random
 * fields sharing the same spatial correlation matrix.
 *
 * At an observed site, r_i = U_i^(2/nu) fixes only the ratio
 * A_i/(A_i+B_i), not the two radii separately.  Conditional on all other
 * sites, write
 *
 *   Z_i = q sqrt(r_i) a,   W_i = q sqrt(1-r_i) b,
 *
 * with ||a||=||b||=1.  The two directions have von Mises-Fisher full
 * conditionals and q has the one-dimensional log-concave density
 *
 *   p(q|.) propto q^(nu+1) exp{-q^2/(2s) + C q/s}, q>0.
 *
 * The routine updates U in place.  A finite number of sweeps gives the usual
 * MCMC approximation to the latent conditional distribution.
 */

static double vec_norm(const double *x, int p)
{
    double ss = 0.0;
    int j;
    for (j = 0; j < p; ++j)
        ss += x[j] * x[j];
    return sqrt(ss);
}

static void uniform_sphere(int p, double *out)
{
    double ss;
    int j;

    if (p <= 0)
        error("Sphere dimension must be positive");

    do {
        ss = 0.0;
        for (j = 0; j < p; ++j) {
            out[j] = norm_rand();
            ss += out[j] * out[j];
        }
    } while (!(ss > 0.0) || !R_FINITE(ss));

    ss = 1.0 / sqrt(ss);
    for (j = 0; j < p; ++j)
        out[j] *= ss;
}

/* Wood's rejection sampler for a vMF distribution on S^(p-1). */
static void sample_vmf(int p, const double *mu, double kappa,
                       double *out, double *work)
{
    int j, accepted;
    double a, prob_plus;
    double pm1, b, x0, c, z, w, log_accept;
    double tangent_scale, vv, denom;

    if (p <= 0)
        error("vMF dimension must be positive");

    if (!R_FINITE(kappa) || kappa < 0.0)
        error("Invalid vMF concentration");

    if (p == 1) {
        /* S^0 = {-1,+1}. */
        a = kappa * mu[0];
        if (a >= 0.0) {
            double e = exp(-2.0 * a);
            prob_plus = 1.0 / (1.0 + e);
        } else {
            double e = exp(2.0 * a);
            prob_plus = e / (1.0 + e);
        }
        out[0] = (unif_rand() < prob_plus) ? 1.0 : -1.0;
        return;
    }

    if (kappa < 1e-10) {
        uniform_sphere(p, out);
        return;
    }

    pm1 = (double) (p - 1);
    /* Rationalized form avoids cancellation when kappa is large. */
    denom = 2.0 * kappa + hypot(2.0 * kappa, pm1);
    b = pm1 / denom;
    x0 = (1.0 - b) / (1.0 + b);
    if (x0 >= 1.0)
        x0 = 1.0 - DBL_EPSILON;
    c = kappa * x0 + pm1 * log1p(-x0 * x0);

    accepted = 0;
    while (!accepted) {
        z = rbeta(0.5 * pm1, 0.5 * pm1);
        w = (1.0 - (1.0 + b) * z) /
            (1.0 - (1.0 - b) * z);
        log_accept = kappa * w + pm1 * log1p(-x0 * w) - c;
        {
            double uu = unif_rand();
            if (uu <= 0.0) uu = DBL_MIN;
            if (log(uu) <= log_accept)
                accepted = 1;
        }
    }

    /* work[0:(p-2)] is a uniform direction on S^(p-2). */
    uniform_sphere(p - 1, work);
    tangent_scale = sqrt(fmax(0.0, 1.0 - w * w));
    out[0] = w;
    for (j = 1; j < p; ++j)
        out[j] = tangent_scale * work[j - 1];

    /* Householder map taking e1 to mu. */
    if (mu[0] < 1.0 - 1e-12) {
        vv = 2.0 * (1.0 - mu[0]);
        if (!(vv > 0.0) || !R_FINITE(vv))
            error("Invalid vMF mean direction");

        /* v = e1 - mu; out <- out - 2 v (v'out)/(v'v). */
        {
            double dot = (1.0 - mu[0]) * out[0];
            for (j = 1; j < p; ++j)
                dot += (-mu[j]) * out[j];
            dot *= 2.0 / vv;
            out[0] -= dot * (1.0 - mu[0]);
            for (j = 1; j < p; ++j)
                out[j] += dot * mu[j];
        }
    }
}

static double log_q_density(double q, double power, double C, double s)
{
    if (!(q > 0.0) || !R_FINITE(q))
        return R_NegInf;
    return power * log(q) - 0.5 * q * q / s + C * q / s;
}

/* Exact one-dimensional slice update for the log-concave q full conditional. */
static double sample_q_slice(double current, double power, double C, double s)
{
    double logy, disc, mode;
    double lo, hi, mid, left, right, step;
    int it;

    if (!(current > 0.0) || !R_FINITE(current))
        current = sqrt((power + 1.0) * s);

    {
        double uu = unif_rand();
        if (uu <= 0.0) uu = DBL_MIN;
        logy = log_q_density(current, power, C, s) + log(uu);
    }

    disc = hypot(C, 2.0 * sqrt(power * s));
    if (C >= 0.0)
        mode = 0.5 * (C + disc);
    else
        mode = (2.0 * power * s) / (disc - C);

    if (!(mode > 0.0) || !R_FINITE(mode))
        mode = sqrt(power * s);

    /* Left root of log f(q) = logy. */
    lo = 0.0;
    hi = fmin(current, mode);
    if (!(hi > 0.0))
        hi = mode;
    for (it = 0; it < 55; ++it) {
        mid = 0.5 * (lo + hi);
        if (log_q_density(mid, power, C, s) >= logy)
            hi = mid;
        else
            lo = mid;
    }
    left = hi;

    /* Right root.  The Gaussian tail guarantees finite bracketing. */
    lo = fmax(current, mode);
    step = fmax(sqrt(s), DBL_EPSILON * fmax(1.0, lo));
    hi = lo + step;
    for (it = 0; it < 100 &&
                 log_q_density(hi, power, C, s) >= logy; ++it) {
        step *= 2.0;
        hi = lo + step;
    }
    if (it == 100 || !R_FINITE(hi))
        error("Unable to bracket the Clayton latent radial full conditional");

    for (it = 0; it < 55; ++it) {
        mid = 0.5 * (lo + hi);
        if (log_q_density(mid, power, C, s) >= logy)
            lo = mid;
        else
            hi = mid;
    }
    right = lo;

    if (!(right >= left) || !R_FINITE(left) || !R_FINITE(right))
        error("Invalid Clayton latent radial slice");

    return left + (right - left) * unif_rand();
}

void clayton_gibbs_sampler(double *SigmaInv,
                           double *ratio,
                           int *n,
                           int *nu,
                           double *U,
                           int *nIter)
{
    int N = *n;
    int V = *nu;
    int D = V + 2;
    int niter = *nIter;
    int i, j, k, iter, pos, swap_index, tmp;
    int *order;
    double qii, s, r, sr, s1r, current_q;
    double mz_norm, mw_norm, kappa, C, power;
    double *m, *a, *b, *mu_a, *mu_b, *work;

    if (N <= 0)
        error("n must be positive");
    if (V <= 0)
        error("Clayton nu must be a positive integer");
    if (niter <= 0)
        error("n_iter must be positive");

    for (i = 0; i < N; ++i) {
        if (!R_FINITE(ratio[i]) || ratio[i] <= 0.0 || ratio[i] >= 1.0)
            error("Clayton latent ratios must lie strictly between 0 and 1");
        qii = SigmaInv[IDX(i, i, N)];
        if (!R_FINITE(qii) || qii <= 0.0)
            error("Non-positive diagonal element in the latent precision matrix");
    }

    order = (int *) R_alloc((size_t) N, sizeof(int));
    m = (double *) R_alloc((size_t) D, sizeof(double));
    a = (double *) R_alloc((size_t) V, sizeof(double));
    b = (double *) R_alloc(2, sizeof(double));
    mu_a = (double *) R_alloc((size_t) V, sizeof(double));
    mu_b = (double *) R_alloc(2, sizeof(double));
    work = (double *) R_alloc((size_t) (V > 2 ? V : 2), sizeof(double));

    power = (double) (V + 1); /* total dimension D=V+2 => q^(D-1) */

    GetRNGstate();

    for (iter = 0; iter < niter; ++iter) {
        for (i = 0; i < N; ++i)
            order[i] = i;
        for (i = N - 1; i > 0; --i) {
            swap_index = (int) floor(unif_rand() * (double) (i + 1));
            if (swap_index > i)
                swap_index = i;
            tmp = order[i];
            order[i] = order[swap_index];
            order[swap_index] = tmp;
        }

        for (pos = 0; pos < N; ++pos) {
            i = order[pos];
            qii = SigmaInv[IDX(i, i, N)];
            s = 1.0 / qii;

            /* Conditional mean for each independent latent Gaussian field. */
            for (j = 0; j < D; ++j) {
                double acc = 0.0;
                for (k = 0; k < N; ++k) {
                    if (k == i)
                        continue;
                    acc += SigmaInv[IDX(i, k, N)] * U[IDX(k, j, N)];
                }
                m[j] = -acc / qii;
            }

            r = ratio[i];
            sr = sqrt(r);
            s1r = sqrt(1.0 - r);

            current_q = 0.0;
            for (j = 0; j < D; ++j)
                current_q += U[IDX(i, j, N)] * U[IDX(i, j, N)];
            current_q = sqrt(current_q);
            if (!(current_q > 0.0) || !R_FINITE(current_q))
                current_q = sqrt((double) D);

            mz_norm = vec_norm(m, V);
            if (mz_norm > 0.0 && R_FINITE(mz_norm)) {
                for (j = 0; j < V; ++j)
                    mu_a[j] = m[j] / mz_norm;
                kappa = current_q * sr * mz_norm / s;
                sample_vmf(V, mu_a, kappa, a, work);
            } else {
                uniform_sphere(V, a);
            }
            /*
             * Numerical guard: a is theoretically on the unit sphere, but
             * very concentrated vMF draws can accumulate a small floating
             * point norm error after the Householder map.  Since the Clayton
             * observation fixes the *ratio* of the two group radii exactly,
             * renormalize the direction before applying that radius.
             */
            {
                double an = vec_norm(a, V);
                if (!(an > 0.0) || !R_FINITE(an))
                    error("Invalid Clayton latent Z direction");
                for (j = 0; j < V; ++j)
                    a[j] /= an;
            }

            mw_norm = hypot(m[V], m[V + 1]);
            if (mw_norm > 0.0 && R_FINITE(mw_norm)) {
                mu_b[0] = m[V] / mw_norm;
                mu_b[1] = m[V + 1] / mw_norm;
                kappa = current_q * s1r * mw_norm / s;
                sample_vmf(2, mu_b, kappa, b, work);
            } else {
                uniform_sphere(2, b);
            }
            {
                double bn = hypot(b[0], b[1]);
                if (!(bn > 0.0) || !R_FINITE(bn))
                    error("Invalid Clayton latent W direction");
                b[0] /= bn;
                b[1] /= bn;
            }

            C = 0.0;
            for (j = 0; j < V; ++j)
                C += sr * a[j] * m[j];
            C += s1r * (b[0] * m[V] + b[1] * m[V + 1]);

            current_q = sample_q_slice(current_q, power, C, s);

            for (j = 0; j < V; ++j)
                U[IDX(i, j, N)] = current_q * sr * a[j];
            U[IDX(i, V, N)] = current_q * s1r * b[0];
            U[IDX(i, V + 1, N)] = current_q * s1r * b[1];
        }
    }

    /*
     * Final projection onto the observed ratio manifold.  In exact arithmetic
     * every update is already on this manifold.  This projection only removes
     * round-off accumulated in long Gibbs runs; it preserves the total radius
     * and both latent directions.
     */
    for (i = 0; i < N; ++i) {
        double zn2 = 0.0, wn2 = 0.0, zn, wn, qtot;
        r = ratio[i];
        sr = sqrt(r);
        s1r = sqrt(1.0 - r);

        for (j = 0; j < V; ++j)
            zn2 += U[IDX(i, j, N)] * U[IDX(i, j, N)];
        wn2 = U[IDX(i, V, N)] * U[IDX(i, V, N)] +
              U[IDX(i, V + 1, N)] * U[IDX(i, V + 1, N)];
        zn = sqrt(zn2);
        wn = sqrt(wn2);
        qtot = hypot(zn, wn);

        if (!(zn > 0.0) || !(wn > 0.0) || !(qtot > 0.0) ||
            !R_FINITE(zn) || !R_FINITE(wn) || !R_FINITE(qtot))
            error("Invalid Clayton latent state before final constraint projection");

        for (j = 0; j < V; ++j)
            U[IDX(i, j, N)] = qtot * sr * U[IDX(i, j, N)] / zn;
        U[IDX(i, V, N)] = qtot * s1r * U[IDX(i, V, N)] / wn;
        U[IDX(i, V + 1, N)] = qtot * s1r * U[IDX(i, V + 1, N)] / wn;
    }

    PutRNGstate();
}
