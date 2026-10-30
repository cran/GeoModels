/* Sparse taper setup for regular space-time and colocated bivariate data.
 * Two scans retain exactly the legacy t,i,v,j CSR order.  Workspace is O(nnz)
 * instead of O(N^2); the pair search itself remains quadratic. */
#include "header.h"

typedef struct {
    SEXP input;
    int started, success, nsite, ntimes, total, biv, metric, tapmodel;
    const double *x, *y, *z, *times;
    double radius, thresholds[2];
    double *out_idx;
    int *out_ja, *out_ia;
} GeoSparseTaperContext;

static SEXP geo_sparse_field(SEXP input, const char *name)
{
    SEXP names = getAttrib(input, R_NamesSymbol);
    if (TYPEOF(input) != VECSXP || TYPEOF(names) != STRSXP)
        error("GeoModels sparse taper setup requires a named list");
    for (R_xlen_t i = 0; i < XLENGTH(input); ++i)
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
            return VECTOR_ELT(input, i);
    error("GeoModels sparse taper setup is missing '%s'", name);
    return R_NilValue;
}

static int geo_sparse_integer(SEXP input, const char *name, int minimum)
{
    SEXP value = geo_sparse_field(input, name);
    if (XLENGTH(value) != 1) error("GeoModels: '%s' must be scalar", name);
    const double v = asReal(value);
    if (!R_FINITE(v) || v < minimum || v > INT_MAX || v != floor(v))
        error("GeoModels: invalid sparse taper '%s'", name);
    return (int)v;
}

static R_xlen_t geo_sparse_scan(GeoSparseTaperContext *ctx, int fill)
{
    R_xlen_t count = 0, visited = 0;
    const int n = ctx->nsite, nt = ctx->ntimes;
    if (fill) ctx->out_ia[0] = 1;
    for (int t = 0; t < nt; ++t) {
        for (int i = 0; i < n; ++i) {
            const int row = t * n + i;
            for (int v = 0; v < nt; ++v) {
                const double temporal = ctx->biv ? 0.0 : fabs(ctx->times[t]-ctx->times[v]);
                for (int j = 0; j < n; ++j) {
                    if ((visited++ & 8191) == 0) R_CheckUserInterrupt();
                    const double spatial = dist(ctx->metric, ctx->x[i], ctx->x[j],
                        ctx->y[i], ctx->y[j], ctx->z[i], ctx->z[j], ctx->radius);
                    int keep;
                    if (ctx->biv) keep = spatial <= dista[t][v];
                    else {
                        double support[2] = {R_NaN, R_NaN};
                        Comp_supp(support, &ctx->tapmodel, spatial, temporal, ctx->thresholds);
                        if (ISNAN(support[0]) || ISNAN(support[1]))
                            error("GeoModels: unsupported space-time taper code %d", ctx->tapmodel);
                        keep = (spatial < support[0] || is_equal(spatial, support[0])) &&
                               (temporal < support[1] || is_equal(temporal, support[1]));
                    }
                    if (!keep) continue;
                    if (count >= INT_MAX-1)
                        error("GeoModels sparse taper has too many nonzero entries for integer CSR indices");
                    if (fill) {
                        const int col = v * n + j;
                        ctx->out_ja[count] = col + 1;
                        /* Legacy row-major flattened order; the matrix is symmetric.
                         * Double indices remain exact under the checked 2^53 bound. */
                        ctx->out_idx[count] = (double)row * ctx->total + col + 1.0;
                        if (ctx->biv) {
                            lags_1[count] = spatial;
                            first_1[count] = j < i ? v : t;
                            second_1[count] = j < i ? t : v;
                        } else {
                            lags[count] = spatial;
                            lagt[count] = temporal;
                        }
                    }
                    ++count;
                }
            }
            if (fill) ctx->out_ia[row+1] = (int)count + 1;
        }
    }
    return count;
}

static SEXP geo_sparse_taper_body(void *data)
{
    GeoSparseTaperContext *ctx = (GeoSparseTaperContext *)data;
    SEXP input = ctx->input;
    int nprotect = 0;
    ctx->nsite = geo_sparse_integer(input, "numcoord", 1);
    ctx->ntimes = geo_sparse_integer(input, "numtime", 1);
    ctx->metric = geo_sparse_integer(input, "distance", 0);
    ctx->tapmodel = geo_sparse_integer(input, "tapmodel", 0);
    SEXP biv = geo_sparse_field(input, "bivariate");
    if (XLENGTH(biv) != 1 || asLogical(biv) == NA_LOGICAL)
        error("GeoModels: invalid sparse taper bivariate flag");
    ctx->biv = asLogical(biv);
    if (ctx->biv && ctx->ntimes != 2)
        error("GeoModels: bivariate sparse taper requires exactly two components");
    if (ctx->metric > 3) error("GeoModels: invalid sparse taper distance code");
    SEXP radius = geo_sparse_field(input, "radius");
    if (XLENGTH(radius) != 1 || !R_FINITE(ctx->radius = asReal(radius)) || ctx->radius <= 0)
        error("GeoModels: invalid sparse taper radius");
    const double total = (double)ctx->nsite * ctx->ntimes;
    if (total > INT_MAX-1 || total * total > 9007199254740991.0)
        error("GeoModels: sparse taper matrix exceeds the supported index range");
    ctx->total = (int)total;
    SEXP xx, yy, zz, tt, sr, tr, sp;
    PROTECT(xx = coerceVector(geo_sparse_field(input, "coordx"), REALSXP)); nprotect++;
    PROTECT(yy = coerceVector(geo_sparse_field(input, "coordy"), REALSXP)); nprotect++;
    PROTECT(zz = coerceVector(geo_sparse_field(input, "coordz"), REALSXP)); nprotect++;
    PROTECT(tt = coerceVector(geo_sparse_field(input, "coordt"), REALSXP)); nprotect++;
    PROTECT(sr = coerceVector(geo_sparse_field(input, "srange"), REALSXP)); nprotect++;
    PROTECT(tr = coerceVector(geo_sparse_field(input, "trange"), REALSXP)); nprotect++;
    PROTECT(sp = coerceVector(geo_sparse_field(input, "tapsep"), REALSXP)); nprotect++;
    if (XLENGTH(xx) < ctx->nsite || XLENGTH(yy) < ctx->nsite || XLENGTH(zz) < ctx->nsite ||
        (!ctx->biv && XLENGTH(tt) < ctx->ntimes) || XLENGTH(sr) < (ctx->biv ? 4 : 2) ||
        XLENGTH(tr) < 2 || XLENGTH(sp) > 5)
        error("GeoModels: inconsistent regular sparse taper dimensions");
    ctx->x = REAL(xx); ctx->y = REAL(yy); ctx->z = REAL(zz); ctx->times = REAL(tt);
    for (int i = 0; i < ctx->nsite; ++i)
        if (!R_FINITE(ctx->x[i]) || !R_FINITE(ctx->y[i]) || !R_FINITE(ctx->z[i]))
            error("GeoModels: sparse taper coordinates must be finite");
    if (!ctx->biv) {
        for (int i = 0; i < ctx->ntimes; ++i)
            if (!R_FINITE(ctx->times[i])) error("GeoModels: sparse taper times must be finite");
    }
    for (R_xlen_t i = 1; i < XLENGTH(sr); ++i)
        if (ISNAN(REAL(sr)[i]) || REAL(sr)[i] < 0)
            error("GeoModels: sparse taper spatial cutoffs must be nonnegative");
    if (ISNAN(REAL(tr)[1]) || REAL(tr)[1] < 0)
        error("GeoModels: sparse taper temporal cutoff must be nonnegative");
    for (R_xlen_t i = 0; i < XLENGTH(sp); ++i)
        if (!R_FINITE(REAL(sp)[i])) error("GeoModels: sparse taper parameters must be finite");

    DeleteGlobalVar();
    ctx->started = 1;
#define GEO_SPARSE_SCALAR(ptr, ctype, value) do { (ptr) = (ctype *)R_Calloc(1, ctype); (ptr)[0] = (value); } while (0)
    GEO_SPARSE_SCALAR(ncoord, int, ctx->nsite);
    GEO_SPARSE_SCALAR(ncoordx, int, ctx->nsite);
    GEO_SPARSE_SCALAR(ncoordy, int, 1);
    GEO_SPARSE_SCALAR(ncoordz, int, 1);
    GEO_SPARSE_SCALAR(ntime, int, ctx->ntimes);
    GEO_SPARSE_SCALAR(npairs, int, 0);
    GEO_SPARSE_SCALAR(isbiv, int, ctx->biv);
    GEO_SPARSE_SCALAR(isst, int, !ctx->biv);
    GEO_SPARSE_SCALAR(ismem, int, 1);
    GEO_SPARSE_SCALAR(istap, int, 1);
    GEO_SPARSE_SCALAR(cdyn, int, 0);
    GEO_SPARSE_SCALAR(type, int, ctx->metric);
    GEO_SPARSE_SCALAR(REARTH, double, ctx->radius);
    GEO_SPARSE_SCALAR(maxtime, double, REAL(tr)[1] ? REAL(tr)[1] : -LOW);
#undef GEO_SPARSE_SCALAR
    maxdist = (double *)R_Calloc(3, double);
    for (int i = 0; i < (ctx->biv ? 3 : 1); ++i)
        maxdist[i] = REAL(sr)[i+1] ? REAL(sr)[i+1] : -LOW;
    tapsep = (double *)R_Calloc(5, double);
    for (R_xlen_t i = 0; i < XLENGTH(sp); ++i) tapsep[i] = REAL(sp)[i];
    if (tapsep[0] == 1.0) tapsep[0] = 0.99999999;
    ctx->thresholds[0] = REAL(sr)[1]; ctx->thresholds[1] = REAL(tr)[1];
    if (ctx->biv) {
        dista = (double **)R_Calloc(2, double *);
        for (int i = 0; i < 2; ++i) dista[i] = (double *)R_Calloc(2, double);
        dista[0][0] = maxdist[0]; dista[0][1] = dista[1][0] = maxdist[1];
        dista[1][1] = maxdist[2];
    }
    const R_xlen_t count = geo_sparse_scan(ctx, 0);
    npairs[0] = (int)count;
    if (ctx->biv) {
        lags_1 = (double *)R_Calloc(count, double);
        first_1 = (int *)R_Calloc(count, int);
        second_1 = (int *)R_Calloc(count, int);
    } else {
        lags = (double *)R_Calloc(count, double);
        lagt = (double *)R_Calloc(count, double);
    }
    SEXP idx, ja, ia, ans, names;
    PROTECT(idx = allocVector(REALSXP, count)); nprotect++;
    PROTECT(ja = allocVector(INTSXP, count)); nprotect++;
    PROTECT(ia = allocVector(INTSXP, (R_xlen_t)ctx->total + 1)); nprotect++;
    ctx->out_idx = REAL(idx); ctx->out_ja = INTEGER(ja); ctx->out_ia = INTEGER(ia);
    if (geo_sparse_scan(ctx, 1) != count) error("GeoModels: inconsistent sparse taper pair count");
    PROTECT(ans = allocVector(VECSXP, 8)); nprotect++;
    PROTECT(names = allocVector(STRSXP, 8)); nprotect++;
    SET_VECTOR_ELT(ans, 0, idx); SET_VECTOR_ELT(ans, 1, ja); SET_VECTOR_ELT(ans, 2, ia);
    SET_VECTOR_ELT(ans, 3, ScalarInteger((int)count));
    SET_VECTOR_ELT(ans, 4, ScalarReal(count / (total * total)));
    SET_VECTOR_ELT(ans, 5, ja); SET_VECTOR_ELT(ans, 6, ia);
    SET_VECTOR_ELT(ans, 7, ScalarInteger(1));
    const char *nms[] = {"idx", "ja", "ia", "numpairs", "nozero", "colidx", "rowidx", "isinit"};
    for (int i = 0; i < 8; ++i) SET_STRING_ELT(names, i, mkChar(nms[i]));
    setAttrib(ans, R_NamesSymbol, names);
    ctx->success = 1;
    UNPROTECT(nprotect);
    return ans;
}

static void geo_sparse_taper_cleanup(void *data)
{
    GeoSparseTaperContext *ctx = (GeoSparseTaperContext *)data;
    if (ctx->started && !ctx->success) DeleteGlobalVar();
}

SEXP GeoSparseTaperSetup(SEXP input)
{
    GeoSparseTaperContext ctx = {0};
    ctx.input = input;
    return R_ExecWithCleanup(geo_sparse_taper_body, &ctx, geo_sparse_taper_cleanup, &ctx);
}
