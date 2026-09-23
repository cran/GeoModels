#include "header.h"
#include <limits.h>
#include <stdlib.h>
#include "correlation_param_count.h"

typedef struct { uint64_t key; double rho; } GeoPairEntry;

typedef struct {
    SEXP neighborhoods, coords, times, cormod, par, type, radius;
    uint64_t *keys, *next_keys;
    double *values, *next_values;
    GeoPairEntry *entries;
    size_t capacity;
    double budget, peak;
} GeoPairBuildContext;

static size_t geo_pair_hash(uint64_t x)
{
    x ^= x >> 33; x *= UINT64_C(0xff51afd7ed558ccd);
    x ^= x >> 33; x *= UINT64_C(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    return (size_t)x;
}

static size_t geo_next_pow2(size_t x)
{
    size_t p = 1;
    while (p < x) {
        if (p > (SIZE_MAX >> 1)) return 0;
        p <<= 1;
    }
    return p;
}

static int geo_pair_compare(const void *aa, const void *bb)
{
    const GeoPairEntry *a = (const GeoPairEntry *)aa, *b = (const GeoPairEntry *)bb;
    return (a->key > b->key) - (a->key < b->key);
}

static void geo_pair_cleanup(void *data)
{
    GeoPairBuildContext *ctx = (GeoPairBuildContext *)data;
    free(ctx->keys); free(ctx->values); free(ctx->entries);
    free(ctx->next_keys); free(ctx->next_values);
    ctx->keys = ctx->next_keys = NULL;
    ctx->values = ctx->next_values = NULL;
    ctx->entries = NULL;
}

/* Budget covers the native workspace and the returned CSR arrays, including
 * simultaneous buffers during rehash. Inputs and other R objects are excluded.
 * Reserve 64 KiB for small objects and allocator bookkeeping.
 */
static int geo_pair_budget_ok(GeoPairBuildContext *ctx, double bytes)
{
    bytes += 65536.0;
    if (!R_FINITE(bytes) || bytes > ctx->budget) return 0;
    if (bytes > ctx->peak) ctx->peak = bytes;
    return 1;
}

static int geo_pair_rehash(GeoPairBuildContext *ctx)
{
    const size_t old = ctx->capacity;
    if (old > SIZE_MAX / 2) return 0;
    const size_t cap = 2 * old;
    if (cap > SIZE_MAX / sizeof(uint64_t) || cap > SIZE_MAX / sizeof(double) ||
        !geo_pair_budget_ok(ctx, (double)(sizeof(uint64_t) + sizeof(double)) *
                                  ((double)old + (double)cap))) return 0;
    ctx->next_keys = (uint64_t *)calloc(cap, sizeof(uint64_t));
    ctx->next_values = (double *)calloc(cap, sizeof(double));
    if (!ctx->next_keys || !ctx->next_values) return 0;
    for (size_t s = 0; s < old; ++s) {
        if ((s & 65535u) == 0) R_CheckUserInterrupt();
        uint64_t key = ctx->keys[s];
        if (key) {
            size_t pos = geo_pair_hash(key) & (cap - 1);
            while (ctx->next_keys[pos]) pos = (pos + 1) & (cap - 1);
            ctx->next_keys[pos] = key;
            ctx->next_values[pos] = ctx->values[s];
        }
    }
    free(ctx->keys); free(ctx->values);
    ctx->keys = ctx->next_keys; ctx->values = ctx->next_values;
    ctx->next_keys = NULL; ctx->next_values = NULL;
    ctx->capacity = cap;
    return 1;
}

static int geo_index_at(SEXP x, R_xlen_t i)
{
    if (TYPEOF(x) == INTSXP) {
        int v = INTEGER(x)[i];
        if (v != NA_INTEGER && v >= 1) return v;
    } else if (TYPEOF(x) == REALSXP) {
        double v = REAL(x)[i];
        if (R_FINITE(v) && v >= 1.0 && v <= INT_MAX && v == floor(v)) return (int)v;
    }
    error("GeoModels pair cache requires finite positive integer indices");
    return 0;
}

static SEXP geo_build_pair_cache_body(void *data)
{
    GeoPairBuildContext *ctx = (GeoPairBuildContext *)data;
    SEXP neighborhoods = ctx->neighborhoods, coords = ctx->coords;
    SEXP dim, coords_r, times_r = R_NilValue, par_r;
    int nprotect = 0;
    if (!isNewList(neighborhoods))
        error("GeoModels pair cache requires a list of neighborhood indices");
    dim = getAttrib(coords, R_DimSymbol);
    if (TYPEOF(dim) != INTSXP || LENGTH(dim) != 2)
        error("GeoModels pair cache requires a coordinate matrix");
    const int nobs = INTEGER(dim)[0], ndim = INTEGER(dim)[1];
    if (nobs < 1 || (ndim != 2 && ndim != 3))
        error("GeoModels pair cache requires N x 2 or N x 3 coordinates");
    PROTECT(coords_r = coerceVector(coords, REALSXP)); nprotect++;
    const double *cc = REAL(coords_r), *ttimes = NULL;
    for (R_xlen_t i = 0; i < XLENGTH(coords_r); ++i) {
        if ((i & 65535) == 0) R_CheckUserInterrupt();
        if (!R_FINITE(cc[i])) error("GeoModels pair cache requires finite coordinates");
    }
    const int use_time = !isNull(ctx->times);
    if (use_time) {
        PROTECT(times_r = coerceVector(ctx->times, REALSXP)); nprotect++;
        if (XLENGTH(times_r) != nobs)
            error("GeoModels pair cache requires one time per observation");
        ttimes = REAL(times_r);
        for (int i = 0; i < nobs; ++i)
            if (!R_FINITE(ttimes[i])) error("GeoModels pair cache requires finite times");
    }
    PROTECT(par_r = coerceVector(ctx->par, REALSXP)); nprotect++;
    int cormod = asInteger(ctx->cormod), type = asInteger(ctx->type);
    double radius = asReal(ctx->radius);
    geo_validate_corr_buffer(cormod, par_r);
    const double *par = REAL(par_r);
    if (!R_FINITE(radius) || radius <= 0.0)
        error("GeoModels pair cache received invalid radius");

    SEXP budget_opt = GetOption1(install("GeoModels.local_pair_cache_max_bytes"));
    ctx->budget = isNull(budget_opt) ? 256.0 * 1024.0 * 1024.0 : asReal(budget_opt);
    if ((!isNull(budget_opt) && XLENGTH(budget_opt) != 1) ||
        !R_FINITE(ctx->budget) || ctx->budget < 0.0)
        error("GeoModels.local_pair_cache_max_bytes must be one finite nonnegative number");

    R_xlen_t ntask = XLENGTH(neighborhoods);
    double total_occ = 0.0;
    for (R_xlen_t tt = 0; tt < ntask; ++tt) {
        if ((tt & 1023) == 0) R_CheckUserInterrupt();
        SEXP idx = VECTOR_ELT(neighborhoods, tt);
        R_xlen_t k = XLENGTH(idx);
        for (R_xlen_t a = 0; a < k; ++a)
            if (geo_index_at(idx, a) > nobs)
                error("GeoModels pair cache found an invalid observation index");
        if (k > 1) total_occ += (double)k * (double)(k - 1) / 2.0;
    }
    if (!R_FINITE(total_occ)) error("GeoModels pair occurrence count overflow");
    if (nobs == INT_MAX) goto no_cache; /* CSR row-vector length must fit int. */
    double est = fmax(16.0, fmin(fmax((double)nobs * 16.0, 4096.0), total_occ));
    if (est > (double)INT_MAX) est = INT_MAX;
    if (est / 0.65 + 1.0 > (double)(SIZE_MAX / 2)) goto no_cache;
    size_t capacity = geo_next_pow2((size_t)(est / 0.65) + 1u);
    if (capacity < 1024u) capacity = 1024u;
    if (capacity > SIZE_MAX / sizeof(uint64_t) || capacity > SIZE_MAX / sizeof(double) ||
        !geo_pair_budget_ok(ctx, (double)(sizeof(uint64_t) + sizeof(double)) * capacity) ||
        !geo_pair_budget_ok(ctx, sizeof(int) * ((double)nobs + 1.0))) goto no_cache;
    ctx->capacity = capacity;
    ctx->keys = (uint64_t *)calloc(capacity, sizeof(uint64_t));
    ctx->values = (double *)calloc(capacity, sizeof(double));
    if (!ctx->keys || !ctx->values) goto no_cache;

    size_t used = 0, visited = 0;
    for (R_xlen_t tt = 0; tt < ntask; ++tt) {
        SEXP idx = VECTOR_ELT(neighborhoods, tt);
        R_xlen_t k = XLENGTH(idx);
        for (R_xlen_t a = 0; a + 1 < k; ++a) {
            int ia = geo_index_at(idx, a) - 1;
            for (R_xlen_t b = a + 1; b < k; ++b) {
                if ((visited++ & 8191u) == 0) R_CheckUserInterrupt();
                int ib = geo_index_at(idx, b) - 1;
                if (ia == ib) continue;
                int i = ia < ib ? ia : ib, j = ia < ib ? ib : ia;
                uint64_t key = (((uint64_t)(uint32_t)i) << 32) | (uint64_t)(uint32_t)j;
                size_t slot = geo_pair_hash(key) & (ctx->capacity - 1u);
                while (ctx->keys[slot] && ctx->keys[slot] != key)
                    slot = (slot + 1u) & (ctx->capacity - 1u);
                if (ctx->keys[slot]) continue;
                if (used >= INT_MAX) goto no_cache; /* CSR offsets are int. */
                if ((double)(used + 1u) > 0.7 * (double)ctx->capacity) {
                    if (!geo_pair_rehash(ctx)) goto no_cache;
                    slot = geo_pair_hash(key) & (ctx->capacity - 1u);
                    while (ctx->keys[slot]) slot = (slot + 1u) & (ctx->capacity - 1u);
                }
                const double zi = ndim == 3 ? cc[i + (R_xlen_t)2 * nobs] : 0.0;
                const double zj = ndim == 3 ? cc[j + (R_xlen_t)2 * nobs] : 0.0;
                double dis = dist(type, cc[i], cc[j], cc[(R_xlen_t)i+nobs], cc[(R_xlen_t)j+nobs], zi, zj, radius);
                double dit = use_time ? fabs(ttimes[i] - ttimes[j]) : 0.0;
                double rho = CorFct(&cormod, dis, dit, (double *)par, 0, 0);
                if (!R_FINITE(rho)) error("GeoModels pair cache produced a non-finite correlation");
                ctx->keys[slot] = key; ctx->values[slot] = rho; ++used;
            }
        }
    }
    if (used > SIZE_MAX / sizeof(GeoPairEntry) ||
        !geo_pair_budget_ok(ctx, (double)(sizeof(uint64_t) + sizeof(double)) * ctx->capacity +
                                  (double)sizeof(GeoPairEntry) * used) ||
        !geo_pair_budget_ok(ctx, (double)(sizeof(GeoPairEntry) + sizeof(int) + sizeof(double)) * used +
                                  sizeof(int) * ((double)nobs + 1.0))) goto no_cache;
    ctx->entries = (GeoPairEntry *)calloc(used ? used : 1u, sizeof(GeoPairEntry));
    if (!ctx->entries) goto no_cache;
    size_t e = 0;
    for (size_t slot = 0; slot < ctx->capacity; ++slot) {
        if ((slot & 65535u) == 0) R_CheckUserInterrupt();
        if (ctx->keys[slot]) {
            ctx->entries[e].key = ctx->keys[slot];
            ctx->entries[e++].rho = ctx->values[slot];
        }
    }
    /* Release the hash before sorting and allocating the CSR output. */
    free(ctx->keys); ctx->keys = NULL;
    free(ctx->values); ctx->values = NULL;
    R_CheckUserInterrupt();
    qsort(ctx->entries, used, sizeof(GeoPairEntry), geo_pair_compare);
    R_CheckUserInterrupt();

    SEXP row_ptr, col_idx, rho_out, ans, names;
    PROTECT(row_ptr = allocVector(INTSXP, (R_xlen_t)nobs + 1)); nprotect++;
    PROTECT(col_idx = allocVector(INTSXP, (R_xlen_t)used)); nprotect++;
    PROTECT(rho_out = allocVector(REALSXP, (R_xlen_t)used)); nprotect++;
    int *rp = INTEGER(row_ptr), *ci = INTEGER(col_idx);
    double *rv = REAL(rho_out);
    memset(rp, 0, ((size_t)nobs + 1u) * sizeof(int));
    for (e = 0; e < used; ++e) {
        if ((e & 65535u) == 0) R_CheckUserInterrupt();
        int i = (int)(ctx->entries[e].key >> 32);
        rp[i + 1]++;
        ci[e] = (int)(ctx->entries[e].key & UINT64_C(0xffffffff));
        rv[e] = ctx->entries[e].rho;
    }
    for (int i = 0; i < nobs; ++i) rp[i + 1] += rp[i];
    free(ctx->entries); ctx->entries = NULL;

    PROTECT(ans = allocVector(VECSXP, 7)); nprotect++;
    PROTECT(names = allocVector(STRSXP, 7)); nprotect++;
    SET_VECTOR_ELT(ans, 0, row_ptr); SET_VECTOR_ELT(ans, 1, col_idx); SET_VECTOR_ELT(ans, 2, rho_out);
    SET_VECTOR_ELT(ans, 3, ScalarReal(total_occ)); SET_VECTOR_ELT(ans, 4, ScalarReal((double)used));
    SET_VECTOR_ELT(ans, 5, ScalarReal(ctx->peak)); SET_VECTOR_ELT(ans, 6, ScalarReal(ctx->budget));
    const char *nm[] = {"row_ptr", "col_idx", "rho", "total_pair_occurrences", "unique_pairs", "peak_bytes", "budget_bytes"};
    for (int i = 0; i < 7; ++i) SET_STRING_ELT(names, i, mkChar(nm[i]));
    setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(nprotect);
    return ans;

no_cache:
    UNPROTECT(nprotect);
    return R_NilValue;
}

static SEXP geo_build_pair_cache(SEXP neighborhoods, SEXP coords, SEXP times,
                                 SEXP cormod, SEXP par, SEXP type, SEXP radius)
{
    GeoPairBuildContext ctx = {0};
    ctx.neighborhoods = neighborhoods; ctx.coords = coords; ctx.times = times;
    ctx.cormod = cormod; ctx.par = par; ctx.type = type; ctx.radius = radius;
    return R_ExecWithCleanup(geo_build_pair_cache_body, &ctx, geo_pair_cleanup, &ctx);
}

SEXP GeoBuildSpatialPairCache(SEXP neighborhoods, SEXP coords,
                              SEXP cormod_s, SEXP par_s,
                              SEXP type_s, SEXP radius_s)
{
    return geo_build_pair_cache(neighborhoods, coords, R_NilValue,
                                cormod_s, par_s, type_s, radius_s);
}

SEXP GeoBuildSpacetimePairCache(SEXP neighborhoods, SEXP coords, SEXP times,
                                SEXP cormod_s, SEXP par_s,
                                SEXP type_s, SEXP radius_s)
{
    return geo_build_pair_cache(neighborhoods, coords, times,
                                cormod_s, par_s, type_s, radius_s);
}

/* Return strict-lower correlations in the same column-major order as
   x[lower.tri(x)].  The cache is CSR with zero-based column indices. */
SEXP GeoPairCachePacked(SEXP indices, SEXP row_ptr, SEXP col_idx, SEXP rho)
{
    SEXP out;
    R_xlen_t k = XLENGTH(indices), npairs, h = 0;
    const int *rp, *ci;
    const double *rv;
    int nrows;
    R_xlen_t jloc, iloc;

    if (TYPEOF(row_ptr) != INTSXP || TYPEOF(col_idx) != INTSXP ||
        TYPEOF(rho) != REALSXP || XLENGTH(col_idx) != XLENGTH(rho) ||
        XLENGTH(row_ptr) < 2 || XLENGTH(row_ptr) > INT_MAX || XLENGTH(rho) > INT_MAX)
        error("GeoModels pair cache has an invalid internal representation");

    rp = INTEGER(row_ptr); ci = INTEGER(col_idx); rv = REAL(rho);
    nrows = LENGTH(row_ptr) - 1;
    if (rp[0] != 0 || rp[nrows] != XLENGTH(rho))
        error("GeoModels pair cache has invalid CSR offsets");
    if ((double)k * (double)(k - 1) / 2.0 > R_XLEN_T_MAX)
        error("GeoModels pair cache lookup is too large");
    npairs = k > 1 ? (k / 2) * (k - 1) + (k % 2) * ((k - 1) / 2) : 0;
    PROTECT(out = allocVector(REALSXP, npairs));

    for (jloc = 0; jloc + 1 < k; ++jloc) {
        int ja = geo_index_at(indices, jloc) - 1;
        if (ja < 0 || ja >= nrows) {
            UNPROTECT(1);
            error("GeoModels pair cache lookup received an invalid index");
        }
        for (iloc = jloc + 1; iloc < k; ++iloc) {
            int ia = geo_index_at(indices, iloc) - 1;
            int i, j, lo, hi, found = -1;
            if (ia < 0 || ia >= nrows) {
                UNPROTECT(1);
                error("GeoModels pair cache lookup received an invalid index");
            }
            i = ia < ja ? ia : ja;
            j = ia < ja ? ja : ia;
            if (rp[i] < 0 || rp[i + 1] < rp[i] || rp[i + 1] > XLENGTH(rho))
                error("GeoModels pair cache has invalid CSR offsets");
            lo = rp[i];
            hi = rp[i + 1] - 1;
            while (lo <= hi) {
                int mid = lo + (hi - lo) / 2;
                if (ci[mid] < j) lo = mid + 1;
                else if (ci[mid] > j) hi = mid - 1;
                else { found = mid; break; }
            }
            if (found < 0) {
                UNPROTECT(1);
                error("GeoModels pair cache lookup missed a neighborhood pair");
            }
            REAL(out)[h++] = rv[found];
        }
    }

    UNPROTECT(1);
    return out;
}

/* Cross-correlation from indexed observations to one prediction point.
   This avoids rebuilding local time blocks in R and supports both fixed and
   dynamic space-time observation layouts through the same flattened index. */
SEXP GeoPairCrossCorrelation(SEXP indices, SEXP coords, SEXP times_s,
                             SEXP loc_s, SEXP pred_time_s,
                             SEXP cormod_s, SEXP par_s,
                             SEXP type_s, SEXP radius_s)
{
    SEXP coords_r, loc_r, par_r, times_r = R_NilValue, out, dim;
    int nobs, ndim, cormod, type, use_time, nprotect = 0;
    double radius, pred_time;
    const double *cc, *ll, *par, *ttimes = NULL;
    R_xlen_t k, a;

    dim = getAttrib(coords, R_DimSymbol);
    if (isNull(dim) || LENGTH(dim) != 2)
        error("GeoModels pair cross-correlation requires a coordinate matrix");
    nobs = INTEGER(dim)[0];
    ndim = INTEGER(dim)[1];
    if (nobs < 1 || (ndim != 2 && ndim != 3))
        error("GeoModels pair cross-correlation requires N x 2 or N x 3 coordinates");

    PROTECT(coords_r = coerceVector(coords, REALSXP)); nprotect++;
    PROTECT(loc_r = coerceVector(loc_s, REALSXP)); nprotect++;
    PROTECT(par_r = coerceVector(par_s, REALSXP)); nprotect++;
    if (XLENGTH(loc_r) != ndim)
        error("GeoModels pair cross-correlation received incompatible prediction coordinates");

    use_time = !isNull(times_s);
    if (use_time) {
        PROTECT(times_r = coerceVector(times_s, REALSXP)); nprotect++;
        if (XLENGTH(times_r) != nobs)
            error("GeoModels pair cross-correlation requires one time per observation");
        ttimes = REAL(times_r);
    }

    cc = REAL(coords_r);
    ll = REAL(loc_r);
    par = REAL(par_r);
    cormod = asInteger(cormod_s);
    geo_validate_corr_buffer(cormod, par_r);
    type = asInteger(type_s);
    radius = asReal(radius_s);
    pred_time = asReal(pred_time_s);
    if (!R_FINITE(radius) || !R_FINITE(pred_time))
        error("GeoModels pair cross-correlation received invalid scalar input");

    k = XLENGTH(indices);
    PROTECT(out = allocVector(REALSXP, k)); nprotect++;
    for (a = 0; a < k; ++a) {
        int i = geo_index_at(indices, a) - 1;
        double zi, zl, dis, dit, val;
        if (i < 0 || i >= nobs)
            error("GeoModels pair cross-correlation received an invalid observation index");
        zi = ndim == 3 ? cc[i + (R_xlen_t) 2 * nobs] : 0.0;
        zl = ndim == 3 ? ll[2] : 0.0;
        dis = dist(type, cc[i], ll[0], cc[(R_xlen_t)i + nobs], ll[1], zi, zl, radius);
        dit = use_time ? fabs(ttimes[i] - pred_time) : 0.0;
        val = CorFct(&cormod, dis, dit, (double *) par, 0, 0);
        if (!R_FINITE(val))
            error("GeoModels pair cross-correlation produced a non-finite value");
        REAL(out)[a] = val;
    }

    UNPROTECT(nprotect);
    return out;
}
