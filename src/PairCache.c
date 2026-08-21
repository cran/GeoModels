#include "header.h"

typedef struct {
    uint64_t key;
    double rho;
} GeoPairEntry;

static size_t geo_pair_hash(uint64_t x)
{
    x ^= x >> 33;
    x *= UINT64_C(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x *= UINT64_C(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    return (size_t) x;
}

static size_t geo_next_pow2(size_t x)
{
    size_t p = 1;
    while (p < x && p <= (SIZE_MAX >> 1)) p <<= 1;
    return p;
}

static int geo_pair_compare(const void *aa, const void *bb)
{
    const GeoPairEntry *a = (const GeoPairEntry *) aa;
    const GeoPairEntry *b = (const GeoPairEntry *) bb;
    if (a->key < b->key) return -1;
    if (a->key > b->key) return 1;
    return 0;
}

static void geo_pair_rehash(uint64_t **keys, double **values,
                            size_t *capacity)
{
    size_t old_capacity = *capacity;
    size_t new_capacity = old_capacity << 1;
    uint64_t *old_keys = *keys;
    double *old_values = *values;
    uint64_t *new_keys;
    double *new_values;
    size_t s;

    if (new_capacity <= old_capacity)
        error("GeoModels pair cache exceeded addressable size");

    new_keys = (uint64_t *) R_Calloc(new_capacity, uint64_t);
    new_values = (double *) R_Calloc(new_capacity, double);

    for (s = 0; s < old_capacity; ++s) {
        uint64_t key = old_keys[s];
        if (key != 0) {
            size_t pos = geo_pair_hash(key) & (new_capacity - 1);
            while (new_keys[pos] != 0)
                pos = (pos + 1) & (new_capacity - 1);
            new_keys[pos] = key;
            new_values[pos] = old_values[s];
        }
    }

    R_Free(old_keys);
    R_Free(old_values);
    *keys = new_keys;
    *values = new_values;
    *capacity = new_capacity;
}

static int geo_index_at(SEXP x, R_xlen_t i)
{
    if (TYPEOF(x) == INTSXP) return INTEGER(x)[i];
    if (TYPEOF(x) == REALSXP) return (int) REAL(x)[i];
    error("GeoModels pair cache requires numeric neighborhood indices");
    return 0;
}

/* Shared builder for spatial and space-time local caches.  When times_s is
   NULL the temporal lag is zero.  Otherwise coords and times are flattened
   in the same global observation order used by neighborhood indices. */
static SEXP geo_build_pair_cache(SEXP neighborhoods, SEXP coords, SEXP times_s,
                                 SEXP cormod_s, SEXP par_s,
                                 SEXP type_s, SEXP radius_s)
{
    SEXP coords_r = R_NilValue, times_r = R_NilValue, par_r = R_NilValue;
    SEXP dim, row_ptr, col_idx, rho_out, ans, names;
    int nobs, ndim, cormod, type, use_time;
    double radius;
    const double *cc, *ttimes = NULL, *par;
    R_xlen_t ntask, tt;
    double total_occ = 0.0;
    size_t capacity, used = 0, estimate, slot;
    uint64_t *keys;
    double *values;
    GeoPairEntry *entries;
    R_xlen_t e;
    int *rp, *ci;
    double *rv;
    int nprotect = 0;

    if (!isNewList(neighborhoods))
        error("GeoModels pair cache requires a list of neighborhood indices");

    dim = getAttrib(coords, R_DimSymbol);
    if (isNull(dim) || LENGTH(dim) != 2)
        error("GeoModels pair cache requires a coordinate matrix");
    nobs = INTEGER(dim)[0];
    ndim = INTEGER(dim)[1];
    PROTECT(coords_r = coerceVector(coords, REALSXP)); nprotect++;
    if (nobs < 1 || (ndim != 2 && ndim != 3))
        error("GeoModels pair cache requires N x 2 or N x 3 coordinates");

    use_time = !isNull(times_s);
    if (use_time) {
        PROTECT(times_r = coerceVector(times_s, REALSXP)); nprotect++;
        if (XLENGTH(times_r) != nobs)
            error("GeoModels pair cache requires one time value per observation");
        ttimes = REAL(times_r);
    }

    PROTECT(par_r = coerceVector(par_s, REALSXP)); nprotect++;
    cc = REAL(coords_r);
    par = REAL(par_r);
    cormod = asInteger(cormod_s);
    type = asInteger(type_s);
    radius = asReal(radius_s);
    if (!R_FINITE(radius)) error("GeoModels pair cache received invalid radius");

    ntask = XLENGTH(neighborhoods);
    for (tt = 0; tt < ntask; ++tt) {
        SEXP idx = VECTOR_ELT(neighborhoods, tt);
        R_xlen_t k = XLENGTH(idx);
        if (k > 1) total_occ += ((double) k * (double) (k - 1)) / 2.0;
    }

    estimate = (size_t) nobs * 16u;
    if (estimate < 4096u) estimate = 4096u;
    if (total_occ > 0.0 && (double) estimate > total_occ)
        estimate = (size_t) total_occ;
    if (estimate < 16u) estimate = 16u;
    capacity = geo_next_pow2((size_t) ((double) estimate / 0.65) + 1u);
    if (capacity < 1024u) capacity = 1024u;

    keys = (uint64_t *) R_Calloc(capacity, uint64_t);
    values = (double *) R_Calloc(capacity, double);

    for (tt = 0; tt < ntask; ++tt) {
        SEXP idx = VECTOR_ELT(neighborhoods, tt);
        R_xlen_t k = XLENGTH(idx), a, b;
        for (a = 0; a + 1 < k; ++a) {
            int ia0 = geo_index_at(idx, a) - 1;
            if (ia0 < 0 || ia0 >= nobs) {
                R_Free(keys); R_Free(values);
                error("GeoModels pair cache found an invalid observation index");
            }
            for (b = a + 1; b < k; ++b) {
                int ib0 = geo_index_at(idx, b) - 1;
                int i, j;
                uint64_t key;
                if (ib0 < 0 || ib0 >= nobs) {
                    R_Free(keys); R_Free(values);
                    error("GeoModels pair cache found an invalid observation index");
                }
                if (ia0 == ib0) continue;
                i = ia0 < ib0 ? ia0 : ib0;
                j = ia0 < ib0 ? ib0 : ia0;
                key = (((uint64_t) (uint32_t) i) << 32) |
                      ((uint64_t) (uint32_t) j);

                if ((used + 1u) * 10u > capacity * 7u)
                    geo_pair_rehash(&keys, &values, &capacity);

                slot = geo_pair_hash(key) & (capacity - 1u);
                while (keys[slot] != 0 && keys[slot] != key)
                    slot = (slot + 1u) & (capacity - 1u);

                if (keys[slot] == 0) {
                    double zi = ndim == 3 ? cc[i + (R_xlen_t) 2 * nobs] : 0.0;
                    double zj = ndim == 3 ? cc[j + (R_xlen_t) 2 * nobs] : 0.0;
                    double dis = dist(type,
                                      cc[i], cc[j],
                                      cc[i + nobs], cc[j + nobs],
                                      zi, zj, radius);
                    double dit = use_time ? fabs(ttimes[i] - ttimes[j]) : 0.0;
                    keys[slot] = key;
                    values[slot] = CorFct(&cormod, dis, dit, (double *) par, 0, 0);
                    if (!R_FINITE(values[slot])) {
                        R_Free(keys); R_Free(values);
                        error("GeoModels pair cache produced a non-finite correlation");
                    }
                    ++used;
                }
            }
        }
    }

    entries = (GeoPairEntry *) R_Calloc(used > 0 ? used : 1u, GeoPairEntry);
    e = 0;
    for (slot = 0; slot < capacity; ++slot) {
        if (keys[slot] != 0) {
            entries[e].key = keys[slot];
            entries[e].rho = values[slot];
            ++e;
        }
    }
    qsort(entries, used, sizeof(GeoPairEntry), geo_pair_compare);

    PROTECT(row_ptr = allocVector(INTSXP, (R_xlen_t) nobs + 1)); nprotect++;
    PROTECT(col_idx = allocVector(INTSXP, (R_xlen_t) used)); nprotect++;
    PROTECT(rho_out = allocVector(REALSXP, (R_xlen_t) used)); nprotect++;
    rp = INTEGER(row_ptr);
    ci = INTEGER(col_idx);
    rv = REAL(rho_out);
    memset(rp, 0, ((size_t) nobs + 1u) * sizeof(int));

    for (e = 0; e < (R_xlen_t) used; ++e) {
        int i = (int) (entries[e].key >> 32);
        rp[i + 1]++;
        ci[e] = (int) (entries[e].key & UINT64_C(0xffffffff));
        rv[e] = entries[e].rho;
    }
    for (int i = 0; i < nobs; ++i) rp[i + 1] += rp[i];

    R_Free(entries);
    R_Free(keys);
    R_Free(values);

    PROTECT(ans = allocVector(VECSXP, 5)); nprotect++;
    PROTECT(names = allocVector(STRSXP, 5)); nprotect++;
    SET_VECTOR_ELT(ans, 0, row_ptr);
    SET_VECTOR_ELT(ans, 1, col_idx);
    SET_VECTOR_ELT(ans, 2, rho_out);
    SET_VECTOR_ELT(ans, 3, ScalarReal(total_occ));
    SET_VECTOR_ELT(ans, 4, ScalarReal((double) used));
    SET_STRING_ELT(names, 0, mkChar("row_ptr"));
    SET_STRING_ELT(names, 1, mkChar("col_idx"));
    SET_STRING_ELT(names, 2, mkChar("rho"));
    SET_STRING_ELT(names, 3, mkChar("total_pair_occurrences"));
    SET_STRING_ELT(names, 4, mkChar("unique_pairs"));
    setAttrib(ans, R_NamesSymbol, names);

    UNPROTECT(nprotect);
    return ans;
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
    const int *rp = INTEGER(row_ptr), *ci = INTEGER(col_idx);
    const double *rv = REAL(rho);
    int nrows = LENGTH(row_ptr) - 1;
    R_xlen_t jloc, iloc;

    if (TYPEOF(row_ptr) != INTSXP || TYPEOF(col_idx) != INTSXP ||
        TYPEOF(rho) != REALSXP || XLENGTH(col_idx) != XLENGTH(rho))
        error("GeoModels pair cache has an invalid internal representation");

    npairs = k > 1 ? k * (k - 1) / 2 : 0;
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
        dis = dist(type, cc[i], ll[0], cc[i + nobs], ll[1], zi, zl, radius);
        dit = use_time ? fabs(ttimes[i] - pred_time) : 0.0;
        val = CorFct(&cormod, dis, dit, (double *) par, 0, 0);
        if (!R_FINITE(val))
            error("GeoModels pair cross-correlation produced a non-finite value");
        REAL(out)[a] = val;
    }

    UNPROTECT(nprotect);
    return out;
}
