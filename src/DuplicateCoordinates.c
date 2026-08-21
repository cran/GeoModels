#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static uint64_t mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

static uint64_t double_bits(double x)
{
    uint64_t u;
    if (x == 0.0) return UINT64_C(0); /* +0 and -0 are the same coordinate */
    memcpy(&u, &x, sizeof(uint64_t));
    return u;
}

static uint64_t row_hash(const double *x, R_xlen_t n, int d, R_xlen_t i)
{
    uint64_t h = UINT64_C(0x9e3779b97f4a7c15);
    int j;
    for (j = 0; j < d; ++j) {
        uint64_t z = mix64(double_bits(x[i + (R_xlen_t)j * n]) +
                           UINT64_C(0x9e3779b97f4a7c15) * (uint64_t)(j + 1));
        h ^= z + UINT64_C(0x9e3779b97f4a7c15) + (h << 6) + (h >> 2);
    }
    return mix64(h);
}

static int rows_equal(const double *x, R_xlen_t n, int d,
                      R_xlen_t i, R_xlen_t k)
{
    int j;
    for (j = 0; j < d; ++j) {
        double a = x[i + (R_xlen_t)j * n];
        double b = x[k + (R_xlen_t)j * n];
        if (a != b) return 0;
    }
    return 1;
}

SEXP GeoDuplicateCoordinates(SEXP coords)
{
    SEXP dim, ans;
    const double *x;
    R_xlen_t n, table_size, mask, i, pos, old;
    int d, j;
    R_xlen_t *table;

    if (!Rf_isReal(coords) || !Rf_isMatrix(coords))
        Rf_error("coordinates must be a numeric matrix");

    dim = Rf_getAttrib(coords, R_DimSymbol);
    n = (R_xlen_t) INTEGER(dim)[0];
    d = INTEGER(dim)[1];
    if (d < 1 || d > 4)
        Rf_error("coordinate matrix must have between one and four columns");

    x = REAL(coords);
    for (j = 0; j < d; ++j) {
        for (i = 0; i < n; ++i) {
            if (!R_FINITE(x[i + (R_xlen_t)j * n]))
                Rf_error("coordinates must be finite");
        }
    }

    PROTECT(ans = Rf_allocVector(REALSXP, 2));
    REAL(ans)[0] = 0.0;
    REAL(ans)[1] = 0.0;
    if (n < 2) {
        UNPROTECT(1);
        return ans;
    }

    table_size = 1;
    while (table_size < 2 * n) {
        table_size <<= 1;
    }
    mask = table_size - 1;
    table = (R_xlen_t *) R_Calloc((size_t)table_size, R_xlen_t);

    for (i = 0; i < n; ++i) {
        pos = (R_xlen_t)(row_hash(x, n, d, i) & (uint64_t)mask);
        for (;;) {
            old = table[pos];
            if (old == 0) {
                table[pos] = i + 1;
                break;
            }
            if (rows_equal(x, n, d, i, old - 1)) {
                REAL(ans)[0] = (double)(old);     /* first occurrence, 1-based */
                REAL(ans)[1] = (double)(i + 1);  /* duplicate, 1-based */
                R_Free(table);
                UNPROTECT(1);
                return ans;
            }
            pos = (pos + 1) & mask;
        }
    }

    R_Free(table);
    UNPROTECT(1);
    return ans;
}
