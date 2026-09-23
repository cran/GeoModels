#ifndef GEOMODELS_CORRELATION_PARAM_COUNT_H
#define GEOMODELS_CORRELATION_PARAM_COUNT_H
/* Keep in sync with R/Utility.r:CorrelationPar; covered by regression tests. */
static int geo_correlation_param_count(int code)
{
    switch (code) {
    case 1: return 2;
    case 2: return 1;
    case 3: return 1;
    case 4: return 1;
    case 5: return 3;
    case 6: return 3;
    case 7: return 3;
    case 8: return 3;
    case 10: return 3;
    case 11: return 2;
    case 12: return 2;
    case 13: return 2;
    case 14: return 2;
    case 15: return 2;
    case 16: return 1;
    case 17: return 2;
    case 18: return 1;
    case 19: return 3;
    case 20: return 2;
    case 21: return 4;
    case 22: return 3;
    case 23: return 3;
    case 24: return 3;
    case 25: return 3;
    case 26: return 4;
    case 27: return 3;
    case 28: return 1;
    case 29: return 4;
    case 30: return 3;
    case 31: return 1;
    case 32: return 1;
    case 42: return 5;
    case 44: return 5;
    case 46: return 5;
    case 48: return 4;
    case 50: return 5;
    case 52: return 5;
    case 54: return 5;
    case 56: return 3;
    case 58: return 4;
    case 60: return 5;
    case 61: return 6;
    case 62: return 6;
    case 63: return 6;
    case 64: return 6;
    case 65: return 6;
    case 66: return 6;
    case 67: return 6;
    case 68: return 6;
    case 69: return 4;
    case 70: return 4;
    case 71: return 4;
    case 72: return 4;
    case 73: return 4;
    case 74: return 4;
    case 75: return 4;
    case 76: return 4;
    case 77: return 4;
    case 78: return 6;
    case 82: return 3;
    case 84: return 2;
    case 85: return 7;
    case 86: return 4;
    case 95: return 5;
    case 87: return 7;
    case 88: return 7;
    case 89: return 5;
    case 94: return 4;
    case 96: return 4;
    case 111: return 7;
    case 112: return 11;
    case 113: return 7;
    case 114: return 11;
    case 115: return 7;
    case 116: return 11;
    case 117: return 11;
    case 118: return 9;
    case 119: return 7;
    case 120: return 9;
    case 121: return 9;
    case 122: return 7;
    case 124: return 7;
    case 126: return 8;
    case 128: return 11;
    case 129: return 9;
    case 130: return 8;
    case 131: return 9;
    case 132: return 14;
    case 134: return 11;
    case 136: return 12;
    case 137: return 14;
    case 138: return 14;
    case 139: return 11;
    default: return -1;
    }
}

static void geo_validate_corr_buffer(int code, SEXP par)
{
    const int n = geo_correlation_param_count(code);
    if (n < 0 || TYPEOF(par) != REALSXP || XLENGTH(par) != n)
        error("GeoModels: invalid correlation model or native parameter length");
    for (int i = 0; i < n; ++i)
        if (!R_FINITE(REAL(par)[i]))
            error("GeoModels: correlation parameters must be finite");
}
#endif
