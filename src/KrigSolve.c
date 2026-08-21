#include "header.h"

/*
 * Dense symmetric-positive-definite solve for local kriging systems.
 *
 * A is an n x n covariance matrix in R column-major storage and B is an
 * n x nrhs right-hand-side matrix.  LAPACK dposv factorizes A and solves
 * A X = B in one native call.  Both A and B are overwritten, exactly as
 * documented by LAPACK; the R wrapper passes private copies.
 */
void GeoDenseSPDSolve(double *A, double *B, int *n, int *nrhs, int *info)
{
    char uplo = 'L';
    int lda = *n;
    int ldb = *n;

    *info = 0;
    if (*n <= 0 || *nrhs <= 0) {
        *info = -1;
        return;
    }

    F77_CALL(dposv)(&uplo, n, nrhs, A, &lda, B, &ldb, info FCONE);
}
