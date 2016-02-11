#include "relapack.h"

static void RELAPACK(dpotrf_rec)(const char *, const int *, double *,
    const int *, int *);


void RELAPACK(dpotrf)(
    const char *uplo, const int *n,
    double *A, const int *ldA,
    int *info
) {

    // Check arguments
    const int lower = LAPACK(lsame)(uplo, "L");
    const int upper = LAPACK(lsame)(uplo, "U");
    *info = 0;
    if (!lower && !upper)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*ldA < MAX(1, *n))
        *info = -4;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("DPOTRF", &minfo);
        return;
    }

    // Clean char * arguments
    const char cleanuplo = lower ? 'L' : 'U';

    RELAPACK(dpotrf_rec)(&cleanuplo, n, A, ldA, info);
}


static void RELAPACK(dpotrf_rec)(
    const char *uplo, const int *n,
    double *A, const int *ldA,
    int *info
){

    if (*n <= MAX(CROSSOVER_DPOTRF, 1)) {
        // Unblocked
        LAPACK(dpotf2)(uplo, n, A, ldA, info);
        return;
    }

    // Constants
    const double ONE[]  = {1};
    const double MONE[] = {-1};

    // Splitting
    const int n1 = REC_SPLIT(*n);
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + *ldA * n1;
    double *const A_BL = A             + n1;
    double *const A_BR = A + *ldA * n1 + n1;

    // recursion(A_TL)
    RELAPACK(dpotrf_rec)(uplo, &n1, A_TL, ldA, info);
    if (*info)
        return;

    if (*uplo == 'L') {
        // A_BL = A_BL / A_TL'
        BLAS(dtrsm)("R", "L", "T", "N", &n2, &n1, ONE, A_TL, ldA, A_BL, ldA);
        // A_BR = A_BR - A_BL * A_BL'
        BLAS(dsyrk)("L", "N", &n2, &n1, MONE, A_BL, ldA, ONE, A_BR, ldA);
    } else {
        // A_TR = A_TL' \ A_TR
        BLAS(dtrsm)("L", "U", "T", "N", &n1, &n2, ONE, A_TL, ldA, A_TR, ldA);
        // A_BR = A_BR - A_TR' * A_TR
        BLAS(dsyrk)("U", "T", &n2, &n1, MONE, A_TR, ldA, ONE, A_BR, ldA);
    }

    // recursion(A_BR)
    RELAPACK(dpotrf_rec)(uplo, &n2, A_BR, ldA, info);
    if (*info)
        *info += n1;
}
