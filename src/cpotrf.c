#include "relapack.h"

static void RELAPACK(cpotrf_rec)(const char *, const int *, float *, 
    const int *, int *);


void RELAPACK(cpotrf)(
    const char *uplo, const int *n,
    float *A, const int *ldA, 
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
        LAPACK(xerbla)("CPOTRF", &minfo);
        return;
    }

    // Clean char * arguments
    const char cleanuplo = lower ? 'L' : 'U';

    RELAPACK(cpotrf_rec)(&cleanuplo, n, A, ldA, info);
}


static void RELAPACK(cpotrf_rec)(
    const char *uplo, const int *n,
    float *A, const int *ldA,
    int *info
){

    if (*n <= CROSSOVER_CPOTRF) {
        // Unblocked
        LAPACK(cpotf2)(uplo, n, A, ldA, info);
        return;
    }

    // Recursive

    // Constants
    // 1, -1
   	const float ONE[] = {1, 0}, MONE[] = {-1, 0};

    // Splitting
    const int n1 = REC_SPLIT(*n);
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + 2 * *ldA * n1;
    float *const A_BL = A                 + 2 * n1;
    float *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // recursion(A_TL)
    RELAPACK(cpotrf_rec)(uplo, &n1, A_TL, ldA, info);
    if (*info)
        return;

    if (*uplo == 'L') {
        // A_BL = A_BL / A_TL'
        BLAS(ctrsm)("R", "L", "C", "N", &n2, &n1, ONE, A_TL, ldA, A_BL, ldA);
        // A_BR = A_BR - A_BL * A_BL'
        BLAS(cherk)("L", "N", &n2, &n1, MONE, A_BL, ldA, ONE, A_BR, ldA);
    } else {
        // A_TR = A_TL' \ A_TR
        BLAS(ctrsm)("L", "U", "C", "N", &n1, &n2, ONE, A_TL, ldA, A_TR, ldA);
        // A_BR = A_BR - A_TR' * A_TR
        BLAS(cherk)("U", "C", &n2, &n1, MONE, A_TR, ldA, ONE, A_BR, ldA);
    }

    // recursion(A_BR)
    RELAPACK(cpotrf_rec)(uplo, &n2, A_BR, ldA, info);
    if (*info)
        *info += n1;
}
