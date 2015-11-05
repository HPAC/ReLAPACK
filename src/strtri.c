#include "larpack.h"

#ifdef LARPACK_SMALL_LAPACK
#define SMALL strtri
#else
#define SMALL strti2
#endif

void LARPACK(strtri)(const char *uplo, const char *diag, const int *n,
        float *A, const int *ldA, int *info) {

    // Check arguments
    const int lower = LAPACK(lsame)(uplo, "L");
    const int upper = LAPACK(lsame)(uplo, "U");
    const int nounit = LAPACK(lsame)(diag, "N");
    const int unit = LAPACK(lsame)(diag, "U");
    *info = 0;
    if (!lower && !upper)
        *info = -1;
    else if (!nounit && !unit)
        *info = -2;
    else if (*n < 0)
        *info = -3;
    else if (*ldA < MAX(1, *n))
        *info = -5;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("STRTRI", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(SMALL)(uplo, diag, n, A, ldA, info);
        return;
    }

    // Recursive

    // Constants
    // 1, -1
    const float ONE[] = {1}, MONE[] = {-1};

    // Splitting
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + *ldA * n1;
    float *const A_BL = A             + n1;
    float *const A_BR = A + *ldA * n1 + n1;

    // recursion(A_TL)
    LARPACK(strtri)(uplo, diag, &n1, A_TL, ldA, info);
    if (*info)
        return;

    if (lower) {
        // A_BL = - A_BL * A_TL
        BLAS(strmm)("R", "L", "N", diag, &n2, &n1, MONE, A_TL, ldA, A_BL, ldA);
        // A_BL = A_BR \ A_BL
        BLAS(strsm)("L", "L", "N", diag, &n2, &n1, ONE, A_BR, ldA, A_BL, ldA);
    } else {
        // A_TR = - A_TL * A_TR
        BLAS(strmm)("L", "U", "N", diag, &n1, &n2, MONE, A_TL, ldA, A_TR, ldA);
        // A_TR = A_TR / A_BR
        BLAS(strsm)("R", "U", "N", diag, &n1, &n2, ONE, A_BR, ldA, A_TR, ldA);
    }

    // recursion(A_BR)
    LARPACK(strtri)(uplo, diag, &n2, A_BR, ldA, info);
    if (*info)
        *info += n1;
}
