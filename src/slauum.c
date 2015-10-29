#include "larpack.h"

void LARPACK(slauum)(const char *uplo, const int *n,
        float *A, const int *ldA, int *info) {
    *info = 0;

    // Check arguments
    int lower = LAPACK(lsame)(uplo, "L");
    int upper = LAPACK(lsame)(uplo, "U");
    if (!lower && !upper)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*ldA < MAX(1, *n))
        *info = -4;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("SLAUUM", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(slauu2)(uplo, n, A, ldA, info);
        return;
    }

    // Recursion
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + *ldA * n1;
    float *const A_BL = A             + n1;
    float *const A_BR = A + *ldA * n1 + n1;

    // 1
   	const float s1[] = {1};

    // recursion(A_TL)
    LARPACK(slauum)(uplo, &n1, A_TL, ldA, info);

    if (lower) {
        // A_TL = A_TL + A_BL' * A_BL
        BLAS(ssyrk)("L", "T", &n1, &n2, s1, A_BL, ldA, s1, A_TL, ldA);
        // A_BL = A_BR' * A_BL
        BLAS(strmm)("L", "L", "T", "N", &n2, &n1, s1, A_BR, ldA, A_BL, ldA);
    } else {
        // A_TL = A_TL + A_TR * A_TR'
        BLAS(ssyrk)("U", "N", &n1, &n2, s1, A_TR, ldA, s1, A_TL, ldA);
        // A_TR = A_TR * A_BR'
        BLAS(strmm)("R", "U", "T", "N", &n1, &n2, s1, A_BR, ldA, A_TR, ldA);
    }

    // recursion(A_BR)
    LARPACK(slauum)(uplo, &n2, A_BR, ldA, info);
}
