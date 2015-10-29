#include "larpack.h"

void LARPACK(clauum)(const char *uplo, const int *n,
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
        LAPACK(xerbla)("CLAUUM", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(clauu2)(uplo, n, A, ldA, info);
        return;
    }

    // Recursion
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + 2 * *ldA * n1;
    float *const A_BL = A                 + 2 * n1;
    float *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // 1
   	const float c1[] = {1, 0};

    // recursion(A_TL)
    LARPACK(clauum)(uplo, &n1, A_TL, ldA, info);

    if (lower) {
        // A_TL = A_TL + A_BL' * A_BL
        BLAS(cherk)("L", "C", &n1, &n2, c1, A_BL, ldA, c1, A_TL, ldA);
        // A_BL = A_BR' * A_BL
        BLAS(ctrmm)("L", "L", "C", "N", &n2, &n1, c1, A_BR, ldA, A_BL, ldA);
    } else {
        // A_TL = A_TL + A_TR * A_TR'
        BLAS(cherk)("U", "N", &n1, &n2, c1, A_TR, ldA, c1, A_TL, ldA);
        // A_TR = A_TR * A_BR'
        BLAS(ctrmm)("R", "U", "C", "N", &n1, &n2, c1, A_BR, ldA, A_TR, ldA);
    }

    // recursion(A_BR)
    LARPACK(clauum)(uplo, &n2, A_BR, ldA, info);
}
