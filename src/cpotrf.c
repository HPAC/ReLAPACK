#include "larpack.h"

void LARPACK(cpotrf)(const char *uplo, const int *n,
        float *A, const int *ldA, int *info) {

    // Check arguments
    int lower = LAPACK(lsame)(uplo, "L");
    int upper = LAPACK(lsame)(uplo, "U");
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

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(cpotf2)(uplo, n, A, ldA, info);
        return;
    }

    // Recursive

    // Constants
    // 1, -1
   	const float c1[] = {1, 0}, cm1[] = {-1, 0};

    // Splitting
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + 2 * *ldA * n1;
    float *const A_BL = A                 + 2 * n1;
    float *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // recursion(A_TL)
    LARPACK(cpotrf)(uplo, &n1, A_TL, ldA, info);
    if (*info)
        return;

    if (lower) {
        // A_BL = A_BL / A_TL'
        BLAS(ctrsm)("R", "L", "C", "N", &n2, &n1, c1, A_TL, ldA, A_BL, ldA);
        // A_BR = A_BR - A_BL * A_BL'
        BLAS(cherk)("L", "N", &n2, &n1, cm1, A_BL, ldA, c1, A_BR, ldA);
    } else {
        // A_TR = A_TL' \ A_TR
        BLAS(ctrsm)("L", "U", "C", "N", &n1, &n2, c1, A_TL, ldA, A_TR, ldA);
        // A_BR = A_BR - A_TR' * A_TR
        BLAS(cherk)("U", "C", &n2, &n1, cm1, A_TR, ldA, c1, A_BR, ldA);
    }

    // recursion(A_BR)
    LARPACK(cpotrf)(uplo, &n2, A_BR, ldA, info);
    if (*info)
        *info += n1;
}
