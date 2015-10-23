#include "larpack.h"

void LARPACK(ctrtri)(const char *uplo, const char *diag, const int *n, float *A, const int *ldA, int *info) {
    *info = 0;

    // Check arguments
    int lower = LAPACK(lsame)(uplo, "L");
    int upper = LAPACK(lsame)(uplo, "U");
    int nounit = LAPACK(lsame)(diag, "N");
    int unit = LAPACK(lsame)(diag, "U");
    if (!upper && !lower)
        *info = -1;
    else if (!nounit && !unit)
        *info = -2;
    else if (*n < 0)
        *info = -3;
    else if (*ldA < MAX(1, *n))
        *info = -5;
    if (*info != 0) {
        int minfo = -*info;
        LAPACK(xerbla)("CTRTRI", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) { 
        // Unblocked
        LAPACK(ctrti2)(uplo, diag, n, A, ldA, info);
        return;
    }

    // Recursive
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A; 
    float *const A_TR = A + 2 * *ldA * n1;
    float *const A_BL = A                 + 2 * n1; 
    float *const A_BR = A + 2 * *ldA * n1 + 2 * n1; 

    // 1, -1
    const float c1[] = {1, 0}, cm1[] = {-1, 0};

    // recursion(A_TL)
    LARPACK(ctrtri)(uplo, diag, &n1, A_TL, ldA, info);
    if (*info)
        return;

    if (lower) {
        // A_BL = - A_BL * A_TL
        BLAS(ctrmm)("R", "L", "N", diag, &n2, &n1, cm1, A_TL, ldA, A_BL, ldA);
        // A_BL = A_BR \ A_BL
        BLAS(ctrsm)("L", "L", "N", diag, &n2, &n1, c1, A_BR, ldA, A_BL, ldA);
    } else {
        // A_TR = - A_TL * A_TR
        BLAS(ctrmm)("L", "U", "N", diag, &n1, &n2, cm1, A_TL, ldA, A_TR, ldA);
        // A_TR = A_TR / A_BR
        BLAS(ctrsm)("R", "U", "N", diag, &n1, &n2, c1, A_BR, ldA, A_TR, ldA);
    }

    // recursion(A_BR)
    LARPACK(ctrtri)(uplo, diag, &n2, A_BR, ldA, info);
    if (*info)
        *info += n1;
}
