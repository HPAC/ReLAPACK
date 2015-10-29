#include "larpack.h"

void LARPACK(zlauum)(const char *uplo, const int *n,
        double *A, const int *ldA, int *info) {
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
        LAPACK(xerbla)("ZLAUUM", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(zlauu2)(uplo, n, A, ldA, info);
        return;
    }

    // Recursion
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + 2 * *ldA * n1;
    double *const A_BL = A                 + 2 * n1;
    double *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // 1
   	const double z1[] = {1, 0};

    // recursion(A_TL)
    LARPACK(zlauum)(uplo, &n1, A_TL, ldA, info);

    if (lower) {
        // A_TL = A_TL + A_BL' * A_BL
        BLAS(zherk)("L", "C", &n1, &n2, z1, A_BL, ldA, z1, A_TL, ldA);
        // A_BL = A_BR' * A_BL
        BLAS(ztrmm)("L", "L", "C", "N", &n2, &n1, z1, A_BR, ldA, A_BL, ldA);
    } else {
        // A_TL = A_TL + A_TR * A_TR'
        BLAS(zherk)("U", "N", &n1, &n2, z1, A_TR, ldA, z1, A_TL, ldA);
        // A_TR = A_TR * A_BR'
        BLAS(ztrmm)("R", "U", "C", "N", &n1, &n2, z1, A_BR, ldA, A_TR, ldA);
    }

    // recursion(A_BR)
    LARPACK(zlauum)(uplo, &n2, A_BR, ldA, info);
}
