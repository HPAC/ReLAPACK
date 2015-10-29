#include "larpack.h"

void LARPACK(dsygst)(const int *itype, const char *uplo, const int *n, double *A, const int *ldA, const double *B, const int *ldB, int *info) {
    *info = 0;

    // Check arguments
    int lower = LAPACK(lsame)(uplo, "L");
    int upper = LAPACK(lsame)(uplo, "U");
    if (*itype < 1 || *itype > 3)
        *info = -1;
    else if (!lower && !upper)
        *info = -2;
    else if (*n < 0)
        *info = -3;
    else if (*ldA < MAX(1, *n))
        *info = -5;
    else if (*ldB < MAX(1, *n))
        *info = -7;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("DSYGST", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(dsygs2)(itype, uplo, n, A, ldA, B, ldB, info);
        return;
    }


    // Recursive
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + *ldA * n1;
    double *const A_BL = A             + n1;
    double *const A_BR = A + *ldA * n1 + n1;

    // B_TL B_TR
    // B_BL B_BR
    const double *const B_TL = B;
    const double *const B_TR = B + *ldB * n1;
    const double *const B_BL = B             + n1;
    const double *const B_BR = B + *ldB * n1 + n1;

    // 1, -1, 1/2, -1/2
   	const double d1[] = {1}, dm1[] = {-1}, dp5[] = {.5}, dmp5[] = {-.5};

    // recursion(A_TL, B_TL)
    LARPACK(dsygst)(itype, uplo, &n1, A_TL, ldA, B_TL, ldB, info);

    if (*itype == 1)
        if (lower) {
            // A_BL = A_BL / B_TL'
            BLAS(dtrsm)("R", "L", "T", "N", &n2, &n1, d1, B_TL, ldB, A_BL, ldA);
            // A_BL = A_BL - 1/2 B_BL * A_TL
            BLAS(dsymm)("R", "L", &n2, &n1, dmp5, A_TL, ldA, B_BL, ldB, d1, A_BL, ldA);
            // A_BR = A_BR - A_BL * B_BL' - B_BL * A_BL'
            BLAS(dsyr2k)("L", "N", &n2, &n1, dm1, A_BL, ldA, B_BL, ldB, d1, A_BR, ldA);
            // A_BL = A_BL - 1/2 B_BL * A_TL
            BLAS(dsymm)("R", "L", &n2, &n1, dmp5, A_TL, ldA, B_BL, ldB, d1, A_BL, ldA);
            // A_BL = B_BR \ A_BL
            BLAS(dtrsm)("L", "L", "N", "N", &n2, &n1, d1, B_BR, ldB, A_BL, ldA);
        } else {
            // A_TR = B_TL' \ A_TR
            BLAS(dtrsm)("L", "U", "T", "N", &n1, &n2, d1, B_TL, ldB, A_TR, ldA);
            // A_TR = A_TR - 1/2 A_TL * B_TR
            BLAS(dsymm)("L", "U", &n1, &n2, dmp5, A_TL, ldA, B_TR, ldB, d1, A_TR, ldA);
            // A_BR = A_BR - A_TR' * B_TR - B_TR' * A_TR
            BLAS(dsyr2k)("U", "T", &n2, &n1, dm1, A_TR, ldA, B_TR, ldB, d1, A_BR, ldA);
            // A_TR = A_TR - 1/2 A_TL * B_TR
            BLAS(dsymm)("L", "U", &n1, &n2, dmp5, A_TL, ldA, B_TR, ldB, d1, A_TR, ldA);
            // A_TR = A_TR / B_BR
            BLAS(dtrsm)("R", "U", "N", "N", &n1, &n2, d1, B_BR, ldB, A_TR, ldA);
        }
    else
        if (lower) {
            // A_BL = A_BL * B_TL
            BLAS(dtrmm)("R", "L", "N", "N", &n2, &n1, d1, B_TL, ldB, A_BL, ldA);
            // A_BL = A_BL + 1/2 A_BR * B_BL
            BLAS(dsymm)("L", "L", &n2, &n1, dp5, A_BR, ldA, B_BL, ldB, d1, A_BL, ldA);
            // A_TL = A_TL + A_BL' * B_BL + B_BL' * A_BL
            BLAS(dsyr2k)("L", "T", &n1, &n2, d1, A_BL, ldA, B_BL, ldB, d1, A_TL, ldA);
            // A_BL = A_BL + 1/2 A_BR * B_BL
            BLAS(dsymm)("L", "L", &n2, &n1, dp5, A_BR, ldA, B_BL, ldB, d1, A_BL, ldA);
            // A_BL = B_BR * A_BL
            BLAS(dtrmm)("L", "L", "T", "N", &n2, &n1, d1, B_BR, ldB, A_BL, ldA);
        } else {
            // A_TR = B_TL * A_TR
            BLAS(dtrmm)("L", "U", "N", "N", &n1, &n2, d1, B_TL, ldB, A_TR, ldA);
            // A_TR = A_TR + 1/2 B_TR A_BR
            BLAS(dsymm)("R", "U", &n1, &n2, dp5, A_BR, ldA, B_TR, ldB, d1, A_TR, ldA);
            // A_TL = A_TL + A_TR * B_TR' + B_TR * A_TR'
            BLAS(dsyr2k)("U", "N", &n1, &n2, d1, A_TR, ldA, B_TR, ldB, d1, A_TL, ldA);
            // A_TR = A_TR + 1/2 B_TR * A_BR
            BLAS(dsymm)("R", "U", &n1, &n2, dp5, A_BR, ldA, B_TR, ldB, d1, A_TR, ldA);
            // A_TR = A_TR * B_BR
            BLAS(dtrmm)("R", "U", "T", "N", &n1, &n2, d1, B_BR, ldB, A_TR, ldA);
        }

    // recursion(A_BR, B_BR)
    LARPACK(dsygst)(itype, uplo, &n2, A_BR, ldA, B_BR, ldB, info);
}
