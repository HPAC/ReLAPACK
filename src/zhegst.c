#include "larpack.h"

void LARPACK(zhegst)(const int *itype, const char *uplo, const int *n,
        double *A, const int *ldA, const double *B, const int *ldB, int *info) {

    // Check arguments
    const int lower = LAPACK(lsame)(uplo, "L");
    const int upper = LAPACK(lsame)(uplo, "U");
    *info = 0;
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
        LAPACK(xerbla)("ZHEGST", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(zhegs2)(itype, uplo, n, A, ldA, B, ldB, info);
        return;
    }

    // Recursive

    // Constants
    // 1, -1, 1/2, -1/2
   	const double z1[] = {1, 0}, zm1[] = {-1, 0}, zp5[] = {.5, 0}, zmp5[] = {-.5, 0};

    // Splitting
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + 2 * *ldA * n1;
    double *const A_BL = A                 + 2 * n1;
    double *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // B_TL B_TR
    // B_BL B_BR
    const double *const B_TL = B;
    const double *const B_TR = B + 2 * *ldB * n1;
    const double *const B_BL = B                 + 2 * n1;
    const double *const B_BR = B + 2 * *ldB * n1 + 2 * n1;

    // recursion(A_TL, B_TL)
    LARPACK(zhegst)(itype, uplo, &n1, A_TL, ldA, B_TL, ldB, info);

    if (*itype == 1)
        if (lower) {
            // A_BL = A_BL / B_TL'
            BLAS(ztrsm)("R", "L", "C", "N", &n2, &n1, z1, B_TL, ldB, A_BL, ldA);
            // A_BL = A_BL - 1/2 B_BL * A_TL
            BLAS(zhemm)("R", "L", &n2, &n1, zmp5, A_TL, ldA, B_BL, ldB, z1, A_BL, ldA);
            // A_BR = A_BR - A_BL * B_BL' - B_BL * A_BL'
            BLAS(zher2k)("L", "N", &n2, &n1, zm1, A_BL, ldA, B_BL, ldB, z1, A_BR, ldA);
            // A_BL = A_BL - 1/2 B_BL * A_TL
            BLAS(zhemm)("R", "L", &n2, &n1, zmp5, A_TL, ldA, B_BL, ldB, z1, A_BL, ldA);
            // A_BL = B_BR \ A_BL
            BLAS(ztrsm)("L", "L", "N", "N", &n2, &n1, z1, B_BR, ldB, A_BL, ldA);
        } else {
            // A_TR = B_TL' \ A_TR
            BLAS(ztrsm)("L", "U", "C", "N", &n1, &n2, z1, B_TL, ldB, A_TR, ldA);
            // A_TR = A_TR - 1/2 A_TL * B_TR
            BLAS(zhemm)("L", "U", &n1, &n2, zmp5, A_TL, ldA, B_TR, ldB, z1, A_TR, ldA);
            // A_BR = A_BR - A_TR' * B_TR - B_TR' * A_TR
            BLAS(zher2k)("U", "C", &n2, &n1, zm1, A_TR, ldA, B_TR, ldB, z1, A_BR, ldA);
            // A_TR = A_TR - 1/2 A_TL * B_TR
            BLAS(zhemm)("L", "U", &n1, &n2, zmp5, A_TL, ldA, B_TR, ldB, z1, A_TR, ldA);
            // A_TR = A_TR / B_BR
            BLAS(ztrsm)("R", "U", "N", "N", &n1, &n2, z1, B_BR, ldB, A_TR, ldA);
        }
    else
        if (lower) {
            // A_BL = A_BL * B_TL
            BLAS(ztrmm)("R", "L", "N", "N", &n2, &n1, z1, B_TL, ldB, A_BL, ldA);
            // A_BL = A_BL + 1/2 A_BR * B_BL
            BLAS(zhemm)("L", "L", &n2, &n1, zp5, A_BR, ldA, B_BL, ldB, z1, A_BL, ldA);
            // A_TL = A_TL + A_BL' * B_BL + B_BL' * A_BL
            BLAS(zher2k)("L", "C", &n1, &n2, z1, A_BL, ldA, B_BL, ldB, z1, A_TL, ldA);
            // A_BL = A_BL + 1/2 A_BR * B_BL
            BLAS(zhemm)("L", "L", &n2, &n1, zp5, A_BR, ldA, B_BL, ldB, z1, A_BL, ldA);
            // A_BL = B_BR * A_BL
            BLAS(ztrmm)("L", "L", "C", "N", &n2, &n1, z1, B_BR, ldB, A_BL, ldA);
        } else {
            // A_TR = B_TL * A_TR
            BLAS(ztrmm)("L", "U", "N", "N", &n1, &n2, z1, B_TL, ldB, A_TR, ldA);
            // A_TR = A_TR + 1/2 B_TR A_BR
            BLAS(zhemm)("R", "U", &n1, &n2, zp5, A_BR, ldA, B_TR, ldB, z1, A_TR, ldA);
            // A_TL = A_TL + A_TR * B_TR' + B_TR * A_TR'
            BLAS(zher2k)("U", "N", &n1, &n2, z1, A_TR, ldA, B_TR, ldB, z1, A_TL, ldA);
            // A_TR = A_TR + 1/2 B_TR * A_BR
            BLAS(zhemm)("R", "U", &n1, &n2, zp5, A_BR, ldA, B_TR, ldB, z1, A_TR, ldA);
            // A_TR = A_TR * B_BR
            BLAS(ztrmm)("R", "U", "C", "N", &n1, &n2, z1, B_BR, ldB, A_TR, ldA);
        }

    // recursion(A_BR, B_BR)
    LARPACK(zhegst)(itype, uplo, &n2, A_BR, ldA, B_BR, ldB, info);
}
