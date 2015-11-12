#include "relapack.h"

void RELAPACK(ztrsyl)(const char *tranA, const char *tranB, const int *isgn,
        const int *m, const int *n,
        const double *A, const int *ldA, const double *B, const int *ldB,
        double *C, const int *ldC, double *scale, int *info) {

    // Check arguments
    const int notransA = LAPACK(lsame)(tranA, "N");
    const int transA = LAPACK(lsame)(tranA, "C");
    const int notransB = LAPACK(lsame)(tranB, "N");
    const int transB = LAPACK(lsame)(tranB, "C");
    *info = 0;
    if (!transA && !notransA)
        *info = -1;
    else if (!transB && !notransB)
        *info = -2;
    else if (*isgn != 1 && *isgn != -1)
        *info = -3;
    else if (*m < 0)
        *info = -4;
    else if (*n < 0)
        *info = -5;
    else if (*ldA < MAX(1, *m))
        *info = -7;
    else if (*ldB < MAX(1, *n))
        *info = -9;
    else if (*ldC < MAX(1, *m))
        *info = -11;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("ZTRSYL", &minfo);
        return;
    }

    if (*m <= RELAPACK_CROSSOVER && *n <= RELAPACK_CROSSOVER) {
        // Unblocked
        LAPACK(ztrsy2)(tranA, tranB, isgn, m, n, A, ldA, B, ldB, C, ldC, scale, info);
        return;
    }

    // Recursive

    // Constants
    // 1, -1, -isgn
   	const double ONE[] = {1, 0}, MONE[] = {-1, 0}, MSGN[] = {-*isgn, 0};
    // 0
    const int iONE[] = {1};

    // Outputs
    double scale1[] = {1, 0}, scale2[] = {1, 0};
    int info1[1], info2[1];

    if (*m > *n) {
        const int m1 = (*m >= 16) ? ((*m + 8) / 16) * 8 : *m / 2;
        const int m2 = *m - m1;

        // A_TL A_TR
        //      A_BR
        const double *const A_TL = A;
        const double *const A_TR = A + 2 * *ldA * m1;
        const double *const A_BR = A + 2 * *ldA * m1 + 2 * m1;

        // C_T
        // C_B
        double *const C_T = C;
        double *const C_B = C + 2 * m1;

        if (notransA) {
            // recusion(A_BR, B, C_B)
            RELAPACK(ztrsyl)(tranA, tranB, isgn, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, scale1, info1);
            // C_T = C_T - A_TR * C_B
            BLAS(zgemm)("N", "N", &m1, n, &m2, MONE, A_TR, ldA, C_B, ldC, scale1, C_T, ldC);
            // recusion(A_TL, B, C_T)
            RELAPACK(ztrsyl)(tranA, tranB, isgn, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(zlascl)("G", iONE, iONE, ONE, scale2, &m2, n, C_B, ldC, info);
        } else {
            // recusion(A_TL, B, C_T)
            RELAPACK(ztrsyl)(tranA, tranB, isgn, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, scale1, info1);
            // C_B = C_B - A_TR' * C_T
            BLAS(zgemm)("C", "N", &m2, n, &m1, MONE, A_TR, ldA, C_T, ldC, scale1, C_B, ldC);
            // recusion(A_BR, B, C_B)
            RELAPACK(ztrsyl)(tranA, tranB, isgn, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(zlascl)("G", iONE, iONE, ONE, scale2, &m1, n, C_B, ldC, info);
        }
    } else {
        const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
        const int n2 = *n - n1;

        // B_TL B_TR
        //      B_BR
        const double *const B_TL = B;
        const double *const B_TR = B + 2 * *ldB * n1;
        const double *const B_BR = B + 2 * *ldB * n1 + 2 * n1;

        // C_L C_R
        double *const C_L = C;
        double *const C_R = C + 2 * *ldC * n1;

        if (notransB) {
            // recusion(A, B_TL, C_L)
            RELAPACK(ztrsyl)(tranA, tranB, isgn, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, scale1, info1);
            // C_R = C_R -/+ C_L * B_TR
            BLAS(zgemm)("N", "N", m, &n2, &n1, MSGN, C_L, ldC, B_TR, ldB, scale1, C_R, ldC);
            // recusion(A, B_BR, C_R)
            RELAPACK(ztrsyl)(tranA, tranB, isgn, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(zlascl)("G", iONE, iONE, ONE, scale2, m, &n1, C_L, ldC, info);
        } else {
            // recusion(A, B_BR, C_R)
            RELAPACK(ztrsyl)(tranA, tranB, isgn, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, scale1, info1);
            // C_L = C_L -/+ C_R * B_TR'
            BLAS(zgemm)("N", "C", m, &n1, &n2, MSGN, C_R, ldC, B_TR, ldB, scale1, C_L, ldC);
            // recusion(A, B_TL, C_L)
            RELAPACK(ztrsyl)(tranA, tranB, isgn, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(zlascl)("G", iONE, iONE, ONE, scale2, m, &n2, C_R, ldC, info);
        }
    }

    *scale = scale1[0] * scale2[0];
    *info = info1[0] || info2[0];
}
