#include "larpack.h"

void LARPACK(ctrsyl)(const char *tranA, const char *tranB, const int *isgn,
        const int *m, const int *n,
        const float *A, const int *ldA, const float *B, const int *ldB,
        float *C, const int *ldC, float *scale, int *info) {
    *info = 0;

    // Check arguments
    int notransA = LAPACK(lsame)(tranA, "N");
    int transA = LAPACK(lsame)(tranA, "T") || LAPACK(lsame)(tranA, "C");
    int notransB = LAPACK(lsame)(tranB, "N");
    int transB = LAPACK(lsame)(tranB, "T") || LAPACK(lsame)(tranB, "C");
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
        LAPACK(xerbla)("CTRSYL", &minfo);
        return;
    }

    if (*m <= LARPACK_CROSSOVER && *n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(ctrsy2)(tranA, tranB, isgn, m, n, A, ldA, B, ldB, C, ldC, scale, info);
        return;
    }

    // Recursive

    // 1, -1, -isgn
   	const float c1[] = {1, 0}, cm1[] = {-1, 0}, cmisgn[] = {-*isgn, 0};
    // 0
    int i0[] = {0};

    float scale1[] = {1, 0}, scale2[] = {1, 0};
    int info1[] = {0}, info2[] = {0};

    if (*m > *n) {
        const int m1 = (*m >= 16) ? ((*m + 8) / 16) * 8 : *m / 2;
        const int m2 = *m - m1;

        // A_TL A_TR
        //      A_BR
        const float *const A_TL = A;
        const float *const A_TR = A + 2 * *ldA * m1;
        const float *const A_BR = A + 2 * *ldA * m1 + 2 * m1;

        // C_T
        // C_B
        float *const C_T = C;
        float *const C_B = C + 2 * m1;

        if (notransA) {
            // recusion(A_BR, B, C_B)
            LARPACK(ctrsyl)(tranA, tranB, isgn, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, scale1, info1);
            // C_T = C_T - A_TR * C_B
            BLAS(cgemm)("N", "N", &m1, n, &m2, cm1, A_TR, ldA, C_B, ldC, scale1, C_T, ldC);
            // recusion(A_TL, B, C_T)
            LARPACK(ctrsyl)(tranA, tranB, isgn, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(clascl)("G", i0, i0, c1, scale2, &m2, n, C_B, ldC, info);
        } else {
            // recusion(A_TL, B, C_T)
            LARPACK(ctrsyl)(tranA, tranB, isgn, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, scale1, info1);
            // C_B = C_B - A_TR' * C_T
            BLAS(cgemm)("C", "N", &m2, n, &m1, cm1, A_TR, ldA, C_T, ldC, scale1, C_B, ldC);
            // recusion(A_BR, B, C_B)
            LARPACK(ctrsyl)(tranA, tranB, isgn, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(clascl)("G", i0, i0, c1, scale2, &m1, n, C_B, ldC, info);
        }
    } else {
        const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
        const int n2 = *n - n1;

        // B_TL B_TR
        //      B_BR
        const float *const B_TL = B;
        const float *const B_TR = B + 2 * *ldB * n1;
        const float *const B_BR = B + 2 * *ldB * n1 + 2 * n1;

        // C_L C_R
        float *const C_L = C;
        float *const C_R = C + 2 * *ldC * n1;

        if (notransB) {
            // recusion(A, B_TL, C_L)
            LARPACK(ctrsyl)(tranA, tranB, isgn, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, scale1, info1);
            // C_R = C_R -/+ C_L * B_TR
            BLAS(cgemm)("N", "N", m, &n2, &n1, cmisgn, C_L, ldC, B_TR, ldB, scale1, C_R, ldC);
            // recusion(A, B_BR, C_R)
            LARPACK(ctrsyl)(tranA, tranB, isgn, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(clascl)("G", i0, i0, c1, scale2, m, &n1, C_L, ldC, info);
        } else {
            // recusion(A, B_BR, C_R)
            LARPACK(ctrsyl)(tranA, tranB, isgn, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, scale1, info1);
            // C_L = C_L -/+ C_R * B_TR'
            BLAS(cgemm)("N", "C", m, &n1, &n2, cmisgn, C_R, ldC, B_TR, ldB, scale1, C_L, ldC);
            // recusion(A, B_TL, C_L)
            LARPACK(ctrsyl)(tranA, tranB, isgn, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(clascl)("G", i0, i0, c1, scale2, m, &n2, C_R, ldC, info);
        }
    }

    *scale = scale1[0] * scale2[0];
    *info = info1[0] || info2[0];
}
