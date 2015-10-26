#include "larpack.h"

void LARPACK(strsyl)(const char *tranA, const char *tranB, const int *isgn, const int *m, const int *n, const float *A, const int *ldA, const float *B, const int *ldB, float *C, const int *ldC, float *scale, int *info) {
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
    if (*info != 0) {
        int minfo = -*info;
        LAPACK(xerbla)("STRSYL", &minfo);
        return;
    }

    *scale = 1;

    if (*m == 1 && *n == 1) { 
        // Single element (fully recursive)
        *C /= (*A + *isgn * *B);
        return;
    }

    // Recursive
    
    // 1, -1, -isgn
   	const float s1[] = {1}, sm1[] = {-1}, smisgn[] = {-*isgn};
   
    if (*m > *n) {
        const int m1 = (*m >= 16) ? ((*m + 8) / 16) * 8 : *m / 2;
        const int m2 = *m - m1;

        // A_TL A_TR
        //      A_BR
        const float *const A_TL = A;
        const float *const A_TR = A + *ldA * m1;
        const float *const A_BR = A + *ldA * m1 + m1;

        // C_T
        // C_B
        float *const C_T = C;
        float *const C_B = C + m1;

        if (notransA) {
            // recusion(A_BR, B, C_B)
            LARPACK(strsyl)(tranA, tranB, isgn, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, scale, info);
            // C_T = C_T - A_TR * C_B
            BLAS(sgemm)("N", "N", &m1, n, &m2, sm1, A_TR, ldA, C_B, ldC, s1, C_T, ldC);
            // recusion(A_TL, B, C_T)
            LARPACK(strsyl)(tranA, tranB, isgn, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, scale, info);
        } else {
            // recusion(A_TL, B, C_T)
            LARPACK(strsyl)(tranA, tranB, isgn, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, scale, info);
            // C_B = C_B - A_TR' * C_T
            BLAS(sgemm)("T", "N", &m2, n, &m1, sm1, A_TR, ldA, C_T, ldC, s1, C_B, ldC);
            // recusion(A_BR, B, C_B)
            LARPACK(strsyl)(tranA, tranB, isgn, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, scale, info);
        }
    } else {
        const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
        const int n2 = *n - n1;

        // B_TL B_TR
        //      B_BR
        const float *const B_TL = B;
        const float *const B_TR = B + *ldB * n1;
        const float *const B_BR = B + *ldB * n1 + n1;

        // C_L C_R
        float *const C_L = C;
        float *const C_R = C + *ldC * n1;

        if (notransB) {
            // recusion(A, B_TL, C_L)
            LARPACK(strsyl)(tranA, tranB, isgn, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, scale, info);
            // C_R = C_R -/+ C_L * B_TR
            BLAS(sgemm)("N", "N", m, &n2, &n1, smisgn, C_L, ldC, B_TR, ldB, s1, C_R, ldC);
            // recusion(A, B_BR, C_R)
            LARPACK(strsyl)(tranA, tranB, isgn, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, scale, info);
        } else {
            // recusion(A, B_BR, C_R)
            LARPACK(strsyl)(tranA, tranB, isgn, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, scale, info);
            // C_L = C_L -/+ C_R * B_TR'
            BLAS(sgemm)("N", "T", m, &n1, &n2, smisgn, C_R, ldC, B_TR, ldB, s1, C_L, ldC);
            // recusion(A, B_TL, C_L)
            LARPACK(strsyl)(tranA, tranB, isgn, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, scale, info);
        }
    }
}
