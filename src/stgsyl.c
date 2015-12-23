#include "relapack.h"
#include <math.h>

void RELAPACK(stgsyl)(const char *trans, const int *ijob, 
        const int *m, const int *n,
        const float *A, const int *ldA, const float *B, const int *ldB,
        float *C, const int *ldC, 
        const float *D, const int *ldD, const float *E, const int *ldE, 
        float *F, const int *ldF,
        float *scale, float *dif, 
        float *Work, const int *lWork, int *iWork, int *info) {
    // Parse arguments
    const int notran = LAPACK(lsame)(trans, "N");
    const int tran = LAPACK(lsame)(trans, "T");

    // Compute work buffer size
    int lwmin = 1;
    if (notran && (*ijob == 1 || *ijob == 2))
        lwmin = MAX(1, 2 * *m * *n);
    *info = 0;

    // Check arguments
    if (!tran && !notran)
        *info = -1;
    else if (notran && (*ijob < 0 || *ijob > 4))
        *info = -2;
    else if (*m <= 0)
        *info = -3;
    else if (*n <= 0)
        *info = -4;
    else if (*ldA < MAX(1, *m))
        *info = -6;
    else if (*ldB < MAX(1, *n))
        *info = -8;
    else if (*ldC < MAX(1, *m))
        *info = -10;
    else if (*ldD < MAX(1, *m))
        *info = -12;
    else if (*ldE < MAX(1, *n))
        *info = -14;
    else if (*ldF < MAX(1, *m))
        *info = -16;
    else if (*lWork < lwmin && *lWork != -1)
        *info = -20;

    // Errors and work space queries
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("STGSYL", &minfo);
        return;
    } else if (*lWork == -1) {
        *Work = lwmin;
        return;
    }
    
    // Constants
    // 0
   	const float ZERO[] = {0};
    
    int isolve = 1;
    int ifunc = 0;
    if (notran) {
        if (*ijob >= 3) {
            ifunc = *ijob - 2;
            LAPACK(slaset)("F", m, n, ZERO, ZERO, C, ldC);
            LAPACK(slaset)("F", m, n, ZERO, ZERO, F, ldF);
        } else if (*ijob >= 1)
            isolve = 2;
    }

    float scale2;
    for (int iround = 1; iround <= isolve; iround++) {
        *scale = 1;
        float dscale = 0;
        float dsum = 1;
        int pq;
        RELAPACK(stgsyl_rec)(trans, &ifunc, m, n, A, ldA, B, ldB, C, ldC, D, ldD, E, ldE, F, ldF, scale, &dsum, &dscale, iWork, &pq, info);
        if (dscale != 0) {
            if (*ijob == 1 || *ijob == 3)
                *dif = sqrt(2 * *m * *n) / (dscale * sqrt(dsum));
            else
                *dif = sqrt(pq) / (dscale * sqrt(dsum));
        }
        if (isolve == 2) {
            if (iround == 1) {
                if (notran)
                    ifunc = *ijob;
                scale2 = *scale;
                LAPACK(slacpy)("F", m, n, C, ldC, Work, m);
                LAPACK(slacpy)("F", m, n, F, ldF, Work + *m * *n, m);
                LAPACK(slaset)("F", m, n, ZERO, ZERO, C, ldC);
                LAPACK(slaset)("F", m, n, ZERO, ZERO, F, ldF);
            } else {
                LAPACK(slacpy)("F", m, n, Work, m, C, ldC);
                LAPACK(slacpy)("F", m, n, Work + *m * *n, m, F, ldF);
                *scale = scale2;
            }
        }
    }
}

void RELAPACK(stgsyl_rec)(const char *trans, const int *ifunc,
        const int *m, const int *n,
        const float *A, const int *ldA, const float *B, const int *ldB,
        float *C, const int *ldC, 
        const float *D, const int *ldD, const float *E, const int *ldE, 
        float *F, const int *ldF,
        float *scale, float *dsum, float *dscale, int *iWork, int *pq, 
        int *info) {

    if (*m <= CROSSOVER_STGSYL && *n <= CROSSOVER_STGSYL) {
        // Unblocked
        LAPACK(stgsy2)(trans, ifunc, m, n, A, ldA, B, ldB, C, ldC, D, ldD, E, ldE, F, ldF, scale, dsum, dscale, iWork, pq, info);
        return;
    }

    // Parse arguments
    const int notran = LAPACK(lsame)(trans, "N");

    // Recursive

    // Constants
    // 1, -1
   	const float ONE[] = {1, 0}, MONE[] = {-1, 0};
    // 0
    const int iONE[] = {1};

    // Outputs
    float scale1[] = {1, 0}, scale2[] = {1, 0};
    int info1[1], info2[1];

    if (*m > *n) {
        int m1 = (*m >= 16) ? ((*m + 8) / 16) * 8 : *m / 2;
        if (A[m1 + *ldA * (m1 - 1)])
            m1++;
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

        // D_TL D_TR
        //      D_BR
        const float *const D_TL = D;
        const float *const D_TR = D + *ldD * m1;
        const float *const D_BR = D + *ldD * m1 + m1;

        // F_T
        // F_B
        float *const F_T = F;
        float *const F_B = F + m1;

        if (notran) {
            // recursion(A_BR, B, C_B, D_BR, E, F_B)
            RELAPACK(stgsyl_rec)(trans, ifunc, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, D_BR, ldD, E, ldE, F_B, ldF, scale1, dsum, dscale, iWork, pq, info1);
            // C_T = C_T - A_TR * C_B
            BLAS(sgemm)("N", "N", &m1, n, &m2, MONE, A_TR, ldA, C_B, ldC, scale1, C_T, ldC);
            // F_T = F_T - D_TR * C_B
            BLAS(sgemm)("N", "N", &m1, n, &m2, MONE, D_TR, ldD, C_B, ldC, scale1, F_T, ldF);
            // recursion(A_TL, B, C_T, D_TL, E, F_T)
            RELAPACK(stgsyl_rec)(trans, ifunc, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, D_TL, ldD, E, ldE, F_T, ldF, scale2, dsum, dscale, iWork, pq, info2);
            // apply scale
            if (scale2[0] != 1) {
                LAPACK(slascl)("G", iONE, iONE, ONE, scale2, &m2, n, C_B, ldC, info);
                LAPACK(slascl)("G", iONE, iONE, ONE, scale2, &m2, n, F_B, ldF, info);
            }
        } else {
            // recursion(A_TL, B, C_T, D_TL, E, F_T)
            RELAPACK(stgsyl_rec)(trans, ifunc, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, D_TL, ldD, E, ldE, F_T, ldF, scale1, dsum, dscale, iWork, pq, info1);
            // apply scale
            if (scale1[0] != 1)
                LAPACK(slascl)("G", iONE, iONE, ONE, scale1, &m2, n, F_B, ldF, info);
            // C_B = C_B - A_TR^H * C_T
            BLAS(sgemm)("T", "N", &m2, n, &m1, MONE, A_TR, ldA, C_T, ldC, scale1, C_B, ldC);
            // C_B = C_B - D_TR^H * F_T
            BLAS(sgemm)("T", "N", &m2, n, &m1, MONE, D_TR, ldD, F_T, ldC, ONE, C_B, ldC);
            // recursion(A_BR, B, C_B, D_BR, E, F_B)
            RELAPACK(stgsyl_rec)(trans, ifunc, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, D_BR, ldD, E, ldE, F_B, ldF, scale2, dsum, dscale, iWork, pq, info2);
            // apply scale
            if (scale2[0] != 1) {
                LAPACK(slascl)("G", iONE, iONE, ONE, scale2, &m1, n, C_T, ldC, info);
                LAPACK(slascl)("G", iONE, iONE, ONE, scale2, &m1, n, F_T, ldF, info);
            }
        }
    } else {
        int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
        if (B[n1 + *ldB * (n1 - 1)])
            n1++;
        const int n2 = *n - n1;

        // B_TL B_TR
        //      B_BR
        const float *const B_TL = B;
        const float *const B_TR = B + *ldB * n1;
        const float *const B_BR = B + *ldB * n1 + n1;

        // C_L C_R
        float *const C_L = C;
        float *const C_R = C + *ldC * n1;

        // E_TL E_TR
        //      E_BR
        const float *const E_TL = E;
        const float *const E_TR = E + *ldE * n1;
        const float *const E_BR = E + *ldE * n1 + n1;

        // F_L F_R
        float *const F_L = F;
        float *const F_R = F + *ldF * n1;

        if (notran) {
            // recursion(A, B_TL, C_L, D, E_TL, F_L)
            RELAPACK(stgsyl_rec)(trans, ifunc, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, D, ldD, E_TL, ldE, F_L, ldF, scale1, dsum, dscale, iWork, pq, info1);
            // C_R = C_R + F_L * B_TR
            BLAS(sgemm)("N", "N", m, &n2, &n1, ONE, F_L, ldF, B_TR, ldB, scale1, C_R, ldC);
            // F_R = F_R + F_L * E_TR
            BLAS(sgemm)("N", "N", m, &n2, &n1, ONE, F_L, ldF, E_TR, ldE, scale1, F_R, ldF);
            // recursion(A, B_BR, C_R, D, E_BR, F_R)
            RELAPACK(stgsyl_rec)(trans, ifunc, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, D, ldD, E_BR, ldE, F_R, ldF, scale2, dsum, dscale, iWork, pq, info2);
            // apply scale
            if (scale2[0] != 1) {
                LAPACK(slascl)("G", iONE, iONE, ONE, scale2, m, &n1, C_L, ldC, info);
                LAPACK(slascl)("G", iONE, iONE, ONE, scale2, m, &n1, F_L, ldF, info);
            }
        } else {
            // recursion(A, B_BR, C_R, D, E_BR, F_R)
            RELAPACK(stgsyl_rec)(trans, ifunc, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, D, ldD, E_BR, ldE, F_R, ldF, scale1, dsum, dscale, iWork, pq, info1);
            // apply scale
            if (scale1[0] != 1)
                LAPACK(slascl)("G", iONE, iONE, ONE, scale1, m, &n1, C_L, ldC, info);
            // F_L = F_L + C_R * B_TR
            BLAS(sgemm)("N", "T", m, &n1, &n2, ONE, C_R, ldC, B_TR, ldB, scale1, F_L, ldF);
            // F_L = F_L + F_R * E_TR
            BLAS(sgemm)("N", "T", m, &n1, &n2, ONE, F_R, ldF, E_TR, ldB, ONE, F_L, ldF);
            // recursion(A, B_TL, C_L, D, E_TL, F_L)
            RELAPACK(stgsyl_rec)(trans, ifunc, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, D, ldD, E_TL, ldE, F_L, ldF, scale2, dsum, dscale, iWork, pq, info2);
            // apply scale
            if (scale2[0] != 1) {
                LAPACK(slascl)("G", iONE, iONE, ONE, scale2, m, &n2, C_R, ldC, info);
                LAPACK(slascl)("G", iONE, iONE, ONE, scale2, m, &n2, F_R, ldF, info);
            }
        }
    }

    *scale = scale1[0] * scale2[0];
    *info = info1[0] || info2[0];
}
