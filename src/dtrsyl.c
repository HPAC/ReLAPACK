#include "relapack.h"

static void RELAPACK(dtrsyl_rec)(const char *, const char *, const int *, 
    const int *, const int *, const double *, const int *, const double *, 
    const int *, double *, const int *, double *, int *);


void RELAPACK(dtrsyl)(
    const char *tranA, const char *tranB, const int *isgn,
    const int *m, const int *n,
    const double *A, const int *ldA, const double *B, const int *ldB,
    double *C, const int *ldC, double *scale, 
    int *info
) {

    // Check arguments
    const int notransA = LAPACK(lsame)(tranA, "N");
    const int transA = LAPACK(lsame)(tranA, "T");
    const int ctransA = LAPACK(lsame)(tranA, "C");
    const int notransB = LAPACK(lsame)(tranB, "N");
    const int transB = LAPACK(lsame)(tranB, "T");
    const int ctransB = LAPACK(lsame)(tranB, "C");
    *info = 0;
    if (!transA && !ctransA && !notransA)
        *info = -1;
    else if (!transB && !ctransB && !notransB)
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
        LAPACK(xerbla)("DTRSYL", &minfo);
        return;
    }

    // Clean char * arguments
    const char cleantranA = notransA ? 'N' : (transA ? 'T' : 'C');
    const char cleantranB = notransB ? 'N' : (transB ? 'T' : 'C');

    RELAPACK(dtrsyl_rec)(&cleantranA, &cleantranB, isgn, m, n, A, ldA, B, ldB, C, ldC, scale, info);
}


static void RELAPACK(dtrsyl_rec)(
    const char *tranA, const char *tranB, const int *isgn,
    const int *m, const int *n,
    const double *A, const int *ldA, const double *B, const int *ldB,
    double *C, const int *ldC, double *scale, 
    int *info
) {

    if (*m <= CROSSOVER_DTRSYL && *n <= CROSSOVER_DTRSYL) {
        // Unblocked
        DTRSY2(tranA, tranB, isgn, m, n, A, ldA, B, ldB, C, ldC, scale, info);
        return;
    }

    // Constants
    // 1, -1, -isgn
   	const double ONE[] = {1}, MONE[] = {-1}, MSGN[] = {-*isgn};
    // 0
    const int iONE[] = {1};

    // Outputs
    double scale1[1], scale2[1];
    int info1[1], info2[1];

    if (*m > *n) {
        int m1 = REC_SPLIT(*m);
        if (A[m1 + *ldA * (m1 - 1)])
            m1++;
        const int m2 = *m - m1;

        // A_TL A_TR
        //      A_BR
        const double *const A_TL = A;
        const double *const A_TR = A + *ldA * m1;
        const double *const A_BR = A + *ldA * m1 + m1;

        // C_T
        // C_B
        double *const C_T = C;
        double *const C_B = C + m1;

        if (*tranA == 'N') {
            // recusion(A_BR, B, C_B)
            RELAPACK(dtrsyl_rec)(tranA, tranB, isgn, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, scale1, info1);
            // C_T = C_T - A_TR * C_B
            BLAS(dgemm)("N", "N", &m1, n, &m2, MONE, A_TR, ldA, C_B, ldC, scale1, C_T, ldC);
            // recusion(A_TL, B, C_T)
            RELAPACK(dtrsyl_rec)(tranA, tranB, isgn, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(dlascl)("G", iONE, iONE, ONE, scale2, &m2, n, C_B, ldC, info);
        } else {
            // recusion(A_TL, B, C_T)
            RELAPACK(dtrsyl_rec)(tranA, tranB, isgn, &m1, n, A_TL, ldA, B, ldB, C_T, ldC, scale1, info1);
            // C_B = C_B - A_TR' * C_T
            BLAS(dgemm)("C", "N", &m2, n, &m1, MONE, A_TR, ldA, C_T, ldC, scale1, C_B, ldC);
            // recusion(A_BR, B, C_B)
            RELAPACK(dtrsyl_rec)(tranA, tranB, isgn, &m2, n, A_BR, ldA, B, ldB, C_B, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(dlascl)("G", iONE, iONE, ONE, scale2, &m1, n, C_B, ldC, info);
        }
    } else {
        int n1 = REC_SPLIT(*n);
        if (B[n1 + *ldB * (n1 - 1)])
            n1++;
        const int n2 = *n - n1;

        // B_TL B_TR
        //      B_BR
        const double *const B_TL = B;
        const double *const B_TR = B + *ldB * n1;
        const double *const B_BR = B + *ldB * n1 + n1;

        // C_L C_R
        double *const C_L = C;
        double *const C_R = C + *ldC * n1;

        if (*tranB == 'N') {
            // recusion(A, B_TL, C_L)
            RELAPACK(dtrsyl_rec)(tranA, tranB, isgn, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, scale1, info1);
            // C_R = C_R -/+ C_L * B_TR
            BLAS(dgemm)("N", "N", m, &n2, &n1, MSGN, C_L, ldC, B_TR, ldB, scale1, C_R, ldC);
            // recusion(A, B_BR, C_R)
            RELAPACK(dtrsyl_rec)(tranA, tranB, isgn, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(dlascl)("G", iONE, iONE, ONE, scale2, m, &n1, C_L, ldC, info);
        } else {
            // recusion(A, B_BR, C_R)
            RELAPACK(dtrsyl_rec)(tranA, tranB, isgn, m, &n2, A, ldA, B_BR, ldB, C_R, ldC, scale1, info1);
            // C_L = C_L -/+ C_R * B_TR'
            BLAS(dgemm)("N", "C", m, &n1, &n2, MSGN, C_R, ldC, B_TR, ldB, scale1, C_L, ldC);
            // recusion(A, B_TL, C_L)
            RELAPACK(dtrsyl_rec)(tranA, tranB, isgn, m, &n1, A, ldA, B_TL, ldB, C_L, ldC, scale2, info2);
            // apply scale
            if (scale2[0] != 1)
                LAPACK(dlascl)("G", iONE, iONE, ONE, scale2, m, &n2, C_R, ldC, info);
        }
    }

    *scale = scale1[0] * scale2[0];
    *info = info1[0] || info2[0];
}
