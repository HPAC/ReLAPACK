#include "relapack.h"
#include <stdlib.h>

void RELAPACK(zhetrf)(const char *uplo, const int *n,
        double *A, const int *ldA, int *ipiv,
        double *Work, const int *lWork, int *info) {

    // Check arguments
    const int lower = LAPACK(lsame)(uplo, "L");
    const int upper = LAPACK(lsame)(uplo, "U");
    *info = 0;
    if (!lower && !upper)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*ldA < MAX(1, *n))
        *info = -4;
    else if (*lWork < 1 && *lWork != -1)
        *info = -7;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("ZHETRF", &minfo);
        return;
    }

    if (*lWork == -1) {
        // Work size query
        *Work = *ldA * *n;
        return;
    }

    double *W = Work;
    if (*lWork < *ldA * *n)
        W = malloc(*ldA * *n * 2 * sizeof(double));

    int nout;
    RELAPACK(zhetrf_rec)(uplo, n, n, &nout, A, ldA, ipiv, W, ldA, info);

    if (W != Work)
        free(W);
}

void RELAPACK(zhetrf_rec)(const char *uplo,
        const int *n_full, const int *n, int *n_out,
        double *A, const int *ldA, int *ipiv,
        double *Work, const int *ldWork, int *info) {

    if (*n <= MAX(CROSSOVER_ZHETRF, 3)) {
        // Unblocked
        if (*n == *n_full) {
            LAPACK(zhetf2)(uplo, n, A, ldA, ipiv, info);
            *n_out = *n;
        } else
            LAPACK(zhetrf_rec2)(uplo, n_full, n, n_out, A, ldA, ipiv, Work, ldWork, info);
        return;
    }

    // Recursive

    const int lower = LAPACK(lsame)(uplo, "L");
    int info1, info2;

    // Constants
    // 1, -1
   	const double ONE[] = {1, 0}, MONE[] = {-1, 0};
    // 1
    const int iONE[] = {1};

    const int n_rest = *n_full - *n;

    if (lower) {
        // Splitting (setup)
        int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
        int n2 = *n - n1;

        // Work_L *
        double *const Work_L = Work;

        // recursion(A_L)
        int n1_out;
        RELAPACK(zhetrf_rec)(uplo, n_full, &n1, &n1_out, A, ldA, ipiv, Work_L, ldWork, &info1);
        n1 = n1_out;

        // Splitting (zontinued)
        n2 = *n - n1;
        const int n_full2 = *n_full - n1;

        // *      *
        // A_BL   A_BR
        // A_BL_B A_BR_B
        double *const A_BL   = A                 + 2 * n1;
        double *const A_BR   = A + 2 * *ldA * n1 + 2 * n1;
        double *const A_BL_B = A                 + 2 * *n;
        double *const A_BR_B = A + 2 * *ldA * n1 + 2 * *n;

        // *        *
        // Work_BL Work_BR
        // *       *
        double *const Work_BL = Work                    + 2 * n1;
        double *const Work_BR = Work + 2 * *ldWork * n1 + 2 * n1;

        // ipiv_T
        // ipiv_B
        int *const ipiv_B = ipiv + n1;

        // A_BR = A_BR - A_BL Work_BL'
        RELAPACK(zgemm_tr)("N", "T", uplo, &n2, &n1, MONE, A_BL, ldA, Work_BL, ldWork, ONE, A_BR, ldA);
        BLAS(zgemm)("N", "T", &n_rest, &n2, &n1, MONE, A_BL_B, ldA, Work_BL, ldWork, ONE, A_BR_B, ldA);

        // recursion(A_BR)
        int n2_out;
        RELAPACK(zhetrf_rec)(uplo, &n_full2, &n2, &n2_out, A_BR, ldA, ipiv_B, Work_BR, ldWork, &info2);

        if (n2_out != n2) {
            // undo 1 column of updates
            const int n_restp1 = n_rest + 1;

            // last column of A_BR
            double *const A_BR_r = A_BR + 2 * *ldA * n2_out + 2 * n2_out;

            // last row of A_BL
            double *const A_BL_b = A_BL + 2 * n2_out;

            // last row of Work_BL
            double *const Work_BL_b = Work_BL + 2 * n2_out;

            // A_BR_r = A_BR_r + A_BL_b Work_BL_b'
            BLAS(zgemv)("N", &n_restp1, &n1, ONE, A_BL_b, ldA, Work_BL_b, ldWork, ONE, A_BR_r, iONE);
        }
        n2 = n2_out;

        // shift pivots
        int i;
        for (i = 0; i < n2; i++)
            if (ipiv_B[i] > 0)
                ipiv_B[i] += n1;
            else
                ipiv_B[i] -= n1;

        *info = info1 || info2;
        *n_out = n1 + n2;
    } else {
        // Splitting (setup)
        int n2 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
        int n1 = *n - n2;

        // * Work_R
        double *const Work_R = Work + 2 * *ldWork * n1;

        // recursion(A_L)
        int n2_out;
        RELAPACK(zhetrf_rec)(uplo, n_full, &n2, &n2_out, A, ldA, ipiv, Work_R, ldWork, &info2);
        n2 = n2_out;

        // Splitting (zontinued)
        n1 = *n - n2;
        const int n_full1 = *n_full - n2;

        // * A_TL_T A_TR_T
        // * A_TL   A_TR
        // * *      *
        double *const A_TL_T = A + 2 * *ldA * n_rest;
        double *const A_TR_T = A + 2 * *ldA * (n_rest + n1);
        double *const A_TL   = A + 2 * *ldA * n_rest        + 2 * n_rest;
        double *const A_TR   = A + 2 * *ldA * (n_rest + n1) + 2 * n_rest;

        // Work_L *
        // *      Work_TR
        // *      *
        double *const Work_L  = Work;
        double *const Work_TR = Work + 2 * *ldWork * n1 + 2 * n_rest;

        // A_TL = A_TL - A_TR Work_TR'
        RELAPACK(zgemm_tr)("N", "T", uplo, &n1, &n2, MONE, A_TR, ldA, Work_TR, ldWork, ONE, A_TL, ldA);
        BLAS(zgemm)("N", "T", &n_rest, &n1, &n2, MONE, A_TR_T, ldA, Work_TR, ldWork, ONE, A_TL_T, ldA);

        // recursion(A_TL)
        int n1_out;
        RELAPACK(zhetrf_rec)(uplo, &n_full1, &n1, &n1_out, A, ldA, ipiv, Work_L, ldWork, &info1);

        if (n1_out != n1) {
            // undo 1 column of updates
            const int n_restp1 = n_rest + 1;

            // A_TL_T(l) = A_TL_T(l) + A_TR_T Work_TR(t)'
            BLAS(zgemv)("N", &n_restp1, &n2, ONE, A_TR_T, ldA, Work_TR, ldWork, ONE, A_TL_T, iONE);
        }
        n1 = n1_out;

        *info = info2 || info1;
        *n_out = n1 + n2;
    }
}
