#include "relapack.h"

void RELAPACK(ssytrf_rec)(const char *uplo,
        const int *n_full, const int *n, int *n_out,
        float *A, const int *ldA, int *ipiv,
        float *Work, const int *ldWork, int *info) {

    if (*n <= MAX(CROSSOVER_SSYTRF, 3)) {
        // Unblocked
        if (*n == *n_full) {
            LAPACK(ssytf2)(uplo, n, A, ldA, ipiv, info);
            *n_out = *n;
        } else
            LAPACK(ssytrf_rec2)(uplo, n_full, n, n_out, A, ldA, ipiv, Work, ldWork, info);
        return;
    }

    // Recursive

    const int lower = LAPACK(lsame)(uplo, "L");
    int info1, info2;

    // Constants
    // 1, -1
   	const float ONE[] = {1}, MONE[] = {-1};
    // 1
    const int iONE[] = {1};

    const int n_rest = *n_full - *n;

    if (lower) {
        // Splitting (setup)
        int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
        int n2 = *n - n1;

        // Work_L *
        float *const Work_L = Work;

        // recursion(A_L)
        int n1_out;
        RELAPACK(ssytrf_rec)(uplo, n_full, &n1, &n1_out, A, ldA, ipiv, Work_L, ldWork, &info1);
        n1 = n1_out;

        // Splitting (continued)
        n2 = *n - n1;
        const int n_full2 = *n_full - n1;

        // *      *
        // A_BL   A_BR
        // A_BL_B A_BR_B
        float *const A_BL = A + n1;
        float *const A_BL_B = A + *n;
        float *const A_BR = A + *ldA * n1 + n1;
        float *const A_BR_B = A + *ldA * n1 + *n;

        // *        *
        // Work_BL Work_BR
        // *       *
        float *const Work_BL = Work + n1;
        float *const Work_BR = Work + *ldWork * n1 + n1;

        // ipiv_T
        // ipiv_B
        int *const ipiv_B = ipiv + n1;

        // A_BR = A_BR - A_BL Work_BL'
        RELAPACK(sgemm_tr)("N", "T", uplo, &n2, &n1, MONE, A_BL, ldA, Work_BL, ldWork, ONE, A_BR, ldA);
        BLAS(sgemm)("N", "T", &n_rest, &n2, &n1, MONE, A_BL_B, ldA, Work_BL, ldWork, ONE, A_BR_B, ldA);

        // recursion(A_BR)
        int n2_out;
        RELAPACK(ssytrf_rec)(uplo, &n_full2, &n2, &n2_out, A_BR, ldA, ipiv_B, Work_BR, ldWork, &info2);

        if (n2_out != n2) {
            // undo 1 column of updates
            const int n_restp1 = n_rest + 1;

            // last column of A_BR
            float *const A_BR_r = A_BR + *ldA * n2_out + n2_out;

            // last row of A_BL
            float *const A_BL_b = A_BL + n2_out;

            // last row of Work_BL
            float *const Work_BL_b = Work_BL + n2_out;

            // A_BR_r = A_BR_r + A_BL_b Work_BL_b'
            BLAS(sgemv)("N", &n_restp1, &n1, ONE, A_BL_b, ldA, Work_BL_b, ldWork, ONE, A_BR_r, iONE);
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
        float *const Work_R = Work + *ldWork * n1;

        // recursion(A_L)
        int n2_out;
        RELAPACK(ssytrf_rec)(uplo, n_full, &n2, &n2_out, A, ldA, ipiv, Work_R, ldWork, &info2);
        n2 = n2_out;

        // Splitting (continued)
        n1 = *n - n2;
        const int n_full1 = *n_full - n2;

        // * A_TL_T A_TR_T
        // * A_TL   A_TR
        // * *      *
        float *const A_TL_T = A + *ldA * n_rest;
        float *const A_TR_T = A + *ldA * (n_rest + n1);
        float *const A_TL   = A + *ldA * n_rest        + n_rest;
        float *const A_TR   = A + *ldA * (n_rest + n1) + n_rest;

        // Work_L *
        // *      Work_TR
        // *      *
        float *const Work_L = Work;
        float *const Work_TR = Work + *ldWork * n1 + n_rest;

        // A_TL = A_TL - A_TR Work_TR'
        RELAPACK(sgemm_tr)("N", "T", uplo, &n1, &n2, MONE, A_TR, ldA, Work_TR, ldWork, ONE, A_TL, ldA);
        BLAS(sgemm)("N", "T", &n_rest, &n1, &n2, MONE, A_TR_T, ldA, Work_TR, ldWork, ONE, A_TL_T, ldA);

        // recursion(A_TL)
        int n1_out;
        RELAPACK(ssytrf_rec)(uplo, &n_full1, &n1, &n1_out, A, ldA, ipiv, Work_L, ldWork, &info1);

        if (n1_out != n1) {
            // undo 1 column of updates
            const int n_restp1 = n_rest + 1;

            // A_TL_T(l) = A_TL_T(l) + A_TR_T Work_TR(t)'
            BLAS(sgemv)("N", &n_restp1, &n2, ONE, A_TR_T, ldA, Work_TR, ldWork, ONE, A_TL_T, iONE);
        }
        n1 = n1_out;

        *info = info2 || info1;
        *n_out = n1 + n2;
    }
}
