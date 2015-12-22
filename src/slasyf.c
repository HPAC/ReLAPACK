#include "relapack.h"

void RELAPACK(slasyf)(const char *uplo, const int *m, const int *n, int *nout,
        float *A, const int *ldA, int *ipiv, 
        float *Work, const int *ldWork, int *info) {

    if (*n <= CROSSOVER_SSYTRF) {
        // Unblocked
        if (*m == *n) {
            LAPACK(ssytf2)(uplo, n, A, ldA, ipiv, info);
            *nout = *n;
        } else
            LAPACK(ssytf3)(uplo, m, n, nout, A, ldA, ipiv, Work, ldWork, info);
        return;
    }

    // Recursive
    
    int info_1, info_2;

    // Constants
    // 1, -1
   	const float ONE[] = {1}, MONE[] = {-1};
    // 1
    const int iONE[] = {1};
    
    // Splitting (setup)
    int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;

    // A_TL *
    // *    *
    float *const A_TL = A;

    // ipiv_T
    // *
    int *const ipiv_T = ipiv;

    // Work_TL *
    // *       *
    float *const Work_TL = Work;

    // recursion(A_L)
    int n1out;
    RELAPACK(slasyf)(uplo, m, &n1, &n1out, A_TL, ldA, ipiv_T, Work_TL, ldWork, &info_1);
    n1 = n1out;

    if (n1 == 0) {
        // Unblocked
        if (*m == *n) {
            LAPACK(ssytf2)(uplo, n, A, ldA, ipiv, info);
            *nout = *n;
        } else
            LAPACK(ssytf3)(uplo, m, n, nout, A, ldA, ipiv, Work, ldWork, info);
        return;
    }

    // Splitting (continued)
    const int n2 = *n - n1;
    const int m2 = *m - n1;
    const int mmn = *m - *n;

    // A_TL  *
    // A_BL  A_BR
    // A_BL2 A_BR2
    float *const A_BL = A + n1;
    float *const A_BL2 = A + *n;
    float *const A_BR = A + *ldA * n1 + n1;
    float *const A_BR2 = A + *ldA * n1 + *n;

    // Work_TL *
    // Work_BL Work_BR
    // *       *
    float *const Work_BL = Work + n1;
    float *const Work_BR = Work + *ldWork * n1 + n1;

    // ipiv_T
    // ipiv_B
    int *const ipiv_B = ipiv + n1;

    // TODO: upper

    // A_BR = A_BR - A_BL Work_BL'
    RELAPACK(sgemm_tr)("N", "T", uplo, &n2, &n1, MONE, A_BL, ldA, Work_BL, ldWork, ONE, A_BR, ldA);
    BLAS(sgemm)("N", "T", &mmn, &n2, &n1, MONE, A_BL2, ldA, Work_BL, ldWork, ONE, A_BR2, ldA);

    // recursion(A_BR)
    int n2out;
    RELAPACK(slasyf)(uplo, &m2, &n2, &n2out, A_BR, ldA, ipiv_B, Work_BR, ldWork, &info_2);

    // shift pivots
    int i;
    for (i = 0; i < n2out; i++)
        if (ipiv_B[i] > 0)
            ipiv_B[i] += n1;
        else
            ipiv_B[i] -= n1;

    *nout = n1 + n2out;

    if (*nout != *n) {
        // undo 1 column of updates
        
        // last column of A_BR
        float *const A_BR_r = A_BR + *ldA * n2out + n2out;

        // last row of A_BL
        float *const A_BL_b = A_BL + n2out;

        // last row of Work_BL
        float *const Work_BL_b = Work_BL + n2out;

        // A_BR_r = A_BR_r + A_BL Work_BL'
        int mmnp1 = mmn + 1;
        BLAS(sgemv)("N", &mmnp1, &n1, ONE, A_BL_b, ldA, Work_BL_b, ldWork, ONE, A_BR_r, iONE);
    }

    // set return values
    *info = info_1 || info_2;
}
