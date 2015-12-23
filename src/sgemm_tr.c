#include "relapack.h"

void RELAPACK(sgemm_tr)(const char *transA, const char *transB, const char *uplo,
        const int *n, const int *k,
        const float *alpha, const float *A, const int *ldA,
        const float *B, const int *ldB,
        const float *beta, float *C, const int *ldC) {

    const int tranA = LAPACK(lsame)(transA, "T");
    const int tranB = LAPACK(lsame)(transB, "T");
    const int lower = LAPACK(lsame)(uplo, "L");

    if (*n <= CROSSOVER_SSYTRF) {
        // Unblocked
        RELAPACK(sgemm_tr2)(transA, transB, uplo, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
        return;
    }

    // Recursion

    // Splitting
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_T
    // A_B
    const float *const A_T = A;
    const float *const A_B = A + (tranA ? *ldA * n1 : n1);

    // B_L B_R
    const float *const B_L = B;
    const float *const B_R = B + (tranB? n1 : *ldB * n1);

    // C_TL C_TR
    // C_BL C_BR
    float *const C_TL = C;
    float *const C_TR = C + *ldC * n1;
    float *const C_BL = C             + n1;
    float *const C_BR = C + *ldC * n1 + n1;

    // recursion(C_TL)
    RELAPACK(sgemm_tr)(transA, transB, uplo, &n1, k, alpha, A_T, ldA, B_L, ldB, beta, C_TL, ldC);

    if (lower)
        // C_BL = alpha A_B B_L + beta C_BL
        BLAS(sgemm)(transA, transB, &n2, &n1, k, alpha, A_B, ldA, B_L, ldB, beta, C_BL, ldC);
    else
        // C_TR = alpha A_T B_R + beta C_TR
        BLAS(sgemm)(transA, transB, &n1, &n2, k, alpha, A_T, ldA, B_R, ldB, beta, C_TR, ldC);

    // recursion(C_BR)
    RELAPACK(sgemm_tr)(transA, transB, uplo, &n2, k, alpha, A_B, ldA, B_R, ldB, beta, C_BR, ldC);
}
