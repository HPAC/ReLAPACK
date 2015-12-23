#include "relapack.h"

void RELAPACK(cgemm_tr)(const char *transA, const char *transB, const char *uplo,
        const int *n, const int *k,
        const float *alpha, const float *A, const int *ldA,
        const float *B, const int *ldB,
        const float *beta, float *C, const int *ldC) {

    const int tranA = LAPACK(lsame)(transA, "T") || LAPACK(lsame)(transA, "C");
    const int tranB = LAPACK(lsame)(transB, "T") || LAPACK(lsame)(transB, "C");
    const int lower = LAPACK(lsame)(uplo, "L");

    if (*n <= CROSSOVER_CGEMM_TR) {
        // Unblocked
        RELAPACK(cgemm_tr2)(transA, transB, uplo, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
        return;
    }

    // Recursion

    // Splitting
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_T
    // A_B
    const float *const A_T = A;
    const float *const A_B = A + 2 * (tranA ? *ldA * n1 : n1);

    // B_L B_R
    const float *const B_L = B;
    const float *const B_R = B + 2 * (tranB? n1 : *ldB * n1);

    // C_TL C_TR
    // C_BL C_BR
    float *const C_TL = C;
    float *const C_TR = C + 2 * *ldC * n1;
    float *const C_BL = C                 + 2 * n1;
    float *const C_BR = C + 2 * *ldC * n1 + 2 * n1;

    // recursion(C_TL)
    RELAPACK(cgemm_tr)(transA, transB, uplo, &n1, k, alpha, A_T, ldA, B_L, ldB, beta, C_TL, ldC);

    if (lower)
        // C_BL = alpha A_B B_L + beta C_BL
        BLAS(cgemm)(transA, transB, &n2, &n1, k, alpha, A_B, ldA, B_L, ldB, beta, C_BL, ldC);
    else
        // C_TR = alpha A_T B_R + beta C_TR
        BLAS(cgemm)(transA, transB, &n1, &n2, k, alpha, A_T, ldA, B_R, ldB, beta, C_TR, ldC);

    // recursion(C_BR)
    RELAPACK(cgemm_tr)(transA, transB, uplo, &n2, k, alpha, A_B, ldA, B_R, ldB, beta, C_BR, ldC);
}

void RELAPACK(cgemm_tr2)(const char *transA, const char *transB, const char *uplo,
        const int *n, const int *k,
        const float *alpha, const float *A, const int *ldA,
        const float *B, const int *ldB,
        const float *beta, float *C, const int *ldC) {

    const int tranA = LAPACK(lsame)(transA, "T") || LAPACK(lsame)(transA, "C");
    const int tranB = LAPACK(lsame)(transB, "T") || LAPACK(lsame)(transB, "C");
    const int lower = LAPACK(lsame)(uplo, "L");

    const int incB = tranB ? *ldB : 1;
    const int incC = 1;

    int i;
    for (i = 0; i < *n; i++) {
        // A_0
        // A_i
        const float *const A_0 = A;
        const float *const A_i = A + 2 * (tranA ? *ldA * i : i);

        // * B_i *
        const float *const B_i = B + 2 * (tranB ? i : *ldB * i);

        // * C_0i *
        // * C_ii *
        float *const C_0i = C + 2 * *ldC * i;
        float *const C_ii = C + 2 * *ldC * i + 2 * i;

        if (lower) {
            const int nmi = *n - i;
            if (tranA)
                BLAS(cgemv)(transA, k, &nmi, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
            else
                BLAS(cgemv)(transA, &nmi, k, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
        } else {
            const int ip1 = i + 1;
            if (tranA)
                BLAS(cgemv)(transA, k, &ip1, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
            else
                BLAS(cgemv)(transA, &ip1, k, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
        }
    }
}
