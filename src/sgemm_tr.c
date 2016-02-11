#include "relapack.h"

static void RELAPACK(sgemm_tr2)(const char *, const char *, const char *,
    const int *, const int *, const float *, const float *, const int *,
    const float *, const int *, const float *, float *, const int *);


void RELAPACK(sgemm_tr_rec)(
    const char *transA, const char *transB, const char *uplo,
    const int *n, const int *k,
    const float *alpha, const float *A, const int *ldA,
    const float *B, const int *ldB,
    const float *beta, float *C, const int *ldC
) {

    if (*n <= MAX(CROSSOVER_SGEMM_TR, 1)) {
        // Unblocked
        RELAPACK(sgemm_tr2)(transA, transB, uplo, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
        return;
    }

    // Splitting
    const int n1 = REC_SPLIT(*n);
    const int n2 = *n - n1;

    // A_T
    // A_B
    const float *const A_T = A;
    const float *const A_B = A + ((*transA == 'T') ? *ldA * n1 : n1);

    // B_L B_R
    const float *const B_L = B;
    const float *const B_R = B + ((*transB == 'T') ? n1 : *ldB * n1);

    // C_TL C_TR
    // C_BL C_BR
    float *const C_TL = C;
    float *const C_TR = C + *ldC * n1;
    float *const C_BL = C             + n1;
    float *const C_BR = C + *ldC * n1 + n1;

    // recursion(C_TL)
    RELAPACK(sgemm_tr_rec)(transA, transB, uplo, &n1, k, alpha, A_T, ldA, B_L, ldB, beta, C_TL, ldC);

    if (*uplo == 'L')
        // C_BL = alpha A_B B_L + beta C_BL
        BLAS(sgemm)(transA, transB, &n2, &n1, k, alpha, A_B, ldA, B_L, ldB, beta, C_BL, ldC);
    else
        // C_TR = alpha A_T B_R + beta C_TR
        BLAS(sgemm)(transA, transB, &n1, &n2, k, alpha, A_T, ldA, B_R, ldB, beta, C_TR, ldC);

    // recursion(C_BR)
    RELAPACK(sgemm_tr_rec)(transA, transB, uplo, &n2, k, alpha, A_B, ldA, B_R, ldB, beta, C_BR, ldC);
}

static void RELAPACK(sgemm_tr2)(
    const char *transA, const char *transB, const char *uplo,
    const int *n, const int *k,
    const float *alpha, const float *A, const int *ldA,
    const float *B, const int *ldB,
    const float *beta, float *C, const int *ldC
) {

    const int incB = (*transB == 'T') ? *ldB : 1;
    const int incC = 1;

    for (int i = 0; i < *n; i++) {
        // A_0
        // A_i
        const float *const A_0 = A;
        const float *const A_i = A + ((*transA == 'T') ? *ldA * i : i);

        // * B_i *
        const float *const B_i = B + ((*transB == 'T') ? i : *ldB * i);

        // * C_0i *
        // * C_ii *
        float *const C_0i = C + *ldC * i;
        float *const C_ii = C + *ldC * i + i;

        if (*uplo == 'L') {
            const int nmi = *n - i;
            if (*transA == 'T')
                BLAS(sgemv)(transA, k, &nmi, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
            else
                BLAS(sgemv)(transA, &nmi, k, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
        } else {
            const int ip1 = i + 1;
            if (*transA == 'T')
                BLAS(sgemv)(transA, k, &ip1, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
            else
                BLAS(sgemv)(transA, &ip1, k, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
        }
    }
}
