#include "relapack.h"

static void RELAPACK(dgemm_tr2)(const char *, const char *, const char *, 
    const int *, const int *, const double *, const double *, const int *, 
    const double *, const int *, const double *, double *, const int *);


void RELAPACK(dgemm_tr_rec)(
    const char *transA, const char *transB, const char *uplo,
    const int *n, const int *k,
    const double *alpha, const double *A, const int *ldA,
    const double *B, const int *ldB,
    const double *beta, double *C, const int *ldC
) {

    if (*n <= CROSSOVER_DGEMM_TR) {
        // Unblocked
        RELAPACK(dgemm_tr2)(transA, transB, uplo, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
        return;
    }

    // Splitting
    const int n1 = REC_SPLIT(*n);
    const int n2 = *n - n1;

    // A_T
    // A_B
    const double *const A_T = A;
    const double *const A_B = A + ((*transA == 'T') ? *ldA * n1 : n1);

    // B_L B_R
    const double *const B_L = B;
    const double *const B_R = B + ((*transB == 'T') ? n1 : *ldB * n1);

    // C_TL C_TR
    // C_BL C_BR
    double *const C_TL = C;
    double *const C_TR = C + *ldC * n1;
    double *const C_BL = C             + n1;
    double *const C_BR = C + *ldC * n1 + n1;

    // recursion(C_TL)
    RELAPACK(dgemm_tr_rec)(transA, transB, uplo, &n1, k, alpha, A_T, ldA, B_L, ldB, beta, C_TL, ldC);

    if (*uplo == 'L')
        // C_BL = alpha A_B B_L + beta C_BL
        BLAS(dgemm)(transA, transB, &n2, &n1, k, alpha, A_B, ldA, B_L, ldB, beta, C_BL, ldC);
    else
        // C_TR = alpha A_T B_R + beta C_TR
        BLAS(dgemm)(transA, transB, &n1, &n2, k, alpha, A_T, ldA, B_R, ldB, beta, C_TR, ldC);

    // recursion(C_BR)
    RELAPACK(dgemm_tr_rec)(transA, transB, uplo, &n2, k, alpha, A_B, ldA, B_R, ldB, beta, C_BR, ldC);
}

static void RELAPACK(dgemm_tr2)(
    const char *transA, const char *transB, const char *uplo, 
    const int *n, const int *k,
    const double *alpha, const double *A, const int *ldA,
    const double *B, const int *ldB,
    const double *beta, double *C, const int *ldC
) {

    const int incB = (*transB == 'T') ? *ldB : 1;
    const int incC = 1;

    int i;
    for (i = 0; i < *n; i++) {
        // A_0
        // A_i
        const double *const A_0 = A;
        const double *const A_i = A + ((*transA == 'T') ? *ldA * i : i);

        // * B_i *
        const double *const B_i = B + ((*transB == 'T') ? i : *ldB * i);

        // * C_0i *
        // * C_ii *
        double *const C_0i = C + *ldC * i;
        double *const C_ii = C + *ldC * i + i;

        if (*uplo == 'L') {
            const int nmi = *n - i;
            if (*transA == 'T')
                BLAS(dgemv)(transA, k, &nmi, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
            else
                BLAS(dgemv)(transA, &nmi, k, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
        } else {
            const int ip1 = i + 1;
            if (*transA == 'T')
                BLAS(dgemv)(transA, k, &ip1, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
            else
                BLAS(dgemv)(transA, &ip1, k, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
        }
    }
}
