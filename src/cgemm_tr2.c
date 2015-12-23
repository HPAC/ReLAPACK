#include "relapack.h"

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
