#include "relapack.h"

void RELAPACK(sgemm_tr2)(const char *transA, const char *transB, const char *uplo,
        const int *n, const int *k,
        const float *alpha, const float *A, const int *ldA,
        const float *B, const int *ldB,
        const float *beta, float *C, const int *ldC) {

    const int tranA = LAPACK(lsame)(transA, "T");
    const int tranB = LAPACK(lsame)(transB, "T");
    const int lower = LAPACK(lsame)(uplo, "L");

    const int incB = tranB ? *ldB : 1;
    const int incC = 1;

    int i;
    for (i = 0; i < *n; i++) {
        // A_0
        // A_i
        const float *const A_0 = A;
        const float *const A_i = A + (tranA ? *ldA * i : i);

        // * B_i *
        const float *const B_i = B + (tranB ? i : *ldB * i);

        // * C_0i *
        // * C_ii *
        float *const C_0i = C + *ldC * i;
        float *const C_ii = C + *ldC * i + i;

        if (lower) {
            const int nmi = *n - i;
            BLAS(sgemv)(transA, &nmi, k, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
        } else {
            const int ip1 = i + 1;
            BLAS(sgemv)(transA, &ip1, k, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
        }
    }
}
