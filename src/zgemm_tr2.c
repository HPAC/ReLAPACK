#include "relapack.h"

void RELAPACK(zgemm_tr2)(const char *transA, const char *transB, const char *uplo,
        const int *n, const int *k,
        const double *alpha, const double *A, const int *ldA,
        const double *B, const int *ldB,
        const double *beta, double *C, const int *ldC) {

    const int tranA = LAPACK(lsame)(transA, "T") || LAPACK(lsame)(transA, "C");
    const int tranB = LAPACK(lsame)(transB, "T") || LAPACK(lsame)(transB, "C");
    const int lower = LAPACK(lsame)(uplo, "L");

    const int incB = tranB ? *ldB : 1;
    const int incC = 1;

    int i;
    for (i = 0; i < *n; i++) {
        // A_0
        // A_i
        const double *const A_0 = A;
        const double *const A_i = A + 2 * (tranA ? *ldA * i : i);

        // * B_i *
        const double *const B_i = B + 2 * (tranB ? i : *ldB * i);

        // * C_0i *
        // * C_ii *
        double *const C_0i = C + 2 * *ldC * i;
        double *const C_ii = C + 2 * *ldC * i + 2 * i;

        if (lower) {
            const int nmi = *n - i;
            if (tranA)
                BLAS(zgemv)(transA, k, &nmi, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
            else
                BLAS(zgemv)(transA, &nmi, k, alpha, A_i, ldA, B_i, &incB, beta, C_ii, &incC);
        } else {
            const int ip1 = i + 1;
            if (tranA)
                BLAS(zgemv)(transA, k, &ip1, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
            else
                BLAS(zgemv)(transA, &ip1, k, alpha, A_0, ldA, B_i, &incB, beta, C_0i, &incC);
        }
    }
}
