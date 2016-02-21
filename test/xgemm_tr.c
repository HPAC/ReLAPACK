#include "test.h"

datatype *A[2], *B[2], *C[2], *Ctmp;

void pre() {
    x2matgen(n, n, A[0], A[1]);
    x2matgen(n, n, B[0], B[1]);
    x2matgen(n, n, C[0], C[1]);
}

void post() {
    error = x2vecerr(n * n, C[0], C[1]);
}

#define ROUTINE XPREF(gemm_tr_rec)

#define xlacpy XPREF(LAPACK(lacpy))

void XLAPACK(ROUTINE)(
    const char *transA, const char *transB, const char *uplo,
    const int *n, const int *k,
    const datatype *alpha, const datatype *A, const int *ldA,
    const datatype *B, const int *ldB,
    const datatype *beta, datatype *C, const int *ldC
) {
    xlacpy(uplo, n, n, C, ldC, Ctmp, n);
    XPREF(BLAS(gemm))(transA, transB, n, n, k, alpha, A, ldA, B, ldB, beta, Ctmp, n);
    xlacpy(uplo, n, n, Ctmp, ldC, C, n);
}

void tests() {
    A[0] = xmalloc(n * n);
    A[1] = xmalloc(n * n);
    B[0] = xmalloc(n * n);
    B[1] = xmalloc(n * n);
    C[0] = xmalloc(n * n);
    C[1] = xmalloc(n * n);
    Ctmp = xmalloc(n * n);

    TEST("N", "N", "L", &n, &n, ONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("N", "N", "L", &n, &n, ONE, A[i], &n, B[i], &n, MONE, C[i], &n);
    TEST("N", "N", "L", &n, &n, MONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("N", "T", "L", &n, &n, ONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("T", "N", "L", &n, &n, ONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("N", "N", "L", &n, &n2, ONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("N", "N", "U", &n, &n, ONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("N", "N", "U", &n, &n, ONE, A[i], &n, B[i], &n, MONE, C[i], &n);
    TEST("N", "N", "U", &n, &n, MONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("N", "T", "U", &n, &n, ONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("T", "N", "U", &n, &n, ONE, A[i], &n, B[i], &n, ONE, C[i], &n);
    TEST("N", "N", "U", &n, &n2, ONE, A[i], &n, B[i], &n, ONE, C[i], &n);

    free(A[0]);
    free(A[1]);
    free(B[0]);
    free(B[1]);
    free(C[0]);
    free(C[1]);
    free(Ctmp);
}
