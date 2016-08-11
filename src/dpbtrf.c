#include "relapack.h"
#include "stdlib.h"

static void RELAPACK_dpbtrf_rec(const char *, const int *, const int *,
    double *, const int *, double *, const int *, int *);


/** DPBTRF computes the Cholesky factorization of a real symmetric positive definite band matrix A.
 *
 * This routine is functionally equivalent to LAPACK's dpbtrf.
 * For details on its interface, see
 * http://www.netlib.org/lapack/explore-html/df/da9/dpbtrf_8f.html
 * */
void RELAPACK_dpbtrf(
    const char *uplo, const int *n, const int *kd,
    double *Ab, const int *ldAb,
    int *info
) {

    // Check arguments
    const int lower = LAPACK(lsame)(uplo, "L");
    const int upper = LAPACK(lsame)(uplo, "U");
    *info = 0;
    if (!lower && !upper)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*kd < 0)
        *info = -3;
    else if (*ldAb < *kd + 1)
        *info = -5;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("DPBTRF", &minfo);
        return;
    }

    // Clean char * arguments
    const char cleanuplo = lower ? 'L' : 'U';

    // Constant
    const double ZERO[] = { 0. };

    // Allocate work space
    const int n1 = REC_SPLIT(*n);
    const int mW = (*kd > n1) ? (lower ? *n - *kd : n1) : *kd;
    const int nW = (*kd > n1) ? (lower ? n1 : *n - *kd) : *kd;
    double *W = malloc(mW * nW * sizeof(double));
    LAPACK(dlaset)("G", &mW, &nW, ZERO, ZERO, W, &mW);

    // Recursive kernel
    RELAPACK_dpbtrf_rec(&cleanuplo, n, kd, Ab, ldAb, W, &mW, info);

    // Free work space
    free(W);
}


/** dpbtrf's recursive compute kernel */
static void RELAPACK_dpbtrf_rec(
    const char *uplo, const int *n, const int *kd,
    double *Ab, const int *ldAb,
    double *W, const int *ldW,
    int *info
){

    if (*n <= MAX(CROSSOVER_DPBTRF, 1)) {
        // Unblocked
        LAPACK(dpbtf2)(uplo, n, kd, Ab, ldAb, info);
        return;
    }

    // Constants
    const double ONE[]  = { 1. };
    const double MONE[] = { -1. };

    // Unskew A
    const int ldA[] = { *ldAb - 1 };
    double *const A = Ab + ((*uplo == 'L') ? 0 : *kd);

    // Splitting
    const int n1 = MIN(REC_SPLIT(*n), *kd);
    const int n2 = *n - n1;

    // * *
    // * Ab_BR
    double *const Ab_BR = Ab + *ldAb * n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + *ldA * n1;
    double *const A_BL = A             + n1;
    double *const A_BR = A + *ldA * n1 + n1;

    // recursion(A_TL)
    RELAPACK_dpotrf(uplo, &n1, A_TL, ldA, info);
    if (*info)
        return;

    // Banded splitting
    const int n21 = MIN(n2, *kd - n1);
    const int n22 = MIN(n2 - n21, n1);

    //     n1    n21    n22
    // n1  *     A_TRl  A_TRr
    // n21 A_BLt A_BRtl A_BRtr
    // n22 A_BLb A_BRbl A_BRbr
    double *const A_TRl  = A_TR;
    double *const A_TRr  = A_TR + *ldA * n21;
    double *const A_BLt  = A_BL;
    double *const A_BLb  = A_BL               + n21;
    double *const A_BRtl = A_BR;
    double *const A_BRtr = A_BR + *ldA * n21;
    double *const A_BRbl = A_BR               + n21;
    double *const A_BRbr = A_BR + *ldA * n21  + n21;

    if (*uplo == 'L') {
        // A_BLt = ABLt / A_TL'
        BLAS(dtrsm)("R", "L", "T", "N", &n21, &n1, ONE, A_TL, ldA, A_BLt, ldA);
        // A_BRtl = A_BRtl - A_BLt * A_BLt'
        BLAS(dsyrk)("L", "N", &n21, &n1, MONE, A_BLt, ldA, ONE, A_BRtl, ldA);
        // W = A_BLb
        LAPACK(dlacpy)("U", &n22, &n1, A_BLb, ldA, W, ldW);
        // W = W / A_TL'
        BLAS(dtrsm)("R", "L", "T", "N", &n22, &n1, ONE, A_TL, ldA, W, ldW);
        // A_BRbl = A_BRbl - W * A_BLt'
        BLAS(dgemm)("N", "T", &n22, &n21, &n1, MONE, W, ldW, A_BLt, ldA, ONE, A_BRbl, ldA);
        // A_BRbr = A_BRbr - W * W'
        BLAS(dsyrk)("L", "N", &n22, &n1, MONE, W, ldW, ONE, A_BRbr, ldA);
        // A_BLb = W
        LAPACK(dlacpy)("U", &n22, &n1, W, ldW, A_BLb, ldA);
    } else {
        // A_TRl = A_TL' \ A_TRl
        BLAS(dtrsm)("L", "U", "T", "N", &n1, &n21, ONE, A_TL, ldA, A_TRl, ldA);
        // A_BRtl = A_BRtl - A_TRl' * A_TRl
        BLAS(dsyrk)("U", "T", &n21, &n1, MONE, A_TRl, ldA, ONE, A_BRtl, ldA);
        // W = A_TRr
        LAPACK(dlacpy)("L", &n1, &n22, A_TRr, ldA, W, ldW);
        // W = A_TL' \ W
        BLAS(dtrsm)("L", "U", "T", "N", &n1, &n22, ONE, A_TL, ldA, W, ldW);
        // A_BRtr = A_BRtr - A_TRl' * W
        BLAS(dgemm)("T", "N", &n21, &n22, &n1, MONE, A_TRl, ldA, W, ldW, ONE, A_BRtr, ldA);
        // A_BRbr = A_BRbr - W' * W
        BLAS(dsyrk)("U", "T", &n22, &n1, MONE, W, ldW, ONE, A_BRbr, ldA);
        // A_TRr = W
        LAPACK(dlacpy)("L", &n1, &n22, W, ldW, A_TRr, ldA);
    }

    // recursion(A_BR)
    if (*kd > n1)
        RELAPACK_dpotrf(uplo, &n2, A_BR, ldA, info);
    else
        RELAPACK_dpbtrf_rec(uplo, &n2, kd, Ab_BR, ldAb, W, ldW, info);
    if (*info)
        *info += n1;
}
