#include "relapack.h"
#include "stdlib.h"

static void RELAPACK_zpbtrf_rec(const char *, const int *, const int *,
    double *, const int *, double *, const int *, int *);


/** ZPBTRF computes the Cholesky factorization of a complex Hermitian positive definite band matrix A.
 *
 * This routine is functionally equivalent to LAPACK's zpbtrf.
 * For details on its interface, see
 * http://www.netlib.org/lapack/explore-html/db/da9/zpbtrf_8f.html
 * */
void RELAPACK_zpbtrf(
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
        LAPACK(xerbla)("ZPBTRF", &minfo);
        return;
    }

    // Clean char * arguments
    const char cleanuplo = lower ? 'L' : 'U';

    // Constant
    const double ZERO[] = { 0., 0. };

    // Allocate work space
    const int n1 = REC_SPLIT(*n);
    const int mW = (*kd > n1) ? (lower ? *n - *kd : n1) : *kd;
    const int nW = (*kd > n1) ? (lower ? n1 : *n - *kd) : *kd;
    double *W = malloc(mW * nW * 2 * sizeof(double));
    LAPACK(zlaset)("G", &mW, &nW, ZERO, ZERO, W, &mW);

    // Recursive kernel
    RELAPACK_zpbtrf_rec(&cleanuplo, n, kd, Ab, ldAb, W, &mW, info);

    // Free work space
    free(W);
}


/** zpbtrf's recursive compute kernel */
static void RELAPACK_zpbtrf_rec(
    const char *uplo, const int *n, const int *kd,
    double *Ab, const int *ldAb,
    double *W, const int *ldW,
    int *info
){

    if (*n <= MAX(CROSSOVER_ZPBTRF, 1)) {
        // Unblocked
        LAPACK(zpbtf2)(uplo, n, kd, Ab, ldAb, info);
        return;
    }

    // Constants
    const double ONE[]  = { 1., 0. };
    const double MONE[] = { -1., 0. };

    // Unskew A
    const int ldA[] = { *ldAb - 1 };
    double *const A = Ab + 2 * ((*uplo == 'L') ? 0 : *kd);

    // Splitting
    const int n1 = MIN(REC_SPLIT(*n), *kd);
    const int n2 = *n - n1;

    // * *
    // * Ab_BR
    double *const Ab_BR = Ab + 2 * *ldAb * n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + 2 * *ldA * n1;
    double *const A_BL = A                 + 2 * n1;
    double *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // recursion(A_TL)
    RELAPACK_zpotrf(uplo, &n1, A_TL, ldA, info);
    if (*info)
        return;

    // Banded splitting
    const int n21 = MIN(n2, *kd - n1);
    const int n22 = MIN(n2 - n21, *kd);

    //     n1    n21    n22
    // n1  *     A_TRl  A_TRr
    // n21 A_BLt A_BRtl A_BRtr
    // n22 A_BLb A_BRbl A_BRbr
    double *const A_TRl  = A_TR;
    double *const A_TRr  = A_TR + 2 * *ldA * n21;
    double *const A_BLt  = A_BL;
    double *const A_BLb  = A_BL                   + 2 * n21;
    double *const A_BRtl = A_BR;
    double *const A_BRtr = A_BR + 2 * *ldA * n21;
    double *const A_BRbl = A_BR                   + 2 * n21;
    double *const A_BRbr = A_BR + 2 * *ldA * n21  + 2 * n21;

    if (*uplo == 'L') {
        // A_BLt = ABLt / A_TL'
        BLAS(ztrsm)("R", "L", "C", "N", &n21, &n1, ONE, A_TL, ldA, A_BLt, ldA);
        // A_BRtl = A_BRtl - A_BLt * A_BLt'
        BLAS(zherk)("L", "N", &n21, &n1, MONE, A_BLt, ldA, ONE, A_BRtl, ldA);
        // W = A_BLb
        LAPACK(zlacpy)("U", &n22, &n1, A_BLb, ldA, W, ldW);
        // W = W / A_TL'
        BLAS(ztrsm)("R", "L", "C", "N", &n22, &n1, ONE, A_TL, ldA, W, ldW);
        // A_BRbl = A_BRbl - W * A_BLt'
        BLAS(zgemm)("N", "C", &n22, &n21, &n1, MONE, W, ldW, A_BLt, ldA, ONE, A_BRbl, ldA);
        // A_BRbr = A_BRbr - W * W'
        BLAS(zherk)("L", "N", &n22, &n1, MONE, W, ldW, ONE, A_BRbr, ldA);
        // A_BLb = W
        LAPACK(zlacpy)("U", &n22, &n1, W, ldW, A_BLb, ldA);
    } else {
        // A_TRl = A_TL' \ A_TRl
        BLAS(ztrsm)("L", "U", "C", "N", &n1, &n21, ONE, A_TL, ldA, A_TRl, ldA);
        // A_BRtl = A_BRtl - A_TRl' * A_TRl
        BLAS(zherk)("U", "C", &n21, &n1, MONE, A_TRl, ldA, ONE, A_BRtl, ldA);
        // W = A_TRr
        LAPACK(zlacpy)("L", &n1, &n22, A_TRr, ldA, W, ldW);
        // W = A_TL' \ W
        BLAS(ztrsm)("L", "U", "C", "N", &n1, &n22, ONE, A_TL, ldA, W, ldW);
        // A_BRtr = A_BRtr - A_TRl' * W
        BLAS(zgemm)("C", "N", &n21, &n22, &n1, MONE, A_TRl, ldA, W, ldW, ONE, A_BRtr, ldA);
        // A_BRbr = A_BRbr - W' * W
        BLAS(zherk)("U", "C", &n22, &n1, MONE, W, ldW, ONE, A_BRbr, ldA);
        // A_TRr = W
        LAPACK(zlacpy)("L", &n1, &n22, W, ldW, A_TRr, ldA);
    }

    // recursion(A_BR)
    if (*kd > n1)
        RELAPACK_zpotrf(uplo, &n2, A_BR, ldA, info);
    else
        RELAPACK_zpbtrf_rec(uplo, &n2, kd, Ab_BR, ldAb, W, ldW, info);
    if (*info)
        *info += n1;
}
