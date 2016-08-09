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

    const double ZERO[] = { 0. };
    const int nW = REC_SPLIT(*n) + 8;
    double *W = malloc(nW * nW * sizeof(double));
    LAPACK(dlaset)("G", &nW, &nW, ZERO, ZERO, W, &nW);

    RELAPACK_dpbtrf_rec(&cleanuplo, n, kd, Ab, ldAb, W, &nW, info);

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

    // Splitting
    const int n1 = REC_SPLIT(*n);
    const int n2 = *n - n1;

    // Ab_TL *
    // *     Ab_BR
    double *const Ab_TL = Ab;
    double *const Ab_BR = Ab + *ldAb * n1;

    // Unskewing A
    const int ldA[] = { *ldAb - 1 };
    double *const A = Ab + ((*uplo == 'L') ? 0 : *kd);

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + *ldA * n1;
    double *const A_BL = A             + n1;
    double *const A_BR = A + *ldA * n1 + n1;

    // recursion(A_TL)
    RELAPACK_dpbtrf_rec(uplo, &n1, kd, Ab_TL, ldAb, W, ldW, info);

    if (*kd > n1) {  // Band is larger than n1
        // Banded splitting
        const int n21 = MIN(n2, *kd - n1);
        const int n22 = n2 - n21;

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
    } else {  // Band is smaller than n1
        // Banded splitting
        const int n11 = n1 - *kd;

        //     n11 kd     kd
        // n11 *   *      0      0
        // kd  *   A_TLbr A_TRbl 0
        // kd  0   A_BLtr A_BRtl *
        //     0   0      *      *
        double *const A_TLbr = A_TL + *ldA * n11 + n11;
        double *const A_TRbl = A_TR              + n11;
        double *const A_BLtr = A_BL + *ldA * n11;
        double *const A_BRtl = A_BR;

        if (*uplo == 'L') {
            // W = A_BLtr
            LAPACK(dlacpy)("U", kd, kd, A_BLtr, ldA, W, ldW);
            // W = W / A_TLbr'
            BLAS(dtrsm)("R", "L", "T", "N", kd, kd, ONE, A_TLbr, ldA, W, ldW);
            // A_BRtl = A_BRtl - W * W'
            BLAS(dsyrk)("L", "N", kd, kd, MONE, W, ldW, ONE, A_BRtl, ldA);
            // A_BLtr = W
            LAPACK(dlacpy)("U", kd, kd, W, ldW, A_BLtr, ldA);
        } else {
            // W = A_TRbl
            LAPACK(dlacpy)("L", kd, kd, A_TRbl, ldA, W, ldW);
            // W = A_TLbr' \ W
            BLAS(dtrsm)("L", "U", "T", "N", kd, kd, ONE, A_TLbr, ldA, W, ldW);
            // A_BRtl = A_BRtl - W' * W
            BLAS(dsyrk)("U", "T", kd, kd, MONE, W, ldW, ONE, A_BRtl, ldA);
            // A_TRbl = W
            LAPACK(dlacpy)("L", kd, kd, W, ldW, A_TRbl, ldA);
        }
    }

    // recursion(A_BR)
    RELAPACK_dpbtrf_rec(uplo, &n2, kd, Ab_BR, ldAb, W, ldW, info);
    if (*info)
        *info += n1;
}
