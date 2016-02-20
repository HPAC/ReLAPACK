#include "relapack.h"

static void RELAPACK(zgetrf_rec)(const int *, const int *, double *,
    const int *, int *, int *);


/** ZGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
 *
 * This routine is functionally equivalent to LAPACK's zgetrf.
 * For details on its interface, see
 * http://www.netlib.org/lapack/explore-html/dd/dd1/zgetrf_8f.html
 * */
void RELAPACK(zgetrf)(
    const int *m, const int *n,
    double *A, const int *ldA, int *ipiv,
    int *info
) {

    // Check arguments
    *info = 0;
    if (*m < 0)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*ldA < MAX(1, *n))
        *info = -4;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("ZGETRF", &minfo);
        return;
    }

    RELAPACK(zgetrf_rec)(m, n, A, ldA, ipiv, info);

    // Right remainder
    if (*m < *n) {
        // Constants
        const double ONE[]  = {1, 0};
        const int   iONE[] = {1};

        // Splitting
        const int rn = *n - *m;

        // A_L A_R
        const double *const A_L = A;
        double *const       A_R = A + 2 * *ldA * *m;

        // A_R = apply(ipiv, A_R)
        LAPACK(zlaswp)(&rn, A_R, ldA, iONE, m, ipiv, iONE);
        // A_R = A_L \ A_R
        BLAS(ztrsm)("L", "L", "N", "U", m, &rn, ONE, A_L, ldA, A_R, ldA);
    }
}


/** zgetrf's recursive compute kernel */
static void RELAPACK(zgetrf_rec)(
    const int *m, const int *n,
    double *A, const int *ldA, int *ipiv,
    int *info
) {

    if (*n <= MAX(CROSSOVER_ZGETRF, 1)) {
        // Unblocked
        LAPACK(zgetf2)(m, n, A, ldA, ipiv, info);
        return;
    }

    // Constants
    const double ONE[]  = {1, 0};
    const double MONE[] = {-1, 0};
    const int   iONE[]  = {1};

    // Splitting
    const int mn = MIN(*m, *n);
    const int n1 = REC_SPLIT(mn);
    const int n2 = mn - n1;
    const int rm = *m - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + 2 * *ldA * n1;
    double *const A_BL = A                 + 2 * n1;
    double *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // ipiv_T
    // ipiv_B
    int *const ipiv_T = ipiv;
    int *const ipiv_B = ipiv + n1;

    // recursion(A_TL, ipiv_T)
    RELAPACK(zgetrf_rec)(m, &n1, A_TL, ldA, ipiv_T, info);
    // apply pivots to A_TR
    LAPACK(zlaswp)(&n2, A_TR, ldA, iONE, &n1, ipiv_T, iONE);

    // A_TR = A_TL \ A_TR
    BLAS(ztrsm)("L", "L", "N", "U", &n1, &n2, ONE, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_BL * A_TR
    BLAS(zgemm)("N", "N", &rm, &n2, &n1, MONE, A_BL, ldA, A_TR, ldA, ONE, A_BR, ldA);

    // recursion(A_BR, ipiv_B)
    RELAPACK(zgetrf_rec)(&rm, &n2, A_BR, ldA, ipiv_B, info);
    if (*info)
        *info += n1;
    // apply pivots to A_BL
    LAPACK(zlaswp)(&n1, A_BL, ldA, iONE, &n2, ipiv_B, iONE);
    // shift pivots
    for (int i = 0; i < n2; i++)
        ipiv_B[i] += n1;
}
