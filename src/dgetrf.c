#include "larpack.h"

void LARPACK(dgetrf)(const int *m, const int *n, double *A, const int *ldA, int *ipiv, int *info) {
    *info = 0;

    // Check arguments
    if (*m < 0)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*ldA < MAX(1, *n))
        *info = -4;
    if (*info != 0) {
        int minfo = -*info;
        LAPACK(xerbla)("DGETRF", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(dgetf2)(m, n, A, ldA, ipiv, info);
        return;
    }

    // Recursion
    const int mn = MIN(*m, *n);
    const int n1 = (mn >= 16) ? ((mn + 8) / 16) * 8 : mn / 2;
    const int n2 = mn - n1;
    const int rm = *m - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + *ldA * n1;
    double *const A_BL = A             + n1;
    double *const A_BR = A + *ldA * n1 + n1;

    // ipiv_T
    // ipiv_B
    int *const ipiv_T = ipiv;
    int *const ipiv_B = ipiv + n1;

    // 1, -1
   	const double d1[] = {1}, dm1[] = {-1};
    // 1
    const int i1[] = {1};

    // recursion(A_TL, ipiv_T)
    LARPACK(dgetrf)(m, &n1, A_TL, ldA, ipiv_T, info);
    // apply pivots to A_TR
    LAPACK(dlaswp)(&n2, A_TR, ldA, i1, &n1, ipiv_T, i1);

    // A_TR = A_TL \ A_TR
    BLAS(dtrsm)("L", "L", "N", "U", &n1, &n2, d1, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_BL * A_TR
    BLAS(dgemm)("N", "N", &rm, &n2, &n1, dm1, A_BL, ldA, A_TR, ldA, d1, A_BR, ldA);

    // recursion(A_BR, ipiv_B)
    LARPACK(dgetrf)(&rm, &n2, A_BR, ldA, ipiv_B, info);
    if (*info)
        *info += n1;
    // apply pivots to A_BL
    LAPACK(dlaswp)(&n1, A_BL, ldA, i1, &n2, ipiv_B, i1);
    // shift pivots
    int i;
    for (i = 0; i < mn; i++)
        ipiv_B[i] += n1;

    if (*n == mn)
        return;

    // Right remainder
    const int rn = *n - mn;

    // A_S A_R
    double *const A_S = A;
    double *const A_R = A + *ldA * mn;

    // A_R = apply(ipiv, A_R)
    LAPACK(dlaswp)(&rn, A_R, ldA, i1, &mn, ipiv, i1);
    // A_R = A_S \ A_R
    BLAS(dtrsm)("L", "L", "N", "U", &mn, &rn, d1, A_S, ldA, A_R, ldA);
}
