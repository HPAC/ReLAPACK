#include "larpack.h"

void LARPACK(sgetrf)(const int *m, const int *n,
        float *A, const int *ldA, int *ipiv, int *info) {

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
        LAPACK(xerbla)("SGETRF", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(sgetf2)(m, n, A, ldA, ipiv, info);
        return;
    }

    // Recursive

    // Constants
    // 1, -1
   	const float s1 = 1, sm1 = -1;
    // 1
    const int i1 = 1;

    // Splitting
    const int mn = MIN(*m, *n);
    const int n1 = (mn >= 16) ? ((mn + 8) / 16) * 8 : mn / 2;
    const int n2 = mn - n1;
    const int rm = *m - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + *ldA * n1;
    float *const A_BL = A             + n1;
    float *const A_BR = A + *ldA * n1 + n1;

    // ipiv_T
    // ipiv_B
    int *const ipiv_T = ipiv;
    int *const ipiv_B = ipiv + n1;

    // recursion(A_TL, ipiv_T)
    LARPACK(sgetrf)(m, &n1, A_TL, ldA, ipiv_T, info);
    // apply pivots to A_TR
    LAPACK(slaswp)(&n2, A_TR, ldA, &i1, &n1, ipiv_T, &i1);

    // A_TR = A_TL \ A_TR
    BLAS(strsm)("L", "L", "N", "U", &n1, &n2, &s1, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_BL * A_TR
    BLAS(sgemm)("N", "N", &rm, &n2, &n1, &sm1, A_BL, ldA, A_TR, ldA, &s1, A_BR, ldA);

    // recursion(A_BR, ipiv_B)
    LARPACK(sgetrf)(&rm, &n2, A_BR, ldA, ipiv_B, info);
    if (*info)
        *info += n1;
    // apply pivots to A_BL
    LAPACK(slaswp)(&n1, A_BL, ldA, &i1, &n2, ipiv_B, &i1);
    // shift pivots
    int i;
    for (i = 0; i < mn; i++)
        ipiv_B[i] += n1;

    if (*n == mn)
        return;

    // Right remainder
    const int rn = *n - mn;

    // A_S A_R
    float *const A_S = A;
    float *const A_R = A + *ldA * mn;

    // A_R = apply(ipiv, A_R)
    LAPACK(slaswp)(&rn, A_R, ldA, &i1, &mn, ipiv, &i1);
    // A_R = A_S \ A_R
    BLAS(strsm)("L", "L", "N", "U", &mn, &rn, &s1, A_S, ldA, A_R, ldA);
}
