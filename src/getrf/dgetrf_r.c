#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"
#include "../util.h"

void dgetrf_r(const int *m, const int *n, double *A, const int *ldA, int *ipiv, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(dgetf2)(m, n, A, ldA, ipiv, info);
        return;
    }

   	const double s1[] = {1}, sm1[] = {-1};
    const int i1[] = {1};

    const int mn = MIN(*m, *n);
    const int n1 = (mn >= 16) ? ((mn + 8) / 16) * 8 : mn / 2;
    const int n2 = mn - n1;
    const int rm = *m - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *A_TL = A;
    double *A_TR = A + *ldA * n1;
    double *A_BL = A             + n1;
    double *A_BR = A + *ldA * n1 + n1;

    // ipiv_T
    // ipiv_B
    int *ipiv_T = ipiv;
    int *ipiv_B = ipiv + n1;

    // A_TL = LU(A_TL)
    dgetrf_r(m, &n1, A_TL, ldA, ipiv_T, info);
    if (*info)
        return;
    // apply pivots to A_TR
    LAPACK(dlaswp)(&n2, A_TR, ldA, i1, &n1, ipiv_T, i1);
    // A_TR = A_TL \ A_TR
    LAPACK(dtrsm)("L", "L", "N", "U", &n1, &n2, s1, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_BL * A_TR
    LAPACK(dgemm)("N", "N", &rm, &n2, &n1, sm1, A_BL, ldA, A_TR, ldA, s1, A_BR, ldA);
    // A_BR = LU(A_BR)
    dgetrf_r(&rm, &n2, A_BR, ldA, ipiv_B, info);
    if (*info) {
        *info += n1;
        return;
    }
    // apply pivots to A_BL
    LAPACK(dlaswp)(&n1, A_BL, ldA, i1, &n2, ipiv_B, i1);
    // shift pivots
    int i;
    for (i = 0; i < mn; i++)
        ipiv_B[i] += n1;

    const int rn = *n - mn;

    if (!rn)
        return;

    // A_S A_R
    double *A_S = A;
    double *A_R = A + *ldA * mn;

    // A_R = apply(ipiv, A_R)
    LAPACK(dlaswp)(&rn, A_R, ldA, i1, &mn, ipiv, i1);
    // A_R = A_S \ A_R
    LAPACK(dtrsm)("L", "L", "N", "U", &mn, &rn, s1, A_S, ldA, A_R, ldA);
}
