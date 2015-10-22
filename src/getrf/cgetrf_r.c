#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"
#include "../util.h"

void cgetrf_r(const int *m, const int *n, float *A, const int *ldA, int *ipiv, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(cgetf2)(m, n, A, ldA, ipiv, info);
        return;
    }

   	const float c1[] = {1, 0}, cm1[] = {-1, 0};
    const int i1[] = {1};

    const int mn = MIN(*m, *n);
    const int n1 = (mn >= 16) ? ((mn + 8) / 16) * 8 : mn / 2;
    const int n2 = mn - n1;
    const int rm = *m - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + 2 * *ldA * n1;
    float *const A_BL = A                 + 2 * n1;
    float *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // ipiv_T
    // ipiv_B
    int *const ipiv_T = ipiv;
    int *const ipiv_B = ipiv + n1;

    // A_TL = LU(A_TL)
    cgetrf_r(m, &n1, A_TL, ldA, ipiv_T, info);
    if (*info)
        return;
    // apply pivots to A_TR
    LAPACK(claswp)(&n2, A_TR, ldA, i1, &n1, ipiv_T, i1);
    // A_TR = A_TL \ A_TR
    LAPACK(ctrsm)("L", "L", "N", "U", &n1, &n2, c1, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_BL * A_TR
    LAPACK(cgemm)("N", "N", &rm, &n2, &n1, cm1, A_BL, ldA, A_TR, ldA, c1, A_BR, ldA);
    // A_BR = LU(A_BR)
    cgetrf_r(&rm, &n2, A_BR, ldA, ipiv_B, info);
    if (*info) {
        *info += n1;
        return;
    }
    // apply pivots to A_BL
    LAPACK(claswp)(&n1, A_BL, ldA, i1, &n2, ipiv_B, i1);
    // shift pivots
    int i;
    for (i = 0; i < mn; i++)
        ipiv_B[i] += n1;

    const int rn = *n - mn;

    if (!rn)
        return;

    // A_S A_R
    float *const A_S = A;
    float *const A_R = A + 2 * *ldA * mn;

    // A_R = apply(ipiv, A_R)
    LAPACK(claswp)(&rn, A_R, ldA, i1, &mn, ipiv, i1);
    // A_R = A_S \ A_R
    LAPACK(ctrsm)("L", "L", "N", "U", &mn, &rn, c1, A_S, ldA, A_R, ldA);
}
