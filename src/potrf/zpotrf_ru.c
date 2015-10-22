#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void zpotrf_ru(const int *n, double *A, const int *ldA, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(zpotf2)("U", n, A, ldA, info);
        return;
    }

   	const double z1[] = {1, 0}, zm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
    double *const A_TL = A;
    double *const A_TR = A + 2 * *ldA * n1;
    double *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // A_TL = Chol(A_TL)
    zpotrf_ru(&n1, A_TL, ldA, info);
    if (*info)
        return;
    // A_TR = A_TL' \ A_TR
    BLAS(ztrsm)("L", "U", "C", "N", &n1, &n2, z1, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_TR' * A_TR
    BLAS(zherk)("U", "C", &n2, &n1, zm1, A_TR, ldA, z1, A_BR, ldA);
    // A_BR = Chol(A_BR)
    zpotrf_ru(&n2, A_BR, ldA, info);
    if (*info) {
        *info += n1;
        return;
    }
}
