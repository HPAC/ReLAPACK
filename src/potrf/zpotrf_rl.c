#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void zpotrf_rl(const int *n, double *A, const int *ldA, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(zpotf2)("L", n, A, ldA, info);
        return;
    }

   	const double z1[] = {1, 0}, zm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
#define A_TL (A)
#define A_BL (A                 + 2 * n1)
#define A_BR (A + 2 * *ldA * n1 + 2 * n1)

    // A_TL = Chol(A_TL)
    zpotrf_rl(&n1, A_TL, ldA, info);
    if (*info)
        return;
    // A_BL = A_BL / A_TL'
    BLAS(ztrsm)("R", "L", "C", "N", &n2, &n1, z1, A_TL, ldA, A_BL, ldA);
    // A_BR = A_BR - A_BL * A_BL'
    BLAS(zherk)("L", "N", &n2, &n1, zm1, A_BL, ldA, z1, A_BR, ldA);
    // A_BR = Chol(A_BR)
    zpotrf_rl(&n2, A_BR, ldA, info);
    if (*info) {
        *info += n1;
        return;
    }
}
