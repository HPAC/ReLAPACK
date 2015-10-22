#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void dpotrf_rl(const int *n, double *A, const int *ldA, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(dpotf2)("L", n, A, ldA, info);
        return;
    }

   	const double d1[] = {1}, dm1[] = {-1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
    double *A_TL = A;
    double *A_BL = A             + n1;
    double *A_BR = A + *ldA * n1 + n1;

    // A_TL = Chol(A_TL)
    dpotrf_rl(&n1, A_TL, ldA, info);
    if (*info)
        return;
    // A_BL = A_BL / A_TL'
    BLAS(dtrsm)("R", "L", "T", "N", &n2, &n1, d1, A_TL, ldA, A_BL, ldA);
    // A_BR = A_BR - A_BL * A_BL'
    BLAS(dsyrk)("L", "N", &n2, &n1, dm1, A_BL, ldA, d1, A_BR, ldA);
    // A_BR = Chol(A_BR)
    dpotrf_rl(&n2, A_BR, ldA, info);
    if (*info) {
        *info += n1;
        return;
    }
}
