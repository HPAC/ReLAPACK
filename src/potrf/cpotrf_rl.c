#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void cpotrf_rl(const int *n, float *A, const int *ldA, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(cpotf2)("L", n, A, ldA, info);
        return;
    }

   	const float c1[] = {1, 0}, cm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
    float *A_TL = A;
    float *A_BL = A                 + 2 * n1;
    float *A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // A_TL = Chol(A_TL)
    cpotrf_rl(&n1, A_TL, ldA, info);
    if (*info)
        return;
    // A_BL = A_BL / A_TL'
    BLAS(ctrsm)("R", "L", "C", "N", &n2, &n1, c1, A_TL, ldA, A_BL, ldA);
    // A_BR = A_BR - A_BL * A_BL'
    BLAS(cherk)("L", "N", &n2, &n1, cm1, A_BL, ldA, c1, A_BR, ldA);
    // A_BR = Chol(A_BR)
    cpotrf_rl(&n2, A_BR, ldA, info);
    if (*info) {
        *info += n1;
        return;
    }
}
