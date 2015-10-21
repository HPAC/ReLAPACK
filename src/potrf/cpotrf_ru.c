#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void cpotrf_ru(const int *n, float *A, const int *ldA, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(cpotf2)("U", n, A, ldA, info);
        return;
    }

   	const float c1[] = {1, 0}, cm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
#define A_TL (A)
#define A_TR (A + 2 * *ldA * n1)
#define A_BR (A + 2 * *ldA * n1 + 2 * n1)

    // A_TL = Chol(A_TL)
    cpotrf_ru(&n1, A_TL, ldA, info);
    if (*info)
        return;
    // A_TR = A_TL' \ A_TR
    BLAS(ctrsm)("L", "U", "C", "N", &n1, &n2, c1, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_TR' * A_TR
    BLAS(cherk)("U", "C", &n2, &n1, cm1, A_TR, ldA, c1, A_BR, ldA);
    // A_BR = Chol(A_BR)
    cpotrf_ru(&n2, A_BR, ldA, info);
    if (*info) {
        *info += n1;
        return;
    }
}
