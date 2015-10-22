#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void spotrf_ru(const int *n, float *A, const int *ldA, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(spotf2)("U", n, A, ldA, info);
        return;
    }

   	const float s1[] = {1}, sm1[] = {-1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
    float *A_TL = A;
    float *A_TR = A + *ldA * n1;
    float *A_BR = A + *ldA * n1 + n1;

    // A_TL = Chol(A_TL)
    spotrf_ru(&n1, A_TL, ldA, info);
    if (*info)
        return;
    // A_TR = A_TL' \ A_TR
    BLAS(strsm)("L", "U", "T", "N", &n1, &n2, s1, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_TR' * A_TR
    BLAS(ssyrk)("U", "T", &n2, &n1, sm1, A_TR, ldA, s1, A_BR, ldA);
    // A_BR = Chol(A_BR)
    spotrf_ru(&n2, A_BR, ldA, info);
    if (*info) {
        *info += n1;
        return;
    }
}
