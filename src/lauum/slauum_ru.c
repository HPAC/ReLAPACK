#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void slauum_ru(const int *n, float *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(slauu2)("U", n, A, ldA, &info);
        return;
    }

   	const float s1[] = {1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
#define A_TL (A)
#define A_TR (A + *ldA * n1)
#define A_BR (A + *ldA * n1 + n1)

    // A_TL = A_TL' * A_TL
    slauum_ru(&n1, A_TL, ldA);
    // A_TL = A_TL + A_TR * A_TR'
    BLAS(ssyrk)("U", "N", &n1, &n2, s1, A_TR, ldA, s1, A_TL, ldA);
    // A_TR = A_TR * A_BR'
    BLAS(strmm)("R", "U", "T", "N", &n1, &n2, s1, A_BR, ldA, A_TR, ldA);
    // A_BR = A_BR' * A_BR
    slauum_ru(&n2, A_BR, ldA);
}
