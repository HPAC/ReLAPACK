#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void slauum_rl(const int *n, float *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(slauu2)("L", n, A, ldA, &info);
        return;
    }

   	const float s1[] = {1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
#define A_TL (A)
#define A_BL (A             + n1)
#define A_BR (A + *ldA * n1 + n1)

    // A_TL = A_TL' * A_TL
    slauum_rl(&n1, A_TL, ldA);
    // A_TL = A_TL + A_BL' * A_BL
    BLAS(ssyrk)("L", "T", &n1, &n2, s1, A_BL, ldA, s1, A_TL, ldA);
    // A_BL = A_BR' * A_BL
    BLAS(strmm)("L", "L", "T", "N", &n2, &n1, s1, A_BR, ldA, A_BL, ldA);
    // A_BR = A_BR' * A_BR
    slauum_rl(&n2, A_BR, ldA);
}
