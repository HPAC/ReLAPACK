#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void clauum_rl(const int *n, float *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(clauu2)("L", n, A, ldA, &info);
        return;
    }

   	const float c1[] = {1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_BL = A                 + 2 * n1;
    float *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // A_TL = A_TL' * A_TL
    clauum_rl(&n1, A_TL, ldA);
    // A_TL = A_TL + A_BL' * A_BL
    BLAS(cherk)("L", "C", &n1, &n2, c1, A_BL, ldA, c1, A_TL, ldA);
    // A_BL = A_BR' * A_BL
    BLAS(ctrmm)("L", "L", "C", "N", &n2, &n1, c1, A_BR, ldA, A_BL, ldA);
    // A_BR = A_BR' * A_BR
    clauum_rl(&n2, A_BR, ldA);
}
