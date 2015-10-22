#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void clauum_ru(const int *n, float *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(clauu2)("U", n, A, ldA, &info);
        return;
    }

   	const float c1[] = {1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
    float *A_TL = A;
    float *A_TR = A + 2 * *ldA * n1;
    float *A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // A_TL = A_TL' * A_TL
    clauum_ru(&n1, A_TL, ldA);
    // A_TL = A_TL + A_TR * A_TR'
    BLAS(cherk)("U", "N", &n1, &n2, c1, A_TR, ldA, c1, A_TL, ldA);
    // A_TR = A_TR * A_BR'
    BLAS(ctrmm)("R", "U", "C", "N", &n1, &n2, c1, A_BR, ldA, A_TR, ldA);
    // A_BR = A_BR' * A_BR
    clauum_ru(&n2, A_BR, ldA);
}
