#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void zlauum_ru(const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(zlauu2)("U", n, A, ldA, &info);
        return;
    }

   	const double z1[] = {1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
#define A_TL (A)
#define A_TR (A + 2 * *ldA * n1)
#define A_BR (A + 2 * *ldA * n1 + 2 * n1)

    // A_TL = A_TL' * A_TL
    zlauum_ru(&n1, A_TL, ldA);
    // A_TL = A_TL + A_TR * A_TR'
    BLAS(zherk)("U", "N", &n1, &n2, z1, A_TR, ldA, z1, A_TL, ldA);
    // A_TR = A_TR * A_BR'
    BLAS(ztrmm)("R", "U", "C", "N", &n1, &n2, z1, A_BR, ldA, A_TR, ldA);
    // A_BR = A_BR' * A_BR
    zlauum_ru(&n2, A_BR, ldA);
}
