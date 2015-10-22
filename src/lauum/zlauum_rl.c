#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void zlauum_rl(const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(zlauu2)("L", n, A, ldA, &info);
        return;
    }

   	const double z1[] = {1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_BL = A                 + 2 * n1;
    double *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // A_TL = A_TL' * A_TL
    zlauum_rl(&n1, A_TL, ldA);
    // A_TL = A_TL + A_BL' * A_BL
    BLAS(zherk)("L", "C", &n1, &n2, z1, A_BL, ldA, z1, A_TL, ldA);
    // A_BL = A_BR' * A_BL
    BLAS(ztrmm)("L", "L", "C", "N", &n2, &n1, z1, A_BR, ldA, A_BL, ldA);
    // A_BR = A_BR' * A_BR
    zlauum_rl(&n2, A_BR, ldA);
}
