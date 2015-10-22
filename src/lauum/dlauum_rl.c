#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void dlauum_rl(const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(dlauu2)("L", n, A, ldA, &info);
        return;
    }

   	const double d1[] = {1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
    double *A_TL = A;
    double *A_BL = A             + n1;
    double *A_BR = A + *ldA * n1 + n1;

    // A_TL = A_TL' * A_TL
    dlauum_rl(&n1, A_TL, ldA);
    // A_TL = A_TL + A_BL' * A_BL
    BLAS(dsyrk)("L", "T", &n1, &n2, d1, A_BL, ldA, d1, A_TL, ldA);
    // A_BL = A_BR' * A_BL
    BLAS(dtrmm)("L", "L", "T", "N", &n2, &n1, d1, A_BR, ldA, A_BL, ldA);
    // A_BR = A_BR' * A_BR
    dlauum_rl(&n2, A_BR, ldA);
}
