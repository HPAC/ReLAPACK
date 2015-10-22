#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void dtrtri_rl(const char *diag, const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(dtrti2)("L", diag, n, A, ldA, &info);
        return;
    }

   	const double d1[] = {1}, dm1[] = {-1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_BL = A             + n1;
    double *const A_BR = A + *ldA * n1 + n1;

    // A_TL = 1 / A_TL
    dtrtri_rl(diag, &n1, A_TL, ldA);
    // A_BL = - A_BL * A_TL
    BLAS(dtrmm)("R", "L", "N", diag, &n2, &n1, dm1, A_TL, ldA, A_BL, ldA);
    // A_BL = A_BR \ A_BL
    BLAS(dtrsm)("L", "L", "N", diag, &n2, &n1, d1, A_BR, ldA, A_BL, ldA);
    // A_BR = 1 / A_BR
    dtrtri_rl(diag, &n2, A_BR, ldA);
}
