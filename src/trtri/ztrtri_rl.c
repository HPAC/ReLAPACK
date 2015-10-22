#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void ztrtri_rl(const char *diag, const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(ztrti2)("L", diag, n, A, ldA, &info);
        return;
    }

    const double z1[] = {1, 0}, zm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
    double *A_TL = A;
    double *A_BL = A                 + 2 * n1;
    double *A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // A_TL = 1 / A_TL
    ztrtri_rl(diag, &n1, A_TL, ldA);
    // A_BL = - A_BL * A_TL
    BLAS(ztrmm)("R", "L", "N", diag, &n2, &n1, zm1, A_TL, ldA, A_BL, ldA);
    // A_BL = A_BR \ A_BL
    BLAS(ztrsm)("L", "L", "N", diag, &n2, &n1, z1, A_BR, ldA, A_BL, ldA);
    // A_BR = 1 / A_BR
    ztrtri_rl(diag, &n2, A_BR, ldA);
}
