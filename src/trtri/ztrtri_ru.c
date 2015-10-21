#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void ztrtri_ru(const char *diag, const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(ztrti2)("U", diag, n, A, ldA, &info);
        return;
    }

    const double z1[] = {1, 0}, zm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
#define A_TL (A)
#define A_TR (A + 2 * *ldA * n1)
#define A_BR (A + 2 * *ldA * n1 + 2 * n1)

    // A_TL = 1 / A_TL
    ztrtri_ru(diag, &n1, A_TL, ldA);
    // A_TR = - A_TL * A_TR
    BLAS(ztrmm)("L", "U", "N", diag, &n1, &n2, zm1, A_TL, ldA, A_TR, ldA);
    // A_TR = A_TR / A_BR
    BLAS(ztrsm)("R", "U", "N", diag, &n1, &n2, z1, A_BR, ldA, A_TR, ldA);
    // A_BR = 1 / A_BR
    ztrtri_ru(diag, &n2, A_BR, ldA);
}
