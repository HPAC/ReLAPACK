#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void ztrtri_r(const char *uplo, const char *diag, const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(ztrti2)(uplo, diag, n, A, ldA, &info);
        return;
    }

    const double z1[] = {1, 0}, zm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + 2 * *ldA * n1;
    double *const A_BL = A                 + 2 * n1;
    double *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // A_TL = 1 / A_TL
    ztrtri_r(uplo, diag, &n1, A_TL, ldA);

    if (uplo[0] == 'L') {
        // A_BL = - A_BL * A_TL
        BLAS(ztrmm)("R", "L", "N", diag, &n2, &n1, zm1, A_TL, ldA, A_BL, ldA);
        // A_BL = A_BR \ A_BL
        BLAS(ztrsm)("L", "L", "N", diag, &n2, &n1, z1, A_BR, ldA, A_BL, ldA);
    } else {
        // A_TR = - A_TL * A_TR
        BLAS(ztrmm)("L", "U", "N", diag, &n1, &n2, zm1, A_TL, ldA, A_TR, ldA);
        // A_TR = A_TR / A_BR
        BLAS(ztrsm)("R", "U", "N", diag, &n1, &n2, z1, A_BR, ldA, A_TR, ldA);
    }

    // A_BR = 1 / A_BR
    ztrtri_r(uplo, diag, &n2, A_BR, ldA);
}
