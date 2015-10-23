#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void dtrtri_r(const char *uplo, const char *diag, const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(dtrti2)(uplo, diag, n, A, ldA, &info);
        return;
    }

   	const double d1[] = {1}, dm1[] = {-1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    double *const A_TL = A;
    double *const A_TR = A + *ldA * n1;
    double *const A_BL = A             + n1;
    double *const A_BR = A + *ldA * n1 + n1;

    // A_TL = 1 / A_TL
    dtrtri_r(uplo, diag, &n1, A_TL, ldA);
    
    if (uplo[0] == 'L') {
        // A_BL = - A_BL * A_TL
        BLAS(dtrmm)("R", "L", "N", diag, &n2, &n1, dm1, A_TL, ldA, A_BL, ldA);
        // A_BL = A_BR \ A_BL
        BLAS(dtrsm)("L", "L", "N", diag, &n2, &n1, d1, A_BR, ldA, A_BL, ldA);
    } else {
        // A_TR = - A_TL * A_TR
        BLAS(dtrmm)("L", "U", "N", diag, &n1, &n2, dm1, A_TL, ldA, A_TR, ldA);
        // A_TR = A_TR / A_BR
        BLAS(dtrsm)("R", "U", "N", diag, &n1, &n2, d1, A_BR, ldA, A_TR, ldA);
    }

    // A_BR = 1 / A_BR
    dtrtri_r(uplo, diag, &n2, A_BR, ldA);
}
