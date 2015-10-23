#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void strtri_r(const char *uplo, const char *diag, const int *n, float *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(strti2)(uplo, diag, n, A, ldA, &info);
        return;
    }

   	const float s1[] = {1}, sm1[] = {-1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + *ldA * n1;
    float *const A_BL = A             + n1;
    float *const A_BR = A + *ldA * n1 + n1;

    // A_TL = 1 / A_TL
    strtri_r(uplo, diag, &n1, A_TL, ldA);

    if (uplo[0] == 'L') {
        // A_BL = - A_BL * A_TL
        BLAS(strmm)("R", "L", "N", diag, &n2, &n1, sm1, A_TL, ldA, A_BL, ldA);
        // A_BL = A_BR \ A_BL
        BLAS(strsm)("L", "L", "N", diag, &n2, &n1, s1, A_BR, ldA, A_BL, ldA);
    } else {
        // A_TR = - A_TL * A_TR
        BLAS(strmm)("L", "U", "N", diag, &n1, &n2, sm1, A_TL, ldA, A_TR, ldA);
        // A_TR = A_TR / A_BR
        BLAS(strsm)("R", "U", "N", diag, &n1, &n2, s1, A_BR, ldA, A_TR, ldA);
    }

    // A_BR = 1 / A_BR
    strtri_r(uplo, diag, &n2, A_BR, ldA);
}
