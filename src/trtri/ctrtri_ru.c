#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void ctrtri_ru(const char *diag, const int *n, float *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(ctrti2)("U", diag, n, A, ldA, &info);
        return;
    }

   	const float c1[] = {1, 0}, cm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
    float *const A_TL = A;
    float *const A_TR = A + 2 * *ldA * n1;
    float *const A_BR = A + 2 * *ldA * n1 + 2 * n1;

    // A_TL = 1 / A_TL
    ctrtri_ru(diag, &n1, A_TL, ldA);
    // A_TR = - A_TL * A_TR
    BLAS(ctrmm)("L", "U", "N", diag, &n1, &n2, cm1, A_TL, ldA, A_TR, ldA);
    // A_TR = A_TR / A_BR
    BLAS(ctrsm)("R", "U", "N", diag, &n1, &n2, c1, A_BR, ldA, A_TR, ldA);
    // A_BR = 1 / A_BR
    ctrtri_ru(diag, &n2, A_BR, ldA);
}
