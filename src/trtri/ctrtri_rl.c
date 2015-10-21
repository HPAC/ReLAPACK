#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void ctrtri_rl(const char *diag, const int *n, float *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(ctrti2)("L", diag, n, A, ldA, &info);
        return;
    }

    const float c1[] = {1, 0}, cm1[] = {-1, 0};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
#define A_TL (A)
#define A_BL (A                 + 2 * n1)
#define A_BR (A + 2 * *ldA * n1 + 2 * n1)

    // A_TL = 1 / A_TL
    ctrtri_rl(diag, &n1, A_TL, ldA);
    // A_BL = - A_BL * A_TL
    BLAS(ctrmm)("R", "L", "N", diag, &n2, &n1, cm1, A_TL, ldA, A_BL, ldA);
    // A_BL = A_BR \ A_BL
    BLAS(ctrsm)("L", "L", "N", diag, &n2, &n1, c1, A_BR, ldA, A_BL, ldA);
    // A_BR = 1 / A_BR
    ctrtri_rl(diag, &n2, A_BR, ldA);
}
