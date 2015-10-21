#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void strtri_ru(const char *diag, const int *n, float *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(strti2)("U", diag, n, A, ldA, &info);
        return;
    }

   	const float s1[] = {1}, sm1[] = {-1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
#define A_TL (A)
#define A_TR (A + *ldA * n1)
#define A_BR (A + *ldA * n1 + n1)

    // A_TL = 1 / A_TL
    strtri_ru(diag, &n1, A_TL, ldA);
    // A_TR = - A_TL * A_TR
    BLAS(strmm)("L", "U", "N", diag, &n1, &n2, sm1, A_TL, ldA, A_TR, ldA);
    // A_TR = A_TR / A_BR
    BLAS(strsm)("R", "U", "N", diag, &n1, &n2, s1, A_BR, ldA, A_TR, ldA);
    // A_BR = 1 / A_BR
    strtri_ru(diag, &n2, A_BR, ldA);
}
