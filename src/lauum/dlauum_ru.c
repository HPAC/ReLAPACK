#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void dlauum_ru(const int *n, double *A, const int *ldA) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        LAPACK(dlauu2)("U", n, A, ldA, &info);
        return;
    }

   	const double d1[] = {1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
    double *A_TL = A;
    double *A_TR = A + *ldA * n1;
    double *A_BR = A + *ldA * n1 + n1;

    // A_TL = A_TL' * A_TL
    dlauum_ru(&n1, A_TL, ldA);
    // A_TL = A_TL + A_TR * A_TR'
    BLAS(dsyrk)("U", "N", &n1, &n2, d1, A_TR, ldA, d1, A_TL, ldA);
    // A_TR = A_TR * A_BR'
    BLAS(dtrmm)("R", "U", "T", "N", &n1, &n2, d1, A_BR, ldA, A_TR, ldA);
    // A_BR = A_BR' * A_BR
    dlauum_ru(&n2, A_BR, ldA);
}
