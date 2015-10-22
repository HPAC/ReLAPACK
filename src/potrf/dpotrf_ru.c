#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void dpotrf_ru(const int *n, double *A, const int *ldA, int *info) {
    if (*n <= LARPACK_CROSSOVER) { 
        LAPACK(dpotf2)("U", n, A, ldA, info);
        return;
    }

   	const double d1[] = {1}, dm1[] = {-1};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    //      A_BR
    double *A_TL = A;
    double *A_TR = A + *ldA * n1;
    double *A_BR = A + *ldA * n1 + n1;

    // A_TL = Chol(A_TL)
    dpotrf_ru(&n1, A_TL, ldA, info);
    if (*info)
        return;
    // A_TR = A_TL' \ A_TR
    BLAS(dtrsm)("L", "U", "T", "N", &n1, &n2, d1, A_TL, ldA, A_TR, ldA);
    // A_BR = A_BR - A_TR' * A_TR
    BLAS(dsyrk)("U", "T", &n2, &n1, dm1, A_TR, ldA, d1, A_BR, ldA);
    // A_BR = Chol(A_BR)
    dpotrf_ru(&n2, A_BR, ldA, info);
    if (*info) {
        *info += n1;
        return;
    }
}
