#include "../../config.h"
#include "../blas.h"
#include "../lapack.h"

void ssygst_ril(const int *n, float *A, const int *ldA, const float *B, const int *ldB) {
    if (*n <= LARPACK_CROSSOVER) { 
        int info;
        int i1 = 1;
        LAPACK(ssygs2)(&i1, "L", n, A, ldA, B, ldB, &info);
        return;
    }

   	const float s1[] = {1}, sm1[] = {-1}, smp5[] = {-.5};

    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_BL = A             + n1;
    float *const A_BR = A + *ldA * n1 + n1;

    // B_TL
    // B_BL B_BR
    const float *const B_TL = B;
    const float *const B_BL = B             + n1;
    const float *const B_BR = B + *ldB * n1 + n1;

    // A_TL = B_TL \ A_TL / B_TL';
    ssygst_ril(&n1, A_TL, ldA, B_TL, ldB);
    // A_BL = A_BL / B_TL'
    BLAS(strsm)("R", "L", "T", "N", &n2, &n1, s1, B_TL, ldB, A_BL, ldA);
    // A_BL = A_BL - 1/2 B_BL * A_TL
    BLAS(ssymm)("R", "L", &n2, &n1, smp5, A_TL, ldA, B_BL, ldB, s1, A_BL, ldA);
    // A_BR = A_BR - A_BL * B_BL' - B_BL * A_BL'
    BLAS(ssyr2k)("L", "N", &n2, &n1, sm1, A_BL, ldA, B_BL, ldB, s1, A_BR, ldA);
    // A_BL = A_BL - 1/2 B_BL * A_TL
    BLAS(ssymm)("R", "L", &n2, &n1, smp5, A_TL, ldA, B_BL, ldB, s1, A_BL, ldA);
    // A_BL = B_BR \ A_BL
    BLAS(strsm)("L", "L", "N", "N", &n2, &n1, s1, B_BR, ldB, A_BL, ldA);
    // A_BR = B_BR \ A_BR / B_BR';
    ssygst_ril(&n2, A_BR, ldA, B_BR, ldB);
}
