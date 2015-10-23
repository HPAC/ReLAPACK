#include "larpack.h"

void LARPACK(spotrf)(const char *uplo, const int *n, float *A, const int *ldA, int *info) {
    *info = 0;

    // Check arguments
    int lower = LAPACK(lsame)(uplo, "L");
    int upper = LAPACK(lsame)(uplo, "U");
    if (!upper && !lower)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*ldA < MAX(1, *n))
        *info = -4;
    if (*info != 0) {
        int minfo = -*info;
        LAPACK(xerbla)("SPOTRF", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) { 
        // Unblocked
        LAPACK(spotf2)(uplo, n, A, ldA, info);
        return;
    }

    // Recursive
    const int n1 = (*n >= 16) ? ((*n + 8) / 16) * 8 : *n / 2;
    const int n2 = *n - n1;

    // A_TL A_TR
    // A_BL A_BR
    float *const A_TL = A;
    float *const A_TR = A + *ldA * n1;
    float *const A_BL = A             + n1;
    float *const A_BR = A + *ldA * n1 + n1;

    // 1, -1
   	const float s1[] = {1}, sm1[] = {-1};

    // recursion(A_TL)
    LARPACK(spotrf)(uplo, &n1, A_TL, ldA, info);
    if (*info)
        return;

    if (lower) {
        // A_BL = A_BL / A_TL'
        BLAS(strsm)("R", "L", "T", "N", &n2, &n1, s1, A_TL, ldA, A_BL, ldA);
        // A_BR = A_BR - A_BL * A_BL'
        BLAS(ssyrk)("L", "N", &n2, &n1, sm1, A_BL, ldA, s1, A_BR, ldA);
    } else {
        // A_TR = A_TL' \ A_TR
        BLAS(strsm)("L", "U", "T", "N", &n1, &n2, s1, A_TL, ldA, A_TR, ldA);
        // A_BR = A_BR - A_TR' * A_TR
        BLAS(ssyrk)("U", "T", &n2, &n1, sm1, A_TR, ldA, s1, A_BR, ldA);
    }
    
    // recursion(A_BR)
    LARPACK(spotrf)(uplo, &n2, A_BR, ldA, info);
    if (*info)
        *info += n1;
}
