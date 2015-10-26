#include "larpack.h"

void LARPACK(ssygst)(const int *itype, const char *uplo, const int *n, float *A, const int *ldA, const float *B, const int *ldB, int *info) {
    *info = 0;

    // Check arguments
    int lower = LAPACK(lsame)(uplo, "L");
    int upper = LAPACK(lsame)(uplo, "U");
    if (*itype < 1 || *itype > 3)
        *info = -1;
    else if (!upper && !lower)
        *info = -2;
    else if (*n < 0)
        *info = -3;
    else if (*ldA < MAX(1, *n))
        *info = -5;
    else if (*ldB < MAX(1, *n))
        *info = -7;
    if (*info != 0) {
        int minfo = -*info;
        LAPACK(xerbla)("SSYGST", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) { 
        // Unblocked
        LAPACK(ssygs2)(itype, uplo, n, A, ldA, B, ldB, info);
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

    // B_TL B_TR
    // B_BL B_BR
    const float *const B_TL = B;
    const float *const B_TR = B + *ldB * n1;
    const float *const B_BL = B             + n1;
    const float *const B_BR = B + *ldB * n1 + n1;

    // 1, -1, 1/2, -1/2
   	const float s1[] = {1}, sm1[] = {-1}, sp5[] = {.5}, smp5[] = {-.5};

    // recursion(A_TL, B_TL)
    LARPACK(ssygst)(itype, uplo, &n1, A_TL, ldA, B_TL, ldB, info);

    if (*itype == 1)
        if (lower) {
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
        } else {
            // A_TR = B_TL' \ A_TR
            BLAS(strsm)("L", "U", "T", "N", &n1, &n2, s1, B_TL, ldB, A_TR, ldA);
            // A_TR = A_TR - 1/2 A_TL * B_TR
            BLAS(ssymm)("L", "U", &n1, &n2, smp5, A_TL, ldA, B_TR, ldB, s1, A_TR, ldA);
            // A_BR = A_BR - A_TR' * B_TR - B_TR' * A_TR
            BLAS(ssyr2k)("U", "T", &n2, &n1, sm1, A_TR, ldA, B_TR, ldB, s1, A_BR, ldA);
            // A_TR = A_TR - 1/2 A_TL * B_TR
            BLAS(ssymm)("L", "U", &n1, &n2, smp5, A_TL, ldA, B_TR, ldB, s1, A_TR, ldA);
            // A_TR = A_TR / B_BR
            BLAS(strsm)("R", "U", "N", "N", &n1, &n2, s1, B_BR, ldB, A_TR, ldA);
        }
    else
        if (lower) {
            BLAS(strmm)("R", "L", "N", "N", &n2, &n1, s1, B_TL, ldB, A_BL, ldA);
            // A_BL = A_BL + 1/2 A_BR * B_BL
            BLAS(ssymm)("L", "L", &n2, &n1, sp5, A_BR, ldA, B_BL, ldB, s1, A_BL, ldA);
            // A_TL = A_TL + A_BL' * B_BL + B_BL' * A_BL
            BLAS(ssyr2k)("L", "T", &n1, &n2, s1, A_BL, ldA, B_BL, ldB, s1, A_TL, ldA);
            // A_BL = A_BL + 1/2 A_BR * B_BL
            BLAS(ssymm)("L", "L", &n2, &n1, sp5, A_BR, ldA, B_BL, ldB, s1, A_BL, ldA);
            // A_BL = B_BR * A_BL
            BLAS(strmm)("L", "L", "T", "N", &n2, &n1, s1, B_BR, ldB, A_BL, ldA);
        } else {
            BLAS(strmm)("L", "U", "N", "N", &n1, &n2, s1, B_TL, ldB, A_TR, ldA);
            // A_TR = A_TR + 1/2 B_TR A_BR
            BLAS(ssymm)("R", "U", &n1, &n2, sp5, A_BR, ldA, B_TR, ldB, s1, A_TR, ldA);
            // A_TL = A_TL + A_TR * B_TR' + B_TR * A_TR'
            BLAS(ssyr2k)("U", "N", &n1, &n2, s1, A_TR, ldA, B_TR, ldB, s1, A_TL, ldA);
            // A_TR = A_TR + 1/2 B_TR * A_BR 
            BLAS(ssymm)("R", "U", &n1, &n2, sp5, A_BR, ldA, B_TR, ldB, s1, A_TR, ldA);
            // A_TR = A_TR * B_BR
            BLAS(strmm)("R", "U", "T", "N", &n1, &n2, s1, B_BR, ldB, A_TR, ldA);
        }

    // recursion(A_BR, B_BR)
    LARPACK(ssygst)(itype, uplo, &n2, A_BR, ldA, B_BR, ldB, info);
}
