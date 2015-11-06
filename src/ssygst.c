#include "larpack.h"
#include "stdlib.h"

void LARPACK(ssygst)(const int *itype, const char *uplo, const int *n,
        float *A, const int *ldA, const float *B, const int *ldB, int *info) {

    // Check arguments
    const int lower = LAPACK(lsame)(uplo, "L");
    const int upper = LAPACK(lsame)(uplo, "U");
    *info = 0;
    if (*itype < 1 || *itype > 3)
        *info = -1;
    else if (!lower && !upper)
        *info = -2;
    else if (*n < 0)
        *info = -3;
    else if (*ldA < MAX(1, *n))
        *info = -5;
    else if (*ldB < MAX(1, *n))
        *info = -7;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("SSYGST", &minfo);
        return;
    }

    if (*n <= LARPACK_CROSSOVER) {
        // Unblocked
        LAPACK(ssygs2)(itype, uplo, n, A, ldA, B, ldB, info);
        return;
    }

    // Recursive

    // Constants
    // 0, 1, -1, 1/2, -1/2
   	const float ZERO[] = {0}, ONE[] = {1}, MONE[] = {-1}, HALF[] = {.5}, MHALF[] = {-.5};
    // 1
    const int IONE[] = {1};

    // Splitting
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

    // recursion(A_TL, B_TL)
    LARPACK(ssygst)(itype, uplo, &n1, A_TL, ldA, B_TL, ldB, info);

#ifdef ALLOW_MALLOC
    float *const T = malloc(n2 * n1 * sizeof(float));
    int i;
#endif

    if (*itype == 1)
        if (lower) {
            // A_BL = A_BL / B_TL'
            BLAS(strsm)("R", "L", "T", "N", &n2, &n1, ONE, B_TL, ldB, A_BL, ldA);
#ifdef ALLOW_MALLOC
            // T = -1/2 * B_BL * A_TL
            BLAS(ssymm)("R", "L", &n2, &n1, MHALF, A_TL, ldA, B_BL, ldB, ZERO, T, &n2);
            // A_BL = A_BL + T
            for (i = 0; i < n1; i++)
                BLAS(saxpy)(&n2, ONE, T + n2 * i, IONE, A_BL + *ldA * i, IONE);
#else
            // A_BL = A_BL - 1/2 B_BL * A_TL
            BLAS(ssymm)("R", "L", &n2, &n1, MHALF, A_TL, ldA, B_BL, ldB, ONE, A_BL, ldA);
#endif
            // A_BR = A_BR - A_BL * B_BL' - B_BL * A_BL'
            BLAS(ssyr2k)("L", "N", &n2, &n1, MONE, A_BL, ldA, B_BL, ldB, ONE, A_BR, ldA);
#ifdef ALLOW_MALLOC
            // A_BL = A_BL + T
            for (i = 0; i < n1; i++)
                BLAS(saxpy)(&n2, ONE, T + n2 * i, IONE, A_BL + *ldA * i, IONE);
#else
            // A_BL = A_BL - 1/2 B_BL * A_TL
            BLAS(ssymm)("R", "L", &n2, &n1, MHALF, A_TL, ldA, B_BL, ldB, ONE, A_BL, ldA);
#endif
            // A_BL = B_BR \ A_BL
            BLAS(strsm)("L", "L", "N", "N", &n2, &n1, ONE, B_BR, ldB, A_BL, ldA);
        } else {
            // A_TR = B_TL' \ A_TR
            BLAS(strsm)("L", "U", "T", "N", &n1, &n2, ONE, B_TL, ldB, A_TR, ldA);
#ifdef ALLOW_MALLOC
            // T = -1/2 * A_TL * B_TR
            BLAS(ssymm)("L", "U", &n1, &n2, MHALF, A_TL, ldA, B_TR, ldB, ZERO, T, &n1);
            // A_TR = A_BL + T
            for (i = 0; i < n2; i++)
                BLAS(saxpy)(&n1, ONE, T + n1 * i, IONE, A_TR + *ldA * i, IONE);
#else
            // A_TR = A_TR - 1/2 A_TL * B_TR
            BLAS(ssymm)("L", "U", &n1, &n2, MHALF, A_TL, ldA, B_TR, ldB, ONE, A_TR, ldA);
#endif
            // A_BR = A_BR - A_TR' * B_TR - B_TR' * A_TR
            BLAS(ssyr2k)("U", "T", &n2, &n1, MONE, A_TR, ldA, B_TR, ldB, ONE, A_BR, ldA);
#ifdef ALLOW_MALLOC
            // A_TR = A_BL + T
            for (i = 0; i < n2; i++)
                BLAS(saxpy)(&n1, ONE, T + n1 * i, IONE, A_TR + *ldA * i, IONE);
#else
            // A_TR = A_TR - 1/2 A_TL * B_TR
            BLAS(ssymm)("L", "U", &n1, &n2, MHALF, A_TL, ldA, B_TR, ldB, ONE, A_TR, ldA);
#endif
            // A_TR = A_TR / B_BR
            BLAS(strsm)("R", "U", "N", "N", &n1, &n2, ONE, B_BR, ldB, A_TR, ldA);
        }
    else
        if (lower) {
            // A_BL = A_BL * B_TL
            BLAS(strmm)("R", "L", "N", "N", &n2, &n1, ONE, B_TL, ldB, A_BL, ldA);
#ifdef ALLOW_MALLOC
            // T = 1/2 * A_BR * B_BL
            BLAS(ssymm)("L", "L", &n2, &n1, HALF, A_BR, ldA, B_BL, ldB, ZERO, T, &n2);
            // A_BL = A_BL + T
            for (i = 0; i < n1; i++)
                BLAS(saxpy)(&n2, ONE, T + n2 * i, IONE, A_BL + *ldA * i, IONE);
#else
            // A_BL = A_BL + 1/2 A_BR * B_BL
            BLAS(ssymm)("L", "L", &n2, &n1, HALF, A_BR, ldA, B_BL, ldB, ONE, A_BL, ldA);
#endif
            // A_TL = A_TL + A_BL' * B_BL + B_BL' * A_BL
            BLAS(ssyr2k)("L", "T", &n1, &n2, ONE, A_BL, ldA, B_BL, ldB, ONE, A_TL, ldA);
#ifdef ALLOW_MALLOC
            // A_BL = A_BL + T
            for (i = 0; i < n1; i++)
                BLAS(saxpy)(&n2, ONE, T + n2 * i, IONE, A_BL + *ldA * i, IONE);
#else
            // A_BL = A_BL + 1/2 A_BR * B_BL
            BLAS(ssymm)("L", "L", &n2, &n1, HALF, A_BR, ldA, B_BL, ldB, ONE, A_BL, ldA);
#endif
            // A_BL = B_BR * A_BL
            BLAS(strmm)("L", "L", "T", "N", &n2, &n1, ONE, B_BR, ldB, A_BL, ldA);
        } else {
            // A_TR = B_TL * A_TR
            BLAS(strmm)("L", "U", "N", "N", &n1, &n2, ONE, B_TL, ldB, A_TR, ldA);
#ifdef ALLOW_MALLOC
            // T = 1/2 * B_TR * A_BR
            BLAS(ssymm)("R", "U", &n1, &n2, HALF, A_BR, ldA, B_TR, ldB, ZERO, T, &n1);
            // A_TR = A_TR + T
            for (i = 0; i < n2; i++)
                BLAS(saxpy)(&n1, ONE, T + n1 * i, IONE, A_TR + *ldA * i, IONE);
#else
            // A_TR = A_TR + 1/2 B_TR A_BR
            BLAS(ssymm)("R", "U", &n1, &n2, HALF, A_BR, ldA, B_TR, ldB, ONE, A_TR, ldA);
#endif
            // A_TL = A_TL + A_TR * B_TR' + B_TR * A_TR'
            BLAS(ssyr2k)("U", "N", &n1, &n2, ONE, A_TR, ldA, B_TR, ldB, ONE, A_TL, ldA);
#ifdef ALLOW_MALLOC
            // A_TR = A_TR + T
            for (i = 0; i < n2; i++)
                BLAS(saxpy)(&n1, ONE, T + n1 * i, IONE, A_TR + *ldA * i, IONE);
#else
            // A_TR = A_TR + 1/2 B_TR * A_BR
            BLAS(ssymm)("R", "U", &n1, &n2, HALF, A_BR, ldA, B_TR, ldB, ONE, A_TR, ldA);
#endif
            // A_TR = A_TR * B_BR
            BLAS(strmm)("R", "U", "T", "N", &n1, &n2, ONE, B_BR, ldB, A_TR, ldA);
        }

#ifdef ALLOW_MALLOC
    free(T);
#endif

    // recursion(A_BR, B_BR)
    LARPACK(ssygst)(itype, uplo, &n2, A_BR, ldA, B_BR, ldB, info);
}
