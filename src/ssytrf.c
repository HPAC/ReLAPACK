#include "relapack.h"
#include <stdlib.h>

void RELAPACK(ssytrf)(const char *uplo, const int *n,
        float *A, const int *ldA, int *ipiv, 
        float *Work, const int *lWork, int *info) {

    // Check arguments
    const int lower = LAPACK(lsame)(uplo, "L");
    const int upper = LAPACK(lsame)(uplo, "U");
    *info = 0;
    if (!lower && !upper)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*ldA < MAX(1, *n))
        *info = -4;
    else if (*lWork < 1 && *lWork != -1)
        *info = -7;
    if (*info) {
        const int minfo = -*info;
        LAPACK(xerbla)("SSYTRF", &minfo);
        return;
    }

    if (*lWork == -1) {
        // Work size query
        *Work = *ldA * *n;
        return;
    }

    float *W = Work;
    if (*lWork < *ldA * *n)
        W = malloc(*ldA * *n * sizeof(float));

    int nout;
    RELAPACK(slasyf)(uplo, n, n, &nout, A, ldA, ipiv, W, ldA, info);

    if (W != Work)
        free(W);
}
