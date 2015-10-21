#include "../../config.h"
#include "../lapack.h"
#include "../util.h"
#include "getrf.h"

void LARPACK(sgetrf)(const int *m, const int *n, float *A, const int *ldA, int *ipiv, int *info) {
    *info = 0;

    // Check arguments
    if (*m < 0)
        *info = -1;
    else if (*n < 0)
        *info = -2;
    else if (*ldA < MAX(1, *n))
        *info = -4;
    if (*info != 0) {
        int minfo = -*info;
        LAPACK(xerbla)("SGETRF", &minfo);
        return;
    }

    // Quick return if possible
    if (*m == 0 || *n == 0)
        return;

    // Call recursive code
    sgetrf_r(m, n, A, ldA, ipiv, info);
}
