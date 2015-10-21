#include "../../config.h"
#include "../lapack.h"
#include "../util.h"
#include "potrf.h"

void LARPACK(dpotrf)(const char *uplo, const int *n, double *A, const int *ldA, int *info) {
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
        LAPACK(xerbla)("DPOTRF", &minfo);
        return;
    }

    // Quick return if possible
    if (*n == 0)
        return;

    // Call recursive code
    if (lower)
        dpotrf_rl(n, A, ldA, info);
    else
        dpotrf_ru(n, A, ldA, info);
}
