#include "../../config.h"
#include "../lapack.h"
#include "../util.h"
#include "lauum.h"

void LARPACK(dlauum)(const char *uplo, const int *n, double *A, const int *ldA, int *info) {
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
        LAPACK(xerbla)("DLAUUM", &minfo);
        return;
    }

    // Quick return if possible
    if (*n == 0)
        return;

    // Call recursive code
    if (lower)
        dlauum_rl(n, A, ldA);
    else
        dlauum_ru(n, A, ldA);
}
