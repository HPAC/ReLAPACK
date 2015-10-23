#include "../../config.h"
#include "../lapack.h"
#include "../util.h"
#include "trtri.h"

void LARPACK(strtri)(const char *uplo, const char *diag, const int *n, float *A, const int *ldA, int *info) {
    *info = 0;

    // Check arguments
    int lower = LAPACK(lsame)(uplo, "L");
    int upper = LAPACK(lsame)(uplo, "U");
    int nounit = LAPACK(lsame)(diag, "N");
    int unit = LAPACK(lsame)(diag, "U");
    if (!upper && !lower)
        *info = -1;
    else if (!nounit && !unit)
        *info = -2;
    else if (*n < 0)
        *info = -3;
    else if (*ldA < MAX(1, *n))
        *info = -5;
    if (*info != 0) {
        int minfo = -*info;
        LAPACK(xerbla)("STRTRI", &minfo);
        return;
    }

    // Quick return if possible
    if (*n == 0)
        return;

    // Check for singularity if non-unit.
    if (nounit) {
        int i;
        for (i = 0; i < *n; i++)
            if (A[i + i * *ldA] == 0) {
                *info = i + 1;
                return;
            }
    }

    // Call recursive code
    strtri_r(uplo, diag, n, A, ldA);
    // if (lower)
    //     strtri_rl(diag, n, A, ldA);
    // else
    //     strtri_ru(diag, n, A, ldA);
}
