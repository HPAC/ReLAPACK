#include "../../config.h"
#include "../lapack.h"
#include "../util.h"
#include "trtri.h"

void LARPACK(ztrtri)(const char *uplo, const char *diag, const int *n, double *A, const int *ldA, int *info) {
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
        LAPACK(xerbla)("ZTRTRI", &minfo);
        return;
    }

    // Quick return if possible
    if (*n == 0)
        return;

    // Check for singularity if non-unit.
    if (nounit) {
        int i;
        for (i = 0; i < *n; i++)
            if (A[2 * (i + i * *ldA)] == 0 && A[2 * (i + i * *ldA) + 1] == 0) {
                *info = i + 1;
                return;
            }
    }

    // Call recursive code
    ztrtri_r(uplo, diag, n, A, ldA);
}
