#include "../../config.h"
#include "../lapack.h"
#include "../util.h"
#include "sygst.h"

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

    // Quick return if possible
    if (*n == 0)
        return;

    // Call recursive code
    if (*itype == 1)
        if (lower)
            ssygst_ril(n, A, ldA, B, ldB);
//        else
//            ssygst_riu(n, A, ldA, B, ldB);
//    else
//        if (lower)
//            ssygst_rnl(n, A, ldA, B, ldB);
//        else
//            ssygst_rnu(n, A, ldA, B, ldB);
}
