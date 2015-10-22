#include "../../config.h"
#include "../lapack.h"
#include "../util.h"
#include "trsyl.h"

void LARPACK(strsyl)(const char *transA, const char *transB, const int *isgn, const int *m, const int *n, const float *A, const int *ldA, const float *B, const int *ldB, float *C, const int *ldC, float *scale, int *info) {
    *info = 0;

    // Check arguments
    int notranA = LAPACK(lsame)(transA, "N");
    int tranA = LAPACK(lsame)(transA, "T") || LAPACK(lsame)(transA, "C");
    int notranB = LAPACK(lsame)(transB, "N");
    int tranB = LAPACK(lsame)(transB, "T") || LAPACK(lsame)(transB, "C");
    if (!tranA && !notranA)
        *info = -1;
    else if (!tranB && !notranB)
        *info = -2;
    else if (*isgn != 1 && *isgn != -1)
        *info = -3;
    else if (*m < 0)
        *info = -4;
    else if (*n < 0)
        *info = -5;
    else if (*ldA < MAX(1, *m))
        *info = -7;
    else if (*ldB < MAX(1, *n))
        *info = -9;
    else if (*ldC < MAX(1, *m))
        *info = -11;
    if (*info != 0) {
        int minfo = -*info;
        LAPACK(xerbla)("STRSYL", &minfo);
        return;
    }

    // Quick return if possible
    *scale = 1;
    if (*m == 0 && *n == 0)
        return;

    // Call recursive code
    if (tranA)
        if (tranB)
            strsyl_rtt(isgn, m, n, A, ldA, B, ldB, C, ldC, scale);
        else
            strsyl_rtn(isgn, m, n, A, ldA, B, ldB, C, ldC, scale);
    else
        if (tranB)
            strsyl_rnt(isgn, m, n, A, ldA, B, ldB, C, ldC, scale);
        else
            strsyl_rnn(isgn, m, n, A, ldA, B, ldB, C, ldC, scale);
}
