#include "test.h"

datatype *A[2];
int info[2];
int ld;

void pre() {
    x2matgen(n + 1, n, A[0], A[1]);
    for (int i = 0; i < n; i++) {
        // set diagonal
        A[0][x1 * (i + ld * i)] =
        A[1][x1 * (i + ld * i)] = (datatype) rand() / RAND_MAX;
        // set first row
        A[0][x1 * (ld * i)] =
        A[1][x1 * (ld * i)] = (datatype) rand() / RAND_MAX + n;
        // set last row
        A[0][x1 * (n + ld * i)] =
        A[1][x1 * (n + ld * i)] = (datatype) rand() / RAND_MAX + n;
    }
}

void post() {
    error = x2vecerr(ld * n, A[0], A[1]);
    printf("info = %d %d\n", info[0], info[1]);
}

void tests() {
    ld = n + 1;
    A[0] = xmalloc(ld * n);
    A[1] = xmalloc(ld * n);

    #define ROUTINE XPREF(pbtrf)

    const int kd1 = n / 4, kd2 = n * 3 / 4;
    TEST("L", &n, &kd1, A[i], &ld, &info[i]);
    TEST("L", &n, &kd2, A[i], &ld, &info[i]);
    TEST("U", &n, &kd1, A[i], &ld, &info[i]);
    TEST("U", &n, &kd2, A[i], &ld, &info[i]);

    free(A[0]);
    free(A[1]);
}
