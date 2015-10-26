#include "../../src/larpack.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

    int info;
    float scale;
    int i1 = 1;

    float A[] = {1, 1000, 2, 3};
    float B[] = {.5, 1000, 1.5, 2.5};
    float C1[] = {4, 3, 2, 1};
    float C2[] = {4, 3, 2, 1};

    int m = 1;
    int n = 2;

    // run
    LARPACK(strsyl)("N", "N", &i1, &m, &n, A, &m, B, &n, C1, &m, &scale, &info);
    LAPACK(strsy2)("N", "N", &i1, &m, &n, A, &m, B, &n, C2, &m, &scale, &info);
    printf("%g\n", scale);

    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf(" %12g", C1[i + j * m]);
        printf("\n");
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf(" %12g", C2[i + j * m]);
        printf("\n");
    }

	return 0;
}
