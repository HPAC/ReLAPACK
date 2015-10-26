#include "../../src/larpack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
    // const int n2 = (n * 3) / 4;
		
	float *A1 = malloc(n * n * sizeof(float));
	float *A2 = malloc(n * n * sizeof(float));
	float *B1 = malloc(n * n * sizeof(float));
	float *B2 = malloc(n * n * sizeof(float));
	float *C1 = malloc(n * n * sizeof(float));
	float *C2 = malloc(n * n * sizeof(float));

    int info;
    float scale;
    int i1 = 1;

    // N N +1 m = n
    {
        // generate matrix
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(n, n, C1, C2);

        // run
        LARPACK(strsyl)("N", "N", &i1, &n, &n, A1, &n, B1, &n, C1, &n, &scale, &info);
        LAPACK(strsy2)("N", "N", &i1, &n, &n, A2, &n, B2, &n, C2, &n, &scale, &info);
        printf("%g\n", scale);

        int i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                printf(" %g", C1[i + j * n]);
            printf("\n");
        }
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                printf(" %g", C2[i + j * n]);
            printf("\n");
        }


        // check error
        float error = s2vecerr(n * n, C1, C2);
        printf("strsyl N N +1 m = n:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);

	return 0;
}
