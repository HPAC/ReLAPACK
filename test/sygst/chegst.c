#include "../../src/larpack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	float *A1 = malloc(2 * n * n * sizeof(float));
	float *A2 = malloc(2 * n * n * sizeof(float));
	float *B1 = malloc(2 * n * n * sizeof(float));
	float *B2 = malloc(2 * n * n * sizeof(float));

    int info;

    int i1[] = {1}, i2[] = {2};

    { // 1 L
        // generate matrix
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);

        // run
        LARPACK(chegst)(i1, "L", &n, A1, &n, B1, &n, &info);
        LAPACK(chegs2)(i1, "L", &n, A2, &n, B2, &n, &info);

        // check error
        float error = c2vecerr(n * n, A1, A2);
        printf("chegst 1 L:\t%g\n", error);
    }

    { // 1 U
        // generate matrix
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);

        // run
        LARPACK(chegst)(i1, "U", &n, A1, &n, B1, &n, &info);
        LAPACK(chegs2)(i1, "U", &n, A2, &n, B2, &n, &info);

        // check error
        float error = c2vecerr(n * n, A1, A2);
        printf("chegst 1 U:\t%g\n", error);
    }

    { // 2 L
        // generate matrix
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);

        // run
        LARPACK(chegst)(i2, "L", &n, A1, &n, B1, &n, &info);
        LAPACK(chegs2)(i2, "L", &n, A2, &n, B2, &n, &info);

        // check error
        float error = c2vecerr(n * n, A1, A2);
        printf("chegst 2 L:\t%g\n", error);
    }

    { // 2 U
        // generate matrix
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);

        // run
        LARPACK(chegst)(i2, "U", &n, A1, &n, B1, &n, &info);
        LAPACK(chegs2)(i2, "U", &n, A2, &n, B2, &n, &info);

        // check error
        float error = c2vecerr(n * n, A1, A2);
        printf("chegst 2 U:\t%g\n", error);
    }

    free(A1); 
    free(A2);
    free(B1); 
    free(B2);

	return 0;
}
