#include "../../src/larpack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	float *A1 = malloc(2 * n * n * sizeof(float));
	float *A2 = malloc(2 * n * n * sizeof(float));

    int info;

    // L
    {
        // generate matrix
        c2matgen(n, n, A1, A2);

        // run
        LARPACK(cpotrf)("L", &n, A1, &n, &info);
        LAPACK(cpotf2)("L", &n, A2, &n, &info);

        // check error
        float error = c2vecerr(n * n, A1, A2);
        printf("cpotrf L:\t%g\n", error);
    }

    // U
    {
        // generate matrix
        c2matgen(n, n, A1, A2);

        // run
        LARPACK(cpotrf)("U", &n, A1, &n, &info);
        LAPACK(cpotf2)("U", &n, A2, &n, &info);

        // check error
        float error = c2vecerr(n * n, A1, A2);
        printf("cpotrf U:\t%g\n", error);
    }

    free(A1);
    free(A2);

	return 0;
}