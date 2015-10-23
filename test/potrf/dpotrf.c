#include "../../src/larpack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	double *A1 = malloc(n * n * sizeof(double));
	double *A2 = malloc(n * n * sizeof(double));

    int info;

    // L
    {
        // generate matrix
        d2matgen(n, n, A1, A2);

        // run
        LAPACK(dpotrf)("L", &n, A1, &n, &info);
        LAPACK(dpotf2)("L", &n, A2, &n, &info);

        // check error
        double error = d2vecerr(n * n, A1, A2);
        printf("dpotrf L:\t%g\n", error);
    }

    // U
    {
        // generate matrix
        d2matgen(n, n, A1, A2);

        // run
        LAPACK(dpotrf)("U", &n, A1, &n, &info);
        LAPACK(dpotf2)("U", &n, A2, &n, &info);

        // check error
        double error = d2vecerr(n * n, A1, A2);
        printf("dpotrf U:\t%g\n", error);
    }

    free(A1);
    free(A2);

	return 0;
}
