#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "LAPACK_ORIG_potrf.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	double *A1 = malloc(n * n * sizeof(double));
	double *A2 = malloc(n * n * sizeof(double));
#define CLEANUP free(A1); free(A2);

    srand(time(NULL));

    int i, j, info;

    // Lower
    {
        // generate matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                A1[i + n * j] = A2[i + n * j] = (double) rand() / RAND_MAX;
            A1[i + n * i] = A2[i + n * i] = (double) rand() / RAND_MAX + n * n;
        }

        // run
        LAPACK(dpotrf)("L", &n, A1, &n, &info);
        LAPACK_ORIG(dpotrf)("L", &n, A2, &n, &info);

        // check error
        double error = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                error += (A1[i + n * j] - A2[i + n * j]) * (A1[i + n * j] - A2[i + n * j]);
        error = sqrtf(error);

        printf("dpotrf Lower:\t%g\n", error);
        if (error > TEST_TOL_DOUBLE) {
            CLEANUP
            return 1;
        }
    }

    // Upper
    {
        // generate matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                A1[i + n * j] = A2[i + n * j] = (double) rand() / RAND_MAX;
            A1[i + n * i] = A2[i + n * i] = (double) rand() / RAND_MAX + n * n;
        }

        // run
        LAPACK(dpotrf)("U", &n, A1, &n, &info);
        LAPACK_ORIG(dpotrf)("U", &n, A2, &n, &info);

        // check error
        double error = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                error += (A1[i + n * j] - A2[i + n * j]) * (A1[i + n * j] - A2[i + n * j]);
        error = sqrt(error);

        printf("dpotrf Upper:\t%g\n", error);
        if (error > TEST_TOL_DOUBLE) {
            CLEANUP
            return 1;
        }
    }

    CLEANUP

	return 0;
}
