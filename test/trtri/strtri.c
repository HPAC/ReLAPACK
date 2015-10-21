#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "LAPACK_ORIG_trtri.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	float *A1 = malloc(n * n * sizeof(float));
	float *A2 = malloc(n * n * sizeof(float));
#define CLEANUP free(A1); free(A2);

    srand(time(NULL));

    int i, j, info;

    // Lower
    {
        // generate matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                A1[i + n * j] = A2[i + n * j] = (float) rand() / RAND_MAX;
            A1[i + n * i] = A2[i + n * i] = (float) rand() / RAND_MAX + n;
        }

        // run
        LAPACK(strtri)("L", "N", &n, A1, &n, &info);
        LAPACK_ORIG(strtri)("L", "N", &n, A2, &n, &info);

        // check error
        float error = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                error += (A1[i + n * j] - A2[i + n * j]) * (A1[i + n * j] - A2[i + n * j]);
        error = sqrtf(error);

        printf("strtri Lower:\t%g\n", error);
        if (error > TEST_TOL_FLOAT) {
            CLEANUP
            return 1;
        }
    }

    // Upper
    {
        // generate matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                A1[i + n * j] = A2[i + n * j] = (float) rand() / RAND_MAX;
            A1[i + n * i] = A2[i + n * i] = (float) rand() / RAND_MAX + n;
        }

        // run
        LAPACK(strtri)("U", "N", &n, A1, &n, &info);
        LAPACK_ORIG(strtri)("U", "N", &n, A2, &n, &info);

        // check error
        float error = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                error += (A1[i + n * j] - A2[i + n * j]) * (A1[i + n * j] - A2[i + n * j]);
        error = sqrtf(error);

        printf("strtri Upper:\t%g\n", error);
        if (error > TEST_TOL_FLOAT) {
            CLEANUP
            return 1;
        }
    }

    CLEANUP

	return 0;
}
