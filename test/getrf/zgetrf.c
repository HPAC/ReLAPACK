#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "LAPACK_ORIG_getrf.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	double *A1 = malloc(2 * n * n * sizeof(double));
	double *A2 = malloc(2 * n * n * sizeof(double));
    int *ipiv1 = malloc(n * sizeof(int));
    int *ipiv2 = malloc(n * sizeof(int));
#define CLEANUP free(A1); free(A2); free(ipiv1); free(ipiv2);

    srand(time(NULL));

    int i, j, info;

    // m = n
    {
        // generate matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                A1[2 * (i + n * j)] = A2[2 * (i + n * j)] = (double) rand() / RAND_MAX;
                A1[2 * (i + n * j) + 1] = A2[2 * (i + n * j) + 1] = (double) rand() / RAND_MAX;
            }
            A1[2 * (i + n * i)] = A2[2 * (i + n * i)] = (double) rand() / RAND_MAX + n;
            A1[2 * (i + n * i) + 1] = A2[2 * (i + n * i) + 1] = (double) rand() / RAND_MAX;
        }

        // run
        LAPACK(zgetrf)(&n, &n, A1, &n, ipiv1, &info);
        LAPACK_ORIG(zgetrf)(&n, &n, A2, &n, ipiv2, &info);

        // check error
        double error = 0;
        for (i = 0; i < n; i++) { 
            for (j = 0; j < n; j++) {
                error += (A1[2 * (i + n * j)] - A2[2 * (i + n * j)]) * (A1[2 * (i + n * j)] - A2[2 * (i + n * j)]);
                error += (A1[2 * (i + n * j) + 1] - A2[2 * (i + n * j) + 1]) * (A1[2 * (i + n * j) + 1] - A2[2 * (i + n * j) + 1]);
            }
            error += (ipiv1[i] - ipiv2[i]) * (ipiv1[i] - ipiv2[i]);
        }
        error = sqrtf(error);

        printf("zgetrf m = n:\t%g\n", error);
        if (error > TEST_TOL_DOUBLE) {
            CLEANUP
            return 1;
        }
    }

    // m > n
    {
        // generate matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                A1[2 * (i + n * j)] = A2[2 * (i + n * j)] = (double) rand() / RAND_MAX;
                A1[2 * (i + n * j) + 1] = A2[2 * (i + n * j) + 1] = (double) rand() / RAND_MAX;
            }
            A1[2 * (i + n * i)] = A2[2 * (i + n * i)] = (double) rand() / RAND_MAX + n;
            A1[2 * (i + n * i) + 1] = A2[2 * (i + n * i) + 1] = (double) rand() / RAND_MAX;
        }

        int n2 = (n * 3) / 4;

        // run
        LAPACK(zgetrf)(&n, &n2, A1, &n, ipiv1, &info);
        LAPACK_ORIG(zgetrf)(&n, &n2, A2, &n, ipiv2, &info);

        // check error
        double error = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                error += (A1[2 * (i + n * j)] - A2[2 * (i + n * j)]) * (A1[2 * (i + n * j)] - A2[2 * (i + n * j)]);
                error += (A1[2 * (i + n * j) + 1] - A2[2 * (i + n * j) + 1]) * (A1[2 * (i + n * j) + 1] - A2[2 * (i + n * j) + 1]);
            }
        for (i = 0; i < n2; i++)
            error += (ipiv1[i] - ipiv2[i]) * (ipiv1[i] - ipiv2[i]);
        error = sqrtf(error);

        printf("zgetrf m > n:\t%g\n", error);
        if (error > TEST_TOL_DOUBLE) {
            CLEANUP
            return 1;
        }
    }

    // m < n
    {
        // generate matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                A1[2 * (i + n * j)] = A2[2 * (i + n * j)] = (double) rand() / RAND_MAX;
                A1[2 * (i + n * j) + 1] = A2[2 * (i + n * j) + 1] = (double) rand() / RAND_MAX;
            }
            A1[2 * (i + n * i)] = A2[2 * (i + n * i)] = (double) rand() / RAND_MAX + n;
            A1[2 * (i + n * i) + 1] = A2[2 * (i + n * i) + 1] = (double) rand() / RAND_MAX;
        }

        int n2 = (n * 3) / 4;

        // run
        LAPACK(zgetrf)(&n2, &n, A1, &n, ipiv1, &info);
        LAPACK_ORIG(zgetrf)(&n2, &n, A2, &n, ipiv2, &info);

        // check error
        double error = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                error += (A1[2 * (i + n * j)] - A2[2 * (i + n * j)]) * (A1[2 * (i + n * j)] - A2[2 * (i + n * j)]);
                error += (A1[2 * (i + n * j) + 1] - A2[2 * (i + n * j) + 1]) * (A1[2 * (i + n * j) + 1] - A2[2 * (i + n * j) + 1]);
            }
        for (i = 0; i < n2; i++)
            error += (ipiv1[i] - ipiv2[i]) * (ipiv1[i] - ipiv2[i]);
        error = sqrtf(error);

        printf("zgetrf m < n:\t%g\n", error);
        if (error > TEST_TOL_DOUBLE) {
            CLEANUP
            return 1;
        }
    }

    CLEANUP

	return 0;
}
