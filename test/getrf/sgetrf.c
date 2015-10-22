#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "LAPACK_ORIG_getrf.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
    const int n2 = (n * 3) / 4;
		
	float *A1 = malloc(n * n * sizeof(float));
	float *A2 = malloc(n * n * sizeof(float));
    int *ipiv1 = malloc(n * sizeof(int));
    int *ipiv2 = malloc(n * sizeof(int));

    int info;

    // m = n
    {
        // generate matrix
        s2matgen(n, n, A1, A2);

        // run
        LAPACK(sgetrf)(&n, &n, A1, &n, ipiv1, &info);
        LAPACK_ORIG(sgetrf)(&n, &n, A2, &n, ipiv2, &info);

        // check error
        float error = s2vecerr(n * n, A1, A2);
        error += i2vecerr(n, ipiv1, ipiv2);
        printf("sgetrf m = n:\t%g\n", error);
    }

    // m > n
    {
        // generate matrix
        s2matgen(n2, n, A1, A2);

        // run
        LAPACK(sgetrf)(&n, &n2, A1, &n, ipiv1, &info);
        LAPACK_ORIG(sgetrf)(&n, &n2, A2, &n, ipiv2, &info);

        // check error
        float error = s2vecerr(n * n2, A1, A2);
        error += i2vecerr(n2, ipiv1, ipiv2);
        printf("sgetrf m > n:\t%g\n", error);
    }

    // m < n
    {
        // generate matrix
        s2matgen(n2, n, A1, A2);

        // run
        LAPACK(sgetrf)(&n2, &n, A1, &n, ipiv1, &info);
        LAPACK_ORIG(sgetrf)(&n2, &n, A2, &n, ipiv2, &info);

        // check error
        float error = s2vecerr(n2 * n, A1, A2);
        error += i2vecerr(n2, ipiv1, ipiv2);
        printf("sgetrf m < n:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(ipiv1);
    free(ipiv2);

	return 0;
}
