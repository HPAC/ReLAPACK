#include "../../src/larpack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
    const int n2 = (n * 3) / 4;
		
	double *A1 = malloc(n * n * sizeof(double));
	double *A2 = malloc(n * n * sizeof(double));
    int *ipiv1 = malloc(n * sizeof(int));
    int *ipiv2 = malloc(n * sizeof(int));

    int info;

    // m = n
    {
        // generate matrix
        d2matgen(n, n, A1, A2);

        // run
        LAPACK(dgetrf)(&n, &n, A1, &n, ipiv1, &info);
        LAPACK(dgetf2)(&n, &n, A2, &n, ipiv2, &info);

        // check error
        double error = d2vecerr(n * n, A1, A2);
        error += i2vecerr(n, ipiv1, ipiv2);
        printf("dgetrf m = n:\t%g\n", error);
    }

    // m > n
    {
        // generate matrix
        d2matgen(n, n2, A1, A2);

        // run
        LAPACK(dgetrf)(&n, &n2, A1, &n, ipiv1, &info);
        LAPACK(dgetf2)(&n, &n2, A2, &n, ipiv2, &info);

        // check error
        double error = d2vecerr(n * n2, A1, A2);
        error += i2vecerr(n2, ipiv1, ipiv2);
        printf("dgetrf m > n:\t%g\n", error);
    }

    // m < n
    {
        // generate matrix
        d2matgen(n2, n, A1, A2);

        // run
        LAPACK(dgetrf)(&n2, &n, A1, &n, ipiv1, &info);
        LAPACK(dgetf2)(&n2, &n, A2, &n, ipiv2, &info);

        // check error
        double error = d2vecerr(n2 * n, A1, A2);
        error += i2vecerr(n2, ipiv1, ipiv2);
        printf("dgetrf m < n:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(ipiv1);
    free(ipiv2);

	return 0;
}
