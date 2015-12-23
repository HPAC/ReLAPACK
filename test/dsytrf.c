#include "../src/relapack.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

    if (argc == 1) {
        fprintf(stderr, "usage: %s n\n", argv[0]);
        return 0;
    }
    const int n = atoi(argv[1]);
		
	double *A1 = malloc(n * n * sizeof(double));
	double *A2 = malloc(n * n * sizeof(double));
	int *ipiv1 = malloc(n * sizeof(int));
	int *ipiv2 = malloc(n * sizeof(int));
    const int lWork = n * n;
	double *Work = malloc(lWork * sizeof(double));

    int info;

    // L
    {
        // generate matrix
        d2matgen(n, n, A1, A2);

        // run
        RELAPACK(dsytrf)("L", &n, A1, &n, ipiv1, Work, &lWork, &info);
        LAPACK(dsytf2)("L", &n, A2, &n, ipiv2, &info);

        // check error
        const double error = d2vecerr(n * n, A1, A2) + i2vecerr(n, ipiv1, ipiv2);
        printf("ssytrf L:\t%g\n", error);
    }

    // U
    {
        // generate matrix
        d2matgen(n, n, A1, A2);

        // run
        RELAPACK(dsytrf)("U", &n, A1, &n, ipiv1, Work, &lWork, &info);
        LAPACK(dsytf2)("U", &n, A2, &n, ipiv2, &info);

        // check error
        const double error = d2vecerr(n * n, A1, A2) + i2vecerr(n, ipiv1, ipiv2);
        printf("ssytrf U:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(ipiv1);
    free(ipiv2);
    free(Work);

	return 0;
}