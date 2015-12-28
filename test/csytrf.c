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
		
	float *A1 = malloc(n * n * 2 * sizeof(float));
	float *A2 = malloc(n * n * 2 * sizeof(float));
	int *ipiv1 = malloc(n * sizeof(int));
	int *ipiv2 = malloc(n * sizeof(int));
    const int lWork = n * n;
	float *Work = malloc(lWork * 2 * sizeof(float));

    int info;

    // L
    {
        // generate matrix
        c2matgen(n, n, A1, A2);

        // run
        RELAPACK(csytrf)("L", &n, A1, &n, ipiv1, Work, &lWork, &info);
        LAPACK(csytf2)("L", &n, A2, &n, ipiv2, &info);

        // check error
        const double error = c2vecerr(n * n, A1, A2) + i2vecerr(n, ipiv1, ipiv2);
        printf("csytrf L:\t%g\n", error);
    }

    // U
    {
        // generate matrix
        c2matgen(n, n, A1, A2);

        // run
        RELAPACK(csytrf)("U", &n, A1, &n, ipiv1, Work, &lWork, &info);
        LAPACK(csytf2)("U", &n, A2, &n, ipiv2, &info);

        // check error
        const double error = c2vecerr(n * n, A1, A2) + i2vecerr(n, ipiv1, ipiv2);
        printf("csytrf U:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(ipiv1);
    free(ipiv2);
    free(Work);

	return 0;
}
