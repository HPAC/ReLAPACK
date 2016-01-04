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
		
    const int lWork = n * n;
	float *A1    = malloc(n * n * 2 * sizeof(float));
	float *A2    = malloc(n * n * 2 * sizeof(float));
	int   *ipiv1 = malloc(n * sizeof(int));
	int   *ipiv2 = malloc(n * sizeof(int));
	float *Work  = malloc(lWork * 2 * sizeof(float));

    // Output
    int info;

    { // L
        // generate matrix
        c2matgen(n, n, A1, A2);

        // run
        RELAPACK(csytrf_rook)("L", &n, A1, &n, ipiv1, Work, &lWork, &info);
        LAPACK(csytf2_rook)("L", &n, A2, &n, ipiv2, &info);

        // check error
        const double error = c2vecerr(n * n, A1, A2) + i2vecerr(n, ipiv1, ipiv2);
        printf("csytrf_rook L:\t%g\n", error);
    }

    { // U
        // generate matrix
        c2matgen(n, n, A1, A2);

        // run
        RELAPACK(csytrf_rook)("U", &n, A1, &n, ipiv1, Work, &lWork, &info);
        LAPACK(csytf2_rook)("U", &n, A2, &n, ipiv2, &info);

        // check error
        const double error = c2vecerr(n * n, A1, A2) + i2vecerr(n, ipiv1, ipiv2);
        printf("csytrf_rook U:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(ipiv1);
    free(ipiv2);
    free(Work);

	return 0;
}
