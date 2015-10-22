#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	double *A1 = malloc(2 * n * n * sizeof(double));
	double *A2 = malloc(2 * n * n * sizeof(double));

    int info;

    // Lower
    {
        // generate matrix
        z2matgen(n, n, A1, A2);

        // run
        LAPACK(ztrtri)("L", "N", &n, A1, &n, &info);
        LAPACK(ztrti2)("L", "N", &n, A2, &n, &info);

        // check error
        double error = z2vecerr(n * n, A1, A2);
        printf("ztrtri Lower:\t%g\n", error);
    }

    // Upper
    {
        // generate matrix
        z2matgen(n, n, A1, A2);

        // run
        LAPACK(ztrtri)("U", "N", &n, A1, &n, &info);
        LAPACK(ztrti2)("U", "N", &n, A2, &n, &info);

        // check error
        double error = z2vecerr(n * n, A1, A2);
        printf("ztrtri Upper:\t%g\n", error);
    }

    free(A1); 
    free(A2);

	return 0;
}
