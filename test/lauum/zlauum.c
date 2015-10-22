#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	double *A1 = malloc(n * n * sizeof(double));
	double *A2 = malloc(n * n * sizeof(double));

    int info;

    // Lower
    {
        // generate matrix
        z2matgen(n, n, A1, A2);

        // run
        LAPACK(zlauum)("L", &n, A1, &n, &info);
        LAPACK(zlauu2)("L", &n, A2, &n, &info);

        // check error
        double error = z2vecerr(n * n, A1, A2);
        printf("zlauum Lower:\t%g\n", error);
    }

    // Upper
    {
        // generate matrix
        z2matgen(n, n, A1, A2);

        // run
        LAPACK(zlauum)("U", &n, A1, &n, &info);
        LAPACK(zlauu2)("U", &n, A2, &n, &info);

        // check error
        double error = z2vecerr(n * n, A1, A2);
        printf("zlauum Upper:\t%g\n", error);
    }

    free(A1); 
    free(A2);

	return 0;
}
