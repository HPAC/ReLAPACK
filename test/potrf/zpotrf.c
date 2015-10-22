#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "LAPACK_ORIG_potrf.h"
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
        LAPACK(zpotrf)("L", &n, A1, &n, &info);
        LAPACK_ORIG(zpotrf)("L", &n, A2, &n, &info);

        // check error
        double error = z2vecerr(n * n, A1, A2);
        printf("zpotrf Lower:\t%g\n", error);
    }

    // Upper
    {
        // generate matrix
        z2matgen(n, n, A1, A2);

        // run
        LAPACK(zpotrf)("U", &n, A1, &n, &info);
        LAPACK_ORIG(zpotrf)("U", &n, A2, &n, &info);

        // check error
        double error = z2vecerr(n * n, A1, A2);
        printf("zpotrf Upper:\t%g\n", error);
    }

    free(A1);
    free(A2);

	return 0;
}
