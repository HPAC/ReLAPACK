#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "LAPACK_ORIG_potrf.h"
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
        d2matgen(n, n, A1, A2);

        // run
        LAPACK(dpotrf)("L", &n, A1, &n, &info);
        LAPACK_ORIG(dpotrf)("L", &n, A2, &n, &info);

        // check error
        double error = d2vecerr(n * n, A1, A2);
        printf("dpotrf Lower:\t%g\n", error);
    }

    // Upper
    {
        // generate matrix
        d2matgen(n, n, A1, A2);

        // run
        LAPACK(dpotrf)("U", &n, A1, &n, &info);
        LAPACK_ORIG(dpotrf)("U", &n, A2, &n, &info);

        // check error
        double error = d2vecerr(n * n, A1, A2);
        printf("dpotrf Upper:\t%g\n", error);
    }

    free(A1);
    free(A2);

	return 0;
}
