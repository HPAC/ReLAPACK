#include "../../config.h"
#include "../../src/lapack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n = TEST_N;
		
	float *A1 = malloc(n * n * sizeof(float));
	float *A2 = malloc(n * n * sizeof(float));

    int info;

    // Lower
    {
        // generate matrix
        s2matgen(n, n, A1, A2);

        // run
        LAPACK(spotrf)("L", &n, A1, &n, &info);
        LAPACK(spotf2)("L", &n, A2, &n, &info);

        // check error
        float error = s2vecerr(n * n, A1, A2);
        printf("spotrf Lower:\t%g\n", error);
    }

    // Upper
    {
        // generate matrix
        s2matgen(n, n, A1, A2);

        // run
        LAPACK(spotrf)("U", &n, A1, &n, &info);
        LAPACK(spotf2)("U", &n, A2, &n, &info);

        // check error
        float error = s2vecerr(n * n, A1, A2);
        printf("spotrf Upper:\t%g\n", error);
    }

    free(A1);
    free(A2);

	return 0;
}
