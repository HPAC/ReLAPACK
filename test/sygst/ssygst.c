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
	float *B1 = malloc(n * n * sizeof(float));
	float *B2 = malloc(n * n * sizeof(float));

    int info;

    int i1[] = {1};

    // Lower
    {
        // generate matrix
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);

        // run
        LAPACK(ssygst)(i1, "L", &n, A1, &n, B1, &n, &info);
        LAPACK(ssygs2)(i1, "L", &n, A2, &n, B1, &n, &info);

        // check error
        float error = s2vecerr(n * n, A1, A2);
        printf("ssygst 1 L:\t%g\n", error);
    }


    free(A1); 
    free(A2);
    free(B1); 
    free(B2);

	return 0;
}
