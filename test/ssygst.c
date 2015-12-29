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
		
	float *A1 = malloc(n * n * sizeof(float));
	float *A2 = malloc(n * n * sizeof(float));
	float *B1 = malloc(n * n * sizeof(float));
	float *B2 = malloc(n * n * sizeof(float));

    // Output
    int info;

    // Constants
    const int iONE[] = {1};
    const int iTWO[] = {2};

    { // 1 L
        // generate matrix
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);

        // run
        RELAPACK(ssygst)(iONE, "L", &n, A1, &n, B1, &n, &info);
        LAPACK(ssygs2)(iONE, "L", &n, A2, &n, B2, &n, &info);

        // check error
        const double error = s2vecerr(n * n, A1, A2);
        printf("ssygst 1 L:\t%g\n", error);
    }

    { // 1 U
        // generate matrix
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);

        // run
        RELAPACK(ssygst)(iONE, "U", &n, A1, &n, B1, &n, &info);
        LAPACK(ssygs2)(iONE, "U", &n, A2, &n, B2, &n, &info);

        // check error
        const double error = s2vecerr(n * n, A1, A2);
        printf("ssygst 1 U:\t%g\n", error);
    }

    { // 2 L
        // generate matrix
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);

        // run
        RELAPACK(ssygst)(iTWO, "L", &n, A1, &n, B1, &n, &info);
        LAPACK(ssygs2)(iTWO, "L", &n, A2, &n, B2, &n, &info);

        // check error
        const double error = s2vecerr(n * n, A1, A2);
        printf("ssygst 2 L:\t%g\n", error);
    }

    { // 2 U
        // generate matrix
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);

        // run
        RELAPACK(ssygst)(iTWO, "U", &n, A1, &n, B1, &n, &info);
        LAPACK(ssygs2)(iTWO, "U", &n, A2, &n, B2, &n, &info);

        // check error
        const double error = s2vecerr(n * n, A1, A2);
        printf("ssygst 2 U:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);

	return 0;
}
