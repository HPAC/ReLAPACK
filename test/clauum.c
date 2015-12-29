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
		
	float *A1 = malloc(2 * n * n * sizeof(float));
	float *A2 = malloc(2 * n * n * sizeof(float));

    // Ouptut
    int info;

    { // L
        // generate matrix
        c2matgen(n, n, A1, A2);

        // run
        RELAPACK(clauum)("L", &n, A1, &n, &info);
        LAPACK(clauu2)("L", &n, A2, &n, &info);

        // check error
        const double error = c2vecerr(n * n, A1, A2);
        printf("clauum L:\t%g\n", error);
    }

    { // U
        // generate matrix
        c2matgen(n, n, A1, A2);

        // run
        RELAPACK(clauum)("U", &n, A1, &n, &info);
        LAPACK(clauu2)("U", &n, A2, &n, &info);

        // check error
        const double error = c2vecerr(n * n, A1, A2);
        printf("clauum U:\t%g\n", error);
    }

    free(A1);
    free(A2);

	return 0;
}
