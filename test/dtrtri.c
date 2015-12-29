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
		
	double *A1 = malloc(n * n * sizeof(double));
	double *A2 = malloc(n * n * sizeof(double));

    // Output
    int info;

    { // L
        // generate matrix
        d2matgen(n, n, A1, A2);

        // run
        RELAPACK(dtrtri)("L", "N", &n, A1, &n, &info);
        LAPACK(dtrti2)("L", "N", &n, A2, &n, &info);

        // check error
        const double error = d2vecerr(n * n, A1, A2);
        printf("dtrtri L:\t%g\n", error);
    }

    { // U
        // generate matrix
        d2matgen(n, n, A1, A2);

        // run
        RELAPACK(dtrtri)("U", "N", &n, A1, &n, &info);
        LAPACK(dtrti2)("U", "N", &n, A2, &n, &info);

        // check error
        const double error = d2vecerr(n * n, A1, A2);
        printf("dtrtri U:\t%g\n", error);
    }

    free(A1);
    free(A2);

	return 0;
}
