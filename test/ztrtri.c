#include "../src/larpack.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

    if (argc == 1) {
        fprintf(stderr, "usage: %s n\n", argv[0]);
        return 0;
    }
    const int n = atoi(argv[1]);
		
	double *A1 = malloc(2 * n * n * sizeof(double));
	double *A2 = malloc(2 * n * n * sizeof(double));

    int info;

    // L
    {
        // generate matrix
        z2matgen(n, n, A1, A2);

        // run
        LARPACK(ztrtri)("L", "N", &n, A1, &n, &info);
        LAPACK(ztrti2)("L", "N", &n, A2, &n, &info);

        // check error
        const double error = z2vecerr(n * n, A1, A2);
        printf("ztrtri L:\t%g\n", error);
    }

    // U
    {
        // generate matrix
        z2matgen(n, n, A1, A2);

        // run
        LARPACK(ztrtri)("U", "N", &n, A1, &n, &info);
        LAPACK(ztrti2)("U", "N", &n, A2, &n, &info);

        // check error
        const double error = z2vecerr(n * n, A1, A2);
        printf("ztrtri U:\t%g\n", error);
    }

    free(A1); 
    free(A2);

	return 0;
}
