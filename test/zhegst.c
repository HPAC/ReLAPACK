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
		
	double *A1 = malloc(n * n * 2 * sizeof(double));
	double *A2 = malloc(n * n * 2 * sizeof(double));
	double *B1 = malloc(n * n * 2 * sizeof(double));
	double *B2 = malloc(n * n * 2 * sizeof(double));

    // Output
    int info;

    // Constants
    const int iONE[] = {1};
    const int iTWO[] = {2};

    { // 1 L
        // generate matrix
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);

        // run
        RELAPACK(zhegst)(iONE, "L", &n, A1, &n, B1, &n, &info);
        LAPACK(zhegs2)(iONE, "L", &n, A2, &n, B2, &n, &info);

        // check error
        const double error = z2vecerr(n * n, A1, A2);
        printf("zhegst 1 L:\t%g\n", error);
    }

    { // 1 U
        // generate matrix
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);

        // run
        RELAPACK(zhegst)(iONE, "U", &n, A1, &n, B1, &n, &info);
        LAPACK(zhegs2)(iONE, "U", &n, A2, &n, B2, &n, &info);

        // check error
        const double error = z2vecerr(n * n, A1, A2);
        printf("zhegst 1 U:\t%g\n", error);
    }

    { // 2 L
        // generate matrix
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);

        // run
        RELAPACK(zhegst)(iTWO, "L", &n, A1, &n, B1, &n, &info);
        LAPACK(zhegs2)(iTWO, "L", &n, A2, &n, B2, &n, &info);

        // check error
        const double error = z2vecerr(n * n, A1, A2);
        printf("zhegst 2 L:\t%g\n", error);
    }

    { // 2 U
        // generate matrix
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);

        // run
        RELAPACK(zhegst)(iTWO, "U", &n, A1, &n, B1, &n, &info);
        LAPACK(zhegs2)(iTWO, "U", &n, A2, &n, B2, &n, &info);

        // check error
        const double error = z2vecerr(n * n, A1, A2);
        printf("zhegst 2 U:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);

	return 0;
}
