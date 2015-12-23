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
    const int n_max = n;
    const int n_min = MAX(1, (n * 3) / 4);
		
	double *A1 = malloc(n * n * sizeof(double));
	double *A2 = malloc(n * n * sizeof(double));
	double *B1 = malloc(n * n * sizeof(double));
	double *B2 = malloc(n * n * sizeof(double));
	double *C1 = malloc(n * n * sizeof(double));
	double *C2 = malloc(n * n * sizeof(double));

    int info;

    // 0, 1, -1
    const double d0[] = {0}, d1[] = {1}, sm1[] = {-1};
    // 0
    const int i0[] = {0};

    { // N N L 1 1
        const int n = n_max;
        // generate matrices
        d2matgen(n, n, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(dgemm_tr)("N", "N", "L", &n, &n, d1, A1, &n, B1, &n, d1, C1, &n);
        BLAS(dgemm)("N", "N", &n, &n, &n, d1, A2, &n, B2, &n, d1, C2, &n);

        // clear upper part of C
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const double error = d2vecerr(n * n, C1, C2);
        printf("cgemm_tr N N L 1 1:\t%g\n", error);
    }

    { // N N L 1 -1
        const int n = n_max;
        // generate matrices
        d2matgen(n, n, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(dgemm_tr)("N", "N", "L", &n, &n, d1, A1, &n, B1, &n, sm1, C1, &n);
        BLAS(dgemm)("N", "N", &n, &n, &n, d1, A2, &n, B2, &n, sm1, C2, &n);

        // clear upper part of C
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const double error = d2vecerr(n * n, C1, C2);
        printf("cgemm_tr N N L 1 -1:\t%g\n", error);
    }

    { // N T L 1 1
        const int n = n_max;
        // generate matrices
        d2matgen(n, n, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(dgemm_tr)("N", "T", "L", &n, &n, d1, A1, &n, B1, &n, d1, C1, &n);
        BLAS(dgemm)("N", "T", &n, &n, &n, d1, A2, &n, B2, &n, d1, C2, &n);

        // clear upper part of C
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const double error = d2vecerr(n * n, C1, C2);
        printf("cgemm_tr N T L 1 1:\t%g\n", error);
    }

    { // T N L 1 1
        const int n = n_max;
        // generate matrices
        d2matgen(n, n, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(dgemm_tr)("T", "N", "L", &n, &n, d1, A1, &n, B1, &n, d1, C1, &n);
        BLAS(dgemm)("T", "N", &n, &n, &n, d1, A2, &n, B2, &n, d1, C2, &n);

        // clear upper part of C
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const double error = d2vecerr(n * n, C1, C2);
        printf("cgemm_tr T N L 1 1:\t%g\n", error);
    }

    { // N N U 1 1
        const int n = n_max;
        // generate matrices
        d2matgen(n, n, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(n, n, C1, C2);

        // clear lower part of C
        const int nm1 = n - 1;
        LAPACK(dlascl)("L", i0, i0, d1, d0, &nm1, &nm1, C1 + 2, &n, &info);

        // run
        RELAPACK(dgemm_tr)("N", "N", "U", &n, &n, d1, A1, &n, B1, &n, d1, C1, &n);
        BLAS(dgemm)("N", "N", &n, &n, &n, d1, A2, &n, B2, &n, d1, C2, &n);

        // clear upper part of C
        LAPACK(dlascl)("L", i0, i0, d1, d0, &nm1, &nm1, C2 + 2, &n, &info);

        // check error
        const double error = d2vecerr(n * n, C1, C2);
        printf("cgemm_tr N N U 1 1:\t%g\n", error);
    }

    { // smallk
        const int n = n_max;
        const int k = n_min;
        // generate matrices
        d2matgen(n, k, A1, A2);
        d2matgen(k, n, B1, B2);
        d2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(dgemm_tr)("N", "N", "L", &n, &k, d1, A1, &n, B1, &n, d1, C1, &n);
        BLAS(dgemm)("N", "N", &n, &n, &k, d1, A2, &n, B2, &n, d1, C2, &n);

        // clear upper part of C
        LAPACK(dlascl)("U", i0, i0, d1, d0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const double error = d2vecerr(n * n, C1, C2);
        printf("cgemm_tr smallk:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);

	return 0;
}
