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
		
	double *A1 = malloc(n * n * 2 * sizeof(double));
	double *A2 = malloc(n * n * 2 * sizeof(double));
	double *B1 = malloc(n * n * 2 * sizeof(double));
	double *B2 = malloc(n * n * 2 * sizeof(double));
	double *C1 = malloc(n * n * 2 * sizeof(double));
	double *C2 = malloc(n * n * 2 * sizeof(double));

    int info;

    // 0, 1, -1
    const double z0[] = {0, 0}, z1[] = {1, 0}, zm1[] = {-1, 0};
    // 0
    const int i0[] = {0};

    { // N N L 1 1
        const int n = n_max;
        // generate matrices
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C1 +  2 * n, &n, &info);

        // run
        RELAPACK(zgemm_tr)("N", "N", "L", &n, &n, z1, A1, &n, B1, &n, z1, C1, &n);
        BLAS(zgemm)("N", "N", &n, &n, &n, z1, A2, &n, B2, &n, z1, C2, &n);

        // clear upper part of C
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = z2vecerr(n * n, C1, C2);
        printf("zgemm_tr N N L 1 1:\t%g\n", error);
    }

    { // N N L 1 -1
        const int n = n_max;
        // generate matrices
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C1 +  2 * n, &n, &info);

        // run
        RELAPACK(zgemm_tr)("N", "N", "L", &n, &n, z1, A1, &n, B1, &n, zm1, C1, &n);
        BLAS(zgemm)("N", "N", &n, &n, &n, z1, A2, &n, B2, &n, zm1, C2, &n);

        // clear upper part of C
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = z2vecerr(n * n, C1, C2);
        printf("zgemm_tr N N L 1 -1:\t%g\n", error);
    }

    { // N T L 1 1
        const int n = n_max;
        // generate matrices
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C1 +  2 * n, &n, &info);

        // run
        RELAPACK(zgemm_tr)("N", "T", "L", &n, &n, z1, A1, &n, B1, &n, z1, C1, &n);
        BLAS(zgemm)("N", "T", &n, &n, &n, z1, A2, &n, B2, &n, z1, C2, &n);

        // clear upper part of C
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = z2vecerr(n * n, C1, C2);
        printf("zgemm_tr N T L 1 1:\t%g\n", error);
    }

    { // T N L 1 1
        const int n = n_max;
        // generate matrices
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C1 +  2 * n, &n, &info);

        // run
        RELAPACK(zgemm_tr)("T", "N", "L", &n, &n, z1, A1, &n, B1, &n, z1, C1, &n);
        BLAS(zgemm)("T", "N", &n, &n, &n, z1, A2, &n, B2, &n, z1, C2, &n);

        // clear upper part of C
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = z2vecerr(n * n, C1, C2);
        printf("zgemm_tr T N L 1 1:\t%g\n", error);
    }

    { // N N U 1 1
        const int n = n_max;
        // generate matrices
        z2matgen(n, n, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(n, n, C1, C2);

        // clear lower part of C
        const int nm1 = n - 1;
        LAPACK(zlascl)("L", i0, i0, z1, z0, &nm1, &nm1, C1 +  2, &n, &info);

        // run
        RELAPACK(zgemm_tr)("N", "N", "U", &n, &n, z1, A1, &n, B1, &n, z1, C1, &n);
        BLAS(zgemm)("N", "N", &n, &n, &n, z1, A2, &n, B2, &n, z1, C2, &n);

        // clear upper part of C
        LAPACK(zlascl)("L", i0, i0, z1, z0, &nm1, &nm1, C2 + 2, &n, &info);

        // check error
        const double error = z2vecerr(n * n, C1, C2);
        printf("zgemm_tr N N U 1 1:\t%g\n", error);
    }

    { // smallk
        const int n = n_max;
        const int k = n_min;
        // generate matrices
        z2matgen(n, k, A1, A2);
        z2matgen(k, n, B1, B2);
        z2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C1 +  2 * n, &n, &info);

        // run
        RELAPACK(zgemm_tr)("N", "N", "L", &n, &k, z1, A1, &n, B1, &n, z1, C1, &n);
        BLAS(zgemm)("N", "N", &n, &n, &k, z1, A2, &n, B2, &n, z1, C2, &n);

        // clear upper part of C
        LAPACK(zlascl)("U", i0, i0, z1, z0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = z2vecerr(n * n, C1, C2);
        printf("zgemm_tr smallk:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);

	return 0;
}
