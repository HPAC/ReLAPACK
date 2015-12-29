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
		
	float *A1 = malloc(n * n * 2 * sizeof(float));
	float *A2 = malloc(n * n * 2 * sizeof(float));
	float *B1 = malloc(n * n * 2 * sizeof(float));
	float *B2 = malloc(n * n * 2 * sizeof(float));
	float *C1 = malloc(n * n * 2 * sizeof(float));
	float *C2 = malloc(n * n * 2 * sizeof(float));

    int info;

    // 0, 1, -1
    const float c0[] = {0, 0}, c1[] = {1, 0}, cm1[] = {-1, 0};
    // 0
    const int i0[] = {0};

    { // N N L 1 1
        const int n = n_max;
        // generate matrices
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C1 + 2 * n, &n, &info);

        // run
        RELAPACK(cgemm_tr_rec)("N", "N", "L", &n, &n, c1, A1, &n, B1, &n, c1, C1, &n);
        BLAS(cgemm)("N", "N", &n, &n, &n, c1, A2, &n, B2, &n, c1, C2, &n);

        // clear upper part of C
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = c2vecerr(n * n, C1, C2);
        printf("cgemm_tr_rec N N L 1 1:\t%g\n", error);
    }

    { // N N L 1 -1
        const int n = n_max;
        // generate matrices
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C1 + 2 * n, &n, &info);

        // run
        RELAPACK(cgemm_tr_rec)("N", "N", "L", &n, &n, c1, A1, &n, B1, &n, cm1, C1, &n);
        BLAS(cgemm)("N", "N", &n, &n, &n, c1, A2, &n, B2, &n, cm1, C2, &n);

        // clear upper part of C
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = c2vecerr(n * n, C1, C2);
        printf("cgemm_tr_rec N N L 1 -1:\t%g\n", error);
    }

    { // N T L 1 1
        const int n = n_max;
        // generate matrices
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C1 + 2 * n, &n, &info);

        // run
        RELAPACK(cgemm_tr_rec)("N", "T", "L", &n, &n, c1, A1, &n, B1, &n, c1, C1, &n);
        BLAS(cgemm)("N", "T", &n, &n, &n, c1, A2, &n, B2, &n, c1, C2, &n);

        // clear upper part of C
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = c2vecerr(n * n, C1, C2);
        printf("cgemm_tr_rec N T L 1 1:\t%g\n", error);
    }

    { // T N L 1 1
        const int n = n_max;
        // generate matrices
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C1 + 2 * n, &n, &info);

        // run
        RELAPACK(cgemm_tr_rec)("T", "N", "L", &n, &n, c1, A1, &n, B1, &n, c1, C1, &n);
        BLAS(cgemm)("T", "N", &n, &n, &n, c1, A2, &n, B2, &n, c1, C2, &n);

        // clear upper part of C
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = c2vecerr(n * n, C1, C2);
        printf("cgemm_tr_rec T N L 1 1:\t%g\n", error);
    }

    { // N N U 1 1
        const int n = n_max;
        // generate matrices
        c2matgen(n, n, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(n, n, C1, C2);

        // clear lower part of C
        const int nm1 = n - 1;
        LAPACK(clascl)("L", i0, i0, c1, c0, &nm1, &nm1, C1 + 2, &n, &info);

        // run
        RELAPACK(cgemm_tr_rec)("N", "N", "U", &n, &n, c1, A1, &n, B1, &n, c1, C1, &n);
        BLAS(cgemm)("N", "N", &n, &n, &n, c1, A2, &n, B2, &n, c1, C2, &n);

        // clear upper part of C
        LAPACK(clascl)("L", i0, i0, c1, c0, &nm1, &nm1, C2 + 2, &n, &info);

        // check error
        const double error = c2vecerr(n * n, C1, C2);
        printf("cgemm_tr_rec N N U 1 1:\t%g\n", error);
    }

    { // smallk
        const int n = n_max;
        const int k = n_min;
        // generate matrices
        c2matgen(n, k, A1, A2);
        c2matgen(k, n, B1, B2);
        c2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C1 + 2 * n, &n, &info);

        // run
        RELAPACK(cgemm_tr_rec)("N", "N", "L", &n, &k, c1, A1, &n, B1, &n, c1, C1, &n);
        BLAS(cgemm)("N", "N", &n, &n, &k, c1, A2, &n, B2, &n, c1, C2, &n);

        // clear upper part of C
        LAPACK(clascl)("U", i0, i0, c1, c0, &nm1, &nm1, C2 + 2 * n, &n, &info);

        // check error
        const double error = c2vecerr(n * n, C1, C2);
        printf("cgemm_tr_rec smallk:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);

	return 0;
}
