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
		
	float *A1 = malloc(n * n * sizeof(float));
	float *A2 = malloc(n * n * sizeof(float));
	float *B1 = malloc(n * n * sizeof(float));
	float *B2 = malloc(n * n * sizeof(float));
	float *C1 = malloc(n * n * sizeof(float));
	float *C2 = malloc(n * n * sizeof(float));

    int info;

    // 0, 1, -1
    const float s0[] = {0}, s1[] = {1}, sm1[] = {-1};
    // 0
    const int i0[] = {0};

    { // N N L 1 1
        const int n = n_max;
        // generate matrices
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(sgemm_tr)("N", "N", "L", &n, &n, s1, A1, &n, B1, &n, s1, C1, &n);
        BLAS(sgemm)("N", "N", &n, &n, &n, s1, A2, &n, B2, &n, s1, C2, &n);

        // clear upper part of C
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const float error = s2vecerr(n * n, C1, C2);
        printf("gemm_tr N N L 1 1:\t%g\n", error);
    }
    
    { // N N L 1 -1
        const int n = n_max;
        // generate matrices
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(sgemm_tr)("N", "N", "L", &n, &n, s1, A1, &n, B1, &n, sm1, C1, &n);
        BLAS(sgemm)("N", "N", &n, &n, &n, s1, A2, &n, B2, &n, sm1, C2, &n);

        // clear upper part of C
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const float error = s2vecerr(n * n, C1, C2);
        printf("gemm_tr N N L 1 -1:\t%g\n", error);
    }

    { // N N L 0 1
        const int n = n_max;
        // generate matrices
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(sgemm_tr)("N", "N", "L", &n, &n, s0, A1, &n, B1, &n, s1, C1, &n);
        BLAS(sgemm)("N", "N", &n, &n, &n, s0, A2, &n, B2, &n, s1, C2, &n);

        // clear upper part of C
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const float error = s2vecerr(n * n, C1, C2);
        printf("gemm_tr N N L 0 1:\t%g\n", error);
    }

    { // N T L 1 1
        const int n = n_max;
        // generate matrices
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(sgemm_tr)("N", "T", "L", &n, &n, s0, A1, &n, B1, &n, s1, C1, &n);
        BLAS(sgemm)("N", "T", &n, &n, &n, s0, A2, &n, B2, &n, s1, C2, &n);

        // clear upper part of C
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const float error = s2vecerr(n * n, C1, C2);
        printf("gemm_tr N T L 1 1:\t%g\n", error);
    }

    { // T N L 1 1
        const int n = n_max;
        // generate matrices
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(sgemm_tr)("T", "N", "L", &n, &n, s0, A1, &n, B1, &n, s1, C1, &n);
        BLAS(sgemm)("T", "N", &n, &n, &n, s0, A2, &n, B2, &n, s1, C2, &n);

        // clear upper part of C
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const float error = s2vecerr(n * n, C1, C2);
        printf("gemm_tr N T L 1 1:\t%g\n", error);
    }

    { // N N U 1 1
        const int n = n_max;
        // generate matrices
        s2matgen(n, n, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(n, n, C1, C2);

        // clear lower part of C
        const int nm1 = n - 1;
        LAPACK(slascl)("L", i0, i0, s1, s0, &nm1, &nm1, C1 + 1, &n, &info);

        // run
        RELAPACK(sgemm_tr)("T", "N", "U", &n, &n, s0, A1, &n, B1, &n, s1, C1, &n);
        BLAS(sgemm)("T", "N", &n, &n, &n, s0, A2, &n, B2, &n, s1, C2, &n);

        // clear upper part of C
        LAPACK(slascl)("L", i0, i0, s1, s0, &nm1, &nm1, C2 + 1, &n, &info);

        // check error
        const float error = s2vecerr(n * n, C1, C2);
        printf("gemm_tr N T U 1 1:\t%g\n", error);
    }

    { // smallk
        const int n = n_max;
        const int k = n_min;
        // generate matrices
        s2matgen(n, k, A1, A2);
        s2matgen(k, n, B1, B2);
        s2matgen(n, n, C1, C2);

        // clear upper part of C
        const int nm1 = n - 1;
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C1 + n, &n, &info);

        // run
        RELAPACK(sgemm_tr)("N", "T", "L", &n, &k, s0, A1, &n, B1, &n, s1, C1, &n);
        BLAS(sgemm)("N", "T", &n, &n, &k, s0, A2, &n, B2, &n, s1, C2, &n);

        // clear upper part of C
        LAPACK(slascl)("U", i0, i0, s1, s0, &nm1, &nm1, C2 + n, &n, &info);

        // check error
        const float error = s2vecerr(n * n, C1, C2);
        printf("gemm_tr smallk:\t%g\n", error);
    }

    free(A1); 
    free(A2);
    free(B1); 
    free(B2);
    free(C1); 
    free(C2);

	return 0;
}
