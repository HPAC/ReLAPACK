#include "../src/relapack.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

void LAPACK(ctgsyl)(const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, float *, float *, const int *, int *, int *); 

int main(int argc, char* argv[]) {

    if (argc == 1) {
        fprintf(stderr, "usage: %s n\n", argv[0]);
        return 0;
    }
    const int n = atoi(argv[1]);
    const int n_max = n;
    const int n_min = MAX(1, (n * 3) / 4);
		
	float *A1 = malloc(2 * n * n * sizeof(float));
	float *A2 = malloc(2 * n * n * sizeof(float));
	float *B1 = malloc(2 * n * n * sizeof(float));
	float *B2 = malloc(2 * n * n * sizeof(float));
	float *C1 = malloc(2 * n * n * sizeof(float));
	float *C2 = malloc(2 * n * n * sizeof(float));
	float *D1 = malloc(2 * n * n * sizeof(float));
	float *D2 = malloc(2 * n * n * sizeof(float));
	float *E1 = malloc(2 * n * n * sizeof(float));
	float *E2 = malloc(2 * n * n * sizeof(float));
	float *F1 = malloc(2 * n * n * sizeof(float));
	float *F2 = malloc(2 * n * n * sizeof(float));
    const int lWork = 2 * n * n;
    float *Work1 = malloc(2 * lWork  * sizeof(float));
    float *Work2 = malloc(2 * lWork  * sizeof(float));
    int *iWork1 = malloc((n + n + 2) * sizeof(int));
    int *iWork2 = malloc((n + n + 2) * sizeof(int));

    int info;
    float scale1, scale2;
    float dif1, dif2;
    // 0, 1, 2, 3, 4
    const int i0[] = {0}, i1[] = {1}, i2[] = {2}, i3[] = {3}, i4[] = {4};

    { // N 0 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);
        c2matgen(m, m, D1, D2);
        c2matgen(n, n, E1, E2);
        c2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const float cmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(cscal)(&m, cmi, A1, &mp1);
        BLAS(cscal)(&m, cmi, A2, &mp1);
        BLAS(cscal)(&m, cmi, D1, &mp1);
        BLAS(cscal)(&m, cmi, D2, &mp1);

        // run
        RELAPACK(ctgsyl)("N", i0, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ctgsyl)("N", i0, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = c2vecerr(m * n, C1, C2) + c2vecerr(m * n, F1, F2);;
        printf("ctgsyl N 0 m = n:\t%g\n", error);
    }

    { // N 0 m < n
        const int m = n_min, n = n_max;
        // generate matrices
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);
        c2matgen(m, m, D1, D2);
        c2matgen(n, n, E1, E2);
        c2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const float cmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(cscal)(&m, cmi, A1, &mp1);
        BLAS(cscal)(&m, cmi, A2, &mp1);
        BLAS(cscal)(&m, cmi, D1, &mp1);
        BLAS(cscal)(&m, cmi, D2, &mp1);

        // run
        RELAPACK(ctgsyl)("N", i0, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ctgsyl)("N", i0, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = c2vecerr(m * n, C1, C2) + c2vecerr(m * n, F1, F2);;
        printf("ctgsyl N 0 m < n:\t%g\n", error);
    }

    { // N 0 m > n
        const int m = n_max, n = n_min;
        // generate matrices
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);
        c2matgen(m, m, D1, D2);
        c2matgen(n, n, E1, E2);
        c2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const float cmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(cscal)(&m, cmi, A1, &mp1);
        BLAS(cscal)(&m, cmi, A2, &mp1);
        BLAS(cscal)(&m, cmi, D1, &mp1);
        BLAS(cscal)(&m, cmi, D2, &mp1);

        // run
        RELAPACK(ctgsyl)("N", i0, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ctgsyl)("N", i0, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = c2vecerr(m * n, C1, C2) + c2vecerr(m * n, F1, F2);;
        printf("ctgsyl N 0 m > n:\t%g\n", error);
    }

    { // N 1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);
        c2matgen(m, m, D1, D2);
        c2matgen(n, n, E1, E2);
        c2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const float cmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(cscal)(&m, cmi, A1, &mp1);
        BLAS(cscal)(&m, cmi, A2, &mp1);
        BLAS(cscal)(&m, cmi, D1, &mp1);
        BLAS(cscal)(&m, cmi, D2, &mp1);

        // run
        RELAPACK(ctgsyl)("N", i1, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ctgsyl)("N", i1, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = c2vecerr(m * n, C1, C2) + c2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("ctgsyl N 1 m = n:\t%g\n", error);
    }

    { // N 2 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);
        c2matgen(m, m, D1, D2);
        c2matgen(n, n, E1, E2);
        c2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const float cmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(cscal)(&m, cmi, A1, &mp1);
        BLAS(cscal)(&m, cmi, A2, &mp1);
        BLAS(cscal)(&m, cmi, D1, &mp1);
        BLAS(cscal)(&m, cmi, D2, &mp1);

        // run
        RELAPACK(ctgsyl)("N", i2, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ctgsyl)("N", i2, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = c2vecerr(m * n, C1, C2) + c2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("ctgsyl N 2 m = n:\t%g\n", error);
    }

    { // N 3 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);
        c2matgen(m, m, D1, D2);
        c2matgen(n, n, E1, E2);
        c2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const float cmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(cscal)(&m, cmi, A1, &mp1);
        BLAS(cscal)(&m, cmi, A2, &mp1);
        BLAS(cscal)(&m, cmi, D1, &mp1);
        BLAS(cscal)(&m, cmi, D2, &mp1);

        // run
        RELAPACK(ctgsyl)("N", i3, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ctgsyl)("N", i3, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = c2vecerr(m * n, C1, C2) + c2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("ctgsyl N 3 m = n:\t%g\n", error);
    }

    { // N 4 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);
        c2matgen(m, m, D1, D2);
        c2matgen(n, n, E1, E2);
        c2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const float cmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(cscal)(&m, cmi, A1, &mp1);
        BLAS(cscal)(&m, cmi, A2, &mp1);
        BLAS(cscal)(&m, cmi, D1, &mp1);
        BLAS(cscal)(&m, cmi, D2, &mp1);

        // run
        RELAPACK(ctgsyl)("N", i4, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ctgsyl)("N", i4, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = c2vecerr(m * n, C1, C2) + c2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("ctgsyl N 4 m = n:\t%g\n", error);
    }

    { // C 0 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);
        c2matgen(m, m, D1, D2);
        c2matgen(n, n, E1, E2);
        c2matgen(m, n, F1, F2);

        // run
        RELAPACK(ctgsyl)("C", i0, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ctgsyl)("C", i0, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = c2vecerr(m * n, C1, C2) + c2vecerr(m * n, F1, F2);
        printf("ctgsyl C 0 m = n:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);
    free(D1);
    free(D2);
    free(E1);
    free(E2);
    free(F1);
    free(F2);
    free(Work1);
    free(Work2);
    free(iWork1);
    free(iWork2);

	return 0;
}
