#include "../src/relapack.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

void LAPACK(ztgsyl)(const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, double *, double *, const int *, int *, int *);

int main(int argc, char* argv[]) {

    if (argc == 1) {
        fprintf(stderr, "usage: %s n\n", argv[0]);
        return 0;
    }
    const int n     = atoi(argv[1]);
    const int n_max = n;
    const int n_min = MAX(1, (n * 3) / 4);
		
    const int lWork = 2 * n * n;
	double *A1     = malloc(n * n * 2 * sizeof(double));
	double *A2     = malloc(n * n * 2 * sizeof(double));
	double *B1     = malloc(n * n * 2 * sizeof(double));
	double *B2     = malloc(n * n * 2 * sizeof(double));
	double *C1     = malloc(n * n * 2 * sizeof(double));
	double *C2     = malloc(n * n * 2 * sizeof(double));
	double *D1     = malloc(n * n * 2 * sizeof(double));
	double *D2     = malloc(n * n * 2 * sizeof(double));
	double *E1     = malloc(n * n * 2 * sizeof(double));
	double *E2     = malloc(n * n * 2 * sizeof(double));
	double *F1     = malloc(n * n * 2 * sizeof(double));
	double *F2     = malloc(n * n * 2 * sizeof(double));
    double *Work1  = malloc(lWork * 2 * sizeof(double));
    double *Work2  = malloc(lWork * 2 * sizeof(double));
    int    *iWork1 = malloc((n + n + 2) * sizeof(int));
    int    *iWork2 = malloc((n + n + 2) * sizeof(int));

    // Outputs
    int info;
    double scale1;
    double scale2;
    double dif1;
    double dif2;

    // Constants
    const int iZERO[]  = {0};
    const int iONE[]   = {1};
    const int iTWO[]   = {2};
    const int iTHREE[] = {3};
    const int iFOUR[]  = {4};

    { // N 0 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);
        z2matgen(m, m, D1, D2);
        z2matgen(n, n, E1, E2);
        z2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const double zmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(zscal)(&m, zmi, A1, &mp1);
        BLAS(zscal)(&m, zmi, A2, &mp1);
        BLAS(zscal)(&m, zmi, D1, &mp1);
        BLAS(zscal)(&m, zmi, D2, &mp1);

        // run
        RELAPACK(ztgsyl)("N", iZERO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ztgsyl)("N", iZERO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = z2vecerr(m * n, C1, C2) + z2vecerr(m * n, F1, F2);;
        printf("ztgsyl N 0 m = n:\t%g\n", error);
    }

    { // N 0 m < n
        const int m = n_min, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);
        z2matgen(m, m, D1, D2);
        z2matgen(n, n, E1, E2);
        z2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const double zmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(zscal)(&m, zmi, A1, &mp1);
        BLAS(zscal)(&m, zmi, A2, &mp1);
        BLAS(zscal)(&m, zmi, D1, &mp1);
        BLAS(zscal)(&m, zmi, D2, &mp1);

        // run
        RELAPACK(ztgsyl)("N", iZERO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ztgsyl)("N", iZERO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = z2vecerr(m * n, C1, C2) + z2vecerr(m * n, F1, F2);;
        printf("ztgsyl N 0 m < n:\t%g\n", error);
    }

    { // N 0 m > n
        const int m = n_max, n = n_min;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);
        z2matgen(m, m, D1, D2);
        z2matgen(n, n, E1, E2);
        z2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const double zmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(zscal)(&m, zmi, A1, &mp1);
        BLAS(zscal)(&m, zmi, A2, &mp1);
        BLAS(zscal)(&m, zmi, D1, &mp1);
        BLAS(zscal)(&m, zmi, D2, &mp1);

        // run
        RELAPACK(ztgsyl)("N", iZERO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ztgsyl)("N", iZERO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = z2vecerr(m * n, C1, C2) + z2vecerr(m * n, F1, F2);;
        printf("ztgsyl N 0 m > n:\t%g\n", error);
    }

    { // N 1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);
        z2matgen(m, m, D1, D2);
        z2matgen(n, n, E1, E2);
        z2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const double zmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(zscal)(&m, zmi, A1, &mp1);
        BLAS(zscal)(&m, zmi, A2, &mp1);
        BLAS(zscal)(&m, zmi, D1, &mp1);
        BLAS(zscal)(&m, zmi, D2, &mp1);

        // run
        RELAPACK(ztgsyl)("N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ztgsyl)("N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = z2vecerr(m * n, C1, C2) + z2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("ztgsyl N 1 m = n:\t%g\n", error);
    }

    { // N 2 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);
        z2matgen(m, m, D1, D2);
        z2matgen(n, n, E1, E2);
        z2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const double zmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(zscal)(&m, zmi, A1, &mp1);
        BLAS(zscal)(&m, zmi, A2, &mp1);
        BLAS(zscal)(&m, zmi, D1, &mp1);
        BLAS(zscal)(&m, zmi, D2, &mp1);

        // run
        RELAPACK(ztgsyl)("N", iTWO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ztgsyl)("N", iTWO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = z2vecerr(m * n, C1, C2) + z2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("ztgsyl N 2 m = n:\t%g\n", error);
    }

    { // N 3 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);
        z2matgen(m, m, D1, D2);
        z2matgen(n, n, E1, E2);
        z2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const double zmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(zscal)(&m, zmi, A1, &mp1);
        BLAS(zscal)(&m, zmi, A2, &mp1);
        BLAS(zscal)(&m, zmi, D1, &mp1);
        BLAS(zscal)(&m, zmi, D2, &mp1);

        // run
        RELAPACK(ztgsyl)("N", iTHREE, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ztgsyl)("N", iTHREE, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = z2vecerr(m * n, C1, C2) + z2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("ztgsyl N 3 m = n:\t%g\n", error);
    }

    { // N 4 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);
        z2matgen(m, m, D1, D2);
        z2matgen(n, n, E1, E2);
        z2matgen(m, n, F1, F2);

        // scale diagonals of A and D
        const double zmi[] = {1. / m, 0};
        const int mp1 = m + 1;
        BLAS(zscal)(&m, zmi, A1, &mp1);
        BLAS(zscal)(&m, zmi, A2, &mp1);
        BLAS(zscal)(&m, zmi, D1, &mp1);
        BLAS(zscal)(&m, zmi, D2, &mp1);

        // run
        RELAPACK(ztgsyl)("N", iFOUR, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ztgsyl)("N", iFOUR, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = z2vecerr(m * n, C1, C2) + z2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("ztgsyl N 4 m = n:\t%g\n", error);
    }

    { // C 0 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);
        z2matgen(m, m, D1, D2);
        z2matgen(n, n, E1, E2);
        z2matgen(m, n, F1, F2);

        // run
        RELAPACK(ztgsyl)("C", iZERO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(ztgsyl)("C", iZERO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = z2vecerr(m * n, C1, C2) + z2vecerr(m * n, F1, F2);
        printf("ztgsyl C 0 m = n:\t%g\n", error);
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
