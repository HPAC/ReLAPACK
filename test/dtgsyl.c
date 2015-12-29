#include "../src/relapack.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

void LAPACK(dtgsyl)(const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, double *, double *, const int *, int *, int *);

int main(int argc, char* argv[]) {

    if (argc == 1) {
        fprintf(stderr, "usage: %s n\n", argv[0]);
        return 0;
    }
    const int n     = atoi(argv[1]);
    const int n_max = n;
    const int n_min = MAX(1, (n * 3) / 4);
		
    const int lWork = 2 * n * n;
	double *A1     = malloc(n * n * sizeof(double));
	double *A2     = malloc(n * n * sizeof(double));
	double *B1     = malloc(n * n * sizeof(double));
	double *B2     = malloc(n * n * sizeof(double));
	double *C1     = malloc(n * n * sizeof(double));
	double *C2     = malloc(n * n * sizeof(double));
	double *D1     = malloc(n * n * sizeof(double));
	double *D2     = malloc(n * n * sizeof(double));
	double *E1     = malloc(n * n * sizeof(double));
	double *E2     = malloc(n * n * sizeof(double));
	double *F1     = malloc(n * n * sizeof(double));
	double *F2     = malloc(n * n * sizeof(double));
    double *Work1  = malloc(lWork * sizeof(double));
    double *Work2  = malloc(lWork * sizeof(double));
    int    *iWork1 = malloc((n + n + 2) * sizeof(int));
    int    *iWork2 = malloc((n + n + 2) * sizeof(int));

    // Outputs
    int info;
    double scale1;
    double scale2;
    double dif1;
    double dif2;

    // Constants
    const double ZERO[]  = {0};
    const int   iZERO[]  = {0};
    const int   iONE[]   = {1};
    const int   iTWO[]   = {2};
    const int   iTHREE[] = {3};
    const int   iFOUR[]  = {4};

    { // N 0 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        // scale diagonals of A and D
        const double dmi = 1. / m;
        BLAS(dscal)(&m, &dmi, A1, &mp1);
        BLAS(dscal)(&m, &dmi, A2, &mp1);
        BLAS(dscal)(&m, &dmi, D1, &mp1);
        BLAS(dscal)(&m, &dmi, D2, &mp1);

        // run
        RELAPACK(dtgsyl)("N", iZERO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("N", iZERO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2);;
        printf("dtgsyl N 0 m = n:\t%g\n", error);
    }

    { // N 0 m < n
        const int m = n_min, n = n_max;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        // scale diagonals of A and D
        const double dmi = 1. / m;
        BLAS(dscal)(&m, &dmi, A1, &mp1);
        BLAS(dscal)(&m, &dmi, A2, &mp1);
        BLAS(dscal)(&m, &dmi, D1, &mp1);
        BLAS(dscal)(&m, &dmi, D2, &mp1);

        // run
        RELAPACK(dtgsyl)("N", iZERO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("N", iZERO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2);;
        printf("dtgsyl N 0 m < n:\t%g\n", error);
    }

    { // N 0 m > n
        const int m = n_max, n = n_min;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        // scale diagonals of A and D
        const double dmi = 1. / m;
        BLAS(dscal)(&m, &dmi, A1, &mp1);
        BLAS(dscal)(&m, &dmi, A2, &mp1);
        BLAS(dscal)(&m, &dmi, D1, &mp1);
        BLAS(dscal)(&m, &dmi, D2, &mp1);

        // run
        RELAPACK(dtgsyl)("N", iZERO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("N", iZERO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2);;
        printf("dtgsyl N 0 m > n:\t%g\n", error);
    }

    { // N 1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        // scale diagonals of A and D
        const double dmi = 1. / m;
        BLAS(dscal)(&m, &dmi, A1, &mp1);
        BLAS(dscal)(&m, &dmi, A2, &mp1);
        BLAS(dscal)(&m, &dmi, D1, &mp1);
        BLAS(dscal)(&m, &dmi, D2, &mp1);

        // run
        RELAPACK(dtgsyl)("N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("dtgsyl N 1 m = n:\t%g\n", error);
    }

    { // N 2 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        // scale diagonals of A and D
        const double dmi = 1. / m;
        BLAS(dscal)(&m, &dmi, A1, &mp1);
        BLAS(dscal)(&m, &dmi, A2, &mp1);
        BLAS(dscal)(&m, &dmi, D1, &mp1);
        BLAS(dscal)(&m, &dmi, D2, &mp1);

        // run
        RELAPACK(dtgsyl)("N", iTWO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("N", iTWO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("dtgsyl N 2 m = n:\t%g\n", error);
    }

    { // N 3 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        // scale diagonals of A and D
        const double dmi = 1. / m;
        BLAS(dscal)(&m, &dmi, A1, &mp1);
        BLAS(dscal)(&m, &dmi, A2, &mp1);
        BLAS(dscal)(&m, &dmi, D1, &mp1);
        BLAS(dscal)(&m, &dmi, D2, &mp1);

        // run
        RELAPACK(dtgsyl)("N", iTHREE, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("N", iTHREE, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("dtgsyl N 3 m = n:\t%g\n", error);
    }

    { // N 4 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        // scale diagonals of A and D
        const double dmi = 1. / m;
        BLAS(dscal)(&m, &dmi, A1, &mp1);
        BLAS(dscal)(&m, &dmi, A2, &mp1);
        BLAS(dscal)(&m, &dmi, D1, &mp1);
        BLAS(dscal)(&m, &dmi, D2, &mp1);

        // run
        RELAPACK(dtgsyl)("N", iFOUR, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("N", iFOUR, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("dtgsyl N 4 m = n:\t%g\n", error);
    }

    { // N 0 offdiag
        const int m = n_max, n = n_max;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        int m1 = (m >= 16) ? ((m + 8) / 16) * 8 : m / 2;
        A1[m1 + m * (m1 - 1)] = A2[m1 + m * (m1 - 1)] = -.5;
        int n1 = (n >= 16) ? ((n + 8) / 16) * 8 : n / 2;
        B1[n1 + n * (n1 - 1)] = B2[n1 + n * (n1 - 1)] = -.3;
        n1 = n / 2;
        B1[n1 + n * (n1 - 1)] = B2[n1 + n * (n1 - 1)] = -.2;

        // scale diagonals of A and D
        const double dmi = 1. / m;
        BLAS(dscal)(&m, &dmi, A1, &mp1);
        BLAS(dscal)(&m, &dmi, A2, &mp1);
        BLAS(dscal)(&m, &dmi, D1, &mp1);
        BLAS(dscal)(&m, &dmi, D2, &mp1);

        // run
        RELAPACK(dtgsyl)("N", iFOUR, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("N", iFOUR, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2) + d2vecerr(1, &dif1, &dif2);
        printf("dtgsyl N 0 offdiag:\t%g\n", error);
    }

    { // T 0 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        d2matgen(m, m, D1, D2);
        d2matgen(n, n, E1, E2);
        d2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(dscal)(&nm1, ZERO, B2 + 1, &np1);

        // run
        RELAPACK(dtgsyl)("T", iZERO, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(dtgsyl)("T", iZERO, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(dlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = d2vecerr(m * n, C1, C2) + d2vecerr(m * n, F1, F2);
        printf("dtgsyl T 0 m = n:\t%g\n", error);
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
