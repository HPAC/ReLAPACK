#include "../src/relapack.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

    if (argc == 1) {
        fprintf(stderr, "usage: %s n\n", argv[0]);
        return 0;
    }
    const int n     = atoi(argv[1]);
    const int n_max = n;
    const int n_min = MAX(1, (n * 3) / 4);
		
	double *A1 = malloc(2 * n * n * sizeof(double));
	double *A2 = malloc(2 * n * n * sizeof(double));
	double *B1 = malloc(2 * n * n * sizeof(double));
	double *B2 = malloc(2 * n * n * sizeof(double));
	double *C1 = malloc(2 * n * n * sizeof(double));
	double *C2 = malloc(2 * n * n * sizeof(double));

    // Outputs
    int info;
    double scale1;
    double scale2;

    // Constants
    const int iZERO[] = {0};
    const int iONE[]  = {1};
    const int iMONE[] = {-1};

    { // N N +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);

        // run
        RELAPACK(ztrsyl)("N", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        ZTRSY2("N", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = z2vecerr(m * n, C1, C2);
        printf("ztrsyl N N +1 m = n:\t%g\n", error);
    }

    { // N N +1 m < n
        const int m = n_min, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);

        // run
        RELAPACK(ztrsyl)("N", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        ZTRSY2("N", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = z2vecerr(m * n, C1, C2);
        printf("ztrsyl N N +1 m < n:\t%g\n", error);
    }

    { // N N +1 m > n
        const int m = n_max, n = n_min;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);

        // run
        RELAPACK(ztrsyl)("N", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        ZTRSY2("N", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = z2vecerr(m * n, C1, C2);
        printf("ztrsyl N N +1 m > n:\t%g\n", error);
    }

    { // C N +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);

        // run
        RELAPACK(ztrsyl)("C", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        ZTRSY2("C", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = z2vecerr(m * n, C1, C2);
        printf("ztrsyl C N +1 m = n:\t%g\n", error);
    }

    { // N C +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);

        // run
        RELAPACK(ztrsyl)("N", "C", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        ZTRSY2("N", "C", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = z2vecerr(m * n, C1, C2);
        printf("ztrsyl N C +1 m = n:\t%g\n", error);
    }

    { // C C +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);

        // run
        RELAPACK(ztrsyl)("C", "C", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        ZTRSY2("C", "C", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = z2vecerr(m * n, C1, C2);
        printf("ztrsyl C C +1 m = n:\t%g\n", error);
    }

    { // N N -1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);

        // scale diagonal of A
        const double smi = 1. / m;
        const int mp1 = m + 1;
        BLAS(zscal)(&m, &smi, A1, &mp1);
        BLAS(zscal)(&m, &smi, A2, &mp1);

        // run
        RELAPACK(ztrsyl)("N", "N", iMONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        ZTRSY2("N", "N", iMONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = z2vecerr(m * n, C1, C2);
        printf("ztrsyl N N -1 m = n:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);

	return 0;
}
