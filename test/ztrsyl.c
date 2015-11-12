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
    const int n_max = n;
    const int n_min = MAX(1, (n * 3) / 4);
		
	double *A1 = malloc(2 * n * n * sizeof(double));
	double *A2 = malloc(2 * n * n * sizeof(double));
	double *B1 = malloc(2 * n * n * sizeof(double));
	double *B2 = malloc(2 * n * n * sizeof(double));
	double *C1 = malloc(2 * n * n * sizeof(double));
	double *C2 = malloc(2 * n * n * sizeof(double));

    int info;
    double scale1, scale2;
    // 0, 1, -1
    const int i0[] = {0}, i1[] = {1}, im1[] = {-1};

    { // N N +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        z2matgen(m, m, A1, A2);
        z2matgen(n, n, B1, B2);
        z2matgen(m, n, C1, C2);

        // run
        LARPACK(ztrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ztrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

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
        LARPACK(ztrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ztrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

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
        LARPACK(ztrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ztrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

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
        LARPACK(ztrsyl)("C", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ztrsy2)("C", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

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
        LARPACK(ztrsyl)("N", "C", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ztrsy2)("N", "C", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

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
        LARPACK(ztrsyl)("C", "C", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ztrsy2)("C", "C", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

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
        LARPACK(ztrsyl)("N", "N", im1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ztrsy2)("N", "N", im1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(zlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

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
