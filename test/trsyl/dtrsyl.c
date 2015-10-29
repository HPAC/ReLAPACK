#include "../../src/larpack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n_max = TEST_N;
    const int n_min = MAX(1, (n_max * 3) / 4);
    const int n = n_max;
		
	double *A1 = malloc(n * n * sizeof(double));
	double *A2 = malloc(n * n * sizeof(double));
	double *B1 = malloc(n * n * sizeof(double));
	double *B2 = malloc(n * n * sizeof(double));
	double *C1 = malloc(n * n * sizeof(double));
	double *C2 = malloc(n * n * sizeof(double));

    int info;
    double scale1, scale2;

    // 0, 1, -1
    const int i0[] = {0}, i1[] = {1}, im1[] = {-1};
    // 0
    const double d0[] = {0};

    // N N +1 m = n
    {
        const int m = n_max, n = n_max;

        // generate matrices
        d2matgen(n, n, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(n, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, d0, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, d0, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, d0, B1 + 1, &np1);
        BLAS(dscal)(&nm1, d0, B2 + 1, &np1);

        // run
        LARPACK(dtrsyl)("N", "N", i1, &n, &n, A1, &n, B1, &n, C1, &n, &scale1, &info);
        LAPACK(dtrsy2)("N", "N", i1, &n, &n, A2, &n, B2, &n, C2, &n, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        double error = d2vecerr(n * n, C1, C2);
        printf("dtrsyl N N +1 m = n:\t%g\n", error);
    }

    { // N N +1 m < n
        const int m = n_min, n = n_max;
        // generate matrix
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, d0, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, d0, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, d0, B1 + 1, &np1);
        BLAS(dscal)(&nm1, d0, B2 + 1, &np1);

        // run
        LARPACK(dtrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(dtrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = d2vecerr(m * n, C1, C2);
        printf("dtrsyl N N +1 m < n:\t%g\n", error);
    }

    { // N N +1 m > n
        const int m = n_max, n = n_min;
        // generate matrix
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, d0, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, d0, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, d0, B1 + 1, &np1);
        BLAS(dscal)(&nm1, d0, B2 + 1, &np1);

        // run
        LARPACK(dtrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(dtrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = d2vecerr(m * n, C1, C2);
        printf("dtrsyl N N +1 m > n:\t%g\n", error);
    }

    { // C N +1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, d0, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, d0, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, d0, B1 + 1, &np1);
        BLAS(dscal)(&nm1, d0, B2 + 1, &np1);

        // run
        LARPACK(dtrsyl)("C", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(dtrsy2)("C", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = d2vecerr(m * n, C1, C2);
        printf("dtrsyl C N +1 m = n:\t%g\n", error);
    }

    { // N C +1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, d0, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, d0, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, d0, B1 + 1, &np1);
        BLAS(dscal)(&nm1, d0, B2 + 1, &np1);

        // run
        LARPACK(dtrsyl)("N", "C", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(dtrsy2)("N", "C", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = d2vecerr(m * n, C1, C2);
        printf("dtrsyl N C +1 m = n:\t%g\n", error);
    }

    { // C C +1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, d0, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, d0, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, d0, B1 + 1, &np1);
        BLAS(dscal)(&nm1, d0, B2 + 1, &np1);

        // run
        LARPACK(dtrsyl)("C", "C", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(dtrsy2)("C", "C", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = d2vecerr(m * n, C1, C2);
        printf("dtrsyl C C +1 m = n:\t%g\n", error);
    }

    { // N N offdiag 
        const int m = n_max, n = n_max;
        // generate matrix
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, d0, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, d0, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, d0, B1 + 1, &np1);
        BLAS(dscal)(&nm1, d0, B2 + 1, &np1);

        int m1 = (m >= 16) ? ((m + 8) / 16) * 8 : m / 2;
        A1[m1 + m * (m1 - 1)] = A2[m1 + m * (m1 - 1)] = -.5;
        int n1 = (n >= 16) ? ((n + 8) / 16) * 8 : n / 2;
        B1[n1 + n * (n1 - 1)] = B2[n1 + n * (n1 - 1)] = -.3;
        n1 = n / 2;
        B1[n1 + n * (n1 - 1)] = B2[n1 + n * (n1 - 1)] = -.2;

        // run
        LARPACK(dtrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(dtrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = d2vecerr(m * n, C1, C2);
        printf("dtrsyl offdiag: \t%g\n", error);
    }

    { // N N full
        const int m = n_max, n = n_max;
        // generate matrix
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);

        // run
        LARPACK(dtrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(dtrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = d2vecerr(m * n, C1, C2);
        printf("dtrsyl full:\t\t%g\n", error);
    }

    { // N N -1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        d2matgen(m, m, A1, A2);
        d2matgen(n, n, B1, B2);
        d2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(dscal)(&mm1, d0, A1 + 1, &mp1);
        BLAS(dscal)(&mm1, d0, A2 + 1, &mp1);
        BLAS(dscal)(&nm1, d0, B1 + 1, &np1);
        BLAS(dscal)(&nm1, d0, B2 + 1, &np1);

        // scale diagonal of A and B
        const double smi = 1. / m, sni = 1. / n;
        BLAS(dscal)(&m, &smi, A1, &mp1);
        BLAS(dscal)(&m, &smi, A2, &mp1);
        BLAS(dscal)(&n, &sni, B1, &np1);
        BLAS(dscal)(&n, &sni, B2, &np1);

        // run
        LARPACK(dtrsyl)("N", "N", im1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(dtrsy2)("N", "N", im1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(dlascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = d2vecerr(m * n, C1, C2);
        printf("dtrsyl N N -1 m = n:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);

	return 0;
}
