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
		
	float *A1 = malloc(n * n * sizeof(float));
	float *A2 = malloc(n * n * sizeof(float));
	float *B1 = malloc(n * n * sizeof(float));
	float *B2 = malloc(n * n * sizeof(float));
	float *C1 = malloc(n * n * sizeof(float));
	float *C2 = malloc(n * n * sizeof(float));

    // Outputs
    int info;
    float scale1;
    float scale2;

    // Constants
    const float ZERO[]  = {0};
    const int   iZERO[] = {0};
    const int   iONE[]  = {1};
    const int   iMONE[] = {-1};

    { // N N +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(sscal)(&nm1, ZERO, B2 + 1, &np1);

        // run
        RELAPACK(strsyl)("N", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        STRSY2("N", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(slascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = s2vecerr(m * n, C1, C2);
        printf("strsyl N N +1 m = n:\t%g\n", error);
    }

    { // N N +1 m < n
        const int m = n_min, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(sscal)(&nm1, ZERO, B2 + 1, &np1);

        // run
        RELAPACK(strsyl)("N", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        STRSY2("N", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(slascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = s2vecerr(m * n, C1, C2);
        printf("strsyl N N +1 m < n:\t%g\n", error);
    }

    { // N N +1 m > n
        const int m = n_max, n = n_min;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(sscal)(&nm1, ZERO, B2 + 1, &np1);

        // run
        RELAPACK(strsyl)("N", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        STRSY2("N", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(slascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = s2vecerr(m * n, C1, C2);
        printf("strsyl N N +1 m > n:\t%g\n", error);
    }

    { // C N +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(sscal)(&nm1, ZERO, B2 + 1, &np1);

        // run
        RELAPACK(strsyl)("C", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        STRSY2("C", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(slascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = s2vecerr(m * n, C1, C2);
        printf("strsyl C N +1 m = n:\t%g\n", error);
    }

    { // N C +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(sscal)(&nm1, ZERO, B2 + 1, &np1);

        // run
        RELAPACK(strsyl)("N", "C", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        STRSY2("N", "C", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(slascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = s2vecerr(m * n, C1, C2);
        printf("strsyl N C +1 m = n:\t%g\n", error);
    }

    { // C C +1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(sscal)(&nm1, ZERO, B2 + 1, &np1);

        // run
        RELAPACK(strsyl)("C", "C", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        STRSY2("C", "C", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(slascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = s2vecerr(m * n, C1, C2);
        printf("strsyl C C +1 m = n:\t%g\n", error);
    }

    { // N N offdiag
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(sscal)(&nm1, ZERO, B2 + 1, &np1);

        int m1 = (m >= 16) ? ((m + 8) / 16) * 8 : m / 2;
        A1[m1 + m * (m1 - 1)] = A2[m1 + m * (m1 - 1)] = -.5;
        int n1 = (n >= 16) ? ((n + 8) / 16) * 8 : n / 2;
        B1[n1 + n * (n1 - 1)] = B2[n1 + n * (n1 - 1)] = -.3;
        n1 = n / 2;
        B1[n1 + n * (n1 - 1)] = B2[n1 + n * (n1 - 1)] = -.2;

        // run
        RELAPACK(strsyl)("N", "N", iONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        STRSY2("N", "N", iONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(slascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = s2vecerr(m * n, C1, C2);
        printf("strsyl offdiag: \t%g\n", error);
    }

    { // N N -1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, ZERO, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, ZERO, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, ZERO, B1 + 1, &np1);
        BLAS(sscal)(&nm1, ZERO, B2 + 1, &np1);

        // scale diagonal of A
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);

        // run
        RELAPACK(strsyl)("N", "N", iMONE, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        STRSY2("N", "N", iMONE, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1)
            LAPACK(slascl)("G", iZERO, iZERO, &scale1, &scale2, &m, &n, C1, &m, &info);

        // check error
        const double error = s2vecerr(m * n, C1, C2);
        printf("strsyl N N -1 m = n:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);

	return 0;
}
