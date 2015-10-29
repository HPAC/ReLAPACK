#include "../../src/larpack.h"
#include "../test_config.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	const int n_max = TEST_N;
    const int n_min = MAX(1, (n_max * 3) / 4);
    const int n = n_max;
		
	float *A1 = malloc(2 * n * n * sizeof(float));
	float *A2 = malloc(2 * n * n * sizeof(float));
	float *B1 = malloc(2 * n * n * sizeof(float));
	float *B2 = malloc(2 * n * n * sizeof(float));
	float *C1 = malloc(2 * n * n * sizeof(float));
	float *C2 = malloc(2 * n * n * sizeof(float));

    int info;
    float scale1, scale2;
    // 0, 1, -1
    const int i0[] = {0}, i1[] = {1}, im1[] = {-1};
    // 1
    const float s1[] = {1};

    { // N N +1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);

        // run
        LARPACK(ctrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ctrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = c2vecerr(m * n, C1, C2);
        printf("ctrsyl N N +1 m = n:\t%g\n", error);
    }

    { // N N +1 m < n
        const int m = n_min, n = n_max;
        // generate matrix
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);

        // run
        LARPACK(ctrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ctrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = c2vecerr(m * n, C1, C2);
        printf("ctrsyl N N +1 m < n:\t%g\n", error);
    }

    { // N N +1 m > n
        const int m = n_max, n = n_min;
        // generate matrix
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);

        // run
        LARPACK(ctrsyl)("N", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ctrsy2)("N", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = c2vecerr(m * n, C1, C2);
        printf("ctrsyl N N +1 m > n:\t%g\n", error);
    }

    { // C N +1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);

        // run
        LARPACK(ctrsyl)("C", "N", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ctrsy2)("C", "N", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = c2vecerr(m * n, C1, C2);
        printf("ctrsyl C N +1 m = n:\t%g\n", error);
    }

    { // N C +1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);

        // run
        LARPACK(ctrsyl)("N", "C", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ctrsy2)("N", "C", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = c2vecerr(m * n, C1, C2);
        printf("ctrsyl N C +1 m = n:\t%g\n", error);
    }

    { // C C +1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);

        // run
        LARPACK(ctrsyl)("C", "C", i1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ctrsy2)("C", "C", i1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        LAPACK(clascl)("G", i0, i0, &scale1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = c2vecerr(m * n, C1, C2);
        printf("ctrsyl C C +1 m = n:\t%g\n", error);
    }

    { // N N -1 m = n
        const int m = n_max, n = n_max;
        // generate matrix
        c2matgen(m, m, A1, A2);
        c2matgen(n, n, B1, B2);
        c2matgen(m, n, C1, C2);

        // scale diagonal of A and B
        const float smi = 1. / m, sni = 1. / n;
        const int imp1 = m + 1, inp1 = n + 1;
        BLAS(cscal)(&m, &smi, A1, &imp1);
        BLAS(cscal)(&m, &smi, A2, &imp1);
        BLAS(cscal)(&n, &sni, B1, &inp1);
        BLAS(cscal)(&n, &sni, B2, &inp1);

        // run
        LARPACK(ctrsyl)("N", "N", im1, &m, &n, A1, &m, B1, &n, C1, &m, &scale1, &info);
        LAPACK(ctrsy2)("N", "N", im1, &m, &n, A2, &m, B2, &n, C2, &m, &scale2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        LAPACK(clascl)("G", i0, i0, s1, &scale1, &m, &n, C2, &m, &info);
        LAPACK(clascl)("G", i0, i0, s1, &scale2, &m, &n, C2, &m, &info);

        // check error
        const double error = c2vecerr(m * n, C1, C2);
        printf("ctrsyl N N -1 m = n:\t%g\n", error);
    }

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(C1);
    free(C2);

	return 0;
}
