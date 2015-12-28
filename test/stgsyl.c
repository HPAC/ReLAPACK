#include "../src/relapack.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

void LAPACK(stgsyl)(const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, float *, float *, const int *, int *, int *);

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
	float *D1 = malloc(n * n * sizeof(float));
	float *D2 = malloc(n * n * sizeof(float));
	float *E1 = malloc(n * n * sizeof(float));
	float *E2 = malloc(n * n * sizeof(float));
	float *F1 = malloc(n * n * sizeof(float));
	float *F2 = malloc(n * n * sizeof(float));
    const int lWork = 2 * n * n;
    float *Work1 = malloc(lWork  * sizeof(float));
    float *Work2 = malloc(lWork  * sizeof(float));
    int *iWork1 = malloc((n + n + 2) * sizeof(int));
    int *iWork2 = malloc((n + n + 2) * sizeof(int));

    int info;
    float scale1, scale2;
    float dif1, dif2;

    // 0, 1, 2, 3, 4
    const int i0[] = {0}, i1[] = {1}, i2[] = {2}, i3[] = {3}, i4[] = {4};
    // 0
    const float s0[] = {0};

    { // N 0 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        // scale diagonals of A and D
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);
        BLAS(sscal)(&m, &smi, D1, &mp1);
        BLAS(sscal)(&m, &smi, D2, &mp1);

        // run
        RELAPACK(stgsyl)("N", i0, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("N", i0, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2);;
        printf("stgsyl N 0 m = n:\t%g\n", error);
    }

    { // N 0 m < n
        const int m = n_min, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        // scale diagonals of A and D
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);
        BLAS(sscal)(&m, &smi, D1, &mp1);
        BLAS(sscal)(&m, &smi, D2, &mp1);

        // run
        RELAPACK(stgsyl)("N", i0, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("N", i0, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2);;
        printf("stgsyl N 0 m < n:\t%g\n", error);
    }

    { // N 0 m > n
        const int m = n_max, n = n_min;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        // scale diagonals of A and D
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);
        BLAS(sscal)(&m, &smi, D1, &mp1);
        BLAS(sscal)(&m, &smi, D2, &mp1);

        // run
        RELAPACK(stgsyl)("N", i0, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("N", i0, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2);;
        printf("stgsyl N 0 m > n:\t%g\n", error);
    }

    { // N 1 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        // scale diagonals of A and D
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);
        BLAS(sscal)(&m, &smi, D1, &mp1);
        BLAS(sscal)(&m, &smi, D2, &mp1);

        // run
        RELAPACK(stgsyl)("N", i1, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("N", i1, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("stgsyl N 1 m = n:\t%g\n", error);
    }

    { // N 2 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        // scale diagonals of A and D
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);
        BLAS(sscal)(&m, &smi, D1, &mp1);
        BLAS(sscal)(&m, &smi, D2, &mp1);

        // run
        RELAPACK(stgsyl)("N", i2, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("N", i2, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("stgsyl N 2 m = n:\t%g\n", error);
    }

    { // N 3 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        // scale diagonals of A and D
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);
        BLAS(sscal)(&m, &smi, D1, &mp1);
        BLAS(sscal)(&m, &smi, D2, &mp1);

        // run
        RELAPACK(stgsyl)("N", i3, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("N", i3, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("stgsyl N 3 m = n:\t%g\n", error);
    }

    { // N 4 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        // scale diagonals of A and D
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);
        BLAS(sscal)(&m, &smi, D1, &mp1);
        BLAS(sscal)(&m, &smi, D2, &mp1);

        // run
        RELAPACK(stgsyl)("N", i4, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("N", i4, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("stgsyl N 4 m = n:\t%g\n", error);
    }

    { // N 0 offdiag
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        int m1 = (m >= 16) ? ((m + 8) / 16) * 8 : m / 2;
        A1[m1 + m * (m1 - 1)] = A2[m1 + m * (m1 - 1)] = -.5;
        int n1 = (n >= 16) ? ((n + 8) / 16) * 8 : n / 2;
        B1[n1 + n * (n1 - 1)] = B2[n1 + n * (n1 - 1)] = -.3;
        n1 = n / 2;
        B1[n1 + n * (n1 - 1)] = B2[n1 + n * (n1 - 1)] = -.2;

        // scale diagonals of A and D
        const float smi = 1. / m;
        BLAS(sscal)(&m, &smi, A1, &mp1);
        BLAS(sscal)(&m, &smi, A2, &mp1);
        BLAS(sscal)(&m, &smi, D1, &mp1);
        BLAS(sscal)(&m, &smi, D2, &mp1);

        // run
        RELAPACK(stgsyl)("N", i4, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("N", i4, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2) + s2vecerr(1, &dif1, &dif2);
        printf("stgsyl N 0 offdiag:\t%g\n", error);
    }

    { // T 0 m = n
        const int m = n_max, n = n_max;
        // generate matrices
        s2matgen(m, m, A1, A2);
        s2matgen(n, n, B1, B2);
        s2matgen(m, n, C1, C2);
        s2matgen(m, m, D1, D2);
        s2matgen(n, n, E1, E2);
        s2matgen(m, n, F1, F2);

        // clear subdiagonal for A and B
        const int mm1 = m - 1, mp1 = m + 1;
        const int nm1 = n - 1, np1 = n + 1;
        BLAS(sscal)(&mm1, s0, A1 + 1, &mp1);
        BLAS(sscal)(&mm1, s0, A2 + 1, &mp1);
        BLAS(sscal)(&nm1, s0, B1 + 1, &np1);
        BLAS(sscal)(&nm1, s0, B2 + 1, &np1);

        // run
        RELAPACK(stgsyl)("T", i0, &m, &n, A1, &m, B1, &n, C1, &m, D1, &m, E1, &n, F1, &m, &scale1, &dif1, Work1, &lWork, iWork1, &info);
        LAPACK(stgsyl)("T", i0, &m, &n, A2, &m, B2, &n, C2, &m, D2, &m, E2, &n, F2, &m, &scale2, &dif2, Work2, &lWork, iWork2, &info);
        if (scale1 != 1 || scale2 != 1)
            printf("scale1 = %12g\tscale2 = %12g\n", scale1, scale2);

        // apply scales
        if (scale1) {
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
            LAPACK(slascl)("G", i0, i0, &scale1, &scale2, &m, &n, C1, &m, &info);
        }

        // check error
        const double error = s2vecerr(m * n, C1, C2) + s2vecerr(m * n, F1, F2);
        printf("stgsyl T 0 m = n:\t%g\n", error);
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
