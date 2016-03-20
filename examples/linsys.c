#include "stdlib.h"
#include "time.h"
#include "stdio.h"
#include "../inc/relapack.h"

// some BLAS and LAPACK prototypes
void dgemv_(const char *, const int *, const int *, const double *, const double
        *, const int *, const double *, const int *, const double *, double *,
        const int *);

void dlaswp_(const int *, double *, const int *, const int *, const int *, const
        int *, const int *);

void dtrsv_(const char *, const char *, const char *, const int *, const double
        *, const int *, const double *, const int *);

void daxpy_(const int *, const double *, const double *, const int *, double *,
        const int *);

double dnrm2_(const int *, const double *, const int *);

const int n = 400;

// This example program generates a random n x n matrix A and a vector x and
// computes y := A x.  It then retrieves x from A and y by solving the linear
// system.  The retireved vector (called x2) is then compared to the original x.
int main(int argc, char *argv[]) {
    // some constants
    const double ONE = 1;
    const double MONE = -1;
    const int iONE = 1;

    // iterators
    int i, j;

    // initialize the random number generator
    srand(time(NULL));

    // generate a random n x n matrix A
    double *A = malloc(n * n * sizeof(double));
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            A[i + j * n] = (double) rand() / RAND_MAX;

    // generate a random vector x of size n
    double *x = malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
        x[i] = (double) rand() / RAND_MAX;

    // compute y := A x
    double *y = malloc(n * sizeof(double));
    dgemv_("N", &n, &n, &ONE, A, &n, x, &iONE, &ONE, y, &iONE);

    // P L U := A (pivoting LU decomposition)
    int *ipiv = malloc(n * sizeof(int)); // pivoting vector representing P
    int info;
    RELAPACK_dgetrf(&n, &n, A, &n, ipiv, &info);

    // compute y := P^-1 y
    dlaswp_(&iONE, y, &n, &iONE, &n, ipiv, &iONE);

    // compute y := L^-1 y
    dtrsv_("L", "N", "U", &n, A, &n, y, &iONE);

    // compute y := U^-1 y
    dtrsv_("U", "N", "N", &n, A, &n, y, &iONE);

    // compute x := -x + y
    daxpy_(&n, &MONE, x, &iONE, y, &iONE);

    // compute norm = ||x||_2
    double norm = dnrm2_(&n, y, &iONE);

    printf("error 2-norm: %f\n", norm);

    // free the matrix and vectors
    free(A);
    free(x);
    free(y);
    free(ipiv);
}
