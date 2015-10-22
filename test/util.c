#include "util.h"

#include <stdlib.h>
#include <time.h>
#include <math.h>

void s2matgen(int m, int n, float *A, float *B) {
    srand(time(NULL));
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            A[i + m * j] = B[i + m * j] = (float) rand() / RAND_MAX + m * (i == j);
}

void d2matgen(int m, int n, double *A, double *B) {
    srand(time(NULL));
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            A[i + m * j] = B[i + m * j] = (double) rand() / RAND_MAX + m * (i == j);
}

void c2matgen(int m, int n, float *A, float *B) {
    srand(time(NULL));
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++) {
            A[2* (i + m * j)] = B[2 * (i + m * j)] = (float) rand() / RAND_MAX + m * (i == j);
            A[2* (i + m * j) + 1] = B[2 * (i + m * j) + 1] = (float) rand() / RAND_MAX;
        }
}

void z2matgen(int m, int n, double *A, double *B) {
    srand(time(NULL));
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++) {
            A[2* (i + m * j)] = B[2 * (i + m * j)] = (double) rand() / RAND_MAX + m * (i == j);
            A[2* (i + m * j) + 1] = B[2 * (i + m * j) + 1] = (double) rand() / RAND_MAX;
        }
}

int i2vecerr(int n, const int *x, const int *y) {
    int i;
    int error = 0;
    for (i = 0; i < n; i++)
        error += (x[i] - y[i]) * (x[i] - y[i]);
    return (int) sqrt((double) error / n);
}

float s2vecerr(int n, const float *x, const float *y) {
    int i;
    float error = 0;
    for (i = 0; i < n; i++)
        error += (x[i] - y[i]) * (x[i] - y[i]);
    return sqrtf(error / n);
}

double d2vecerr(int n, const double *x, const double *y) {
    int i;
    double error = 0;
    for (i = 0; i < n; i++)
        error += (x[i] - y[i]) * (x[i] - y[i]);
    return sqrt(error / n);
}

float c2vecerr(int n, const float *x, const float *y) {
    int i;
    float error = 0;
    for (i = 0; i < n; i++) {
        error += (x[2 * i] - y[2 * i]) * (x[2 * i] - y[2 * i]);
        error += (x[2 * i + 1] - y[2 * i + 1]) * (x[2 * i + 1] - y[2 * i + 1]);
    }
    return sqrtf(error / n);
}

double z2vecerr(int n, const double *x, const double *y) {
    int i;
    double error = 0;
    for (i = 0; i < n; i++) {
        error += (x[2 * i] - y[2 * i]) * (x[2 * i] - y[2 * i]);
        error += (x[2 * i + 1] - y[2 * i + 1]) * (x[2 * i + 1] - y[2 * i + 1]);
    }
    return sqrtf(error / n);
}
