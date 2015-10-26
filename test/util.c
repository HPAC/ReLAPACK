#include "util.h"

#include <stdlib.h>
#include <time.h>
#include <math.h>

void s2matgen(int m, int n, float *A, float *B) {
    srand(time(NULL) + (size_t) A);
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            A[i + m * j] = B[i + m * j] = (float) rand() / RAND_MAX + m * (i == j);
}

void d2matgen(int m, int n, double *A, double *B) {
    srand(time(NULL) + (size_t) A);
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            A[i + m * j] = B[i + m * j] = (double) rand() / RAND_MAX + m * (i == j);
}

void c2matgen(int m, int n, float *A, float *B) {
    srand(time(NULL) + (size_t) A);
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++) {
            A[2* (i + m * j)] = B[2 * (i + m * j)] = (float) rand() / RAND_MAX + m * (i == j);
            A[2* (i + m * j) + 1] = B[2 * (i + m * j) + 1] = ((float) rand() / RAND_MAX) * (i != j);
        }
}

void z2matgen(int m, int n, double *A, double *B) {
    srand(time(NULL) + (size_t) A);
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++) {
            A[2* (i + m * j)] = B[2 * (i + m * j)] = (double) rand() / RAND_MAX + m * (i == j);
            A[2* (i + m * j) + 1] = B[2 * (i + m * j) + 1] = ((double) rand() / RAND_MAX) * (i != j);
        }
}

int i2vecerr(int n, const int *x, const int *y) {
    int i;
    int error = 0;
    for (i = 0; i < n; i++) {
        int nom = abs(x[i] - y[i]);
        int den = abs(y[i]);
        error = MAX(error, MIN(nom, nom / (float) den));
    }
    return error;
}

float s2vecerr(int n, const float *x, const float *y) {
    int i;
    float error = 0;
    for (i = 0; i < n; i++) {
        float nom = fabsf(x[i] - y[i]);
        float den = fabsf(y[i]);
        error = MAX(error, MIN(nom, nom / den));
    }
    return error;
}

double d2vecerr(int n, const double *x, const double *y) {
    int i;
    double error = 0;
    for (i = 0; i < n; i++) {
        double nom = fabs(x[i] - y[i]);
        double den = fabs(y[i]);
        error = MAX(error, MIN(nom, nom / den));
    }
    return error;
}

float c2vecerr(int n, const float *x, const float *y) {
    int i;
    float error = 0;
    for (i = 0; i < n; i++) {
        float nom = sqrtf((x[2 * i] - y[2 * i]) * (x[2 * i] - y[2 * i]) + (x[2 * i + 1] - y[2 * i + 1]) * (x[2 * i + 1] - y[2 * i + 1]));
        float den = sqrtf(y[2 * i] * y[2 * i] + y[2 * i + 1] * y[2 * i + 1]);
        error = MAX(error, MIN(nom, nom / den));
    }
    return error;
}

double z2vecerr(int n, const double *x, const double *y) {
    int i;
    double error = 0;
    for (i = 0; i < n; i++) {
        double nom = sqrt((x[2 * i] - y[2 * i]) * (x[2 * i] - y[2 * i]) + (x[2 * i + 1] - y[2 * i + 1]) * (x[2 * i + 1] - y[2 * i + 1]));
        double den = sqrt(y[2 * i] * y[2 * i] + y[2 * i + 1] * y[2 * i + 1]);
        error = MAX(error, MIN(nom, nom / den));
    }
    return error;
}
