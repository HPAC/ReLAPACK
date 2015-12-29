#include "util.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

void s2matgen(const int m, const int n, float *A, float *B) {
    srand(time(NULL) + (size_t) A);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            A[i + m * j] = B[i + m * j] = (float) rand() / RAND_MAX + m * (i == j);
}

void d2matgen(const int m, const int n, double *A, double *B) {
    srand(time(NULL) + (size_t) A);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            A[i + m * j] = B[i + m * j] = (double) rand() / RAND_MAX + m * (i == j);
}

void c2matgen(const int m, const int n, float *A, float *B) {
    srand(time(NULL) + (size_t) A);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            A[2* (i + m * j)]     = B[2 * (i + m * j)]     = (float) rand() / RAND_MAX + m * (i == j);
            A[2* (i + m * j) + 1] = B[2 * (i + m * j) + 1] = ((float) rand() / RAND_MAX) * (i != j);
        }
}

void z2matgen(const int m, const int n, double *A, double *B) {
    srand(time(NULL) + (size_t) A);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            A[2* (i + m * j)]     = B[2 * (i + m * j)]     = (double) rand() / RAND_MAX + m * (i == j);
            A[2* (i + m * j) + 1] = B[2 * (i + m * j) + 1] = ((double) rand() / RAND_MAX) * (i != j);
        }
}

double i2vecerr(const int n, const int *x, const int *y) {
    double error = 0;
    for (int i = 0; i < n; i++) {
        double nom = abs(x[i] - y[i]);
        double den = abs(y[i]);
        error = MAX(error, (den > 0) ? MIN(nom, nom / den) : nom);
    }
    return error;
}

double s2vecerr(const int n, const float *x, const float *y) {
    float error = 0;
    for (int i = 0; i < n; i++) {
        double nom = fabs((double) x[i] - y[i]);
        double den = fabs(y[i]);
        error = MAX(error, (den > 0) ? MIN(nom, nom / den) : nom);
    }
    return error;
}

double d2vecerr(const int n, const double *x, const double *y) {
    double error = 0;
    for (int i = 0; i < n; i++) {
        double nom = fabs(x[i] - y[i]);
        double den = fabs(y[i]);
        error = MAX(error, (den > 0) ? MIN(nom, nom / den) : nom);
    }
    return error;
}

double c2vecerr(const int n, const float *x, const float *y) {
    double error = 0;
    for (int i = 0; i < n; i++) {
        double nom = sqrt(((double) x[2 * i] - y[2 * i]) * ((double) x[2 * i] - y[2 * i]) + ((double) x[2 * i + 1] - y[2 * i + 1]) * ((double) x[2 * i + 1] - y[2 * i + 1]));
        double den = sqrt((double) y[2 * i] * y[2 * i] + (double) y[2 * i + 1] * y[2 * i + 1]);
        error = MAX(error, (den > 0) ? MIN(nom, nom / den) : nom);
    }
    return error;
}

double z2vecerr(const int n, const double *x, const double *y) {
    double error = 0;
    for (int i = 0; i < n; i++) {
        double nom = sqrt((x[2 * i] - y[2 * i]) * (x[2 * i] - y[2 * i]) + (x[2 * i + 1] - y[2 * i + 1]) * (x[2 * i + 1] - y[2 * i + 1]));
        double den = sqrt(y[2 * i] * y[2 * i] + y[2 * i + 1] * y[2 * i + 1]);
        error = MAX(error, (den > 0) ? MIN(nom, nom / den) : nom);
    }
    return error;
}
