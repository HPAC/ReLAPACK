#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include "../src/util.h"

void s2matgen(int, int, float *, float *);
void d2matgen(int, int, double *, double *);
void c2matgen(int, int, float *, float *);
void z2matgen(int, int, double *, double *);

int i2vecerr(int, const int *, const int *);
float s2vecerr(int, const float *, const float *);
double d2vecerr(int, const double *, const double *);
float c2vecerr(int, const float *, const float *);
double z2vecerr(int, const double *, const double *);

#endif /* TEST_UTIL_H */
