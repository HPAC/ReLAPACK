#ifndef LAPACK_OROG_GETRF_H
#define LAPACK_OROG_GETRF_H

#include "../test_config.h"

extern void LAPACK_ORIG(sgetrf)(const int*, const int *, float *, const int *, int *, int *);
extern void LAPACK_ORIG(dgetrf)(const int*, const int *, double *, const int *, int *, int *);
extern void LAPACK_ORIG(cgetrf)(const int*, const int *, float *, const int *, int *, int *);
extern void LAPACK_ORIG(zgetrf)(const int*, const int *, double *, const int *, int *, int *);

#endif /* LAPACK_OROG_GETRF_H */
