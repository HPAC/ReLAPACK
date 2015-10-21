#ifndef LAPACK_OROG_POTRF_H
#define LAPACK_OROG_POTRF_H

#include "../test_config.h"

extern void LAPACK_ORIG(spotrf)(const char*, const int *, float *, const int *, int *);
extern void LAPACK_ORIG(dpotrf)(const char*, const int *, double *, const int *, int *);
extern void LAPACK_ORIG(cpotrf)(const char*, const int *, float *, const int *, int *);
extern void LAPACK_ORIG(zpotrf)(const char*, const int *, double *, const int *, int *);

#endif /* LAPACK_OROG_POTRF_H */
