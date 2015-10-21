#ifndef LAPACK_OROG_TRTRI_H
#define LAPACK_OROG_TRTRI_H

#include "../test_config.h"

extern void LAPACK_ORIG(strtri)(const char*, const char *, const int *, float *, const int *, int *);
extern void LAPACK_ORIG(dtrtri)(const char*, const char *, const int *, double *, const int *, int *);
extern void LAPACK_ORIG(ctrtri)(const char*, const char *, const int *, float *, const int *, int *);
extern void LAPACK_ORIG(ztrtri)(const char*, const char *, const int *, double *, const int *, int *);

#endif /* LAPACK_OROG_TRTRI_H */
