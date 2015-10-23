#ifndef LARPACK_H
#define LARPACK_H

#include "../config.h"
#include "lapack.h"
#include "blas.h"
#include "util.h"

extern void LARPACK(slauum)(const char *, const int *, float *, const int *, int *);
extern void LARPACK(dlauum)(const char *, const int *, double *, const int *, int *);
extern void LARPACK(clauum)(const char *, const int *, float *, const int *, int *);
extern void LARPACK(zlauum)(const char *, const int *, double *, const int *, int *);

extern void LARPACK(strtri)(const char *, const char *, const int *, float *, const int *, int *);
extern void LARPACK(dtrtri)(const char *, const char *, const int *, double *, const int *, int *);
extern void LARPACK(ctrtri)(const char *, const char *, const int *, float *, const int *, int *);
extern void LARPACK(ztrtri)(const char *, const char *, const int *, double *, const int *, int *);

extern void LARPACK(spotrf)(const char *, const int *, float *, const int *, int *);
extern void LARPACK(dpotrf)(const char *, const int *, double *, const int *, int *);
extern void LARPACK(cpotrf)(const char *, const int *, float *, const int *, int *);
extern void LARPACK(zpotrf)(const char *, const int *, double *, const int *, int *);

extern void LARPACK(sgetrf)(const int *, const int *, float *, const int *, int *, int *);
extern void LARPACK(dgetrf)(const int *, const int *, double *, const int *, int *, int *);
extern void LARPACK(cgetrf)(const int *, const int *, float *, const int *, int *, int *);
extern void LARPACK(zgetrf)(const int *, const int *, double *, const int *, int *, int *);

extern void LARPACK(ssygst)(const int *, const char *, const int *, float *, const int *, const float *, const int *, int *);
extern void LARPACK(dsygst)(const int *, const char *, const int *, double *, const int *, const double *, const int *, int *);

extern void LARPACK(chegst)(const int *, const char *, const int *, float *, const int *, const float *, const int *, int *);
extern void LARPACK(zhegst)(const int *, const char *, const int *, double *, const int *, const double *, const int *, int *);

#endif /*  LARPACK_H */
