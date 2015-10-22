#ifndef LAPACK_H
#define LAPACK_H

#include "fortran_interface.h"

extern int LAPACK(lsame)(const char *, const char *);
extern int LAPACK(xerbla)(const char *, const int *);

extern void LAPACK(slaswp)(const int *, float *, const int *, const int *, const int *, const int *, const int *);
extern void LAPACK(dlaswp)(const int *, double *, const int *, const int *, const int *, const int *, const int *);
extern void LAPACK(claswp)(const int *, float *, const int *, const int *, const int *, const int *, const int *);
extern void LAPACK(zlaswp)(const int *, double *, const int *, const int *, const int *, const int *, const int *);

extern void LAPACK(slauum)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(dlauum)(const char *, const int *, double *, const int *, int *);
extern void LAPACK(clauum)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(zlauum)(const char *, const int *, double *, const int *, int *);

extern void LAPACK(slauu2)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(dlauu2)(const char *, const int *, double *, const int *, int *);
extern void LAPACK(clauu2)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(zlauu2)(const char *, const int *, double *, const int *, int *);

extern void LAPACK(strtri)(const char *, const char *, const int *, float *, const int *, int *);
extern void LAPACK(dtrtri)(const char *, const char *, const int *, double *, const int *, int *);
extern void LAPACK(ctrtri)(const char *, const char *, const int *, float *, const int *, int *);
extern void LAPACK(ztrtri)(const char *, const char *, const int *, double *, const int *, int *);

extern void LAPACK(strti2)(const char *, const char *, const int *, float *, const int *, int *);
extern void LAPACK(dtrti2)(const char *, const char *, const int *, double *, const int *, int *);
extern void LAPACK(ctrti2)(const char *, const char *, const int *, float *, const int *, int *);
extern void LAPACK(ztrti2)(const char *, const char *, const int *, double *, const int *, int *);

extern void LAPACK(spotrf)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(dpotrf)(const char *, const int *, double *, const int *, int *);
extern void LAPACK(cpotrf)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(zpotrf)(const char *, const int *, double *, const int *, int *);

extern void LAPACK(spotf2)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(dpotf2)(const char *, const int *, double *, const int *, int *);
extern void LAPACK(cpotf2)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(zpotf2)(const char *, const int *, double *, const int *, int *);

extern void LAPACK(sgetrf)(const int *, const int *, float *, const int *, int *, int *);
extern void LAPACK(dgetrf)(const int *, const int *, double *, const int *, int *, int *);
extern void LAPACK(cgetrf)(const int *, const int *, float *, const int *, int *, int *);
extern void LAPACK(zgetrf)(const int *, const int *, double *, const int *, int *, int *);

extern void LAPACK(sgetf2)(const int *, const int *, float *, const int *, int *, int *);
extern void LAPACK(dgetf2)(const int *, const int *, double *, const int *, int *, int *);
extern void LAPACK(cgetf2)(const int *, const int *, float *, const int *, int *, int *);
extern void LAPACK(zgetf2)(const int *, const int *, double *, const int *, int *, int *);

#endif /* LAPACK_H */
