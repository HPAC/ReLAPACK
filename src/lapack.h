#ifndef LAPACK_H
#define LAPACK_H

#include "fortran_interface.h"

extern int LAPACK(lsame)(const char *, const char *);
extern int LAPACK(xerbla)(const char *, const int *);

extern void LAPACK(slaswp)(const int *, float *, const int *, const int *, const int *, const int *, const int *);
extern void LAPACK(dlaswp)(const int *, double *, const int *, const int *, const int *, const int *, const int *);
extern void LAPACK(claswp)(const int *, float *, const int *, const int *, const int *, const int *, const int *);
extern void LAPACK(zlaswp)(const int *, double *, const int *, const int *, const int *, const int *, const int *);

extern void LAPACK(slascl)(const char *, const int *, const int *, const float *, const float *, const int *, const int *, float *, const int *, int *);
extern void LAPACK(dlascl)(const char *, const int *, const int *, const double *, const double *, const int *, const int *, double *, const int *, int *);
extern void LAPACK(clascl)(const char *, const int *, const int *, const float *, const float *, const int *, const int *, float *, const int *, int *);
extern void LAPACK(zlascl)(const char *, const int *, const int *, const double *, const double *, const int *, const int *, double *, const int *, int *);

extern void LAPACK(slauu2)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(dlauu2)(const char *, const int *, double *, const int *, int *);
extern void LAPACK(clauu2)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(zlauu2)(const char *, const int *, double *, const int *, int *);

extern void LAPACK(strti2)(const char *, const char *, const int *, float *, const int *, int *);
extern void LAPACK(dtrti2)(const char *, const char *, const int *, double *, const int *, int *);
extern void LAPACK(ctrti2)(const char *, const char *, const int *, float *, const int *, int *);
extern void LAPACK(ztrti2)(const char *, const char *, const int *, double *, const int *, int *);

extern void LAPACK(spotf2)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(dpotf2)(const char *, const int *, double *, const int *, int *);
extern void LAPACK(cpotf2)(const char *, const int *, float *, const int *, int *);
extern void LAPACK(zpotf2)(const char *, const int *, double *, const int *, int *);

extern void LAPACK(sgetf2)(const int *, const int *, float *, const int *, int *, int *);
extern void LAPACK(dgetf2)(const int *, const int *, double *, const int *, int *, int *);
extern void LAPACK(cgetf2)(const int *, const int *, float *, const int *, int *, int *);
extern void LAPACK(zgetf2)(const int *, const int *, double *, const int *, int *, int *);

extern void LAPACK(ssygs2)(const int *, const char *, const int *, float *, const int *, const float *, const int *, int *);
extern void LAPACK(dsygs2)(const int *, const char *, const int *, double *, const int *, const double *, const int *, int *);
extern void LAPACK(chegs2)(const int *, const char *, const int *, float *, const int *, const float *, const int *, int *);
extern void LAPACK(zhegs2)(const int *, const char *, const int *, double *, const int *, const double *, const int *, int *);

#endif /* LAPACK_H */
