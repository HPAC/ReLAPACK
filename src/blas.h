#ifndef BLAS_H
#define BLAS_H

#include "fortran_interface.h"

extern void BLAS(sgemm)(const char *, const char *, const int *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, const float *, const int*);
extern void BLAS(dgemm)(const char *, const char *, const int *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, const double *, const int*);
extern void BLAS(cgemm)(const char *, const char *, const int *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, const float *, const int*);
extern void BLAS(zgemm)(const char *, const char *, const int *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, const double *, const int*);

extern void BLAS(strsm)(const char *, const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, float *, const int *);
extern void BLAS(dtrsm)(const char *, const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, double *, const int *);
extern void BLAS(ctrsm)(const char *, const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, float *, const int *);
extern void BLAS(ztrsm)(const char *, const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, double *, const int *);

extern void BLAS(strmm)(const char *, const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, float *, const int *);
extern void BLAS(dtrmm)(const char *, const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, double *, const int *);
extern void BLAS(ctrmm)(const char *, const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, float *, const int *);
extern void BLAS(ztrmm)(const char *, const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, double *, const int *);

extern void BLAS(ssyrk)(const char *, const char *, const int *, const int *, const float *, float *, const int *, const float *, float *, const int *);
extern void BLAS(dsyrk)(const char *, const char *, const int *, const int *, const double *, double *, const int *, const double *, double *, const int *);

extern void BLAS(cherk)(const char *, const char *, const int *, const int *, const float *, float *, const int *, const float *, float *, const int *);
extern void BLAS(zherk)(const char *, const char *, const int *, const int *, const double *, double *, const int *, const double *, double *, const int *);

extern void BLAS(ssymm)(const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
extern void BLAS(dsymm)(const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
extern void BLAS(csymm)(const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
extern void BLAS(zsymm)(const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);

extern void BLAS(ssyr2k)(const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
extern void BLAS(dsyr2k)(const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
extern void BLAS(csyr2k)(const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
extern void BLAS(zsyr2k)(const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);

#endif /* BLAS_H */
