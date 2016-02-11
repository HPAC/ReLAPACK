#ifndef RELAPACK_INT_H
#define RELAPACK_INT_H

#include "../inc/relapack.h"

#if BLAS_UNDERSCORE
#define BLAS(routine) routine ## _
#else
#define BLAS(routine) routine
#endif

#include "lapack.h"
#include "blas.h"
#include "util.h"


// sytrf helper routines
void RELAPACK(sgemm_tr_rec)(const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
void RELAPACK(dgemm_tr_rec)(const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
void RELAPACK(cgemm_tr_rec)(const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
void RELAPACK(zgemm_tr_rec)(const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);

void LAPACK(ssytrf_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(dsytrf_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void LAPACK(csytrf_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(chetrf_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(zsytrf_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void LAPACK(zhetrf_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void LAPACK(ssytrf_rook_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(dsytrf_rook_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void LAPACK(csytrf_rook_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(chetrf_rook_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(zsytrf_rook_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void LAPACK(zhetrf_rook_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);

// trsyl helper routines
#ifdef RELAPACK_AS_LAPACK
#define STRSY2 LAPACK(strsy2)
#define DTRSY2 LAPACK(dtrsy2)
#define CTRSY2 LAPACK(ctrsy2)
#define ZTRSY2 LAPACK(ztrsy2)
#else
#define STRSY2 LAPACK(strsyl)
#define DTRSY2 LAPACK(dtrsyl)
#define CTRSY2 LAPACK(ctrsyl)
#define ZTRSY2 LAPACK(ztrsyl)
#endif
void STRSY2(const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, int *);
void DTRSY2(const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, int *);
void CTRSY2(const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, int *);
void ZTRSY2(const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, int *);

#endif /*  RELAPACK_INT_H */
