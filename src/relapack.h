#ifndef RELAPACK_H
#define RELAPACK_H

#include "../config.h"

#ifdef BLAS_UNDERSCORE
#define BLAS(routine) routine ## _
#else
#define BLAS(routine) routine
#endif

#ifdef LAPACK_UNDERSCORE
#define LAPACK(routine) routine ## _
#else
#define LAPACK(routine) routine
#endif

#ifdef RELAPACK_AS_LAPACK
#define RELEPACK(routine) LAPACK(routine)
#else
#define RELAPACK(routine) RELAPACK_ ## routine
#endif

#include "lapack.h"
#include "blas.h"
#include "util.h"

void RELAPACK(slauum)(const char *, const int *, float *, const int *, int *);
void RELAPACK(dlauum)(const char *, const int *, double *, const int *, int *);
void RELAPACK(clauum)(const char *, const int *, float *, const int *, int *);
void RELAPACK(zlauum)(const char *, const int *, double *, const int *, int *);

void RELAPACK(strtri)(const char *, const char *, const int *, float *, const int *, int *);
void RELAPACK(dtrtri)(const char *, const char *, const int *, double *, const int *, int *);
void RELAPACK(ctrtri)(const char *, const char *, const int *, float *, const int *, int *);
void RELAPACK(ztrtri)(const char *, const char *, const int *, double *, const int *, int *);

void RELAPACK(spotrf)(const char *, const int *, float *, const int *, int *);
void RELAPACK(dpotrf)(const char *, const int *, double *, const int *, int *);
void RELAPACK(cpotrf)(const char *, const int *, float *, const int *, int *);
void RELAPACK(zpotrf)(const char *, const int *, double *, const int *, int *);

void RELAPACK(ssytrf)(const char *, const int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(dsytrf)(const char *, const int *, double *, const int *, int *, double *, const int *, int *);
void RELAPACK(chetrf)(const char *, const int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(zhetrf)(const char *, const int *, double *, const int *, int *, double *, const int *, int *);
void RELAPACK(csytrf)(const char *, const int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(zsytrf)(const char *, const int *, double *, const int *, int *, double *, const int *, int *);

void RELAPACK(sgetrf)(const int *, const int *, float *, const int *, int *, int *);
void RELAPACK(dgetrf)(const int *, const int *, double *, const int *, int *, int *);
void RELAPACK(cgetrf)(const int *, const int *, float *, const int *, int *, int *);
void RELAPACK(zgetrf)(const int *, const int *, double *, const int *, int *, int *);

void RELAPACK(ssygst)(const int *, const char *, const int *, float *, const int *, const float *, const int *, int *);
void RELAPACK(dsygst)(const int *, const char *, const int *, double *, const int *, const double *, const int *, int *);
void RELAPACK(chegst)(const int *, const char *, const int *, float *, const int *, const float *, const int *, int *);
void RELAPACK(zhegst)(const int *, const char *, const int *, double *, const int *, const double *, const int *, int *);

void RELAPACK(strsyl)(const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, int *);
void RELAPACK(dtrsyl)(const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, int *);
void RELAPACK(ctrsyl)(const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, int *);
void RELAPACK(ztrsyl)(const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, int *);

void RELAPACK(stgsyl)(const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, float *, float *, const int *, int *, int *); 
void RELAPACK(dtgsyl)(const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, double *, double *, const int *, int *, int *); 
void RELAPACK(ctgsyl)(const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, float *, float *, const int *, int *, int *); 
void RELAPACK(ztgsyl)(const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, double *, double *, const int *, int *, int *); 

// sytrf helper routines
void RELAPACK(sgemm_tr)(const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
void RELAPACK(dgemm_tr)(const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
void RELAPACK(cgemm_tr)(const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
void RELAPACK(zgemm_tr)(const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);

void RELAPACK(sgemm_tr2)(const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
void RELAPACK(dgemm_tr2)(const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
void RELAPACK(cgemm_tr2)(const char *, const char *, const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *);
void RELAPACK(zgemm_tr2)(const char *, const char *, const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);

void RELAPACK(ssytrf_rec)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(dsytrf_rec)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void RELAPACK(csytrf_rec)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(chetrf_rec)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(zsytrf_rec)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void RELAPACK(zhetrf_rec)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);

void LAPACK(ssytrf_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(dsytrf_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void LAPACK(csytrf_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(chetrf_rec2)(const char *, const int *, const int *, int *, float *, const int *, int *, float *, const int *, int *);
void LAPACK(zsytrf_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);
void LAPACK(zhetrf_rec2)(const char *, const int *, const int *, int *, double *, const int *, int *, double *, const int *, int *);

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

// tgsyl helper routines
void RELAPACK(stgsyl_rec)(const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, float *, float *, int *, int *, int *);
void RELAPACK(dtgsyl_rec)(const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, double *, double *, int *, int *, int *);
void RELAPACK(ctgsyl_rec)(const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, float *, float *, int *);
void RELAPACK(ztgsyl_rec)(const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, double *, double *, int *);

#endif /*  RELAPACK_H */
