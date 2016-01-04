#ifndef RELAPACK_H
#define RELAPACK_H

#include "../config.h"

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
void RELAPACK(csytrf)(const char *, const int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(chetrf)(const char *, const int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(zsytrf)(const char *, const int *, double *, const int *, int *, double *, const int *, int *);
void RELAPACK(zhetrf)(const char *, const int *, double *, const int *, int *, double *, const int *, int *);
void RELAPACK(ssytrf_rook)(const char *, const int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(dsytrf_rook)(const char *, const int *, double *, const int *, int *, double *, const int *, int *);
void RELAPACK(csytrf_rook)(const char *, const int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(chetrf_rook)(const char *, const int *, float *, const int *, int *, float *, const int *, int *);
void RELAPACK(zsytrf_rook)(const char *, const int *, double *, const int *, int *, double *, const int *, int *);
void RELAPACK(zhetrf_rook)(const char *, const int *, double *, const int *, int *, double *, const int *, int *);

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

#endif /*  RELAPACK_H */
