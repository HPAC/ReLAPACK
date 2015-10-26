#ifndef LARPACK_H
#define LARPACK_H

#include "../config.h"
#include "lapack.h"
#include "blas.h"
#include "util.h"

void LARPACK(slauum)(const char *, const int *, float *, const int *, int *);
void LARPACK(dlauum)(const char *, const int *, double *, const int *, int *);
void LARPACK(clauum)(const char *, const int *, float *, const int *, int *);
void LARPACK(zlauum)(const char *, const int *, double *, const int *, int *);

void LARPACK(strtri)(const char *, const char *, const int *, float *, const int *, int *);
void LARPACK(dtrtri)(const char *, const char *, const int *, double *, const int *, int *);
void LARPACK(ctrtri)(const char *, const char *, const int *, float *, const int *, int *);
void LARPACK(ztrtri)(const char *, const char *, const int *, double *, const int *, int *);

void LARPACK(spotrf)(const char *, const int *, float *, const int *, int *);
void LARPACK(dpotrf)(const char *, const int *, double *, const int *, int *);
void LARPACK(cpotrf)(const char *, const int *, float *, const int *, int *);
void LARPACK(zpotrf)(const char *, const int *, double *, const int *, int *);

void LARPACK(sgetrf)(const int *, const int *, float *, const int *, int *, int *);
void LARPACK(dgetrf)(const int *, const int *, double *, const int *, int *, int *);
void LARPACK(cgetrf)(const int *, const int *, float *, const int *, int *, int *);
void LARPACK(zgetrf)(const int *, const int *, double *, const int *, int *, int *);

void LARPACK(ssygst)(const int *, const char *, const int *, float *, const int *, const float *, const int *, int *);
void LARPACK(dsygst)(const int *, const char *, const int *, double *, const int *, const double *, const int *, int *);

void LARPACK(chegst)(const int *, const char *, const int *, float *, const int *, const float *, const int *, int *);
void LARPACK(zhegst)(const int *, const char *, const int *, double *, const int *, const double *, const int *, int *);

void LAPACK(strsy2)(const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, int *);
void LAPACK(dtrsy2)(const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, int *);
void LAPACK(ctrsy2)(const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, int *);
void LAPACK(ztrsy2)(const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, int *);

void LARPACK(strsyl)(const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, int *);
void LARPACK(dtrsyl)(const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, int *);
void LARPACK(ctrsyl)(const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, float *, int *);
void LARPACK(ztrsyl)(const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, double *, int *);

#endif /*  LARPACK_H */
