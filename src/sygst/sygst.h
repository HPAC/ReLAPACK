#ifndef SYGST_H
#define SYGST_H

void ssygst_ril(const int *, float *, const int *, const float *, const int *);
void dsygst_ril(const int *, double *, const int *, const double *, const int *);
void csygst_ril(const int *, float *, const int *, const float *, const int *);
void zsygst_ril(const int *, double *, const int *, const double *, const int *);

void ssygst_riu(const int *, float *, const int *, const float *, const int *);
void dsygst_riu(const int *, double *, const int *, const double *, const int *);
void csygst_riu(const int *, float *, const int *, const float *, const int *);
void zsygst_riu(const int *, double *, const int *, const double *, const int *);

void ssygst_rnl(const int *, float *, const int *, const float *, const int *);
void dsygst_rnl(const int *, double *, const int *, const double *, const int *);
void csygst_rnl(const int *, float *, const int *, const float *, const int *);
void zsygst_rnl(const int *, double *, const int *, const double *, const int *);

void ssygst_rnu(const int *, float *, const int *, const float *, const int *);
void dsygst_rnu(const int *, double *, const int *, const double *, const int *);
void csygst_rnu(const int *, float *, const int *, const float *, const int *);
void zsygst_rnu(const int *, double *, const int *, const double *, const int *);

#endif /* SYGST_H */
