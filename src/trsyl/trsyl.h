#ifndef TRSYL_H
#define TRSYL_H

void strsyl_rnn(const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *);
void dtrsyl_rnn(const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *);
void ctrsyl_rnn(const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *);
void ztrsyl_rnn(const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *);

void strsyl_rnt(const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *);
void dtrsyl_rnt(const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *);
void ctrsyl_rnt(const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *);
void ztrsyl_rnt(const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *);

void strsyl_rtn(const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *);
void dtrsyl_rtn(const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *);
void ctrsyl_rtn(const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *);
void ztrsyl_rtn(const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *);

void strsyl_rtt(const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *);
void dtrsyl_rtt(const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *);
void ctrsyl_rtt(const int *, const int *, const int *, const float *, const int *, const float *, const int *, float *, const int *, const float *);
void ztrsyl_rtt(const int *, const int *, const int *, const double *, const int *, const double *, const int *, double *, const int *, const double *);

#endif /* TRSYL_H */
