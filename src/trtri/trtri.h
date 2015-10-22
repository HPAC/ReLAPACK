#ifndef TRTRI_H
#define TRTRI_H

void strtri_rl(const char *, const int *, float *, const int *);
void dtrtri_rl(const char *, const int *, double *, const int *);
void ctrtri_rl(const char *, const int *, float *, const int *);
void ztrtri_rl(const char *, const int *, double *, const int *);

void strtri_ru(const char *, const int *, float *, const int *);
void dtrtri_ru(const char *, const int *, double *, const int *);
void ctrtri_ru(const char *, const int *, float *, const int *);
void ztrtri_ru(const char *, const int *, double *, const int *);

#endif /* TRTRI_H */
