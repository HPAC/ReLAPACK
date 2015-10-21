#ifndef POTRF_H
#define POTRF_H

void spotrf_rl(const int *, float *, const int *, int *);
void dpotrf_rl(const int *, double *, const int *, int *);
void cpotrf_rl(const int *, float *, const int *, int *);
void zpotrf_rl(const int *, double *, const int *, int *);

void dpotrf_ru(const int *, double *, const int *, int *);
void spotrf_ru(const int *, float *, const int *, int *);
void cpotrf_ru(const int *, float *, const int *, int *);
void zpotrf_ru(const int *, double *, const int *, int *);

#endif /* POTRF_H */
