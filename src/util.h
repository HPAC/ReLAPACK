#ifndef UTIL_H
#define UTIL_H

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define REC_SPLIT(n) (n >= 16) ? MIN(((n + 8) / 16) * 8, n - ((n + 8) / 16) * 8) : n / 2 

#endif /* UTIL_H */
