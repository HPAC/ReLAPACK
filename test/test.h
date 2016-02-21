#ifndef TEST_H
#define TEST_H

#include "config.h"
#include "../src/relapack.h"
#include "lapack.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int i, n, n2, err_bound, fail, info;
double error;

#define CAT(A, B) A ## B
#define XCAT(A, B) CAT(A, B)

#define XLAPACK(X) LAPACK(X)
#define STR(X) #X
#define XSTR(X) STR(X)
#define PRE pre
#define POST post

#define TEST(...) \
    PRE(); \
    i = 0; \
    XCAT(RELAPACK_, ROUTINE)(__VA_ARGS__); \
    i = 1; \
    XLAPACK(ROUTINE)(__VA_ARGS__); \
    POST(); \
    fail |= error > ERR_BOUND; \
    printf("%s(%s)\t%g\n", XSTR(ROUTINE), #__VA_ARGS__, error);

// data type stuff
#define XPREF(A) XCAT(DT_PREFIX, A)
#define xmalloc  XPREF(malloc)
#define x2matgen XPREF(2matgen)
#define x2vecerr XPREF(2vecerr)
#define datatype XPREF(datatype_)
#define ERR_BOUND XPREF(ERR_BOUND_)
#define MONE XPREF(MONE)
#define ZERO XPREF(ZERO)
#define ONE XPREF(ONE)
#define x1 XPREF(DT_MULT)
#define xCTRANS XPREF(CTRANS)

#define sdatatype_ float
#define ddatatype_ double
#define cdatatype_ float
#define zdatatype_ double

#define sERR_BOUND_ SINGLE_ERR_BOUND
#define dERR_BOUND_ DOUBLE_ERR_BOUND
#define cERR_BOUND_ SINGLE_ERR_BOUND
#define zERR_BOUND_ DOUBLE_ERR_BOUND

#define imalloc(S) malloc((S) * sizeof(int))
#define smalloc(S) malloc((S) * sizeof(float))
#define dmalloc(S) malloc((S) * sizeof(double))
#define cmalloc(S) malloc((S) * 2 * sizeof(float))
#define zmalloc(S) malloc((S) * 2 * sizeof(double))

#define sDT_MULT 1
#define dDT_MULT 1
#define cDT_MULT 2
#define zDT_MULT 2

#define sCTRANS "T"
#define dCTRANS "T"
#define cCTRANS "C"
#define zCTRANS "C"

// some constants
const float  sMONE[]  = {-1};
const double dMONE[]  = {-1};
const float  cMONE[]  = {-1, 0};
const double zMONE[]  = {-1, 0};
const float  sZERO[]  = {0};
const double dZERO[]  = {0};
const float  cZERO[]  = {0, 0};
const double zZERO[]  = {0, 0};
const float  sONE[]   = {1};
const double dONE[]   = {1};
const float  cONE[]   = {1, 0};
const double zONE[]   = {1, 0};

const int iMONE[]  = {-1};
const int iZERO[]  = {0};
const int iONE[]   = {1};
const int iTWO[]   = {2};
const int iTHREE[] = {3};
const int iFOUR[]  = {4};

void tests();

int main(int argc, char* argv[]) {
    n = TEST_SIZE;
    n2 = (3 * n) / 4;
    fail = 0;

    tests();

	return fail;
}

#endif /* TEST_H */
