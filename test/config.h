#ifndef TEST_CONFIG_H
#define TEST_CONFIG_H

// error bound for single and single complex routines
#ifndef SINGLE_ERR_BOUND
#define SINGLE_ERR_BOUND 1e-5
#endif

// error bound for double an double complex routines
#ifndef DOUBLE_ERR_BOUND
#define DOUBLE_ERR_BOUND 1e-13
#endif

// size of test matrices
#ifndef TEST_SIZE
#define TEST_SIZE 100
#endif

#endif /* TEST_CONFIG_H */
