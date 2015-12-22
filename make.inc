# library options
# select routines to be included in the library.
#
# Shorthands:
# all:		xlauum	xsygst	xtrtri	xpotrf	xgetrf	xtrsyl
# xlauum:	slauum	dlauum	clauum	zlauum
# xsygst:  	ssygst	dsygst	chegst	zhegst
# xpotrf:  	spotrf	dpotrf	cpotrf	zpotrf
# xgetrf:  	sgetrf	dgetrf	cgetrf	zgetrf
# xtrsyl:  	strsyl	dtrsyl	ctrsyl	ztrsyl
ROUTINES = all

# compiler options
CC = gcc
CFLAGS = -g -O3 -Wall
FC = gfortran
FFLAGS = -g -O3 -Wall
AR = ar

# TESTING CONFIGURATION

# BLAS / LAPACK link lines (for
LINK_BLAS = $(HOME)/phd/software/blas/blas_LINUX.a -lgfortran
LINK_LAPACK = $(HOME)/phd/software/lapack/liblapack.a
LINK_LAPACK = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -ldl -lpthread -lm -fopenmp

# test matrix size
TEST_SIZE = 100

# vim: ft=make
