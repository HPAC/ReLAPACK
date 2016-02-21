# library options
# select routines to be included in the library.
#
# Shorthands:
# all:		xlauum	xsygst	xtrtri	xpotrf	xsytrf	xgetrf	xtrsyl	xtgsyl
# xlauum:	slauum	dlauum	clauum	zlauum
# xsygst:	ssygst	dsygst	chegst	zhegst
# xpotrf:	spotrf	dpotrf	cpotrf	zpotrf
# xsytrf:	ssytrf	dsytrf	csytrf	chetrf zsytrf zhetrf
# 			ssytrf_rook	dsytrf_rook	csytrf_rook	chetrf_rook	zsytrf_rook	zhetrf_rook
# xgetrf:	sgetrf	dgetrf	cgetrf	zgetrf
# xtrsyl:	strsyl	dtrsyl	ctrsyl	ztrsyl
# xtgsyl:	stgsyl	dtgsyl	ctgsyl	zgrsyl
ROUTINES = all

# compiler options
CC = gcc
CFLAGS = -O3 -Wall
FC = gfortran
FFLAGS = -O3 -Wall
AR = ar

# TESTING CONFIGURATION

# BLAS / LAPACK link line (for test)
LINK_TEST = -L$(HOME)/phd/software/lapack/ -L$(HOME)/phd/software/blas/ -llapack -lblas -lgfortran

# test matrix size
TEST_SIZE = 100

# vim: ft=make
