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
CFLAGS = -O3 -Wall
FC = gfortran
FFLAGS = -O3 -Wall
AR = ar

# BLAS / LAPACK link lines (for testing)
LINK_BLAS = $(HOME)/phd/software/blas/blas_LINUX.a -lgfortran
LINK_LAPACK = $(HOME)/phd/software/lapack/liblapack.a

# vim: ft=make
