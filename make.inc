# library nam
LIBNAME = librelapack.a

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
