# library nam
LIBNAME = librelapack.a

# compiler options
CC = gcc
CFLAGS = -O3 -std=c99 -Wpedantic
FC = gfortran
FFLAGS = -O3 -Wpedantic
AR = ar

# BLAS / LAPACK link line for test (relative to test/)
LINK_TEST = -L../lapack/ -L../blas/ -llapack -lblas -lgfortran

# vim: ft=make
