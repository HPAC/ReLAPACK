# ReLAPACK library name
LIBNAME = librelapack.a

# compiler options
CC = gcc
CFLAGS = -O3 -std=c99 -pedantic

# archive creation
AR = ar

# BLAS / LAPACK link line for test
LINK_TEST = -Llapack -Lblas -llapack -lblas -lgfortran -lm

# vim: ft=make
