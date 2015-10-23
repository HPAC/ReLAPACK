CC = gcc
CFLAGS = -O3 -Wall
FC = gfortran
FFLAGS = -O3 -Wall

LINK_BLAS = $(HOME)/phd/software/blas/blas_LINUX.a -lgfortran
LINK_LAPACK = $(HOME)/phd/software/lapack/liblapack.a

