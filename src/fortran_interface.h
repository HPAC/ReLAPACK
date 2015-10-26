#ifndef FORTRAN_INTERFACE_H
#define FORTRAN_INTERFACE_H

#ifdef BLAS_UNDERSCORE
#define BLAS(routine) routine ## _
#else
#define BLAS(routine) routine
#endif

#ifdef LAPACK_UNDERSCORE
#define LAPACK(routine) routine ## _
#else
#define LAPACK(routine) routine
#endif

#endif /* FORTRAN_INTERFACE_H */
