ReLAPACK Examples
=================

Since ReLAPACK replaces parts of LAPACK, any LAPACK example involving the
covered routines applies directly to ReLAPACK.


LAPACK drivers
--------------
Since ReLAPACK only covers compute kernels, many of LAPACK's "driver" routines
are good examples of how LAPACK, and thus ReLAPACK are used.  For instance the
driver `chegv` uses the reduction `chegst` as part of the solution of a
generalized Hermitian-definite eigenvalue problem.

The FORTRAN source code for this driver (copied from
[LAPACK](http://www.netlib.org/lapack/explore-html/d0/db9/chegv_8f_source.html))
is provided in `chegv.f`.


Using the ReLAPACK interface in C
---------------------------------
To demonstrate the use of ReLAPACK's interface (i.e., not its LAPACK wrappers),
`linsys.c` provides an example that solves a linear system with a general matrix
directly using `ReLAPACK_dgetrf`.

The example uses the FORTRAN interface for BLAS and LAPACK (trailing underscore)
and needs to be linked to ReLAPACK, LAPACK and BLAS.  E.g., on Mac OS X: `gcc
linsys.c -L.. -lrelapack -framework Accelerate`.  Note that in this example, we
can replace `RELAPACK_dgetrf` with `dgetrf_` in line 55 to use the LAPACK
interface to ReLAPACK; in that case, we must link with ReLAPACK before linking
with LAPACK.
