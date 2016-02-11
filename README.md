ReLAPACK
========

Recursive LAPACK Collection

ReLAPACK offers a collection of recursive algorithms for many of LAPACK's
compute kernels.  Since it preserves LAPACK's established interfaces, ReLAPACK
integrates effortlessly into existing application codes.  ReLAPACK's routines
not only to outperform reference LAPACK but also improve upon the performance of
tuned implementations, such as OpenBLAS and MKL.

Coverage
--------
For a detailed list of covered operations and an overview of operations to which
recursion is not efficiently applicable, see [coverage.md](coverage.md).

Installation
------------
To compile with the default configuration, simply run `make` to create the
library `librelapack.a`.  

Configuration
-------------
ReLAPACK has two configuration files: `make.inc`, which is included by the
Makefile, and `config.h` which is included in the source files.

### Build and Testing Environment
The build environment (compilers and flags) and the test configuration (linker
flags for BLAS and LAPACK, and the test matrix size) are specified in `make.inc`

### Routine Selection
Which of the available ReLAPACK routines are included in `librelapack.a` is
configured in `make.inc`.  The default (`all`) is to include all routines;
alternatively routines can be included or excluded individually or by operation.

## Routine Names
By default, ReLAPACK's routine names coincide with the functionally equivalent
LAPACK routines. By setting `RELAPACK_AS_LAPACK` to 0 in `config.h`, they will
receive the prefix `RELAPACK_`; e.g. the LU decomposition `dgetrf` would become
`RELAPACK_dgetrf`.

ReLAPACK uses the FORTRAN interfaces of BLAS and unblocked LAPACK routines,
i.e., all routines have a suffix `_`.  Setting `BLAS_UNDERSCORE` or
`LAPACK_UNDERSCORE` to 0 in `config.h` removes this suffix for calls to the
respective library.

### Crossover Size
The crossover size determines below which matrix sizes ReLAPACK's recursive
algorithms switch to LAPACK's unblocked routines to avoid tiny BLAS Level 3
routines.  The crossover size is set in `config.h` and can be chosen either
globally for the entire library, by operation, or individually by routine.

### Allowing Temporary Buffers
Two of ReLAPACK's routines make use of temporary buffers, which are allocated
and freed within ReLAPACK.  Setting `ALLOW_MALLOC` (or one of the routine
specific counterparts) to 0 in `config.h` will disable these buffers.  The
affected routines are:

 * `xsytrf`: The LDL decomposition requires a buffer of size n^2 / 2.  As in
   LAPACK, this size can be queried by setting `lWork = -1` and the passed
   buffer will be used if it is large enough; only if it is not, a local buffer
   will be allocated.  
   
   The advantage of this mechanism is that ReLAPACK will seamlessly work even
   with codes that statically provide too little memory instead of breaking
   them.

 * `xsygst`: The reduction of a real symmetric-definite generalized eigenproblem
   to standard form can use an auxiliary buffer of size n^2 / 2 to avoid
   redundant computations.  It thereby performs about 30% less FLOPs than
   LAPACK.
