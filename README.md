ReLAPACK
========

Recursive LAPACK Collection

ReLAPACK offers a collection of recursive algorithms for many of LAPACK's
compute kernels.  Since it preserves LAPACK's established interfaces, ReLAPACK
integrates effortlessly into existing application codes.  ReLAPACK's routines
not only outperform the reference LAPACK but also improve upon the performance
of tuned implementations, such as OpenBLAS and MKL.


Coverage
--------
For a detailed list of covered operations and an overview of operations to which
recursion is not efficiently applicable, see [coverage.md](coverage.md).


Installation
------------
To compile with the default configuration, simply run `make` to create the
library `librelapack.a`.  

### Dependencies
ReLAPACK builds on top of [BLAS](http://www.netlib.org/blas/) and unblocked
kernels from [LAPACK](http://www.netlib.org/lapack/).  There are many optimized
and machine specific implementations of these libraries, which are commonly
provided by hardware vendors or available as open source (e.g.,
[OpenBLAS](http://www.openblas.net/)).


Configuration
-------------
For an overview of the configuration options see [config.md](config.md).


Testing
-------
ReLAPACK's test suite compares its routines numerically with LAPACK's
counterparts.  To set up the test located int `test/` you need to specify link
flags for BLAS and LAPACK (version 3.5.0 or newer) in `make.inc`; then `make
test` runs the tests.  For details, see [test/README.md](test/README.md).


Examples
--------
Since ReLAPACK replaces parts of LAPACK, any LAPACK example involving the
covered routines applies directly to ReLAPACK.  A few separate examples are
given in `examples/`. For details, see [examples/README.md](examples/README.md).
