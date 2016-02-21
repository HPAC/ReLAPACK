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
For an overview of the configuration options see [config.md](config.md).

Testing
-------
ReLAPACK's test suite compares its routines numerically with LAPACK's
counterparts.  The test are located int `test/` and started by running `make
test`.
