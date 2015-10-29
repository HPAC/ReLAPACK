LARPACK
=======

Linear Algebra Recursive Package

LARPACK provides recursive implementations of many blocked algorithms within
LAPACK

Coverage
--------

* xlauum (A = L * L^T)
* xtrtri (L = L^-1)
* xpotrf (L L^T = A)
* xgetrf (L U = A)
* xhegst (A = L A L^T)

Not covered
-----------
The following routines are not covered because recursive variants would require
considerably more FLOPs:

* QR decomposition (and related)
* symmetric to tridiagonal reduction
* full to upper Hessenberg reduction
