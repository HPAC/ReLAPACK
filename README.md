LARPACK
=======

Linear Algebra Recursive Package

LARPACK provides recursive implementations of many blocked algorithms within
LAPACK

Coverage
--------

# `slauum`, `dlauum`, `clauum`, `zlauum`
Symmetric triangular matrix squaring
* A = L^T L
* A = U U^T

# `strtri`, `dtrtri`, `ctrtri`, `ztrtri`
Triangular matrix inversion
* L = inv(L)
* U = inv(U)

# `spotrf`, `dpotrf`, `cpotrf`, `zpotrf`
Cholesky decomposition
* L L^T = A
* U^T U = A

# `sgetrf`, `dgetrf`, `cgetrf`, `zgetrf`
LU decomposition with pivoting
* P L U = A

# `ssygst`, `dsygst`, `chegst`, `zhegst`
Two-sided symmetric `trmm` or `trsm`
* A = inv(L) A inv(L^T)
* A = inv(U^T) A inv(U)
* A = L^T A L
* A = U A U^T

# `strsyl`, `dtrsyl`, `ctrsyl`, `ztrsyl`
Sylvester equation solver
* A X + B Y = C -> X
* A^T X + B Y = C -> X
* A X + B^T Y = C -> X
* A^T X + B^T Y = C -> X
* A X - B Y = C -> X
* A^T X - B Y = C -> X
* A X - B^T Y = C -> X
* A^T X - B^T Y = C -> X

Not covered yet
---------------
The following routines can potentially be translated to recursive algorithms but
are not covered yet:

# `stgsyl`, `dtgsyl`, `ctgsyl`, `ztgsyl`
Generalized Sylvester solver
* A R - L B = C, D R - L E = F
* A^T R + D^T L = C, R B^T - L E^T = -F

# `ssytrf`, `dsytrf`, `chetrf`, `zhetrf`
LDL symmetric decomposition
* L D L^T = A
* U^T D U = A

Recursion not applicable
------------------------
The following routines are not covered because recursive variants would require
considerably more FLOPs:

# QR decomposition (and related)
* `sgeqrf`, `dgeqrf`, `cgeqrf`, `zgeqrf`
* `sgerqf`, `dgerqf`, `cgerqf`, `zgerqf`
* `sgeqlf`, `dgeqlf`, `cgeqlf`, `zgeqlf`
* `sgelqf`, `dgelqf`, `cgelqf`, `zgelqf`

# symmetric reduction to tridiagonal
* `ssytrd`, `dsytrd`, `csytrd`, `zsytrd`

# reduction to upper Hessenberg
* `sgehrd`, `dgehrd`, `cgehrd`, `zgehrd`
