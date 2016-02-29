/*  -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/** RELAPACK_CTRSYL_REC2 solves the complex Sylvester matrix equation (unblocked algorithm)
 *
 * This routine is an exact copy of LAPACK's ctrsyl.
 * It serves as an unblocked kernel in the recursive algorithms. 
 * */
/* Subroutine */ void RELAPACK_ctrsyl_rec2(char *trana, char *tranb, integer 
	*isgn, integer *m, integer *n, complex *a, integer *lda, complex *b, 
	integer *ldb, complex *c__, integer *ldc, real *scale, integer *info, 
	ftnlen trana_len, ftnlen tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer j, k, l;
    static complex a11;
    static real db;
    static complex x11;
    static real da11;
    static complex vec;
    static real dum[1], eps, sgn, smin;
    static complex suml, sumr;
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Complex */ VOID cdotu_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int slabad_(real *, real *);
    extern doublereal clange_(char *, integer *, integer *, complex *, 
	    integer *, real *, ftnlen);
    extern /* Complex */ VOID cladiv_(complex *, complex *, complex *);
    static real scaloc;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *), xerbla_(char *, integer *, ftnlen);
    static real bignum;
    static logical notrna, notrnb;
    static real smlnum;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
    notrnb = lsame_(tranb, "N", (ftnlen)1, (ftnlen)1);
    *info = 0;
    if (! notrna && ! lsame_(trana, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! notrnb && ! lsame_(tranb, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*isgn != 1 && *isgn != -1) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < max(1,*m)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*m)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTRSY2", &i__1, (ftnlen)6);
	return;
    }
    *scale = 1.f;
    if (*m == 0 || *n == 0) {
	return;
    }
    eps = slamch_("P", (ftnlen)1);
    smlnum = slamch_("S", (ftnlen)1);
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    smlnum = smlnum * (real) (*m * *n) / eps;
    bignum = 1.f / smlnum;
/* Computing MAX */
    r__1 = smlnum, r__2 = eps * clange_("M", m, m, &a[a_offset], lda, dum, (
	    ftnlen)1), r__1 = max(r__1,r__2), r__2 = eps * clange_("M", n, n, 
	    &b[b_offset], ldb, dum, (ftnlen)1);
    smin = dmax(r__1,r__2);
    sgn = (real) (*isgn);
    if (notrna && notrnb) {
	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    for (k = *m; k >= 1; --k) {
		i__2 = *m - k;
/* Computing MIN */
		i__3 = k + 1;
/* Computing MIN */
		i__4 = k + 1;
		cdotu_(&q__1, &i__2, &a[k + min(i__3,*m) * a_dim1], lda, &c__[
			min(i__4,*m) + l * c_dim1], &c__1);
		suml.r = q__1.r, suml.i = q__1.i;
		i__2 = l - 1;
		cdotu_(&q__1, &i__2, &c__[k + c_dim1], ldc, &b[l * b_dim1 + 1]
			, &c__1);
		sumr.r = q__1.r, sumr.i = q__1.i;
		i__2 = k + l * c_dim1;
		q__3.r = sgn * sumr.r, q__3.i = sgn * sumr.i;
		q__2.r = suml.r + q__3.r, q__2.i = suml.i + q__3.i;
		q__1.r = c__[i__2].r - q__2.r, q__1.i = c__[i__2].i - q__2.i;
		vec.r = q__1.r, vec.i = q__1.i;
		scaloc = 1.f;
		i__2 = k + k * a_dim1;
		i__3 = l + l * b_dim1;
		q__2.r = sgn * b[i__3].r, q__2.i = sgn * b[i__3].i;
		q__1.r = a[i__2].r + q__2.r, q__1.i = a[i__2].i + q__2.i;
		a11.r = q__1.r, a11.i = q__1.i;
		da11 = (r__1 = a11.r, dabs(r__1)) + (r__2 = r_imag(&a11), 
			dabs(r__2));
		if (da11 <= smin) {
		    a11.r = smin, a11.i = 0.f;
		    da11 = smin;
		    *info = 1;
		}
		db = (r__1 = vec.r, dabs(r__1)) + (r__2 = r_imag(&vec), dabs(
			r__2));
		if (da11 < 1.f && db > 1.f) {
		    if (db > bignum * da11) {
			scaloc = 1.f / db;
		    }
		}
		q__3.r = scaloc, q__3.i = 0.f;
		q__2.r = vec.r * q__3.r - vec.i * q__3.i, q__2.i = vec.r * 
			q__3.i + vec.i * q__3.r;
		cladiv_(&q__1, &q__2, &a11);
		x11.r = q__1.r, x11.i = q__1.i;
		if (scaloc != 1.f) {
		    i__2 = *n;
		    for (j = 1; j <= i__2; ++j) {
			csscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L10: */
		    }
		    *scale *= scaloc;
		}
		i__2 = k + l * c_dim1;
		c__[i__2].r = x11.r, c__[i__2].i = x11.i;
/* L20: */
	    }
/* L30: */
	}
    } else if (! notrna && notrnb) {
	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    i__2 = *m;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = k - 1;
		cdotc_(&q__1, &i__3, &a[k * a_dim1 + 1], &c__1, &c__[l * 
			c_dim1 + 1], &c__1);
		suml.r = q__1.r, suml.i = q__1.i;
		i__3 = l - 1;
		cdotu_(&q__1, &i__3, &c__[k + c_dim1], ldc, &b[l * b_dim1 + 1]
			, &c__1);
		sumr.r = q__1.r, sumr.i = q__1.i;
		i__3 = k + l * c_dim1;
		q__3.r = sgn * sumr.r, q__3.i = sgn * sumr.i;
		q__2.r = suml.r + q__3.r, q__2.i = suml.i + q__3.i;
		q__1.r = c__[i__3].r - q__2.r, q__1.i = c__[i__3].i - q__2.i;
		vec.r = q__1.r, vec.i = q__1.i;
		scaloc = 1.f;
		r_cnjg(&q__2, &a[k + k * a_dim1]);
		i__3 = l + l * b_dim1;
		q__3.r = sgn * b[i__3].r, q__3.i = sgn * b[i__3].i;
		q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		a11.r = q__1.r, a11.i = q__1.i;
		da11 = (r__1 = a11.r, dabs(r__1)) + (r__2 = r_imag(&a11), 
			dabs(r__2));
		if (da11 <= smin) {
		    a11.r = smin, a11.i = 0.f;
		    da11 = smin;
		    *info = 1;
		}
		db = (r__1 = vec.r, dabs(r__1)) + (r__2 = r_imag(&vec), dabs(
			r__2));
		if (da11 < 1.f && db > 1.f) {
		    if (db > bignum * da11) {
			scaloc = 1.f / db;
		    }
		}
		q__3.r = scaloc, q__3.i = 0.f;
		q__2.r = vec.r * q__3.r - vec.i * q__3.i, q__2.i = vec.r * 
			q__3.i + vec.i * q__3.r;
		cladiv_(&q__1, &q__2, &a11);
		x11.r = q__1.r, x11.i = q__1.i;
		if (scaloc != 1.f) {
		    i__3 = *n;
		    for (j = 1; j <= i__3; ++j) {
			csscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L40: */
		    }
		    *scale *= scaloc;
		}
		i__3 = k + l * c_dim1;
		c__[i__3].r = x11.r, c__[i__3].i = x11.i;
/* L50: */
	    }
/* L60: */
	}
    } else if (! notrna && ! notrnb) {
	for (l = *n; l >= 1; --l) {
	    i__1 = *m;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = k - 1;
		cdotc_(&q__1, &i__2, &a[k * a_dim1 + 1], &c__1, &c__[l * 
			c_dim1 + 1], &c__1);
		suml.r = q__1.r, suml.i = q__1.i;
		i__2 = *n - l;
/* Computing MIN */
		i__3 = l + 1;
/* Computing MIN */
		i__4 = l + 1;
		cdotc_(&q__1, &i__2, &c__[k + min(i__3,*n) * c_dim1], ldc, &b[
			l + min(i__4,*n) * b_dim1], ldb);
		sumr.r = q__1.r, sumr.i = q__1.i;
		i__2 = k + l * c_dim1;
		r_cnjg(&q__4, &sumr);
		q__3.r = sgn * q__4.r, q__3.i = sgn * q__4.i;
		q__2.r = suml.r + q__3.r, q__2.i = suml.i + q__3.i;
		q__1.r = c__[i__2].r - q__2.r, q__1.i = c__[i__2].i - q__2.i;
		vec.r = q__1.r, vec.i = q__1.i;
		scaloc = 1.f;
		i__2 = k + k * a_dim1;
		i__3 = l + l * b_dim1;
		q__3.r = sgn * b[i__3].r, q__3.i = sgn * b[i__3].i;
		q__2.r = a[i__2].r + q__3.r, q__2.i = a[i__2].i + q__3.i;
		r_cnjg(&q__1, &q__2);
		a11.r = q__1.r, a11.i = q__1.i;
		da11 = (r__1 = a11.r, dabs(r__1)) + (r__2 = r_imag(&a11), 
			dabs(r__2));
		if (da11 <= smin) {
		    a11.r = smin, a11.i = 0.f;
		    da11 = smin;
		    *info = 1;
		}
		db = (r__1 = vec.r, dabs(r__1)) + (r__2 = r_imag(&vec), dabs(
			r__2));
		if (da11 < 1.f && db > 1.f) {
		    if (db > bignum * da11) {
			scaloc = 1.f / db;
		    }
		}
		q__3.r = scaloc, q__3.i = 0.f;
		q__2.r = vec.r * q__3.r - vec.i * q__3.i, q__2.i = vec.r * 
			q__3.i + vec.i * q__3.r;
		cladiv_(&q__1, &q__2, &a11);
		x11.r = q__1.r, x11.i = q__1.i;
		if (scaloc != 1.f) {
		    i__2 = *n;
		    for (j = 1; j <= i__2; ++j) {
			csscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L70: */
		    }
		    *scale *= scaloc;
		}
		i__2 = k + l * c_dim1;
		c__[i__2].r = x11.r, c__[i__2].i = x11.i;
/* L80: */
	    }
/* L90: */
	}
    } else if (notrna && ! notrnb) {
	for (l = *n; l >= 1; --l) {
	    for (k = *m; k >= 1; --k) {
		i__1 = *m - k;
/* Computing MIN */
		i__2 = k + 1;
/* Computing MIN */
		i__3 = k + 1;
		cdotu_(&q__1, &i__1, &a[k + min(i__2,*m) * a_dim1], lda, &c__[
			min(i__3,*m) + l * c_dim1], &c__1);
		suml.r = q__1.r, suml.i = q__1.i;
		i__1 = *n - l;
/* Computing MIN */
		i__2 = l + 1;
/* Computing MIN */
		i__3 = l + 1;
		cdotc_(&q__1, &i__1, &c__[k + min(i__2,*n) * c_dim1], ldc, &b[
			l + min(i__3,*n) * b_dim1], ldb);
		sumr.r = q__1.r, sumr.i = q__1.i;
		i__1 = k + l * c_dim1;
		r_cnjg(&q__4, &sumr);
		q__3.r = sgn * q__4.r, q__3.i = sgn * q__4.i;
		q__2.r = suml.r + q__3.r, q__2.i = suml.i + q__3.i;
		q__1.r = c__[i__1].r - q__2.r, q__1.i = c__[i__1].i - q__2.i;
		vec.r = q__1.r, vec.i = q__1.i;
		scaloc = 1.f;
		i__1 = k + k * a_dim1;
		r_cnjg(&q__3, &b[l + l * b_dim1]);
		q__2.r = sgn * q__3.r, q__2.i = sgn * q__3.i;
		q__1.r = a[i__1].r + q__2.r, q__1.i = a[i__1].i + q__2.i;
		a11.r = q__1.r, a11.i = q__1.i;
		da11 = (r__1 = a11.r, dabs(r__1)) + (r__2 = r_imag(&a11), 
			dabs(r__2));
		if (da11 <= smin) {
		    a11.r = smin, a11.i = 0.f;
		    da11 = smin;
		    *info = 1;
		}
		db = (r__1 = vec.r, dabs(r__1)) + (r__2 = r_imag(&vec), dabs(
			r__2));
		if (da11 < 1.f && db > 1.f) {
		    if (db > bignum * da11) {
			scaloc = 1.f / db;
		    }
		}
		q__3.r = scaloc, q__3.i = 0.f;
		q__2.r = vec.r * q__3.r - vec.i * q__3.i, q__2.i = vec.r * 
			q__3.i + vec.i * q__3.r;
		cladiv_(&q__1, &q__2, &a11);
		x11.r = q__1.r, x11.i = q__1.i;
		if (scaloc != 1.f) {
		    i__1 = *n;
		    for (j = 1; j <= i__1; ++j) {
			csscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L100: */
		    }
		    *scale *= scaloc;
		}
		i__1 = k + l * c_dim1;
		c__[i__1].r = x11.r, c__[i__1].i = x11.i;
/* L110: */
	    }
/* L120: */
	}
    }
    return;
} /* relapack_ctrsyl_rec2__ */

