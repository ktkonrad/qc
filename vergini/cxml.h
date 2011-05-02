/* cxml.h : BLAS/LAPACK/other CXML headers, from Michael Haggerty.

   Barnett 99/10/28
*/

#ifndef ALEX_CXML_H
#define ALEX_CXML_H 1

// BLAS C++ headers, hacked myself. Uses nice pass-by-reference interface.
// See M. Haggerty for original headers.

extern "C" void dgemm_(
	const char &transa,
	const char &transb,
	const int &m,
	const int &n,
	const int &k,
	const double &alpha,
	const double *a,
	const int &lda,
	const double *b,
	const int &ldb,
	const double &beta,
	double *c,
	const int &ldc
	);

extern "C" void dgemt_(
	const char &trans,
	const int &m,
	const int &n,
	const double &alpha,
	const double *a,
	const int &lda,
	const double *b,
	const int &ldb
	);

extern "C" void dgema_(
	const char &transa,
	const char &transb,
	const int &m,
	const int &n,
	const double &alpha,
	const double *a,
	const int &lda,
	const double &beta,
	const double *b,
	const int &ldb,
	double *c,
	const int &ldc
	);

extern "C" void dgemv_(
	const char &trans,
	const int &m,
	const int &n,
	const double &alpha,
	const double *a,
	const int &lda,
	const double *x,
	const int &incx,
	const double &beta,
	double *y,
	const int &incy
	);

extern "C" double dnrm2_(
	const int &n,
	const double *x,
	const int &incx
	);

extern "C" double ddot_(
	const int &n,
	const double *x,
	const int &incx,
	const double *y,
	const int &incy
	);

extern "C" double dscal_(
	const int &n,
	const double &alpha,
	const double *x,
	const int &incx
	);

extern "C" double dcopy_(
	const int &n,
	const double *x,
	const int &incx,
	const double *y,
	const int &incy
	);

// Signal Processing headers, for bounce/ fft routines:

extern "C" int dfft_(
	const char &informat,
	const char &outformat,
	const char &direction,
	const double *in,
	const double *out,
	const int &size,
	const int &stride
	);


// Sort routine headers, for wigner.cc:

extern "C" int dsortqx_(
			const char &order,
			const int &size,
			const double *x,
			const int &incx,
			const int *index
			);

// LAPACK headers:
extern "C" void dposv_(
		      const char &uplo,
		      const int &n,
		      const int &nrhs,
		      const double *A,
		      const int &lda,
		      const double *B,
		      const int &ldb,
		      const int &info
		      );

extern "C" void dgelss_(
	const int &m,			// (input)
	const int &n,			// (input)
	const int &nrhs,		// (input)
	double *a,			// a[n][lda] (input/output)
	const int &lda,			// (input)
	double *b,			// b[nrhs][ldb] (input/output)
	const int &ldb,			// (input)
	double *s,			// s[min(m,n)] (output)
	const double &rcond,		// (input)
	int &rank,			// (output)
	double *work,			// work[lwork] (workspace/output)
	const int &lwork,		// (input)
	int &info			// (output)
	);

extern "C" void dsyev_(
        const char &jobz,               /* (input) */
        const char &uplo,               /* (input) */
        const int &n,                   /* (input) */
        double *a,                      /* a[n][lda] (input/output) */
        const int &lda,                 /* (input) */
        double *w,                      /* w[n] (output) */
        double *work,                   /* work[lwork] (workspace/output) */
        const int &lwork,               /* (input) */
        int &info                       /* (output) */
        );
		

// LAPACK headers from M. Haggerty's (not-C++ part of) header files:

extern "C" void dgesvd_(
        const char *jobu,               /* (input) */
        const char *jobvt,              /* (input) */
        const int *m,                   /* (input) */
        const int *n,                   /* (input) */
        double *a,                      /* a[n][lda] (input/output) */
        const int *lda,                 /* (input) */
        double *s,                      /* s[min(m,n)] (output) */
        double *u,                      /* u[ucol][ldu] (output) */
        const int *ldu,                 /* (input) */
        double *vt,                     /* vt[n][ldvt] (output) */
        const int *ldvt,                /* (input) */
        double *work,                   /* work[lwork] (workspace/output) */
        const int *lwork,               /* (input) */
        int *info                       /* (output) */
        );

extern "C" void dsygv_(
        const int *itype,               /* (input) */
        const char *jobz,               /* (input) */
        const char *uplo,               /* (input) */
        const int *n,                   /* (input) */
        double *a,                      /* a[n][lda] (input/output) */
        const int *lda,                 /* (input) */
        double *b,                      /* b[n][ldb] (input/output) */
        const int *ldb,                 /* (input) */
        double *w,                      /* w[n] (output) */
        double *work,                   /* work[lwork] (workspace/output) */
        const int *lwork,               /* (input) */
        int *info                       /* (output) */
        );

extern "C" void dsygvx_(
        const int *itype,               /* (input) */
        const char *jobz,               /* (input) */
        const char *range,              /* (input) */
        const char *uplo,               /* (input) */
        const int *n,                   /* (input) */
        double *a,                      /* a[n][lda] (input) */
        const int *lda,                 /* (input) */
        double *b,                      /* b[n][ldb] (input/output) */
        const int *ldb,                 /* (input) */
        double *vl,                     /* vl (input) */
        double *vu,                     /* vu (input) */
        int *il,                        /* (input) */
        int *iu,                        /* (input) */
        double *abstol,                 /* (input) */
        int *m,                         /* m (output) */
        double *w,                      /* w[n] (output) */
        double *z,                      /* z[n][ldz] (output) */
        const int *ldz,                 /* (input) */
        double *work,                   /* work[lwork] (workspace/output) */
        const int *lwork,               /* (input) */
        int *iwork,                     /* (work) */
        int *ifail,                     /* (work/output) */
        int *info                       /* (output) */
        );

#endif /* ALEX_CXML_H */
