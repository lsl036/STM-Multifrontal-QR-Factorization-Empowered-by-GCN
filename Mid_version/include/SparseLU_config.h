/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:   	SpareseLU_config.h
 * BRIEF:   叙述LU配置
 *****************************************************************************/




#ifndef NRECIPROCAL
#if defined (MATHWORKS)
#define NRECIPROCAL
#endif
#endif





#if defined (LU_WINDOWS) || defined (LU_MINGW)

#define NPOSIX
#endif







#ifndef CHOLMOD_BLAS_H
#define CHOLMOD_BLAS_H





#if defined (__sun) || defined (MSOL2) || defined (ARCH_SOL2)
#define CHOLMOD_SOL2
#define CHOLMOD_ARCHITECTURE "Sun Solaris"

#elif defined (__sgi) || defined (MSGI) || defined (ARCH_SGI)
#define CHOLMOD_SGI
#define CHOLMOD_ARCHITECTURE "SGI Irix"

#elif defined (__linux) || defined (MGLNX86) || defined (ARCH_GLNX86)
#define CHOLMOD_LINUX
#define CHOLMOD_ARCHITECTURE "Linux"

#elif defined (__APPLE__)
#define CHOLMOD_MAC
#define CHOLMOD_ARCHITECTURE "Mac"

#elif defined (_AIX) || defined (MIBM_RS) || defined (ARCH_IBM_RS)
#define CHOLMOD_AIX
#define CHOLMOD_ARCHITECTURE "IBM AIX"



#elif defined (__alpha) || defined (MALPHA) || defined (ARCH_ALPHA)
#define CHOLMOD_ALPHA
#define CHOLMOD_ARCHITECTURE "Compaq Alpha"

#elif defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#if defined (__MINGW32__) || defined (__MINGW32__)
#define CHOLMOD_MINGW
#elif defined (__CYGWIN32__) || defined (__CYGWIN32__)
#define CHOLMOD_CYGWIN
#else
#define CHOLMOD_WINDOWS
#define BLAS_NO_UNDERSCORE
#endif
#define CHOLMOD_ARCHITECTURE "Microsoft Windows"

#elif defined (__hppa) || defined (__hpux) || defined (MHPUX) || defined (ARCH_HPUX)
#define CHOLMOD_HP
#define CHOLMOD_ARCHITECTURE "HP Unix"
#define BLAS_NO_UNDERSCORE

#elif defined (__hp700) || defined (MHP700) || defined (ARCH_HP700)
#define CHOLMOD_HP
#define CHOLMOD_ARCHITECTURE "HP 700 Unix"
#define BLAS_NO_UNDERSCORE

#else


#define CHOLMOD_ARCHITECTURE "unknown"
#endif








#if defined(CHOLMOD_SOL2) && !defined(NSUNPERF) && defined(BLAS64)
#define SUN64
#endif

#ifdef SUN64

#define BLAS_DTRSV dtrsv_64_
#define BLAS_DGEMV dgemv_64_
#define BLAS_DTRSM dtrsm_64_
#define BLAS_DGEMM dgemm_64_
#define BLAS_DSYRK dsyrk_64_
#define BLAS_DGER  dger_64_
#define BLAS_DSCAL dscal_64_
#define LAPACK_DPOTRF dpotrf_64_

#elif defined (BLAS_NO_UNDERSCORE)

#define BLAS_DTRSV dtrsv
#define BLAS_DGEMV dgemv
#define BLAS_DTRSM dtrsm
#define BLAS_DGEMM dgemm
#define BLAS_DSYRK dsyrk
#define BLAS_DGER  dger
#define BLAS_DSCAL dscal
#define LAPACK_DPOTRF dpotrf

#else

#define BLAS_DTRSV dtrsv_
#define BLAS_DGEMV dgemv_
#define BLAS_DTRSM dtrsm_
#define BLAS_DGEMM dgemm_
#define BLAS_DSYRK dsyrk_
#define BLAS_DGER  dger_
#define BLAS_DSCAL dscal_
#define LAPACK_DPOTRF dpotrf_

#endif







#if defined (LONGBLAS) || defined (BLAS64)
#define BLAS_INT Sparse_long
#else
#define BLAS_INT int
#endif



#define CHECK_BLAS_INT (sizeof (BLAS_INT) < sizeof (Int))
#define EQ(K,k) (((BLAS_INT) K) == ((Int) k))





void BLAS_DGEMV (char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
	double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
	double *Y, BLAS_INT *incy) ;

#define BLAS_dgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
	BLAS_DGEMV (trans, &M, &N, alpha, A, &LDA, X, &INCX, beta, Y, &INCY) ; \
}

void BLAS_DTRSV (char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
	BLAS_INT *lda, double *X, BLAS_INT *incx) ;

#define BLAS_dtrsv(uplo,trans,diag,n,A,lda,X,incx) \
{ \
    BLAS_INT N = n, LDA = lda, INCX = incx ; \
	BLAS_DTRSV (uplo, trans, diag, &N, A, &LDA, X, &INCX) ; \
}

void BLAS_DTRSM (char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
	BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
	BLAS_INT *ldb) ;

#define BLAS_dtrsm(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, LDB = ldb ; \
	BLAS_DTRSM (side, uplo, transa, diag, &M, &N, alpha, A, &LDA, B, &LDB);\
}

void BLAS_DGEMM (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
	BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
	BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc) ;

#define BLAS_dgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
{ \
    BLAS_INT M = m, N = n, K = k, LDA = lda, LDB = ldb, LDC = ldc ; \
	BLAS_DGEMM (transa, transb, &M, &N, &K, alpha, A, &LDA, B, &LDB, beta, \
	    C, &LDC) ; \
}

void BLAS_DSYRK (char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
	double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
	BLAS_INT *ldc) ;

#define BLAS_dsyrk(uplo,trans,n,k,alpha,A,lda,beta,C,ldc) \
{ \
    BLAS_INT N = n, K = k, LDA = lda, LDC = ldc ; \
	BLAS_DSYRK (uplo, trans, &N, &K, alpha, A, &LDA, beta, C, &LDC) ; \
} \

void LAPACK_DPOTRF (char *uplo, BLAS_INT *n, double *A, BLAS_INT *lda,
	BLAS_INT *info) ;

#define LAPACK_dpotrf(uplo,n,A,lda,info) \
{ \
    BLAS_INT N = n, LDA = lda, INFO = 1 ; \
	LAPACK_DPOTRF (uplo, &N, A, &LDA, &INFO) ; \
    info = INFO ; \
}



void BLAS_DSCAL (BLAS_INT *n, double *alpha, double *Y, BLAS_INT *incy) ;

#define BLAS_dscal(n,alpha,Y,incy) \
{ \
    BLAS_INT N = n, INCY = incy ; \
	BLAS_DSCAL (&N, alpha, Y, &INCY) ; \
}

void BLAS_DGER (BLAS_INT *m, BLAS_INT *n, double *alpha,
	double *X, BLAS_INT *incx, double *Y, BLAS_INT *incy,
	double *A, BLAS_INT *lda) ;

#define BLAS_dger(m,n,alpha,X,incx,Y,incy,A,lda) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
	BLAS_DGER (&M, &N, alpha, X, &INCX, Y, &INCY, A, &LDA) ; \
}

#endif








#define BLAS_GEMM(m,n,k,A,B,ldb,C,ldac) \
{ \
    double alpha = -1, beta = 1 ; \
    BLAS_dgemm ("N", "T", m, n, k, &alpha, A, ldac, B, ldb, &beta, C, ldac) ; \
}







#define BLAS_GER(m,n,x,y,A,d) \
{ \
    double alpha = -1 ; \
    BLAS_dger (m, n, &alpha, x, 1, y, 1, A, d) ; \
}







#define BLAS_GEMV(m,n,A,x,y,d) \
{ \
    double alpha = -1, beta = 1 ; \
    BLAS_dgemv ("N", m, n, &alpha, A, d, x, 1, &beta, y, 1) ; \
}







#define BLAS_TRSV(m,A,b,d) \
{ \
    BLAS_dtrsv ("L", "N", "U", m, A, d, b, 1) ; \
}







#define BLAS_TRSM_RIGHT(m,n,A,lda,B,ldb) \
{ \
    double alpha = 1 ; \
    BLAS_dtrsm ("R", "L", "T", "U", m, n, &alpha, A, lda, B, ldb) ; \
}







#define BLAS_SCAL(n,s,x) \
{ \
    double alpha = REAL_COMPONENT (s) ; \
    BLAS_dscal (n, &alpha, (double *) x, 1) ; \
}
