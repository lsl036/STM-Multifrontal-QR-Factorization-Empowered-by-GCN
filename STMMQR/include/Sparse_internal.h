/* ========================================================================== */
/* === Include/SparseCore_internal.h =========================================== */
/* ========================================================================== */

#ifndef SPARSE_INTERNAL_H
#define SPARSE_INTERNAL_H

/* 假设支持大文件。如果出现问题，则用-DNLARGEFILE 编译 */
#ifndef NLARGEFILE

#undef  _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#undef  _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64

#endif

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>

#undef TRUE
#undef FALSE
#define TRUE 1
#define FALSE 0
#define BOOLEAN(x) ((x) ? TRUE : FALSE)

#ifndef NULL
#define NULL ((void *) 0)
#endif

#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) (((i) < EMPTY) ? FLIP (i) : (i))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MAX3(a,b,c) (((a) > (b)) ? (MAX (a,c)) : (MAX (b,c)))
#define MAX4(a,b,c,d) (((a) > (b)) ? (MAX3 (a,c,d)) : (MAX3 (b,c,d)))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define IMPLIES(p,q) (!(p) || (q))

#define SIGN(x) (((x) < 0) ? (-1) : (((x) > 0) ? 1 : 0))

#define ROUNDUP(x,s) ((s) * (((x) + ((s) - 1)) / (s)))

#define ERROR(status,msg) \
    SparseCore_error (status, __FILE__, __LINE__, msg, Common)

#define RETURN_IF_NULL(A,result) \
{ \
    if ((A) == NULL) \
    { \
	if (Common->status != SPARSE_OUT_OF_MEMORY) \
	{ \
	    ERROR (SPARSE_INVALID, "argument missing") ; \
	} \
	return (result) ; \
    } \
}

/* Return if Common is NULL or invalid */
#define RETURN_IF_NULL_COMMON(result) \
{ \
    if (Common == NULL) \
    { \
	return (result) ; \
    } \
    if (Common->itype != ITYPE || Common->dtype != DTYPE) \
    { \
	Common->status = SPARSE_INVALID ; \
	return (result) ; \
    } \
}

#define IS_NAN(x)	SPARSE_IS_NAN(x)
#define IS_ZERO(x)	SPARSE_IS_ZERO(x)
#define IS_NONZERO(x)	SPARSE_IS_NONZERO(x)
#define IS_LT_ZERO(x)	SPARSE_IS_LT_ZERO(x)
#define IS_GT_ZERO(x)	SPARSE_IS_GT_ZERO(x)
#define IS_LE_ZERO(x)	SPARSE_IS_LE_ZERO(x)

#define HUGE_DOUBLE 1e308

#include "SparseBase_config.h"

#define Size_max ((size_t) (-1))

/* routines for doing arithmetic on size_t, and checking for overflow */
size_t SparseCore_add_size_t (size_t a, size_t b, int *ok) ;
size_t SparseCore_mult_size_t (size_t a, size_t k, int *ok) ;

/* -------------------------------------------------------------------------- */
/* double , Sparse_long */
/* -------------------------------------------------------------------------- */

#define Real double
#define Int Sparse_long
#define Int_max Sparse_long_max
#define LONG
#define DOUBLE
#define ITYPE SPARSE_LONG
#define DTYPE SPARSE_DOUBLE
#define ID Sparse_long_id

/* ========================================================================== */
/* === real arithmetic ============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* pattern 二进制 */
/* -------------------------------------------------------------------------- */

#define P_TEMPLATE(name)		p_ ## name
#define P_ASSIGN2(x,z,p,ax,az,q)	x [p] = 1

/* -------------------------------------------------------------------------- */
/* real 实数 */
/* -------------------------------------------------------------------------- */

#define R_TEMPLATE(name)			r_ ## name
#define R_ASSEMBLE(x,z,p,ax,az,q)		x [p] += ax [q]
#define R_ASSIGN(x,z,p,ax,az,q)			x [p]  = ax [q]
#define R_ASSIGN_CONJ(x,z,p,ax,az,q)		x [p]  = ax [q]
#define R_ASSIGN_REAL(x,p,ax,q)			x [p]  = ax [q]
#define R_XTYPE_OK(type)			((type) == SPARSE_REAL)
#define R_IS_NONZERO(ax,az,q)			IS_NONZERO (ax [q])
#define R_IS_ZERO(ax,az,q)			IS_ZERO (ax [q])
#define R_IS_ONE(ax,az,q)			(ax [q] == 1)
#define R_MULT(x,z,p, ax,az,q, bx,bz,r)		x [p]  = ax [q] * bx [r]
#define R_MULTADD(x,z,p, ax,az,q, bx,bz,r)	x [p] += ax [q] * bx [r]
#define R_MULTSUB(x,z,p, ax,az,q, bx,bz,r)	x [p] -= ax [q] * bx [r]
#define R_MULTADDCONJ(x,z,p, ax,az,q, bx,bz,r)	x [p] += ax [q] * bx [r]
#define R_MULTSUBCONJ(x,z,p, ax,az,q, bx,bz,r)	x [p] -= ax [q] * bx [r]
#define R_ADD(x,z,p, ax,az,q, bx,bz,r)		x [p]  = ax [q] + bx [r]
#define R_ADD_REAL(x,p, ax,q, bx,r)		x [p]  = ax [q] + bx [r]
#define R_CLEAR(x,z,p)				x [p]  = 0
#define R_CLEAR_IMAG(x,z,p)
#define R_DIV(x,z,p,ax,az,q)			x [p] /= ax [q]
#define R_LLDOT(x,p, ax,az,q)			x [p] -= ax [q] * ax [q]

#define R_DIV_REAL(x,z,p, ax,az,q, bx,r)	x [p] = ax [q] / bx [r]
#define R_MULT_REAL(x,z,p, ax,az,q, bx,r)	x [p] = ax [q] * bx [r]

#define R_LDLDOT(x,p, ax,az,q, bx,r)		x [p] -=(ax[q] * ax[q])/ bx[r]


/* 检查 A->xtype 和两个数组 A->x , A->z 是否有效. 如果状态为 SPARSE_OUT_OF_MEMORY 
 * 则 报错.  A 可以是稀疏矩阵, 稠密矩阵, 分解因子factor 或者为 三元组. */

#define RETURN_IF_XTYPE_INVALID(A,xtype1,xtype2,result) \
{ \
    if ((A)->xtype < (xtype1) || (A)->xtype > (xtype2) || \
        ((A)->xtype != SPARSE_PATTERN && ((A)->x) == NULL)) \
    { \
	if (Common->status != SPARSE_OUT_OF_MEMORY) \
	{ \
	    ERROR (SPARSE_INVALID, "invalid xtype") ; \
	} \
	return (result) ; \
    } \
}

/* ========================================================================== */
/* === Architecture and BLAS ================================================ */
/* ========================================================================== */

#define BLAS_OK Common->blas_ok

#if defined (__sun) || defined (MSOL2) || defined (ARCH_SOL2)
#define SPARSE_SOL2
#define SPARSE_ARCHITECTURE "Sun Solaris"

#elif defined (__sgi) || defined (MSGI) || defined (ARCH_SGI)
#define SPARSE_SGI
#define SPARSE_ARCHITECTURE "SGI Irix"

#elif defined (__linux) || defined (MGLNX86) || defined (ARCH_GLNX86)
#define SPARSE_LINUX
#define SPARSE_ARCHITECTURE "Linux"

#elif defined (__APPLE__)
#define SPARSE_MAC
#define SPARSE_ARCHITECTURE "Mac"

#elif defined (_AIX) || defined (MIBM_RS) || defined (ARCH_IBM_RS)
#define SPARSE_AIX
#define SPARSE_ARCHITECTURE "IBM AIX"

#elif defined (__alpha) || defined (MALPHA) || defined (ARCH_ALPHA)
#define SPARSE_ALPHA
#define SPARSE_ARCHITECTURE "Compaq Alpha"

#elif defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#if defined (__MINGW32__) || defined (__MINGW32__)
#define SPARSE_MINGW
#elif defined (__CYGWIN32__) || defined (__CYGWIN32__)
#define SPARSE_CYGWIN
#else
#define SPARSE_WINDOWS
#define BLAS_NO_UNDERSCORE
#endif
#define SPARSE_ARCHITECTURE "Microsoft Windows"

#elif defined (__hppa) || defined (__hpux) || defined (MHPUX) || defined (ARCH_HPUX)
#define SPARSE_HP
#define SPARSE_ARCHITECTURE "HP Unix"
#define BLAS_NO_UNDERSCORE

#elif defined (__hp700) || defined (MHP700) || defined (ARCH_HP700)
#define SPARSE_HP
#define SPARSE_ARCHITECTURE "HP 700 Unix"
#define BLAS_NO_UNDERSCORE

#else
// 如果架构未知 BLAS_BY_VALUE, BLAS_NO_UNDERSCORE 和BLAS_CHAR_ARG 需要自己定义
#define SPARSE_ARCHITECTURE "unknown"
#endif

/* ========================================================================== */
/* === BLAS and LAPACK names ================================================ */
/* ========================================================================== */

/* BLA各种版本的原型  */

/* 确定是否使用64位Sun BLAS */
#if defined(SPARSE_SOL2) && !defined(NSUNPERF) && defined(BLAS64)
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

/* ========================================================================== */
/* === BLAS and LAPACK integer arguments ==================================== */
/* ========================================================================== */

/* 如果BLAS使用64位整数，则使用-DBLAS64 编译CHOL、LU和QR */

#if defined (LONGBLAS) || defined (BLAS64)
#define BLAS_INT Sparse_long
#else
#define BLAS_INT int
#endif

/* 如果BLAS整数小于基本的SPARSE整数，那么在从Int转换到BLAS_INT时，
 * 我们需要检查整数是否溢出。如果有任何整数溢出，外部定义的BLAS_OK变量将被设置为FALSE。
 * 在调用BLAS_*宏之前，应该将BLAS_OK设置为TRUE。
 */

#define CHECK_BLAS_INT (sizeof (BLAS_INT) < sizeof (Int))
#define EQ(K,k) (((BLAS_INT) K) == ((Int) k))

/* ================================================================ */
/* === BLAS、 LAPACK 的 原型 和 宏 ================================ */
/* ================================================================ */

void BLAS_DGEMV (char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
	double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
	double *Y, BLAS_INT *incy) ;

#define BLAS_dgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGEMV (trans, &M, &N, alpha, A, &LDA, X, &INCX, beta, Y, &INCY) ; \
    } \
}

void BLAS_DTRSV (char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
	BLAS_INT *lda, double *X, BLAS_INT *incx) ;

#define BLAS_dtrsv(uplo,trans,diag,n,A,lda,X,incx) \
{ \
    BLAS_INT N = n, LDA = lda, INCX = incx ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda) && EQ (INCX,incx))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DTRSV (uplo, trans, diag, &N, A, &LDA, X, &INCX) ; \
    } \
}

void BLAS_DTRSM (char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
	BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
	BLAS_INT *ldb) ;

#define BLAS_dtrsm(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, LDB = ldb ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (LDB,ldb))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DTRSM (side, uplo, transa, diag, &M, &N, alpha, A, &LDA, B, &LDB);\
    } \
}

void BLAS_DGEMM (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
	BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
	BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc) ;

#define BLAS_dgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
{ \
    BLAS_INT M = m, N = n, K = k, LDA = lda, LDB = ldb, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (K,k) && \
        EQ (LDA,lda) && EQ (LDB,ldb) && EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGEMM (transa, transb, &M, &N, &K, alpha, A, &LDA, B, &LDB, beta, \
	    C, &LDC) ; \
    } \
}

void BLAS_DSYRK (char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
	double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
	BLAS_INT *ldc) ;

#define BLAS_dsyrk(uplo,trans,n,k,alpha,A,lda,beta,C,ldc) \
{ \
    BLAS_INT N = n, K = k, LDA = lda, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (K,k) && EQ (LDA,lda) && \
        EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DSYRK (uplo, trans, &N, &K, alpha, A, &LDA, beta, C, &LDC) ; \
    } \
} \

void LAPACK_DPOTRF (char *uplo, BLAS_INT *n, double *A, BLAS_INT *lda,
	BLAS_INT *info) ;

#define LAPACK_dpotrf(uplo,n,A,lda,info) \
{ \
    BLAS_INT N = n, LDA = lda, INFO = 1 ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	LAPACK_DPOTRF (uplo, &N, A, &LDA, &INFO) ; \
    } \
    info = INFO ; \
}

/* ========================================================================== */

void BLAS_DSCAL (BLAS_INT *n, double *alpha, double *Y, BLAS_INT *incy) ;

#define BLAS_dscal(n,alpha,Y,incy) \
{ \
    BLAS_INT N = n, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DSCAL (&N, alpha, Y, &INCY) ; \
    } \
}

void BLAS_DGER (BLAS_INT *m, BLAS_INT *n, double *alpha,
	double *X, BLAS_INT *incx, double *Y, BLAS_INT *incy,
	double *A, BLAS_INT *lda) ;

#define BLAS_dger(m,n,alpha,X,incx,Y,incy,A,lda) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
          EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGER (&M, &N, alpha, X, &INCX, Y, &INCY, A, &LDA) ; \
    } \
}


#endif
