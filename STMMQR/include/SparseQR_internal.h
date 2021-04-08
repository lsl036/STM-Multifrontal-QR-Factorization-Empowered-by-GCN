// =============================================================================
// === SparseQR_internal.h =====================================================
// =============================================================================

#ifndef SparseQR_INTERNAL_H
#define SparseQR_INTERNAL_H

#include "SparseQR_struct.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#define Long Sparse_long

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define EMPTY (-1)
#define TRUE 1
#define FALSE 0 
#define IMPLIES(p,q) (!(p) || (q))

#ifndef NULL
#define NULL ((void *) 0)
#endif

// 矩阵A的元素 对应下标转换
#define INDEX(i,j,lda) ((i) + ((j)*(lda)))

// FLIP 是 "关于-1的否定", 用于标记一个通常是非负的整数 i
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) (((i) < EMPTY) ? FLIP (i) : (i))

#define ITYPE SPARSE_LONG
#define DTYPE SPARSE_DOUBLE
#define ID Sparse_long_id

// -----------------------------------------------------------------------------
// 统计浮点计算量； 如果使用了TBB则禁用
// -----------------------------------------------------------------------------

#define FLOP_COUNT(f) { if (cc->SPQR_grain <= 1) cc->SPQR_flopcount += (f) ; }

// 数值分解中 factorize 和 kernel 的工作空间申请。
typedef struct qr_work_struct
{
    Long *Stair1 ;         //
    Long *Cmap ;           // 大小为 maxfn
    Long *Fmap ;           // 大小为 n 
    double *WTwork ;       // 大小为 (fchunk + (keepH ? 0:1)) * maxfn

    double *Stack_head ;     // 栈底
    double *Stack_top ;      // 栈顶

    Long sumfrank ;         //  波前阵f 的秩和
    Long maxfrank ;         //  波前阵f 秩的最大值

    // 用于计算w的二范数， dead column 的范数向量
    double wscale ;         //  范数标量因子
    double wssq ;           //  范数平方和

    Long numa_node_rank;
} qr_work ;

// qr_blob 是 qr_kernel 运算需要的对象集合

typedef struct qr_blob_struct
{
    double tol ;
    qr_symbolic *QRsym ; // 符号对象, 定义在 HnuSparseQR.h 中
    qr_numeric  *QRnum ; // 数值对象, 定义在 HnuSparseQR.h 中
    qr_work  *Work ; 
    Long *Cm ;
    double **Cblock ;
    double *Sx ;
    Long ntol ;
    Long fchunk ;
    sparse_common *cc ;
} qr_blob;

//并行时传输属性的任务结构体。
typedef struct task_struct
{
    Long id;
    qr_blob *Blob;
    //TPSM_t *tpool;
    int lock_id;
}tasktodo;

inline double qr_abs (double x) 
{
    return (fabs (x)) ;
}

// 把两个非负 Long 相加，并检查是否越界

inline Long qr_add (Long a, Long b, int *ok)
{
    Long c = a + b ;
    if (c < 0)
    {
        (*ok) = FALSE ;
        return (EMPTY) ;
    }
    return (c) ;
}

// 把两个非负 Long 相乘，并检查是否越界

inline Long qr_mult (Long a, Long b, int *ok)
{
    Long c = a * b ;
    if (((double) c) != ((double) a) * ((double) b))
    {
        (*ok) = FALSE ;
        return (EMPTY) ;
    }
    return (c) ;
}


// =============================================================================
// === BLAS interface ==========================================================
// =============================================================================

#include "Sparse_internal.h"

#undef CHECK_BLAS_INT
#undef EQ
#define CHECK_BLAS_INT (sizeof (BLAS_INT) < sizeof (Long))
#define EQ(K,k) (((BLAS_INT) K) == ((Long) k))

#ifdef SUN64

#define BLAS_DNRM2    dnrm2_64_
#define LAPACK_DLARF  dlarf_64_
#define LAPACK_DLARFG dlarfg_64_
#define LAPACK_DLARFT dlarft_64_
#define LAPACK_DLARFB dlarfb_64_

#elif defined (BLAS_NO_UNDERSCORE)

#define BLAS_DNRM2    dnrm2
#define LAPACK_DLARF  dlarf
#define LAPACK_DLARFG dlarfg
#define LAPACK_DLARFT dlarft
#define LAPACK_DLARFB dlarfb

#else

#define BLAS_DNRM2    dnrm2_
#define LAPACK_DLARF  dlarf_
#define LAPACK_DLARFG dlarfg_
#define LAPACK_DLARFT dlarft_
#define LAPACK_DLARFB dlarfb_

#endif

// =============================================================================
// === BLAS and LAPACK prototypes ==============================================
// =============================================================================

void LAPACK_DLARFT (char *direct, char *storev, BLAS_INT *n, BLAS_INT *k,
    double *V, BLAS_INT *ldv, double *Tau, double *T, BLAS_INT *ldt) ;

void LAPACK_DLARFB (char *side, char *trans, char *direct, char *storev,
    BLAS_INT *m, BLAS_INT *n, BLAS_INT *k, double *V, BLAS_INT *ldv,
    double *T, BLAS_INT *ldt, double *C, BLAS_INT *ldc, double *Work,
    BLAS_INT *ldwork) ;

double BLAS_DNRM2 (BLAS_INT *n, double *X, BLAS_INT *incx) ;

void LAPACK_DLARFG (BLAS_INT *n, double *alpha, double *X, BLAS_INT *incx,
    double *tau) ;

void LAPACK_DLARF (char *side, BLAS_INT *m, BLAS_INT *n, double *V,
    BLAS_INT *incv, double *tau, double *C, BLAS_INT *ldc, double *Work) ;

#endif
