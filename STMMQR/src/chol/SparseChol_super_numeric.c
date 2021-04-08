/**
 * @file SparseChol_super_numeric.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-22
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef NGPL
#ifndef NSUPERNODAL

#include "Sparse_internal.h"
#include "SparseChol.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include "tpsm.h"
#include <sys/time.h>


/* ========================================================================== */
/* === 用于常规数值分解的模板代码 ============= */
/* ========================================================================== */

#define REAL

/* SparseChol_super_numeric的模板例程 */

#include "Sparse_template.h"
#include "tpsm_base.h"

#undef L_ENTRY
#undef L_CLEAR
#undef L_ASSIGN
#undef L_MULTADD
#undef L_ASSEMBLE
#undef L_ASSEMBLESUB

#ifdef REAL

/* -------------------------------------------------------------------------- */
/* A, F, L都是实数型 */
/* -------------------------------------------------------------------------- */

#define L_ENTRY 1
#define L_CLEAR(Lx,p)               Lx [p] = 0
#define L_ASSIGN(Lx,q, Ax,p)     Lx [q] = Ax [p]
#define L_MULTADD(Lx,q, Ax,p, f) Lx [q] += Ax [p] * f [0]
#define L_ASSEMBLE(Lx,q,b)          Lx [q] += b [0]
#define L_ASSEMBLESUB(Lx,q,C,p)     Lx [q] -= C [p]

#endif

// 分解所需的基本参数
typedef struct chol_base_args_struct
{
    Int *Super;
    Int *Lpi;
    Int *Lpx;
    double *Lx;
    Int *Ls;
    Int *Map;
    int stype;
    Int *Ap;
    int Apacked;
    Int *Anz;
    Int *Ai;
    double *Ax;
    Int *Fp;
    int Fpacked;
    Int *Fnz;
    Int *Fi;
    double *Fx;
    double *beta;
    int repeat_supernode;
    Int *Head;
    Int *Next;
    Int *Lpos_save;
    Int *Lpos;
    Int *Next_save;
    // double *C;
    size_t Lmaxcsize;
    double *one;
    double *zero;
    // Int *RelativeMap;
    Int *SuperMap;
    Int nsuper;
    sparse_factor *L;
    sparse_common *Common;
}chol_base_args;

// 传入线程池的参数，包括超节点ID和基本参数
typedef struct chol_snode_args_struct
{
    int snode_id;
    int nscol, psx, nsrow, k1, psi, nsrow2 ;
    chol_base_args *base_args;
}chol_snode_args;

void* numeric_update(chol_snode_args *snode_args)
{
    int s = snode_args->snode_id;
    Int *Super = snode_args->base_args->Super;
    Int *Lpi = snode_args->base_args->Lpi;
    Int *Lpx = snode_args->base_args->Lpx;
    double *Lx = snode_args->base_args->Lx;
    Int *Ls = snode_args->base_args->Ls;
    Int *Map = snode_args->base_args->Map;
    int stype = snode_args->base_args->stype;
    Int *Ap = snode_args->base_args->Ap;
    int Apacked = snode_args->base_args->Apacked;
    Int *Anz = snode_args->base_args->Anz;
    Int *Ai = snode_args->base_args->Ai;
    double *Ax = snode_args->base_args->Ax;
    Int *Fp = snode_args->base_args->Fp;
    int Fpacked = snode_args->base_args->Fpacked;
    Int *Fnz = snode_args->base_args->Fnz;
    Int *Fi = snode_args->base_args->Fi;
    double *Fx = snode_args->base_args->Fx;
    double *beta = snode_args->base_args->beta;
    int repeat_supernode = snode_args->base_args->repeat_supernode;
    Int *Head = snode_args->base_args->Head;
    Int *Next = snode_args->base_args->Next;
    Int *Lpos_save = snode_args->base_args->Lpos_save;
    Int *Lpos = snode_args->base_args->Lpos;
    Int *Next_save = snode_args->base_args->Next_save;
    // double *C = snode_args->base_args->C;   
    size_t Lmaxcsize = snode_args->base_args->Lmaxcsize;
    double *one = snode_args->base_args->one;
    double *zero = snode_args->base_args->zero;
    // Int *RelativeMap = snode_args->base_args->RelativeMap;
    Int *SuperMap = snode_args->base_args->SuperMap;
    Int nsuper = snode_args->base_args->nsuper;
    sparse_factor *L = snode_args->base_args->L;
    sparse_common *Common = snode_args->base_args->Common;

    int i, j, p, k, imap, pf, pfend, d, px, q, dancestor, ss, sparent;
    /* ------------------------------------------------------------------ */
    /* 获取超节点s的大小 */
    /* ------------------------------------------------------------------ */
    int k1 = Super [s] ;            /* s包含L的k1到k2-1的列 */
    int k2 = Super [s+1] ;
    int nscol = k2 - k1 ;           /* 所有s的列的# */
    int psi = Lpi [s] ;             /* 指向Ls中s的第一行的指针 */
    int psx = Lpx [s] ;             /* 指向Lx中s的第一行的指针 */
    int psend = Lpi [s+1] ;         /* 指针刚经过Ls中s的最后一行 */
    int nsrow = psend - psi ;       /* 所有s中的行数 */

    /* ------------------------------------------------------------------ */
    /* 超节点s置0 */
    /* ------------------------------------------------------------------ */
    int pend = psx + nsrow * nscol ;        /* s 大小为 nsrow*nscol */

// #pragma omp parallel for num_threads(CHOLESKY_OMP_NUM_THREADS)   \
//     schedule (static) if ( pend - psx > 1024 )

    for (p = psx ; p < pend ; p++) {
        L_CLEAR (Lx,p);
    }

// #pragma omp parallel for num_threads(CHOLESKY_OMP_NUM_THREADS)   \
//     if ( nsrow > 128 )

    for (k = 0 ; k < nsrow ; k++)
    {
        Map [Ls [psi + k]] = k ;
    }
    
    int pk = psx ;
    if (stype != 0)
    {
// #pragma omp parallel for private ( p, pend, pfend, pf, i, j, imap, q )  \
//     num_threads(CHOLESKY_OMP_NUM_THREADS) if ( k2-k1 > 64 )
        for (k = k1 ; k < k2 ; ++k)
        {
            /* 将A的第k列拷贝到超节点中 */
            p = Ap [k] ;
            pend = (Apacked) ? (Ap [k+1]) : (p + Anz [k]) ;
            for ( ; p < pend ; p++)
            {
                /* L的第i行在s的行映射Map[i]中 */
                i = Ai [p] ;
                if (i >= k)
                {
                    /* 这个测试在这里只是为了避免段错误。
                     * 如果测试为假，则A的数字因数分解为未定义。
                     * 它不会检测所有无效条目，只检测其中的一部分。*/
                    imap = Map [i] ;
                    if (imap >= 0 && imap < nsrow)
                    {
                        // Lx无重复写入，安全
                        /* Lx [Map [i] + pk] = Ax [p] ; */
                        L_ASSIGN (Lx,(imap+(psx+(k-k1)*nsrow)), Ax,p) ;
                    }
                }
            }               
        }
    }
    else
    {
        for (k = k1 ; k < k2 ; ++k)
        {
            double fjk[2];
            /* 拷贝A*F的第k列到超节点中 */
            pf = Fp [k] ;
            pfend = (Fpacked) ? (Fp [k+1]) : (p + Fnz [k]) ;
            for ( ; pf < pfend ; pf++)
            {
                j = Fi [pf] ;
                /* fjk = Fx [pf] ; */
                L_ASSIGN (fjk,0, Fx,pf) ;
                p = Ap [j] ;
                pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
                for ( ; p < pend ; p++)
                {
                    i = Ai [p] ;
                    if (i >= k)
                    {
                        /* 参见上面关于imap的讨论 */
                        imap = Map [i] ;
                        if (imap >= 0 && imap < nsrow)
                        {
                            /* Lx [Map [i] + pk] += Ax [p] * fjk ; */
                            L_MULTADD (Lx,(imap+(psx+(k-k1)*nsrow)),
                                    Ax,p, fjk) ;
                        }
                    }
                }
            }
        }
    }

    // 默认跳过
    /* 如果非零，在超节点的对角线上加上beta */
    if (beta [0] != 0.0)
    {
        /* 注意，只使用了beta的实部 */
        pk = psx ;
        for (k = k1 ; k < k2 ; k++)
        {
            /* Lx [pk] += beta [0] ; */
            L_ASSEMBLE (Lx,pk, beta) ;
            pk += nsrow + 1 ;       /* 前进到下一个对角线条目 */
        }
    }
    /* ------------------------------------------------------------------ */
    /* 保存/恢复超级节点列表 */
    /* ------------------------------------------------------------------ */
    // d无重复写入，安全
    if (!repeat_supernode)
    {
        /* 在矩阵非正定的情况下，保存挂起的后代列表。
         * 还要为每个后代d保存Lpos，以便我们可以找到d的哪一部分用于更新s。 */
        for (d = Head [s] ; d != EMPTY ; d = Next [d])
        {
            Lpos_save [d] = Lpos [d] ;
            Next_save [d] = Next [d] ;
        }
    }
    else
    {
        for (d = Head [s] ; d != EMPTY ; d = Next [d])
        {
            Lpos [d] = Lpos_save [d] ;
            Next [d] = Next_save [d] ;
        }
    }
    /* ------------------------------------------------------------------ */
    /* 用每个挂起的后代d更新超节点s */
    /* ------------------------------------------------------------------ */
    int dnext = Head[s];
    int kd1, kd2, ndcol, pdi, pdx, pdend, ndrow, pdi1, pdx1, pdi2, ndrow1, ndrow2, ndrow3;
    while ( dnext != EMPTY )
    {
        d = dnext ;

        kd1 = Super [d] ;       /* d包含L的kd1到kd2-1的列 */
        kd2 = Super [d+1] ;
        ndcol = kd2 - kd1 ;     /* d中所有列的# */
        pdi = Lpi [d] ;         /* 在Ls中指向d的第一行的指针 */
        pdx = Lpx [d] ;         /* 指向Lx中d的第一行的指针 */
        pdend = Lpi [d+1] ;     /* 指针刚好经过Ls中d的最后一行 */
        ndrow = pdend - pdi ;   /* d中所有的行数 */

        p = Lpos [d] ;          /* d的第一行偏移量影响s */
        pdi1 = pdi + p ;        /* 在Ls中，ptr到d的第一行影响s */
        pdx1 = pdx + p ;        /* ptr到d的第一行影响Lx中的s */

        /* 要更新s, d中必须至少还有一行 */
        for (pdi2 = pdi1 ; pdi2 < pdend && Ls [pdi2] < k2 ; pdi2++) ;
        ndrow1 = pdi2 - pdi1 ;      /* d的第一部分中的#行 */
        ndrow2 = pdend - pdi1 ;     /* 剩余d中的行数 */
        ndrow3 = ndrow2 - ndrow1 ;  /* C2的行数 */

	    double *C = (double*)malloc(Lmaxcsize*sizeof(double));
        // 执行其中一个对称的秩k运算
        // C = alpha*AA' + beta*C
        BLAS_dsyrk ("L", "N",
            ndrow1, ndcol,              /* N, K: L1大小为ndrow1*ndcol*/
            one,                        /* ALPHA:  1 */
            Lx + L_ENTRY*pdx1, ndrow,   /* A, LDA: L1, ndrow */
            zero,                       /* BETA:   0 */
            C, ndrow2) ;                /* C, LDC: C1 */

        /* 计算C的剩余部分(ndrow2-ndrow1)*ndrow1大小的块
        * C2 = L2*L1' */
        if (ndrow3 > 0)
        {
            // 执行矩阵矩阵乘
            // C = alpha*op(A)*op(B) + beta*C
            // op(X) = X or X'
            BLAS_dgemm ("N", "C",
                ndrow3, ndrow1, ndcol,          /* M, N, K */
                one,                            /* ALPHA:  1 */
                Lx + L_ENTRY*(pdx1 + ndrow1),   /* A, LDA: L2 */
                ndrow,                          /* ndrow */
                Lx + L_ENTRY*pdx1,              /* B, LDB: L1 */
                ndrow,                          /* ndrow */
                zero,                           /* BETA:   0 */
                C + L_ENTRY*ndrow1,             /* C, LDC: C2 */
                ndrow2) ;
        }

        int *RelativeMap_local;
		RelativeMap_local = (int*)malloc(ndrow2*sizeof(int));

// #pragma omp parallel for num_threads(CHOLESKY_OMP_NUM_THREADS)   \
//     if ( ndrow2 > 64 )
        for (i = 0 ; i < ndrow2 ; ++i)
        {
            // RelativeMap [i] = Map [Ls [pdi1 + i]] ;
            RelativeMap_local [i] = Map [Ls [pdi1 + i]] ;
        }

        /* ---------------------------------------------------------- */
        /* 使用相对映射将C组装成超节点s */
        /* ---------------------------------------------------------- */
// #pragma omp parallel for private ( j, i, px, q )                \
//     num_threads(CHOLESKY_OMP_NUM_THREADS) if (ndrow1 > 64 )

        for (j = 0 ; j < ndrow1 ; ++j)              /* cols k1:k2-1 */
        {
            // px = psx + RelativeMap [j] * nsrow ;
            px = psx + RelativeMap_local [j] * nsrow ;
            for (i = j ; i < ndrow2 ; ++i)          /* rows k1:n-1 */
            {
                /* Lx [px + RelativeMap [i]] -= C [i + pj] ; */
                // q = px + RelativeMap [i] ;
                q = px + RelativeMap_local [i] ;
                L_ASSEMBLESUB (Lx,q, C, i+ndrow2*j) ;
            }
        }


        /* -------------------------------------------------------------- */
        /* 为它的下一个祖先准备这个超节点d */
        /* -------------------------------------------------------------- */

        dnext = Next [d] ;
        // printf("node %d dnext = %d\n", d, dnext);
        if (!repeat_supernode)
        {
            // printf("pdi2 = %d Ls [pdi2] = %d\n", pdi2, Ls [pdi2]);
            Lpos [d] = pdi2 - pdi ;
            if (Lpos [d] < ndrow)
            {
                dancestor = SuperMap [Ls [pdi2]] ;
                /* 将d放在它下一个祖先的链接列表中 */
                Next [d] = Head [dancestor] ;
                Head [dancestor] = d ;
            }
        }
        free(RelativeMap_local);
        free(C);
    }  /* 子代超节点循环的结束 */

    snode_args->nsrow = nsrow ;
    snode_args->nscol = nscol ;
    snode_args->psx = psx ;
    snode_args->k1 = k1 ;
    snode_args->psi = psi ;
}

// 数值分解的线程池并行函数
void* parallel_numeric_factor(void * arg)
{   
    chol_snode_args *snode_args = (chol_snode_args *) arg;
    int s = snode_args->snode_id;
    int nsuper = snode_args->base_args->nsuper;
    int repeat_supernode = snode_args->base_args->repeat_supernode;
    Int *Head = snode_args->base_args->Head;
    Int *Next = snode_args->base_args->Next;
    double *Lx = snode_args->base_args->Lx;
    double *one = snode_args->base_args->one;
    sparse_factor *L = snode_args->base_args->L;
    sparse_common *Common = snode_args->base_args->Common;
    int nsrow = snode_args->nsrow;
    int nscol = snode_args->nscol;
    int psx = snode_args->psx;
    int k1 = snode_args->k1;

    int nscol_new = 0;
    int nscol2 = (repeat_supernode) ? (nscol_new) : (nscol) ;
    int info, p, pend, ss, nsrow2;   
    // A = LL'
    LAPACK_dpotrf ("L",
        nscol2,                     /* N: nscol2 */
        Lx + L_ENTRY*psx, nsrow,    /* A, LDA: S1, nsrow */
        info) ;                     /* INFO */

    /* ------------------------------------------------------------------ */
    /* 检查矩阵是否是正定的 */
    /* ------------------------------------------------------------------ */
    if (repeat_supernode)
    {
        /* 主导部分已重构;它一定成功了 */
        info = 0 ;
        /* 将剩下的超结点归零 */
        p = psx + nsrow * nscol_new ;
        pend = psx + nsrow * nscol ;            /* s大小为nsrow*nscol */
        for ( ; p < pend ; p++)
        {
            /* Lx [p] = 0 ; */
            L_CLEAR (Lx,p) ;
        }
    }
    /* 如果blas_ok为假，info在LAPACK_*potrf中设置为1。如果分解成功，则在dpotrf/zpotrf中将其设置为零。 */
    if (CHECK_BLAS_INT && !Common->blas_ok)
    {
        ERROR (SPARSE_TOO_LARGE, "problem too large for the BLAS") ;
    }
    if (info != 0)
    {
        /* 矩阵不是正定的。如果L的对角线上有NaN, dpotrf/zpotrf只在对角线上有0时才报告错误。 */
        if (Common->status == SPARSE_OK)
        {
            ERROR (SPARSE_NOT_POSDEF, "matrix not positive definite") ;
        }

        /* L->minor 是L的列，它包含一个零或负的对角项 */
        L->minor = k1 + info - 1 ;

        /* 清除所有后续超级节点的链接列表 */
        for (ss = s+1 ; ss < nsuper ; ++ss)
        {
            Head [ss] = EMPTY ;
        }

        /* 将这个超节点，以及所有剩余的超节点置为0 */
        pend = L->xsize ;
        for (p = psx ; p < pend ; p++)
        {
            /* Lx [p] = 0. ; */
            L_CLEAR (Lx,p) ;
        }

        if (info == 1 || Common->quick_return_if_not_posdef)
        {
            Head [s] = EMPTY ;
            // return (Common->status >= SPARSE_OK) ;
        }
        else
        {
            repeat_supernode = TRUE ;
            s-- ;
            nscol_new = info - 1 ;
            return 0;
        }
    }

    /* ------------------------------------------------------------------ */
    /* 计算次对角块并为父块准备超节点 */
    /* ------------------------------------------------------------------ */
    nsrow2 = nsrow - nscol2 ;
    if (nsrow2 > 0)
    {
        // 解其中一个矩阵方程：op(A)*X = alpha*B, or X*op(A) = alpha*B
        BLAS_dtrsm ("R", "L", "C", "N",
            nsrow2, nscol2,                 /* M, N */
            one,                            /* ALPHA: 1 */
            Lx + L_ENTRY*psx, nsrow,        /* A, LDA: L1, nsrow */
            Lx + L_ENTRY*(psx + nscol2),    /* B, LDB, L2, nsrow */
            nsrow) ;

        if (CHECK_BLAS_INT && !Common->blas_ok)
        {
            ERROR (SPARSE_TOO_LARGE, "problem too large for the BLAS") ;
        }
    }
    snode_args->nsrow2 = nsrow2 ;
}

int numeric_finish(chol_snode_args *snode_args)
{
    int s = snode_args->snode_id;
    int nsuper = snode_args->base_args->nsuper;
    int repeat_supernode = snode_args->base_args->repeat_supernode;
    Int *Head = snode_args->base_args->Head;
    Int *Next = snode_args->base_args->Next;
    Int *SuperMap = snode_args->base_args->SuperMap;
    Int *Lpos = snode_args->base_args->Lpos;
    Int *Ls = snode_args->base_args->Ls;
    double *Lx = snode_args->base_args->Lx;
    double *one = snode_args->base_args->one;
    sparse_factor *L = snode_args->base_args->L;
    sparse_common *Common = snode_args->base_args->Common;
    int nsrow = snode_args->nsrow;
    int nscol = snode_args->nscol;
    int psx = snode_args->psx;
    int k1 = snode_args->k1;
    int psi = snode_args->psi;
    int nsrow2 = snode_args->nsrow2;

    int sparent ;
    if (nsrow2 > 0)
    {
        if (!repeat_supernode)
        {
            /* Lpos [s]是影响s父节点的第一行的偏移量 */
            Lpos [s] = nscol ;
            sparent = SuperMap [Ls [psi + nscol]] ;
            /* 将s放在其父节点的链接列表中 */
            Next [s] = Head [sparent] ;
            Head [sparent] = s ;
        }
    }

    Head [s] = EMPTY ;  /* 链接列表的超级节点s不再需要 */

    if (repeat_supernode)
    {
        /* 矩阵非正定;完成了对包含负对角的超节点的清理 */
        return (Common->status >= SPARSE_OK) ;
    }
}

/**
 * @brief   仅当BLAS中发生整型溢出时，此函数返回FALSE。
 *          否则，不管矩阵是否正，它都返回真。
 * 
 */
static int TEMPLATE (SparseChol_super_numeric)
(
    /* ---- input ---- */
    //TPSM_t *tpool,
    sparse_csc *A,          /* 分解的稀疏矩阵 */
    sparse_csc *F,          /* F = A' or A(:,f)' */
    double beta [2],            /* beta*I 加到被分解的稀疏矩阵的对角线上 */
    /* ---- in/out --- */
    sparse_factor *L,          /* 分解 */
    /* -- workspace -- */
    dense_array *Cwork,       /* 大小为(L->maxcsize)-by-1 */
    /* --------------- */
    sparse_common *Common
)
{
    int max_depth, *Snode_num, **Snode_distribution;
    double one [2], zero [2], tstart ;
    double *Lx, *Ax, *Fx, *C ;
    Int *Super, *Head, *Ls, *Lpi, *Lpx, *Map, *SuperMap, *RelativeMap, *Next,
        *Lpos, *Fp, *Fi, *Fnz, *Ap, *Ai, *Anz, *Iwork, *Next_save, *Lpos_save,
        *Previous;
    Int nsuper, n, j, i, k, s, p, pend, k1, k2, nscol, psi, psx, psend, nsrow,
        pj, d, kd1, kd2, info, ndcol, ndrow, pdi, pdx, pdend, pdi1, pdi2, pdx1,
        ndrow1, ndrow2, px, dancestor, sparent, dnext, nsrow2, ndrow3, pk, pf,
        pfend, stype, Apacked, Fpacked, q, imap, repeat_supernode, nscol2, ss,
        tail, nscol_new = 0;

    /* ---------------------------------------------------------------------- */
    /* 防止BLAS中的整型溢出 */
    /* ---------------------------------------------------------------------- */

    /* 如果BLAS中发生整数溢出，Common->status设置为SPARSE_TOO_LARGE，并且Lx的内容未定义。 */
    Common->blas_ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    nsuper = L->nsuper ;
    n = L->n ;
    C = Cwork->x ;      /* 大小为L->maxcsize的工作空间 */

    one [0] =  1.0 ;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
    one [1] =  0. ;
    zero [0] = 0. ;     /* BETA for *syrk, *herk, and *gemm */
    zero [1] = 0. ;

    /* Iwork必须是大小2n + 5*nsuper，调用SparseChol_super_numeric分配。
     * 这里不能分配内存，因为SparseChol_super_numeric初始化SuperMap，
     * 并且如果需要增加空间大小，SparseCore_allocate_work不会保留现有的工作空间。 */

    /* 分配的整数工作区 */
    Iwork = Common->Iwork ;
    SuperMap    = Iwork ;                                   /* size n (i/i/l) */
    RelativeMap = Iwork + n ;                               /* size n (i/i/l) */
    Next        = Iwork + 2*((size_t) n) ;                  /* size nsuper*/
    Lpos        = Iwork + 2*((size_t) n) + nsuper ;         /* size nsuper*/
    Next_save   = Iwork + 2*((size_t) n) + 2*((size_t) nsuper) ;/* size nsuper*/
    Lpos_save   = Iwork + 2*((size_t) n) + 3*((size_t) nsuper) ;/* size nsuper*/
    Previous    = Iwork + 2*((size_t) n) + 4*((size_t) nsuper) ;/* size nsuper*/

    Map  = Common->Flag ;   /* 大小为n，使用Flag作为Map数组的工作空间 */
    Head = Common->Head ;   /* 大小为n+1,仅有Head [0..nsuper-1]被使用 */

    Ls = L->s ;
    Lpi = L->pi ;
    Lpx = L->px ;

    Super = L->super ;

    Lx = L->x ;

    max_depth = L->max_depth;
    Snode_num = L->Snode_num;
    Snode_distribution = L->Snode_distribution;

    stype = A->stype ;

    if (stype != 0)
    {
        /* F无法访问 */
        Fp = NULL ;
        Fi = NULL ;
        Fx = NULL ;

        Fnz = NULL ;
        Fpacked = TRUE ;
    }
    else
    {
        Fp = F->p ;
        Fi = F->i ;
        Fx = F->x ;
        Fnz = F->nz ;
        Fpacked = F->packed ;
    }

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Anz = A->nz ;
    Apacked = A->packed ;

    /* 清除Map，以便可以检测到A的pattern的变化 */
// #pragma omp parallel for num_threads(CHOLESKY_OMP_NUM_THREADS) \
//     if ( n > 128 ) schedule (static)
    for (i = 0 ; i < n ; ++i)
    {
        Map [i] = EMPTY ;
    }

    /* 如果矩阵不是正定的，则包含L的第一个零或负对角项的超节点s是
     * 重复的(但仅分解到有问题的对角项之前)。目的是为MATLAB提供
     * [R,p]=chol(A);L=R'的第1到p-1列是必需的，其中L(p,p)是有问题的对角项。
     * repeat_supernode标志告诉我们这是否是重复的超节点。
     * 一旦超节点s被重复，分解就终止了。 */
    repeat_supernode = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* 超节点数值分解 */
    /* ---------------------------------------------------------------------- */

    chol_base_args base_args;
    {
        base_args.Super = Super;
        base_args.Lpi = Lpi;
        base_args.Lpx = Lpx;
        base_args.Lx = Lx;
        base_args.Ls = Ls;
        base_args.Map = Map;
        base_args.stype = stype;
        base_args.Ap = Ap;
        base_args.Apacked = Apacked;
        base_args.Anz = Anz;
        base_args.Ai = Ai;
        base_args.Ax = Ax;
        base_args.Fp = Fp;
        base_args.Fpacked = Fpacked;
        base_args.Fnz = Fnz;
        base_args.Fi = Fi;
        base_args.Fx = Fx;
        base_args.beta = beta;
        base_args.repeat_supernode = repeat_supernode;
        base_args.Head = Head;
        base_args.Next = Next;
        base_args.Lpos_save = Lpos_save;
        base_args.Lpos = Lpos;
        base_args.Next_save = Next_save;
        // base_args.C = C;
        base_args.Lmaxcsize = L->maxcsize;
        base_args.one = one;
        base_args.zero = zero;
        // base_args.RelativeMap = RelativeMap;
        base_args.SuperMap = SuperMap;
        base_args.nsuper = nsuper;
        base_args.L = L;
        base_args.Common = Common;
    }

    for(i = max_depth ; i >= 0; --i)
    {
        int snode_num = Snode_num[i];
        chol_snode_args snode_args[snode_num];
        for(j = 0; j < snode_num; ++j)
        {
            int s = Snode_distribution[i][j];
            snode_args[j].snode_id = s;
            snode_args[j].base_args = &base_args;     
        }
        for(j = 0; j < snode_num; ++j)
        {
            numeric_update(&snode_args[j]);
        }
        for(j = 0; j < snode_num; ++j)
        {
        #ifdef PTHREAD_POOL
            // if(TPSM_addTask( parallel_numeric_factor, &(snode_args[j]), 1)==1)
            // {
            //     printf("TASK_BUFFER_SIZE is too small\n");
            // } 
            if( TPSM_addTask( parallel_numeric_factor, &(snode_args[j]), 1)==1)
            {
                printf("TASK_BUFFER_SIZE is too small\n");
            }
        #else
            parallel_numeric_factor(&snode_args[j]);
        #endif       
        }
        #ifdef PTHREAD_POOL
        TPSM_barrier_tag(1);
        #endif       
        for(j = 0; j < snode_num; ++j)
        {
            numeric_finish(&snode_args[j]);
        }
    }

    /* success; matrix is positive definite */
    L->minor = n ;
    free(Snode_num);
    free(Snode_distribution);
    return (Common->status >= SPARSE_OK) ;
}

#undef PATTERN
#undef REAL

/**
 * @brief   如果成功，或矩阵不正定，则返回TRUE。
 *          如果内存不足、输入无效或发生其他致命错误，则返回FALSE。
 * 
 */
int SparseChol_super_numeric
(
    /* ---- input ---- */
    //TPSM_t *tpool,
    sparse_csc *A,	/* 分解的稀疏矩阵 */
    sparse_csc *F,	/* F = A' or A(:,f)' */
    double beta [2],	/* beta*I加到被分解的稀疏矩阵的对角线上 */
    /* ---- in/out --- */
    sparse_factor *L,	/* 分解 */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *C ;
    Int *Super, *Map, *SuperMap ;
    size_t maxcsize ;
    Int nsuper, n, i, k, s, stype, nrow ;
    int ok = TRUE, symbolic ;
    size_t t, w ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    stype = A->stype ;
    if (stype < 0)
    {
        if (A->nrow != A->ncol || A->nrow != L->n)
        {
            ERROR (SPARSE_INVALID, "invalid dimensions") ;
            return (FALSE) ;
        }
    }
    else if (stype == 0)
    {
        if (A->nrow != L->n)
        {
            ERROR (SPARSE_INVALID, "invalid dimensions") ;
            return (FALSE) ;
        }
        RETURN_IF_NULL (F, FALSE) ;
        RETURN_IF_XTYPE_INVALID (F, SPARSE_REAL, SPARSE_REAL, FALSE) ;
        if (A->nrow != F->ncol || A->ncol != F->nrow || F->stype != 0)
        {
            ERROR (SPARSE_INVALID, "F invalid") ;
            return (FALSE) ;
        }
        if (A->xtype != F->xtype)
        {
            ERROR (SPARSE_INVALID, "A and F must have same xtype") ;
            return (FALSE) ;
        }
    }
    else
    {
        /* 对称上三角情况不支持 */
        ERROR (SPARSE_INVALID, "symmetric upper case not supported") ;
        return (FALSE) ;
    }
    if (!(L->is_super))
    {
        ERROR (SPARSE_INVALID, "L not supernodal") ;
        return (FALSE) ;
    }
    if (L->xtype != SPARSE_PATTERN)
    {
        if (! ((A->xtype == SPARSE_REAL && L->xtype == SPARSE_REAL)))
        {
            ERROR (SPARSE_INVALID, "complex type mismatch") ;
            return (FALSE) ;
        }
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 共同分配工作空间 */
    /* ---------------------------------------------------------------------- */

    nsuper = L->nsuper ;
    maxcsize = L->maxcsize ;
    nrow = A->nrow ;
    n = nrow ;

    /* w = 2*n + 5*nsuper */
    w = SparseCore_mult_size_t (n, 2, &ok) ;
    t = SparseCore_mult_size_t (nsuper, 5, &ok) ;
    w = SparseCore_add_size_t (w, t, &ok) ;
    if (!ok)
    {
        ERROR (SPARSE_TOO_LARGE, "problem too large") ;
        return (FALSE) ;
    }

    SparseCore_allocate_work (n, w, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	    return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 如果需要，获取当前的分解因子L并且指定数值部分 */
    /* ---------------------------------------------------------------------- */

    Super = L->super ;
    symbolic = (L->xtype == SPARSE_PATTERN) ;
    if (symbolic)
    {
        /* 通过分配L->x转换到超节点数值 */
        SparseCore_change_factor ( SPARSE_REAL, TRUE, TRUE, TRUE, TRUE, L, Common) ;
        if (Common->status < SPARSE_OK)
        {
            /* 因子L保持符号超节点形式 */
            return (FALSE) ;
        }
    }

    /* 不支持超节点LDL' */
    L->is_ll = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 获取更多的工作空间 */
    /* ---------------------------------------------------------------------- */

    C = SparseCore_allocate_dense (maxcsize, 1, maxcsize, L->xtype, Common) ;
    if (Common->status < SPARSE_OK)
    {
        int status = Common->status ;
        if (symbolic)
        {
            /* 将L更改回符号，因为数值没有初始化。这不会失败。 */
            SparseCore_change_factor (SPARSE_PATTERN, TRUE, TRUE, TRUE, TRUE,
                L, Common) ;
        }
        /* 因子L现在回到了它在输入时的形式 */
        Common->status = status ;
        return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 获取工作空间 */
    /* ---------------------------------------------------------------------- */

    SuperMap = Common->Iwork ;		/* 大小为n (i/i/l) */
    Map = Common->Flag ;            /* 大小为n, 使用Flag作为Map数组的工作空间 */

// #pragma omp parallel for num_threads(CHOLESKY_OMP_NUM_THREADS) \
//     if ( n > 128 ) schedule (static)

    for (i = 0 ; i < n ; i++)
    {
	    Map [i] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* 找出节点到松弛超节点的映射 */
    /* ---------------------------------------------------------------------- */
    // 找出节点到松弛超节点的映射

    /* 如果超节点s包含第k列，则SuperMap [k] = s */
    for (s = 0 ; s < nsuper ; s++)
    {
        for (k = Super [s] ; k < Super [s+1] ; k++)
        {
            SuperMap [k] = s ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* 超节点数值分解，使用模板例程 */
    /* ---------------------------------------------------------------------- */

    if (A->xtype == SPARSE_REAL)
    {
        ok = r_SparseChol_super_numeric ( A, F, beta, L, C, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 清除公共工作区，释放临时工作区C，并返回 */
    /* ---------------------------------------------------------------------- */

    /* 标记数组被用作工作区，清除它 */
    Common->mark = EMPTY ;
    /* SparseCore_clear_flag (Common) ; */
    SPARSE_CLEAR_FLAG (Common) ;
    SparseCore_free_dense (&C, Common) ;
    return (ok) ;
}

#endif
#endif
