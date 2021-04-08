/**
 * @file SparseCore_metis.c
 * @author your name (you@domain.com)
 * @brief 连接METIS的一个接口，引用了METIS自己的nested dissection 算法。
 *        通常比chol_nested_dissection快，主要是因为它只对分隔树的叶子
 *        使用最小度，而不是对整个矩阵使用最小度。
 * 
 *  METIS 在内存溢出时 不报错，而是会中断程序。这个接口试图避免这个
 *  问题，方法是在调用前先分配足够大的空间，然后再释放空间。
 *  以便在METIS内进行任何内存分配。但这么做也并不能确保工程总是成功。
 *  当发生这个问题的时候，增大 Common->metis_memory 。
 *  如果您不介意终止您的程序，请将Common->metis_memory设置为零(2.0通常是安全的)。
 *  在这个文件的例程中有其他几个METIS解决方案。
 *  有关更多细节，请参见下面对metis_memory_ok的描述。
 * 
 * @version 0.1
 * @date 2020-11-18
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef NMETIS

#include "Sparse_internal.h"
#include "metis.h"
#include "Sparse_partition.h"
#include "SparseChol.h"

/* ========================================================================== */
/* === metis_memory_ok ====================================================== */
/* ========================================================================== */

/* 这个例程将会分配一块内存空间，大小等同于METIS使用的内存上限，如果
 * 内存开辟失败，将不会调用 METIS. 经典的上界为 10*nz + 50*n + 4096
 */
#define GUESS(nz,n) (10 * (nz) + 50 * (n) + 4096)

static int metis_memory_ok
(
    Int n,
    Int nz,
    sparse_common *Common
)
{
    double s ;
    void *p ;
    size_t metis_guard ;

    if (Common->metis_memory <= 0)
    {
	return (TRUE) ;
    }

    n  = MAX (1, n) ;
    nz = MAX (0, nz) ;

    /* 以双精度计算，以避免整数溢出 */
    s = GUESS ((double) nz, (double) n) ;
    s *= Common->metis_memory ;

    if (s * sizeof (idx_t) >= ((double) Size_max))
    {
	/* 不使用 malloc 开启那么大的块 */
	return (FALSE) ;
    }

    /* 用 size_t 重新计算 */
    metis_guard = GUESS ((size_t) nz, (size_t) n) ;
    metis_guard *= Common->metis_memory ;

    /* 尝试开空间 */
    p = SparseCore_malloc (metis_guard, sizeof (idx_t), Common) ;
    if (p == NULL)
    {
	/* 失败-返回内存不足的情况 */
	return (FALSE) ;
    }

    /* 成功了，证明空间足够，释放掉 */
    SparseCore_free (metis_guard, sizeof (idx_t), p, Common) ;
    return (TRUE) ;
}
/* =========================================================== */
/* === SparseCore_metis_bisector ============================= */
/* =========================================================== */

/* 找到一组节点平分 A 或 AA' 的图形 ( METIS_ComputeVertexSeparator 直接接口 )
 *
 * 输入矩阵 A 必须为对称矩阵 且 没有对角线项目. 不做检查
 */

Sparse_long SparseCore_metis_bisector	/* 返回 分割的 尺寸 */
(
    /* ---- input ---- */
    sparse_csc *A,	/* 待分割矩阵 */
    Int *Anw,		/* 大小为 A->nrow, 结点权重, 可以是 NULL, */
    Int *Aew,		/* 大小为 nz, 边权重 (silently ignored). */
                        /* 这个选项在METIS 4 可用, 而在 5 不可用 */
                        /* 这个参数现在没有使用，但是为了向后兼容 */
                        /*它保留了下来，以便不更改 API 。 */
    /* ---- output --- */
    Int *Partition,	/* 大小为 A->nrow */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Ap, *Ai ;
    idx_t *Mp, *Mi, *Mnw, *Mpart ;
    Int n, nleft, nright, j, p, csep, total_weight, lightest, nz ;
    idx_t nn, csp ;
    size_t n1 ;
    int ok ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_NULL (Partition, EMPTY) ;
    if (A->stype || A->nrow != A->ncol)
    {
	/* A必须为方形矩阵, 上面和下面都有 */
	ERROR (SPARSE_INVALID, "matrix must be square, symmetric,"
		" and with both upper/lower parts present") ;
	return (EMPTY) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 快速返回 */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    if (n == 0)
    {
	return (0) ;
    }
    n1 = ((size_t) n) + 1 ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    nz = Ap [n] ;

    /* ---------------------------------------------------------------------- */
    /* 如果需要，复制Int到METIS idx_t */
    /* ---------------------------------------------------------------------- */

    if (sizeof (Int) == sizeof (idx_t))
    {
	/* 这是一个典型的例子 */
	Mi    = (idx_t *) Ai ;
	Mp    = (idx_t *) Ap ;
	Mnw   = (idx_t *) Anw ;
	Mpart = (idx_t *) Partition ;
    }
    else
    {
	/* idx_t 与 Int 不同; 将图复制到 METIS idx_t */
	Mi    = SparseCore_malloc (nz, sizeof (idx_t), Common) ;
	Mp    = SparseCore_malloc (n1, sizeof (idx_t), Common) ;
	Mnw   = Anw ? (SparseCore_malloc (n,  sizeof (idx_t), Common)) : NULL ;
	Mpart = SparseCore_malloc (n,  sizeof (idx_t), Common) ;
	if (Common->status < SPARSE_OK)
	{
	    SparseCore_free (nz, sizeof (idx_t), Mi,    Common) ;
	    SparseCore_free (n1, sizeof (idx_t), Mp,    Common) ;
	    SparseCore_free (n,  sizeof (idx_t), Mnw,   Common) ;
	    SparseCore_free (n,  sizeof (idx_t), Mpart, Common) ;
	    return (EMPTY) ;
	}
	for (p = 0 ; p < nz ; p++)
	{
	    Mi [p] = Ai [p] ;
	}
	for (j = 0 ; j <= n ; j++)
	{
	    Mp [j] = Ap [j] ;
	}
        if (Anw != NULL)
        {
            for (j = 0 ; j <  n ; j++)
            {
                Mnw [j] = Anw [j] ;
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* METIS解决方案:尽量确保METIS不会耗尽内存 */
    /* ---------------------------------------------------------------------- */

    if (!metis_memory_ok (n, nz, Common))
    {
        printf ("Metis memory error \n");
	/* METIS可能会要求太多的内存，从而终止程序 */
	if (sizeof (Int) != sizeof (idx_t))
	{
	    SparseCore_free (nz, sizeof (idx_t), Mi,    Common) ;
	    SparseCore_free (n1, sizeof (idx_t), Mp,    Common) ;
	    SparseCore_free (n,  sizeof (idx_t), Mnw,   Common) ;
	    SparseCore_free (n,  sizeof (idx_t), Mpart, Common) ;
	}
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 切分 图形 */
    /* ---------------------------------------------------------------------- */
    /* 
    METIS_ComputeVertexSeparator(
            idx_t *nvtxs,       number of nodes
            idx_t *xadj,        column pointers
            idx_t *adjncy,      row indices
            idx_t *vwgt,        vertex weights (NULL means unweighted)
            idx_t *options,     options (NULL means defaults)
            idx_t *sepsize,     separator size
            idx_t *part);       partition.  part [i] = 0,1,2, where:
                                0:left, 1:right, 2:separator
    */
    nn = n ;
    ok = METIS_ComputeVertexSeparator (&nn, Mp, Mi, Mnw, NULL, &csp, Mpart) ;
    csep = csp ;

    /* ---------------------------------------------------------------------- */
    /* 从 idx_t 复制回结果 */
    /* ---------------------------------------------------------------------- */

    if (ok == METIS_OK && (sizeof (Int) != sizeof (idx_t)))
    {
	for (j = 0 ; j < n ; j++)
	{
	    Partition [j] = Mpart [j] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 给 METIS 释放空间 */
    /* ---------------------------------------------------------------------- */

    if (sizeof (Int) != sizeof (idx_t))
    {
	SparseCore_free (nz, sizeof (idx_t), Mi,    Common) ;
	SparseCore_free (n1, sizeof (idx_t), Mp,    Common) ;
	SparseCore_free (n,  sizeof (idx_t), Mnw,   Common) ;
	SparseCore_free (n,  sizeof (idx_t), Mpart, Common) ;
    }

    if (ok == METIS_ERROR_MEMORY)
    {
        ERROR (SPARSE_OUT_OF_MEMORY, "out of memory in METIS") ;
	return (EMPTY) ;
    }
    else if (ok == METIS_ERROR_INPUT)
    {
        ERROR (SPARSE_INVALID, "invalid input to METIS") ;
	return (EMPTY) ;
    }
    else if (ok == METIS_ERROR)
    {
        ERROR (SPARSE_INVALID, "unspecified METIS error") ;
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 确保合理的分割 */
    /* ---------------------------------------------------------------------- */

    if (csep == 0)
    {
	/* 分隔符为空，选择最轻的节点作为分隔符。如果符合条件，选择编号最高的节点。 */
        if (Anw == NULL)
        {
	    lightest = n-1 ;
        }
        else
        {
	    lightest = 0 ;
            for (j = 0 ; j < n ; j++)
            {
                if (Anw [j] <= Anw [lightest])
                {
                    lightest = j ;
                }
            }
        }
	Partition [lightest] = 2 ;
	csep = (Anw ? (Anw [lightest]) : 1) ;
    }

    /* 确定图的左右部分的节点权重 */
    nleft = 0 ;
    nright = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	if (Partition [j] == 0)
	{
	    nleft += (Anw ? (Anw [j]) : 1) ;
	}
	else if (Partition [j] == 1)
	{
	    nright += (Anw ? (Anw [j]) : 1) ;
	}
    }


    total_weight = nleft + nright + csep ;

    if (csep < total_weight)
    {
	/* 分隔符小于整个图。 确保左右都空 或者 都为非空 */

	if ((nleft == 0 && nright > 0) || (nleft > 0 && nright == 0))
	{
	    /* 左或右为空，将所有的结点放到分割器中 */
	    csep = total_weight ;
	    for (j = 0 ; j < n ; j++)
	    {
		Partition [j] = 2 ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 返回分隔符中节点权重的总和 */
    /* ---------------------------------------------------------------------- */

    return (csep) ;
}


/* ========================================================================== */
/* === SparseCore_metis ======================================================== */
/* ========================================================================== */

int SparseCore_metis
(
    /* ---- input ---- */
    sparse_csc *A,	/* 排序矩阵 */
    Int *fset,		
    size_t fsize,	
    int postorder,	/* 如果为 TRUE, 使用 etree 或者 coletree 后序排序 */
    /* ---- output --- */
    Int *Perm,		/* 大小为 A->nrow, 输出重排序矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    double d ;
    Int *Iperm, *Iwork, *Bp, *Bi ;
    idx_t *Mp, *Mi, *Mperm, *Miperm ;
    sparse_csc *B ;
    Int i, j, n, nz, p, identity, uncol ;
    idx_t nn, zero = 0 ;
    size_t n1, s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Perm, FALSE) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* quick return */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    if (n == 0)
    {
	return (TRUE) ;
    }
    n1 = ((size_t) n) + 1 ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    /* s = 4*n + uncol */
    uncol = (A->stype == 0) ? A->ncol : 0 ;
    s = SparseCore_mult_size_t (n, 4, &ok) ;
    s = SparseCore_add_size_t (s, uncol, &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (n, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 把矩阵转换为邻接表形式 */
    /* ---------------------------------------------------------------------- */

    /* METIS的输入图必须是对称的，上和下三角都有，并且没有对角线条目。列不需要排序。
     * B = A+A', A*A', or A(:,f)*A(:,f)' */
    if (A->stype)
    {
	/* 通过转换为非对称模式，将上/下部分添加到对称上/下矩阵中 */
	/* 工作区: Iwork (nrow) */
	B = SparseCore_copy (A, 0, -1, Common) ;
    }
    else
    {
	/* B = A*A' or A(:,f)*A(:,f)', 没有对角项 */
	/* 工作区: Flag (nrow), Iwork (max (nrow,ncol)) */
	B = SparseCore_aat (A, fset, fsize, -1, Common) ;
    }

    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }


    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Iperm = Iwork ;		/* 大小 n (i/i/l) */

    Bp = B->p ;
    Bi = B->i ;
    nz = Bp [n] ;

    /* B不包括对角线，以及上下部分。Common->anz包括对角线部分，以及B部分的下半部分*/
    Common->anz = nz / 2 + n ;

    /* ---------------------------------------------------------------------- */
    /* 如果需要，分配METIS输入数组 */
    /* ---------------------------------------------------------------------- */

    if (sizeof (Int) == sizeof (idx_t))
    {
	Miperm = (idx_t *) Iperm ;
	Mperm  = (idx_t *) Perm ;
	Mp     = (idx_t *) Bp ;
	Mi     = (idx_t *) Bi ;
    }
    else
    {
	/* 仅当Int和idx_t不同时，为METIS分配图形 */
	Miperm = SparseCore_malloc (n,  sizeof (idx_t), Common) ;
	Mperm  = SparseCore_malloc (n,  sizeof (idx_t), Common) ;
	Mp     = SparseCore_malloc (n1, sizeof (idx_t), Common) ;
	Mi     = SparseCore_malloc (nz, sizeof (idx_t), Common) ;
	if (Common->status < SPARSE_OK)
	{
	    /* 内存越界 */
	    SparseCore_free_sparse (&B, Common) ;
	    SparseCore_free (n,  sizeof (idx_t), Miperm, Common) ;
	    SparseCore_free (n,  sizeof (idx_t), Mperm, Common) ;
	    SparseCore_free (n1, sizeof (idx_t), Mp, Common) ;
	    SparseCore_free (nz, sizeof (idx_t), Mi, Common) ;
	    return (FALSE) ;
	}
	for (j = 0 ; j <= n ; j++)
	{
	    Mp [j] = Bp [j] ;
	}
	for (p = 0 ; p < nz ; p++)
	{
	    Mi [p] = Bi [p] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* METIS 工作区 */
    /* ---------------------------------------------------------------------- */

    identity = FALSE ;
    if (nz == 0)
    {
	/* 矩阵没有非对角元素。在这种情况下METIS_NodeND失败，所以避免使用它。*/
    /* 最好的排列方式是恒等，所以这很容易解决。 */
	identity = TRUE ;
    }
    else if (Common->metis_nswitch > 0)
    {
	/* METIS 4.0.1中的METIS_NodeND给出了一个段故障，其阶矩阵n = 3005, nz = 6,036,025，包括对角项。
     * 解决方法是返回单位排列而不是使用METIS矩阵维度3000或更多的和密度的66%或更多。
	 */
	d = ((double) nz) / (((double) n) * ((double) n)) ;
	if (n > (Int) (Common->metis_nswitch) && d > Common->metis_dswitch)
	{
	    identity = TRUE ;
	}
    }

    if (!identity && !metis_memory_ok (n, nz, Common))
    {
	/* METIS 可能会申请太多内存空间 然后终止程序 */
	identity = TRUE ;
    }

    /* ---------------------------------------------------------------------- */
    /* 找寻到重排序矩阵 */
    /* ---------------------------------------------------------------------- */

    if (identity)
    {
	/* 不需要后序排序 */
	postorder = FALSE ;
	for (i = 0 ; i < n ; i++)
	{
	    Mperm [i] = i ;
	}
    }
    else
    {
        /*
        int METIS_NodeND(
            idx_t *nvtxs,       number of nodes
            idx_t *xadj,        column pointers
            idx_t *adjncy,      row indices
            idx_t *vwgt,        vertex weights (NULL means unweighted)
            idx_t *options,     options (NULL means defaults)
            idx_t *perm,        fill-reducing ordering
            idx_t *iperm);      inverse of perm
        */

	nn = n ;
	METIS_NodeND (&nn, Mp, Mi, NULL, NULL, Mperm, Miperm) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 释放METIS输入数组 */
    /* ---------------------------------------------------------------------- */

    if (sizeof (Int) != sizeof (idx_t))
    {
	for (i = 0 ; i < n ; i++)
	{
	    Perm [i] = (Int) (Mperm [i]) ;
	}
	SparseCore_free (n,   sizeof (idx_t), Miperm, Common) ;
	SparseCore_free (n,   sizeof (idx_t), Mperm, Common) ;
	SparseCore_free (n+1, sizeof (idx_t), Mp, Common) ;
	SparseCore_free (nz,  sizeof (idx_t), Mi, Common) ;
    }

    SparseCore_free_sparse (&B, Common) ;

    /* ---------------------------------------------------------------------- */
    /* etree 或者 column-etree 后序排序, 利用了 Cholesky 模块 */
    /* ---------------------------------------------------------------------- */

    if (postorder)
    {
	Int *Parent, *Post, *NewPerm ;
	Int k ;

	Parent = Iwork + 2*((size_t) n) + uncol ;   
	Post   = Parent + n ;			 

	/* 工作区: Iwork (2*nrow+uncol), Flag (nrow), Head (nrow+1) */
	SparseChol_analyze_ordering (A, SPARSE_METIS, Perm, fset, fsize,
		Parent, Post, NULL, NULL, NULL, Common) ;
	if (Common->status == SPARSE_OK)
	{
	    /* 利用后序遍历结合 METIS 转置矩阵 */
	    NewPerm = Parent ;	  
	    for (k = 0 ; k < n ; k++)
	    {
		NewPerm [k] = Perm [Post [k]] ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		Perm [k] = NewPerm [k] ;
	    }
	}
    }
    return (Common->status == SPARSE_OK) ;
}
#endif