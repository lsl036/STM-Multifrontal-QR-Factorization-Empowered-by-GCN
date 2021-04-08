/**
 * @file SparseCore_matrixops.c
 * @author your name (you@domain.com)
 * @brief C = A*A' or C = A(:,f)*A(:,f)'
 * @version 0.1
 * @date 2020-09-20
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_internal.h"
#include "SparseCore.h"

/**
 * @brief 
 * 
 */
sparse_csc *SparseCore_aat
(
    /* ---- input ---- */
    sparse_csc *A,	/* 输入矩阵，用于构造C=A*A' */
    Int *fset,			/* 0:(A->ncol)-1的子集 */
    size_t fsize,		/* 子集的大小 */
    int mode,			/* >0: numerical, 0: pattern, <0: pattern (no diag)
			 			 * -2: pattern only, no diagonal, add 50% + n extra
			 			 * space to C */
    /* --------------- */
    sparse_common *Common
)
{
    double fjt ;
    double *Ax, *Fx, *Cx, *W ;
    Int *Ap, *Anz, *Ai, *Fp, *Fi, *Cp, *Ci, *Flag ;
    sparse_csc *C, *F ;
    Int packed, j, i, pa, paend, pf, pfend, n, mark, cnz, t, p, values, diag,
	extra ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    values = (mode > 0) && (A->xtype != SPARSE_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    if (A->stype)
    {
	ERROR (SPARSE_INVALID, "matrix cannot be symmetric") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    diag = (mode >= 0) ;
    n = A->nrow ;
    SparseCore_allocate_work (n, MAX (A->ncol, A->nrow), values ? n : 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    /* 得到矩阵A */
    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    packed = A->packed ;

    /* 获取工作空间 */
    W = Common->Xwork ;		/* 大小为n, 如果值为FALSE则不使用 */
    Flag = Common->Flag ;	/* 大小为n, Flag [0..n-1] < 输入标记 */

    /* ---------------------------------------------------------------------- */
    /* F = A' or A(:,f)' */
    /* ---------------------------------------------------------------------- */

    /* 工作空间:Iwork(如果没有fset，则为nrow;如果有则为MAX (nrow,ncol)) */
    F = SparseCore_ptranspose (A, values, NULL, fset, fsize, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    Fp = F->p ;
    Fi = F->i ;
    Fx = F->x ;

    /* ---------------------------------------------------------------------- */
    /* 计算结果C中的条目数 */
    /* ---------------------------------------------------------------------- */

    cnz = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	/* 清除标记数组 */
	/* mark = SparseCore_clear_flag (Common) ; */
	SPARSE_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	/* 如果需要，排除对角线 */
	if (!diag)
	{
	    Flag [j] = mark ;
	}

	/* 对于第j列中每个非零元F(t,j)，做 */
	pfend = Fp [j+1] ;
	for (pf = Fp [j] ; pf < pfend ; pf++)
	{
	    /* F(t,j)是非零元 */
	    t = Fi [pf] ;

	    /* 将A(:，t)的非零pattern加入到C(:，j)的pattern中 */
	    pa = Ap [t] ;
	    paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
	    for ( ; pa < paend ; pa++)
	    {
		i = Ai [pa] ;
		if (Flag [i] != mark)
		{
		    Flag [i] = mark ;
		    cnz++ ;
		}
	    }
	}
	if (cnz < 0)
	{
	    break ;	    /* 整型溢出 */
	}
    }

    extra = (mode == -2) ? (cnz/2 + n) : 0 ;

    mark = SparseCore_clear_flag (Common) ;

    /* ---------------------------------------------------------------------- */
    /* 检查整型溢出 */
    /* ---------------------------------------------------------------------- */

    if (cnz < 0 || (cnz + extra) < 0)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	SparseCore_clear_flag (Common) ;
	SparseCore_free_sparse (&F, Common) ;
	return (NULL) ;	    /* 问题太大 */
    }

    /* ---------------------------------------------------------------------- */
    /* 给C分配空间 */
    /* ---------------------------------------------------------------------- */

    C = SparseCore_allocate_sparse (n, n, cnz + extra, FALSE, TRUE, 0,
	    values ? A->xtype : SPARSE_PATTERN, Common) ;
    if (Common->status < SPARSE_OK)
    {
	SparseCore_free_sparse (&F, Common) ;
	return (NULL) ;	    /* 内存溢出 */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* C = A*A' */
    /* ---------------------------------------------------------------------- */

    cnz = 0 ;

    if (values)
    {

	/* pattern and values */
	for (j = 0 ; j < n ; j++)
	{
	    /* 清楚标记数组 */
	    mark = SparseCore_clear_flag (Common) ;

	    /* 开始做C的j列 */
	    Cp [j] = cnz ;

	    /* 对于第j列中每个非零F(t,j)，做: */
	    pfend = Fp [j+1] ;
	    for (pf = Fp [j] ; pf < pfend ; pf++)
	    {
		/* F(t,j) is nonzero */
		t = Fi [pf] ;
		fjt = Fx [pf] ;

		/* 将A(:，t)的非零pattern添加到C(:，j)的pattern中，并将值分散到W中 */
		pa = Ap [t] ;
		paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
		for ( ; pa < paend ; pa++)
		{
		    i = Ai [pa] ;
		    if (Flag [i] != mark)
		    {
			Flag [i] = mark ;
			Ci [cnz++] = i ;
		    }
		    W [i] += Ax [pa] * fjt ;
		}
	    }

	    /* 将值收集到C(:，j)中 */
	    for (p = Cp [j] ; p < cnz ; p++)
	    {
		i = Ci [p] ;
		Cx [p] = W [i] ;
		W [i] = 0 ;
	    }
	}

    }
    else
    {

	/* pattern only */
	for (j = 0 ; j < n ; j++)
	{
	    /* 清楚标记数组 */
	    mark = SparseCore_clear_flag (Common) ;

	    /* 如果需要，排除对角线 */
	    if (!diag)
	    {
		Flag [j] = mark ;
	    }

	    /* 开始做C的第j列 */
	    Cp [j] = cnz ;

	    /* 对于第j列中每个非零F(t,j)，做: */
	    pfend = Fp [j+1] ;
	    for (pf = Fp [j] ; pf < pfend ; pf++)
	    {
		/* F(t,j)是非零元 */
		t = Fi [pf] ;

		/* 将A(:，t)的非零pattern加入到C(:，j)的pattern中 */
		pa = Ap [t] ;
		paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
		for ( ; pa < paend ; pa++)
		{
		    i = Ai [pa] ;
		    if (Flag [i] != mark)
		    {
			Flag [i] = mark ;
			Ci [cnz++] = i ;
		    }
		}
	    }
	}
    }

    Cp [n] = cnz ;

    /* ---------------------------------------------------------------------- */
    /* 清除工作区和释放临时矩阵并返回结果 */
    /* ---------------------------------------------------------------------- */

    SparseCore_free_sparse (&F, Common) ;
    SparseCore_clear_flag (Common) ;
    return (C) ;
}


/**
 * @brief C = alpha*A + beta*B, or spones(A+B)
 * 
 */
sparse_csc *SparseCore_add
(
    /* ---- input ---- */
    sparse_csc *A,	    /* 相加的矩阵A */
    sparse_csc *B,	    /* 相加的矩阵B  */
    double alpha [2],	    /* A的标量因子 */
    double beta [2],	    /* B的标量因子 */
    int values,		    	/* 如果为真，则计算C的数值 */
    int sorted,		    	/* 如果为真，对C的列进行排序 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Ax, *Bx, *Cx, *W ;
    Int apacked, up, lo, nrow, ncol, bpacked, nzmax, pa, paend, pb, pbend, i,
	j, p, mark, nz ;
    Int *Ap, *Ai, *Anz, *Bp, *Bi, *Bnz, *Flag, *Cp, *Ci ;
    sparse_csc *A2, *B2, *C ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_NULL (B, NULL) ;
    values = values &&
	(A->xtype != SPARSE_PATTERN) && (B->xtype != SPARSE_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    RETURN_IF_XTYPE_INVALID (B, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    if (A->nrow != B->nrow || A->ncol != B->ncol)
    {
	/* A和B必须有相同的维度 */
	ERROR (SPARSE_INVALID, "A and B dimesions do not match") ;
	return (NULL) ;
    }
    /* 如果值为真，A和B必须具有相同的数值类型(都必须是SPARSE_REAL，这在上面已经隐式检查过了) */

    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;
    SparseCore_allocate_work (nrow, MAX (nrow,ncol), values ? nrow : 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    if (nrow <= 1)
    {
	/* C是隐式排序的，所以这里不需要排序 */
	sorted = FALSE ;
    }

    /* 如有必要，将A或B转换为不对称 */
    A2 = NULL ;
    B2 = NULL ;

    if (A->stype != B->stype)
    {
	if (A->stype)
	{
	    /* 工作空间: Iwork (max (nrow,ncol)) */
	    A2 = SparseCore_copy (A, 0, values, Common) ;
	    if (Common->status < SPARSE_OK)
	    {
		return (NULL) ;	    /* 内存溢出 */
	    }
	    A = A2 ;
	}
	if (B->stype)
	{
	    /* 工作空间: Iwork (max (nrow,ncol)) */
	    B2 = SparseCore_copy (B, 0, values, Common) ;
	    if (Common->status < SPARSE_OK)
	    {
		SparseCore_free_sparse (&A2, Common) ;
		return (NULL) ;	    /* 内存溢出 */
	    }
	    B = B2 ;
	}
    }

    /* 得到矩阵A */
    up = (A->stype > 0) ;
    lo = (A->stype < 0) ;

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    apacked = A->packed ;

    /* 得到矩阵B */
    Bp  = B->p ;
    Bnz = B->nz ;
    Bi  = B->i ;
    Bx  = B->x ;
    bpacked = B->packed ;

    /* 得到工作空间 */
    W = Common->Xwork ;	    /* 值为真时使用，大小为nrow */
    Flag = Common->Flag ;   /* Flag [0..nrow-1] < 输入标记，大小为nrow */

    /* ---------------------------------------------------------------------- */
    /* 给结果矩阵C分配空间 */
    /* ---------------------------------------------------------------------- */

    /* 如果发生整数溢出，nzmax < 0和分配失败正确(同样在大多数其他矩阵操作例程). */

    nzmax = SparseCore_nnz (A, Common) + SparseCore_nnz (B, Common) ;

    C = SparseCore_allocate_sparse (nrow, ncol, nzmax, FALSE, TRUE,
	    SIGN (A->stype), values ? A->xtype : SPARSE_PATTERN, Common) ;
    if (Common->status < SPARSE_OK)
    {
	SparseCore_free_sparse (&A2, Common) ;
	SparseCore_free_sparse (&B2, Common) ;
	return (NULL) ;	    /* 内存溢出 */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* 计算C = alpha*A + beta*B */
    /* ---------------------------------------------------------------------- */

    nz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	Cp [j] = nz ;

	/* 清楚标记 */
	/* mark = SparseCore_clear_flag (Common) ; */
	SPARSE_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	/* 将B分散到W中 */
	pb = Bp [j] ;
	pbend = (bpacked) ? (Bp [j+1]) : (pb + Bnz [j]) ;
	for (p = pb ; p < pbend ; p++)
	{
	    i = Bi [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    Flag [i] = mark ;
	    if (values)
	    {
		W [i] = beta [0] * Bx [p] ;
	    }
	}

	/* 加上A，从W收集到C(:，j) */
	pa = Ap [j] ;
	paend = (apacked) ? (Ap [j+1]) : (pa + Anz [j]) ;
	for (p = pa ; p < paend ; p++)
	{
	    i = Ai [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    Flag [i] = EMPTY ;
	    Ci [nz] = i ;
	    if (values)
	    {
		Cx [nz] = W [i] + alpha [0] * Ax [p] ;
		W [i] = 0 ;
	    }
	    nz++ ;
	}

	/* 使用B的pattern将剩余的条目收集到C(:，j)中 */
	for (p = pb ; p < pbend ; p++)
	{
	    i = Bi [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    if (Flag [i] == mark)
	    {
		Ci [nz] = i ;
		if (values)
		{
		    Cx [nz] = W [i] ;
		    W [i] = 0 ;
		}
		nz++ ;
	    }
	}
    }

    Cp [ncol] = nz ;

    /* ---------------------------------------------------------------------- */
    /* 减少C的大小和自由临时矩阵 */
    /* ---------------------------------------------------------------------- */

    SparseCore_reallocate_sparse (nz, C, Common) ;

    /* 清楚标记数组 */
    mark = SparseCore_clear_flag (Common) ;

    SparseCore_free_sparse (&A2, Common) ;
    SparseCore_free_sparse (&B2, Common) ;

    /* ---------------------------------------------------------------------- */
    /* 给C排序 */
    /* ---------------------------------------------------------------------- */

    if (sorted)
    {
	/* 工作空间: Iwork (max (nrow,ncol)) */
	if (!SparseCore_sort (C, Common))
	{
	    SparseCore_free_sparse (&C, Common) ;
	    if (Common->status < SPARSE_OK)
	    {
		return (NULL) ;		/* 内存溢出 */
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */
    return (C) ;
}


/**
 * @brief C = tril (triu (A,k1), k2)
 * 
 */
static sparse_csc *band		/* 返回C，如果失败返回NULL */
(
    /* ---- input，inplace=TRUE时为in/out --- */
    sparse_csc *A,
    /* ---- input ---- */
    Sparse_long k1,    	/* 忽略第k1条对角线以下的条目 */
    Sparse_long k2,    	/* 忽略第k2条对角线上的条目 */
    int mode,	    		/* >0: numerical, 0: pattern, <0: pattern (no diagonal) */
    int inplace,    		/* TRUE时：在适当的位置转换A */
    /* --------------- */
    sparse_common *Common
)
{
    double *Ax, *Cx ;
    Int packed, nz, j, p, pend, i, ncol, nrow, jlo, jhi, ilo, ihi, sorted,
	values, diag ;
    Int *Ap, *Anz, *Ai, *Cp, *Ci ;
    sparse_csc *C ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入*/
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    values = (mode > 0) && (A->xtype != SPARSE_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    packed = A->packed ;
    diag = (mode >= 0) ;
    if (inplace && !packed)
    {
	/* 不能在适当的位置对未填充矩阵进行操作 */
	ERROR (SPARSE_INVALID, "cannot operate on unpacked matrix in-place") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    sorted = A->sorted ;


    if (A->stype > 0)
    {
	/* 忽略A的严格下三角形部分中的任何元素 */
	k1 = MAX (k1, 0) ;
    }
    if (A->stype < 0)
    {
	/* 忽略A严格上三角形部分中的任何元素 */
	k2 = MIN (k2, 0) ;
    }
    ncol = A->ncol ;
    nrow = A->nrow ;

    /* 确保k1和k2在-nrow到+ncol的范围内，以避免k1和k2很大时可能出现的整数溢出 */
    k1 = MAX (-nrow, k1) ;
    k1 = MIN (k1, ncol) ;
    k2 = MAX (-nrow, k2) ;
    k2 = MIN (k2, ncol) ;

    /* 考虑从jlo到jhi列。此范围之外的列为空 */
    jlo = MAX (k1, 0) ;
    jhi = MIN (k2+nrow, ncol) ;

    if (k1 > k2)
    {
	jlo = ncol ;
	jhi = ncol ;
    }

    /* ---------------------------------------------------------------------- */
    /* 为C分配空间，或在适当的位置操作A */
    /* ---------------------------------------------------------------------- */

    if (inplace)
    {
	/* 原矩阵上转换A */
	C = A ;
    }
    else
    {
	/* 计算结果C中的条目数 */
	nz = 0 ;
	if (sorted)
	{
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i > ihi)
		    {
			break ;
		    }
		    if (i >= ilo && (diag || i != j))
		    {
			nz++ ;
		    }
		}
	    }
	}
	else
	{
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i >= ilo && i <= ihi && (diag || i != j))
		    {
			nz++ ;
		    }
		}
	    }
	}
	/* 给C分配空间，A不会被修改，如果A是有序的，C也是有序的 */
	C = SparseCore_allocate_sparse (A->nrow, ncol, nz, sorted, TRUE,
		A->stype, values ? A->xtype : SPARSE_PATTERN, Common) ;
	if (Common->status < SPARSE_OK)
	{
	    return (NULL) ;	/* 内存溢出 */
	}
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* 构造C */
    /* ---------------------------------------------------------------------- */

    /* 0到jlo-1列为空 */
    for (j = 0 ; j < jlo ; j++)
    {
	Cp [j] = 0 ;
    }

    nz = 0 ;
    if (sorted)
    {
	if (values)
	{
	    /* pattern and values */
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Cp [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i > ihi)
		    {
			break ;
		    }
		    if (i >= ilo)
		    {
			Ci [nz] = i ;
			Cx [nz] = Ax [p] ;
			nz++ ;
		    }
		}
	    }
	}
	else
	{
	    /* pattern only, 也许没有对角线 */
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Cp [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i > ihi)
		    {
			break ;
		    }
		    if (i >= ilo && (diag || i != j))
		    {
			Ci [nz++] = i ;
		    }
		}
	    }
	}
    }
    else
    {
	if (values)
	{
	    /* pattern and values */
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Cp [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i >= ilo && i <= ihi)
		    {
			Ci [nz] = i ;
			Cx [nz] = Ax [p] ;
			nz++ ;
		    }
		}
	    }
	}
	else
	{
	    /* pattern only,也许没有对角线 */
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Cp [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i >= ilo && i <= ihi && (diag || i != j))
		    {
			Ci [nz++] = i ;
		    }
		}
	    }
	}
    }

    /* jhi到ncol-1列为空 */
    for (j = jhi ; j <= ncol ; j++)
    {
	Cp [j] = nz ;
    }

    /* ---------------------------------------------------------------------- */
    /* 如果在适当的地方减小A的大小 */
    /* ---------------------------------------------------------------------- */

    if (inplace)
    {
	/* 释放A中未使用的部分，并减小A->i和A->x的大小 */
	SparseCore_reallocate_sparse (nz, A, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果矩阵C */
    /* ---------------------------------------------------------------------- */
    return (C) ;
}

/**
 * @brief band
 * 
 */
sparse_csc *SparseCore_band
(
    /* ---- input ---- */
    sparse_csc *A,		/* 从矩阵中提取带矩阵 */
    Sparse_long k1,    	/* 忽略第k1条对角线以下的条目 */
    Sparse_long k2,    	/* 忽略第k2条对角线上的条目 */
    int mode,				/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    /* --------------- */
    sparse_common *Common
)
{
    return (band (A, k1, k2, mode, FALSE, Common)) ;
}

/**
 * @brief band_inplace
 * 
 */
int SparseCore_band_inplace
(
    /* ---- input ---- */
    Sparse_long k1,    /* 忽略第k1条对角线以下的条目*/
    Sparse_long k2,    /* 忽略第k2条对角线上的条目 */
    int mode,			  /* >0: numerical, 0: pattern, <0: pattern (no diag) */
    /* ---- in/out --- */
    sparse_csc *A,	  /* 去掉不属于带的项的矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    return (band (A, k1, k2, mode, TRUE, Common) != NULL) ;
}


/**
 * @brief 	构造对称稀疏矩阵的非对称副本。当A是对称的时候，
 * 			它执行C = SparseCore_copy (A, 0, mode, Common)的工作。
 * 			在这种情况下，可以向C中添加额外的空间。
 * 
 */
static sparse_csc *copy_sym_to_unsym
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要拷贝的矩阵 */
    int mode,			/* >0: numerical, 0: pattern, <0: pattern (no diag)
			 			 * -2: pattern only, 没有对角线，加50% + n 的空间给C */
    /* --------------- */
    sparse_common *Common
)
{
    double aij ;
    double *Ax, *Cx ;
    Int *Ap, *Ai, *Anz, *Cp, *Ci, *Wj, *Iwork ;
    sparse_csc *C ;
    Int nrow, ncol, nz, packed, j, p, pend, i, pc, up, lo, values, diag,
	astype, extra ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;
    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    packed = A->packed ;
    values = (mode > 0) && (A->xtype != SPARSE_PATTERN) ;
    diag = (mode >= 0) ;

    astype = SIGN (A->stype) ;
    up = (astype > 0) ;
    lo = (astype < 0) ;

    /* ---------------------------------------------------------------------- */
    /* 创建一个对称矩阵的非对称副本 */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wj = Iwork ;		    /* 大小为ncol (i/i/l) */

    /* 
     * 用于将一个对称/下矩阵转换为不对称矩阵:
     *	L = tril (A) ;
     *	U = triu (L',1) ;
     *	C = L + U ;
     */

    /* 计算C的每一列的条目数 */
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = 0 ;
    }
    for (j = 0 ; j < ncol ; j++)
    {
	p = Ap [j] ;
	pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	for ( ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if (i == j)
	    {
		/* 对角线上的条目A(i,i)将只出现一次(除非它被模式< 0排除) */
		if (diag)
		{
		    Wj [j]++ ;
		}
	    }
	    else if ((up && i < j) || (lo && i > j))
	    {
		/* 大写:A(i,j)在严格的上部;A(j,i)将被加到c的严格较低的部分，小写则相反。 */
		Wj [j]++ ;
		Wj [i]++ ;
	    }
	}
    }
    nz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	nz += Wj [j] ;
    }

    extra = (mode == -2) ? (nz/2 + ncol) : 0 ;

    /* 给C分配空间.只有A有序，C才有序  */
    C = SparseCore_allocate_sparse (nrow, ncol, nz + extra, A->sorted, TRUE, 0,
	    values ? A->xtype : SPARSE_PATTERN, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* 为C构造列指针 */
    p = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	Cp [j] = p ;
	p += Wj [j] ;
    }
    Cp [ncol] = p ;
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = Cp [j] ;
    }

    /* 构造 C */
    if (values)
    {

	/* pattern and values */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		aij = Ax [p] ;
		if (i == j)
		{
		    /* 将对角项A(i,i)加到第i列 */
		    pc = Wj [i]++ ;
		    Ci [pc] = i ;
		    Cx [pc] = aij ;
		}
		else if ((up && i < j) || (lo && i > j))
		{
		    /* 将A(i,j)加到第j列 */
		    pc = Wj [j]++ ;
		    Ci [pc] = i ;
		    Cx [pc] = aij ;
		    /* 将A(j,i)加到第i列 */
		    pc = Wj [i]++ ;
		    Ci [pc] = j ;
		    Cx [pc] = aij ;
		}
	    }
	}

    }
    else
    {

	/* pattern only, 可能不包括对角线 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i == j)
		{
		    /* 将对角项A(i,i)添加到第i列(除非模式< 0排除) */
		    if (diag)
		    {
			Ci [Wj [i]++] = i ;
		    }
		}
		else if ((up && i < j) || (lo && i > j))
		{
		    /* 将A(i,j)加到第j列 */
		    Ci [Wj [j]++] = i ;
		    /* 将A(j,i)加到第i列 */
		    Ci [Wj [i]++] = j ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */
    return (C) ;
}

/**
 * @brief copy
 * 
 */
sparse_csc *SparseCore_copy
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要拷贝的矩阵 */
    int stype,			/* C的数据类型 */
    int mode,			/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *C ;
    Int nrow, ncol, up, lo, values, diag, astype ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    values = (mode > 0) && (A->xtype != SPARSE_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if ((stype || A->stype) && nrow != ncol)
    {
	/* 非法输入 */
	ERROR (SPARSE_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配内存空间 */
    /* ---------------------------------------------------------------------- */

    SparseCore_allocate_work (0, MAX (nrow,ncol), 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	/* 内存溢出 */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    diag = (mode >= 0) ;
    astype = SIGN (A->stype) ;
    stype = SIGN (stype) ;
    up = (astype > 0) ;
    lo = (astype < 0) ;

    /* ---------------------------------------------------------------------- */
    /* 拷贝矩阵 */
    /* ---------------------------------------------------------------------- */

    if (astype == stype)
    {

	/* ------------------------------------------------------------------ */
	/* A和C的对称性相同 */
	/* ------------------------------------------------------------------ */

	/* 将A复制到C中，保持相同的对称。如果A是对称的，则A中被忽略部分的项不会被复制到C中 */
	C = SparseCore_band (A, -nrow, ncol, mode, Common) ;

    }
    else if (!astype)
    {

	/* ------------------------------------------------------------------ */
	/* 将非对称矩阵A转换为对称矩阵C */
	/* ------------------------------------------------------------------ */

	if (stype > 0)
	{
	    /* C = triu (A) */
	    C = SparseCore_band (A, 0, ncol, mode, Common) ;
	}
	else
	{
	    /* C = tril (A) */
	    C = SparseCore_band (A, -nrow, 0, mode, Common) ;
	}
	if (Common->status < SPARSE_OK)
	{
	    /* 内存溢出 */
	    return (NULL) ;
	}
	C->stype = stype ;

    }
    else if (astype == -stype)
    {

	/* ------------------------------------------------------------------ */
	/* 对称矩阵的转置 */
	/* ------------------------------------------------------------------ */

	/* 从上到下或从下到上的转换 */
	/* 工作空间: Iwork (nrow) */
	C = SparseCore_transpose (A, values, Common) ;
	if (!diag)
	{
	    /* 如果需要，去掉对角线 */
	    SparseCore_band_inplace (-nrow, ncol, -1, C, Common) ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* 创建一个对称矩阵的非对称副本 */
	/* ------------------------------------------------------------------ */

	C = copy_sym_to_unsym (A, mode, Common) ;
    }

    if (Common->status < SPARSE_OK)
    {
	/* 内存溢出 */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */
    return (C) ;
}

/* ========================================================================== */
/* === 模板 ============================================================= */
/* ========================================================================== */

#define PATTERN
#include "SparseCore_t_transpose.c"
#define REAL
#include "SparseCore_t_transpose.c"


/**
 * @brief 	计算F = A'， A (:， F)'，或A (p, F)'，其中A不对称且F已经分配。
 * 			参见SparseCore_transpose获得一个更简单的例程。
 * 
 */
int SparseCore_transpose_unsym
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要转置的矩阵 */
    int values,			/* 1: 数组转置, 0: 不转置数值 */
    Int *Perm,			/* 大小nrow，如果存在(可以为NULL) */
    Int *fset,			/* 0:(A->ncol)-1的子集 */
    size_t fsize,		/* 子集的大小 */
    /* ---- output --- */
    sparse_csc *F,	/* F = A', A(:,f)', or A(p,f)' */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Fp, *Fnz, *Ap, *Ai, *Anz, *Wi ;
    Int nrow, ncol, permute, use_fset, Apacked, Fpacked, p, pend,
	i, j, k, Fsorted, nf, jj, jlast ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (F, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (F, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    if (A->nrow != F->ncol || A->ncol != F->nrow)
    {
	ERROR (SPARSE_INVALID, "F has the wrong dimensions") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    nf = fsize ;
    use_fset = (fset != NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;

    Ap = A->p ;		/* A的列指针，大小为A->ncol+1 */
    Ai = A->i ;		/* A的行索引，大小为nz = Ap [A->ncol] */
    Anz = A->nz ;
    Apacked = A->packed ;
    permute = (Perm != NULL) ;

    Fp = F->p ;		/* F的行指针，大小为A->nrow+1 */
    Fnz = F->nz ;
    Fpacked = F->packed ;

    nf = (use_fset) ? nf : ncol ;

    /* ---------------------------------------------------------------------- */
    /* 分配内存空间 */
    /* ---------------------------------------------------------------------- */

    /* s = nrow + ((fset != NULL) ? ncol : 0) */
    s = SparseCore_add_size_t (nrow, ((fset != NULL) ? ncol : 0), &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (0, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;	/* 内存溢出 */
    }

    Wi = Common->Iwork ;	/* 大小为nrow (i/l/l) */

    /* ---------------------------------------------------------------------- */
    /* 检查Perm和fset */
    /* ---------------------------------------------------------------------- */

    if (permute)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = 1 ;
	}
	for (k = 0 ; k < nrow ; k++)
	{
	    i = Perm [k] ;
	    if (i < 0 || i > nrow || Wi [i] == 0)
	    {
		ERROR (SPARSE_INVALID, "invalid permutation") ;
		return (FALSE) ;
	    }
	    Wi [i] = 0 ;
	}
    }

    if (use_fset)
    {
	for (j = 0 ; j < ncol ; j++)
	{
	    Wi [j] = 1 ;
	}
	for (k = 0 ; k < nf ; k++)
	{
	    j = fset [k] ;
	    if (j < 0 || j > ncol || Wi [j] == 0)
	    {
		ERROR (SPARSE_INVALID, "invalid fset") ;
		return (FALSE) ;
	    }
	    Wi [j] = 0 ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 计算A或A(:，f)每一行的条目 */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nrow ; i++)
    {
	Wi [i] = 0 ;
    }

    jlast = EMPTY ;
    Fsorted = TRUE ;

    if (use_fset)
    {
	/* 计算A(:，f)每行中的条目 */
	for (jj = 0 ; jj < nf ; jj++)
	{
	    j = fset [jj] ;
	    if (j <= jlast)
	    {
		Fsorted = FALSE ;
	    }
	    p = Ap [j] ;
	    pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Wi [Ai [p]]++ ;
	    }
	    jlast = j ;
	}

	/* 如果F是拆开的，则保存nz的值,并重新计算所有的A */
	if (!Fpacked)
	{
	    if (permute)
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [Perm [i]] ;
		}
	    }
	    else
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [i] ;
		}
	    }
	    for (i = 0 ; i < nrow ; i++)
	    {
		Wi [i] = 0 ;
	    }

	    /* 计算A的每一行中的条目 */
	    for (j = 0 ; j < ncol ; j++)
	    {
		p = Ap [j] ;
		pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    Wi [Ai [p]]++ ;
		}
	    }
	}

    }
    else
    {

	/* 计算A的每一行中的条目 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Wi [Ai [p]]++ ;
	    }
	}

	/* 如果F是拆开的，则保存nz的值 */
	if (!Fpacked)
	{
	    if (permute)
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [Perm [i]] ;
		}
	    }
	    else
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [i] ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 计算行指针 */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    if (permute)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Fp [i] = p ;
	    p += Wi [Perm [i]] ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [Perm [i]] = Fp [i] ;
	}
    }
    else
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Fp [i] = p ;
	    p += Wi [i] ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = Fp [i] ;
	}
    }
    Fp [nrow] = p ;

    if (p > (Int) (F->nzmax))
    {
	ERROR (SPARSE_INVALID, "F is too small") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 转置矩阵，使用模板例程 */
    /* ---------------------------------------------------------------------- */

    ok = FALSE ;
    if (values == 0 || F->xtype == SPARSE_PATTERN)
    {
	ok = p_SparseCore_transpose_unsym (A, Perm, fset, nf, F, Common) ;
    }
    else if (F->xtype == SPARSE_REAL)
    {
	ok = r_SparseCore_transpose_unsym (A, Perm, fset, nf, F, Common) ;
    }
    /* ---------------------------------------------------------------------- */
    /* 完成结果F */
    /* ---------------------------------------------------------------------- */

    if (ok)
    {
	F->sorted = Fsorted ;
    }
    return (ok) ;
}


/**
 * @brief 计算F = A'或A (p,p)'，其中A是对称的，F已经分配。参见SparseCore_transpose获得一个更简单的例程。
 * 
 */
int SparseCore_transpose_sym
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要转置的矩阵 */
    int values,			/* 2: 1: 数组转置, 0: 数值不转置 */
    Int *Perm,			/* 大小nrow，如果存在(可以为NULL) */
    /* ---- output --- */
    sparse_csc *F,	/* F = A' or A(p,p)' */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Ap, *Anz, *Ai, *Fp, *Wi, *Pinv, *Iwork ;
    Int p, pend, packed, upper, permute, jold, n, i, j, k, iold ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (F, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (F, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    if (A->nrow != A->ncol || A->stype == 0)
    {
		/* 这个例程只处理对称方阵 */
		ERROR (SPARSE_INVALID, "matrix must be symmetric") ;
		return (FALSE) ;
    }
    if (A->nrow != F->ncol || A->ncol != F->nrow)
    {
		ERROR (SPARSE_INVALID, "F has the wrong dimensions") ;
		return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    permute = (Perm != NULL) ;
    n = A->nrow ;
    Ap = A->p ;		/* A的列指针，大小为A->ncol+1 */
    Ai = A->i ;		/* A的行索引，大小为nz = Ap [A->ncol] */
    Anz = A->nz ;
    packed = A->packed ;
    upper = (A->stype > 0) ;

    Fp = F->p ;		/* F的行指针，大小为A->nrow+1 */

    /* ---------------------------------------------------------------------- */
    /* 分配内存空间 */
    /* ---------------------------------------------------------------------- */

    /* s = (Perm != NULL) ? 2*n : n */
    s = SparseCore_add_size_t (n, ((Perm != NULL) ? n : 0), &ok) ;
    if (!ok)
    {
		ERROR (SPARSE_TOO_LARGE, "problem too large") ;
		return (FALSE) ;
    }

    SparseCore_allocate_work (0, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
		return (FALSE) ;	/* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wi   = Iwork ;	    /* 大小为n (i/l/l) */
    Pinv = Iwork + n ;	    /* 大小为n (i/i/l) , 如果Perm为空，则未使用 */

    /* ---------------------------------------------------------------------- */
    /* 检查Perm并构造逆置换 */
    /* ---------------------------------------------------------------------- */

    if (permute)
    {
	for (i = 0 ; i < n ; i++)
	{
	    Pinv [i] = EMPTY ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    i = Perm [k] ;
	    if (i < 0 || i > n || Pinv [i] != EMPTY)
	    {
			ERROR (SPARSE_INVALID, "invalid permutation") ;
			return (FALSE) ;
	    }
	    Pinv [i] = k ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 计算F中每一行的条目数 */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
		Wi [i] = 0 ;
    }

    if (packed)
    {
		if (permute)
		{
			if (upper)
			{
				/* packed, permuted, upper */
				for (j = 0 ; j < n ; j++)
				{
					jold = Perm [j] ;
					pend = Ap [jold+1] ;
					for (p = Ap [jold] ; p < pend ; p++)
					{
						iold = Ai [p] ;
						if (iold <= jold)
						{
							i = Pinv [iold] ;
							Wi [MIN (i, j)]++ ;
						}
					}
				}
			}
			else
			{
				/* packed, permuted, lower */
				for (j = 0 ; j < n ; j++)
				{
					jold = Perm [j] ;
					pend = Ap [jold+1] ;
					for (p = Ap [jold] ; p < pend ; p++)
					{
						iold = Ai [p] ;
						if (iold >= jold)
						{
							i = Pinv [iold] ;
							Wi [MAX (i, j)]++ ;
						}
					}
				}
			}
		}
		else
		{
			if (upper)
			{
				/* packed, unpermuted, upper */
				for (j = 0 ; j < n ; j++)
				{
					pend = Ap [j+1] ;
					for (p = Ap [j] ; p < pend ; p++)
					{
						i = Ai [p] ;
						if (i <= j)
						{
							Wi [i]++ ;
						}
					}
				}
			}
			else
			{
				/* packed, unpermuted, lower */
				for (j = 0 ; j < n ; j++)
				{
					pend = Ap [j+1] ;
					for (p = Ap [j] ; p < pend ; p++)
					{
						i = Ai [p] ;
						if (i >= j)
						{
							Wi [i]++ ;
						}
					}
				}
			}
		}	
	}
    else
    {
	if (permute)
	{
	    if (upper)
	    {
		/* unpacked, permuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    p = Ap [jold] ;
		    pend = p + Anz [jold] ;
		    for ( ; p < pend ; p++)
		    {
			iold = Ai [p] ;
			if (iold <= jold)
			{
			    i = Pinv [iold] ;
			    Wi [MIN (i, j)]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* unpacked, permuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    p = Ap [jold] ;
		    pend = p + Anz [jold] ;
		    for ( ; p < pend ; p++)
		    {
			iold = Ai [p] ;
			if (iold >= jold)
			{
			    i = Pinv [iold] ;
			    Wi [MAX (i, j)]++ ;
			}
		    }
		}
	    }
	}
	else
	{
	    if (upper)
	    {
		/* unpacked, unpermuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			if (i <= j)
			{
			    Wi [i]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* unpacked, unpermuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			if (i >= j)
			{
			    Wi [i]++ ;
			}
		    }
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 计算行指针 */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (i = 0 ; i < n ; i++)
    {
		Fp [i] = p ;
		p += Wi [i] ;
    }
    Fp [n] = p ;
    for (i = 0 ; i < n ; i++)
    {
		Wi [i] = Fp [i] ;
    }

    if (p > (Int) (F->nzmax))
    {
		ERROR (SPARSE_INVALID, "F is too small") ;
		return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 转置矩阵，使用模板例程 */
    /* ---------------------------------------------------------------------- */

    ok = FALSE ;
    if (values == 0 || F->xtype == SPARSE_PATTERN)
    {
	ok = p_SparseCore_transpose_sym (A, Perm, F, Common) ;
    }
    else if (F->xtype == SPARSE_REAL)
    {
	ok = r_SparseCore_transpose_sym (A, Perm, F, Common) ;
    }
    /* ---------------------------------------------------------------------- */
    /* 完成结果F */
    /* ---------------------------------------------------------------------- */

    /* 如果没有置换向量F就是有序的 */
    if (ok)
    {
		F->sorted = !permute ;
		F->packed = TRUE ;
		F->stype = - SIGN (A->stype) ;	/* 转换数据类型 */
    }
    return (ok) ;
}

/* Returns A'.  See also SparseCore_ptranspose below. */

/**
 * @brief 返回A' 可看SparseCore_ptranspose
 * 
 */
sparse_csc *SparseCore_transpose
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要转置的矩阵 */
    int values,			/* 1: 数组转置, 0: 不转置数值(以HNUCHOL模式返回结果) */
    /* --------------- */
    sparse_common *Common
)
{
    return (SparseCore_ptranspose (A, values, NULL, NULL, 0, Common)) ;
}

/**
 * @brief 如果A是对称的，返回A'或A(p,p)'。如果A不对称，返回A'， A(:，f)'，或A(p,f)'。
 * 
 */
sparse_csc *SparseCore_ptranspose
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要转置的矩阵 */
    int values,			/* 1: 数组转置, 0: 不转置数值 */
    Int *Perm,			/* 如果非空 F = A(p,f) or A(p,p) */
    Int *fset,			/* 0:(A->ncol)-1的子集 */
    size_t fsize,		/* 子集的大小 */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Ap, *Anz ;
    sparse_csc *F ;
    Int nrow, ncol, use_fset, j, jj, fnz, packed, stype, nf, xtype ;
    size_t ineed ;
    int ok = TRUE ;

    nf = fsize ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    stype = A->stype ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配内存空间 */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;

    if (stype != 0)
    {
		use_fset = FALSE ;
		if (Perm != NULL)
		{
			ineed = SparseCore_mult_size_t (A->nrow, 2, &ok) ;
		}
		else
		{
			ineed = A->nrow ;
		}
    }
    else
    {
		use_fset = (fset != NULL) ;
		if (use_fset)
		{
			ineed = MAX (A->nrow, A->ncol) ;
		}
		else
		{
			ineed = A->nrow ;
		}
    }

    if (!ok)
    {
		ERROR (SPARSE_TOO_LARGE, "problem too large") ;
		return (NULL) ;
    }

    SparseCore_allocate_work (0, ineed, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
		return (NULL) ;	    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Anz = A->nz ;
    packed = A->packed ;
    xtype = values ? A->xtype : SPARSE_PATTERN ;

    /* ---------------------------------------------------------------------- */
    /* 给F分配空间 */
    /* ---------------------------------------------------------------------- */

    /* 确定F中非零的# */
    if (stype != 0)
    {
		/* F=A' or F=A(p,p)' 不管fset */
		fnz = SparseCore_nnz (A, Common) ;
    }
    else
    {
		nf = (use_fset) ? nf : ncol ;
		if (use_fset)
		{
			fnz = 0 ;
			/* F=A(:,f)' or F=A(p,f)' */
			for (jj = 0 ; jj < nf ; jj++)
			{
				/* 
				 * fset尚未被检查;它将在SparseCore_transpose_unsym中被彻底检查。
				 * 现在，只要确保我们没有越界访问Ap和Anz。
				 */
				j = fset [jj] ;
				if (j >= 0 && j < ncol)
				{
					fnz += packed ? (Ap [j+1] - Ap [j]) : MAX (0, Anz [j]) ;
				}
			}
		}
		else
		{
			/* F=A' or F=A(p,:)' */
			fnz = SparseCore_nnz (A, Common) ;
		}
    }

	/*
	 * F是ncole*nrow大小的, fnz个非零元，只有F存在、
	 * 未经过排序、填充、与A相反的样式、有/没有数值才对它进行排序
	 */
    F = SparseCore_allocate_sparse (ncol, nrow, fnz, TRUE, TRUE, -SIGN(stype),
	    xtype, Common) ;
    if (Common->status < SPARSE_OK)
    {
		return (NULL) ;	    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 转置和选择性地置换矩阵A */
    /* ---------------------------------------------------------------------- */

    if (stype != 0)
    {
		/* F = A (p,p)', 只是用A的上三角或者下三角部分 */
		ok = SparseCore_transpose_sym (A, values, Perm, F, Common) ;
    }
    else
    {
		/* F = A (p,f)' */
		ok = SparseCore_transpose_unsym (A, values, Perm, fset, nf, F, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 返回矩阵F，如果出现错误则返回NULL */
    /* ---------------------------------------------------------------------- */

    if (!ok)
    {
		SparseCore_free_sparse (&F, Common) ;
    }
    return (F) ;
}

/**
 * @brief 对A的列进行排序。返回已打包的A，即使它开始时是未填充的。删除对称矩阵中被忽略的部分中的项。
 * 
 */
int SparseCore_sort
(
    /* ---- in/out --- */
    sparse_csc *A,	/* 待排序稀疏矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Ap ;
    sparse_csc *F ;
    Int anz, ncol, nrow, stype ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    Common->status = SPARSE_OK ;
    nrow = A->nrow ;
    if (nrow <= 1)
    {
	/* 1*n的稀疏矩阵必须被排序 */
	A->sorted = TRUE ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 分配内存空间 */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    SparseCore_allocate_work (0, MAX (nrow, ncol), 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;	/* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    anz = SparseCore_nnz (A, Common) ;
    stype = A->stype ;

    /* ---------------------------------------------------------------------- */
    /* 对矩阵进行列排序 */
    /* ---------------------------------------------------------------------- */

    /* 为转置分配工作区:ncole*nrow，与A相同的非零的#，排序，填充，与A相同的样式，与A相同的数值类型。 */
    F = SparseCore_allocate_sparse (ncol, nrow, anz, TRUE, TRUE, stype,
	    A->xtype, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;	/* 内存溢出 */
    }

    if (stype != 0)
    {
	/* F = A', 只用上三角或者下三角 */
	SparseCore_transpose_sym (A, 1, NULL, F, Common) ;
	A->packed = TRUE ;
	/* A = F' */
	SparseCore_transpose_sym (F, 1, NULL, A, Common) ;
    }
    else
    {
	/* F = A' */
	SparseCore_transpose_unsym (A, 1, NULL, NULL, 0, F, Common) ;
	A->packed = TRUE ;
	/* A = F' */
	SparseCore_transpose_unsym (F, 1, NULL, NULL, 0, A, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 如果需要，缩小A的尺寸。这必须成功。*/
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    anz = Ap [ncol] ;
    SparseCore_reallocate_sparse (anz, A, Common) ;

    /* ---------------------------------------------------------------------- */
    /* 释放工作空间 */
    /* ---------------------------------------------------------------------- */

    SparseCore_free_sparse (&F, Common) ;
    return (TRUE) ;
}
