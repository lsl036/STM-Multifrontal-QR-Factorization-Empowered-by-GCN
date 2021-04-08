/**
 * @file SparseChol_factorize_change.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-22
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef NCHOLESKY

#include "Sparse_internal.h"
#include "SparseChol.h"
#include "tpsm.h"

/**
 * @brief 	分解PAP'，L既是输入，也是输出。矩阵A必须与通过分析的矩阵有相同的非零模式
 * 
 */
int SparseChol_factorize
(
    /* ---- input ---- */
	//TPSM_t *tpool,
    sparse_csc *A,	/* 分解的稀疏矩阵 */
    /* ---- in/out --- */
    sparse_factor *L,	/* 结果因子 */
    /* --------------- */
    sparse_common *Common
)
{
    double zero [2] ;
    zero [0] = 0 ;
    zero [1] = 0 ;
    return (SparseChol_factorize_p ( A, zero, NULL, 0, L, Common)) ;
}

/**
 * @brief 	和SparseChol_factorize相同，但是有更多的操作
 * 
 */
int SparseChol_factorize_p
(
    /* ---- input ---- */
	//TPSM_t *tpool,
    sparse_csc *A,	/* 分解的稀疏矩阵 */
    double beta [2],	/* beta*I+A或者beta*I+A'*A */
    Int *fset,			/* 0:(A->ncol)-1的子集 */
    size_t fsize,		/* fset的大小 */
    /* ---- in/out --- */
    sparse_factor *L,	/* 结果因子 */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *S, *F, *A1, *A2 ;
    Int nrow, ncol, stype, convert, n, nsuper, grow2, status ;
    size_t s, t, uncol ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    n = L->n ;
    stype = A->stype ;
    if (L->n != A->nrow)
    {
		ERROR (SPARSE_INVALID, "A and L dimensions do not match") ;
		return (FALSE) ;
    }
    if (stype != 0 && nrow != ncol)
    {
		ERROR (SPARSE_INVALID, "matrix invalid") ;
		return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    nsuper = (L->is_super ? L->nsuper : 0) ;
    uncol = ((stype != 0) ? 0 : ncol) ;

    /* s = 2*nrow + MAX (uncol, 2*nsuper) */
    s = SparseCore_mult_size_t (nsuper, 2, &ok) ;
    s = MAX (uncol, s) ;
    t = SparseCore_mult_size_t (nrow, 2, &ok) ;
    s = SparseCore_add_size_t (s, t, &ok) ;
    if (!ok)
    {
		ERROR (SPARSE_TOO_LARGE, "problem too large") ;
		return (FALSE) ;
    }

    SparseCore_allocate_work (nrow, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
		return (FALSE) ;
    }

    S  = NULL ;
    F  = NULL ;
    A1 = NULL ;
    A2 = NULL ;

    /* 完成后，如果需要，转换为另一个表单 */
    convert = !(Common->final_asis) ;
    /* ---------------------------------------------------------------------- */
    /* 执行超节点LL'或者简单LDL'分解 */
    /* ---------------------------------------------------------------------- */
	// 执行超节点LL'分解or简单LDL'分解
    if (L->is_super)
    {

#ifndef NSUPERNODAL

	/* ------------------------------------------------------------------ */
	/* 超节点分解 */
	/* ------------------------------------------------------------------ */

	if (L->ordering == SPARSE_NATURAL)
	{

	    /* -------------------------------------------------------------- */
	    /* 自然排序 */
	    /* -------------------------------------------------------------- */

	    if (stype > 0)
	    {
		/* S = tril (A'), 不需要F */
		/* 工作空间: Iwork (nrow) */
		A1 = SparseCore_ptranspose (A, 2, NULL, NULL, 0, Common) ;
		S = A1 ;
	    }
	    else if (stype < 0)
	    {
		/* 这是自然排序最快的选择 */
		/* S = A; 不需要F */
		S = A ;
	    }
	    else
	    {
		/* F = A(:,f)' */
		/* 工作空间: Iwork (nrow) */
		/* 工作空间: Iwork (如果不是fset--nrow; 如果是fset--MAX (nrow,ncol) */
		A1 = SparseCore_ptranspose (A, 2, NULL, fset, fsize, Common) ;
		F = A1 ;
		/* S = A */
		S = A ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* 分解之前置换输入矩阵 */
	    /* -------------------------------------------------------------- */
		// 分解之前置换输入矩阵
	    if (stype > 0)
	    {
			/* 这是分解一个组合矩阵的最快选择 */
			/* S = tril (PAP'); 不需要F */
			/* 工作空间： Iwork (2*nrow) */
			A1 = SparseCore_ptranspose (A, 2, L->Perm, NULL, 0, Common) ;
			S = A1 ;
	    }
	    else if (stype < 0)
	    {
		/* A2 = triu (PAP') */
		/* 工作空间： Iwork (2*nrow) */
		A2 = SparseCore_ptranspose (A, 2, L->Perm, NULL, 0, Common) ;
		/* S = tril (A2'); 不需要F */
		/* 工作空间： Iwork (nrow) */
		A1 = SparseCore_ptranspose (A2, 2, NULL, NULL, 0, Common) ;
		S = A1 ;
		SparseCore_free_sparse (&A2, Common) ;
	    }
	    else
	    {
		/* F = A(p,f)' */
		/* 工作空间： Iwork (如果无fset--nrow; 如果有fset--MAX (nrow,ncol)) */
		A1 = SparseCore_ptranspose (A, 2, L->Perm, fset, fsize, Common) ;
		F = A1 ;
		/* S = F' */
		/* 工作空间： Iwork (nrow) */
		A2 = SparseCore_ptranspose (F, 2, NULL, NULL, 0, Common) ;
		S = A2 ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* 超节点分解 */
	/* ------------------------------------------------------------------ */

	/* 工作空间： Flag (nrow), Head (nrow+1), Iwork (2*nrow+2*nsuper) */
	if (Common->status == SPARSE_OK)
	{
		// F, beta 为空
	    SparseChol_super_numeric ( S, F, beta, L, Common) ;
	}
	status = Common->status ;

	/* ------------------------------------------------------------------ */
	/* 需要的话，转换成最终形式 */
	/* ------------------------------------------------------------------ */

	if (Common->status >= SPARSE_OK && convert)
	{
	    /* 工作空间： none */
	    ok = SparseCore_change_factor (L->xtype, Common->final_ll,
		    Common->final_super, Common->final_pack,
		    Common->final_monotonic, L, Common) ;
	    if (ok && Common->final_resymbol && !(L->is_super))
	    {
		/* 工作空间： Flag (nrow), Head (nrow+1),
		 *	对称:   Iwork (2*nrow)
		 *	非对称: Iwork (2*nrow+ncol) */
		SparseChol_resymbol_noperm (S, fset, fsize, Common->final_pack,
		    L, Common) ;
	    }
	}

#else

	/* ------------------------------------------------------------------ */
	/* HNUCHOL超节点模式未安装 */
	/* ------------------------------------------------------------------ */

	status = SPARSE_NOT_INSTALLED ;
	ERROR (SPARSE_NOT_INSTALLED,"Supernodal module not installed") ;

#endif

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* 简单LDL'分解 */
	/* ------------------------------------------------------------------ */

	/* 如有必要，对输入矩阵A进行排列。SparseChol_rowfac要求对称情况下的列形式的triu(A)，
	 * 非对称情况下的列形式的A，非对称情况下的列形式的A，或者等价的列形式的A'(矩阵F).
	 */

	if (L->ordering == SPARSE_NATURAL)
	{

	    /* -------------------------------------------------------------- */
	    /* 自然排序 */
	    /* -------------------------------------------------------------- */

	    if (stype > 0)
	    {
		/* F不需要, S = A */
		S = A ;
	    }
	    else if (stype < 0)
	    {
		/* F不需要, S = A' */
		/* 工作空间： Iwork (nrow) */
		A2 = SparseCore_ptranspose (A, 2, NULL, NULL, 0, Common) ;
		S = A2 ;
	    }
	    else
	    {
		/* F = A (:,f)' */
		/* 工作空间： Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
		A1 = SparseCore_ptranspose (A, 2, NULL, fset, fsize, Common) ;
		F = A1 ;
		S = A ;
	    }

	}
	else
	{
	    /* -------------------------------------------------------------- */
	    /* 在分解之前对输入矩阵进行置换 */
	    /* -------------------------------------------------------------- */
	    if (stype > 0)
	    {
			/* F = tril (A (p,p)') */
			/* 工作空间： Iwork (2*nrow) */
			A1 = SparseCore_ptranspose (A, 2, L->Perm, NULL, 0, Common) ;
			/* A2 = triu (F') */
			/* 工作空间： Iwork (nrow) */
			A2 = SparseCore_ptranspose (A1, 2, NULL, NULL, 0, Common) ;
			/* 对称情况不需要F，释放它并设为NULL*/
			SparseCore_free_sparse (&A1, Common) ;
	    }
	    else if (stype < 0)
	    {
		/* A2 = triu (A (p,p)'), 不需要F.  这是使用简单例程(SparseChol_rowfac)分解矩阵的最快方法。 */
		/* 工作空间： Iwork (2*nrow) */
		A2 = SparseCore_ptranspose (A, 2, L->Perm, NULL, 0, Common) ;
	    }
	    else
	    {
			/* F = A (p,f)' */
			/* 工作空间： Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
			A1 = SparseCore_ptranspose (A, 2, L->Perm, fset, fsize, Common) ;
			F = A1 ;
			/* A2 = F' */
			/* 工作空间： Iwork (nrow) */
			A2 = SparseCore_ptranspose (F, 2, NULL, NULL, 0, Common) ;
	    }
	    S = A2 ;
	}

	/* ------------------------------------------------------------------ */
	/* 简单LDL'或者LL'分解 */
	/* ------------------------------------------------------------------ */

	/* 分解beta*I+S (对称)或者beta*I+F*F' (非对称) */
	/* 工作空间： Flag (nrow), W (nrow), Iwork (2*nrow) */
	if (Common->status == SPARSE_OK)
	{
	    grow2 = Common->grow2 ;
	    L->is_ll = BOOLEAN (Common->final_ll) ;
	    if (L->xtype == SPARSE_PATTERN && Common->final_pack)
	    {
		/* 精确分配所需的空间 */
		Common->grow2 = 0 ;
	    }
	    SparseChol_rowfac (S, F, beta, 0, nrow, L, Common) ;
	    Common->grow2 = grow2 ;
	}
	status = Common->status ;

	/* ------------------------------------------------------------------ */
	/* 转换成最终形式，如果需要 */
	/* ------------------------------------------------------------------ */

	if (Common->status >= SPARSE_OK && convert)
	{
	    /* 工作空间： none */
	    SparseCore_change_factor (L->xtype, L->is_ll, FALSE,
		    Common->final_pack, Common->final_monotonic, L, Common) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 释放A1 A2 */
    /* ---------------------------------------------------------------------- */

    SparseCore_free_sparse (&A1, Common) ;
    SparseCore_free_sparse (&A2, Common) ;
    Common->status = MAX (Common->status, status) ;
    return (Common->status >= SPARSE_OK) ;
}

/* Remove entries from L that are not in the factorization of P*A*P', P*A*A'*P',
 * or P*F*F'*P' (depending on A->stype and whether fset is NULL or not). */

/**
 * @brief  	从L中删除不属于P*A*P'、P*A*A'*P'或P*F*F'*P'
 * 			(取决于A->stype以及fset是否为空)因数分解的项。
 * 
 */
int SparseChol_resymbol
(
    /* ---- input ---- */
    sparse_csc *A,	/* 要分析的矩阵 */
    Int *fset,			/* 0:(A->ncol)-1的子集 */
    size_t fsize,		/* fset的大小 */
    int pack,			/* 如果为真，填充L的列 */
    /* ---- in/out --- */
    sparse_factor *L,	/* 因式分解，在输出中删除条目 */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *H, *F, *G ;
    Int stype, nrow, ncol ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    Common->status = SPARSE_OK ;
    if (L->is_super)
    {
	/* 不能作用于超节点因数分解 */
	ERROR (SPARSE_INVALID, "cannot operate on supernodal L") ;
	return (FALSE) ;
    }
    if (L->n != A->nrow)
    {
	ERROR (SPARSE_INVALID, "A and L dimensions do not match") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    stype = A->stype ;
    nrow = A->nrow ;
    ncol = A->ncol ;

    /* s = 2*nrow + (stype ? 0 : ncol) */
    s = SparseCore_mult_size_t (nrow, 2, &ok) ;
    s = SparseCore_add_size_t (s, (stype ? 0 : ncol), &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (nrow, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 如有必要，对输入矩阵进行排列 */
    /* ---------------------------------------------------------------------- */

    H = NULL ;
    G = NULL ;

    if (stype > 0)
    {
	if (L->ordering == SPARSE_NATURAL)
	{
	    /* F = triu(A)' */
	    /* 工作空间： Iwork (nrow) */
	    G = SparseCore_ptranspose (A, 0, NULL, NULL, 0, Common) ;
	}
	else
	{
	    /* F = triu(A(p,p))' */
	    /* 工作空间： Iwork (2*nrow) */
	    G = SparseCore_ptranspose (A, 0, L->Perm, NULL, 0, Common) ;
	}
	F = G ;
    }
    else if (stype < 0)
    {
	if (L->ordering == SPARSE_NATURAL)
	{
	    F = A ;
	}
	else
	{
	    /* G = triu(A(p,p))' */
	    /* 工作空间： Iwork (2*nrow) */
	    G = SparseCore_ptranspose (A, 0, L->Perm, NULL, 0, Common) ;
	    /* H = G' */
	    /* 工作空间： Iwork (nrow) */
	    H = SparseCore_ptranspose (G, 0, NULL, NULL, 0, Common) ;
	    F = H ;
	}
    }
    else
    {
	if (L->ordering == SPARSE_NATURAL)
	{
	    F = A ;
	}
	else
	{
	    /* G = A(p,f)' */
	    /* 工作空间： Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
	    G = SparseCore_ptranspose (A, 0, L->Perm, fset, fsize, Common) ;
	    /* H = G' */
	    /* 工作空间： Iwork (ncol) */
	    H = SparseCore_ptranspose (G, 0, NULL, NULL, 0, Common) ;
	    F = H ;
	}
    }

    /* 这里不需要检查失败。如果F为空，SparseCore_resymbol_noperm将返回FALSE。 */

    /* ---------------------------------------------------------------------- */
    /* resymbol */
    /* ---------------------------------------------------------------------- */

    ok = SparseChol_resymbol_noperm (F, fset, fsize, pack, L, Common) ;

    /* ---------------------------------------------------------------------- */
    /* 释放临时矩阵，如果它们存在的话 */
    /* ---------------------------------------------------------------------- */

    SparseCore_free_sparse (&H, Common) ;
    SparseCore_free_sparse (&G, Common) ;
    return (ok) ;
}

/**
 * @brief 	重做I+ F*F'或I+A的符号LDL'或' LL'分解，其中F=A(:， F)。
 * 
 */
int SparseChol_resymbol_noperm
(
    /* ---- input ---- */
    sparse_csc *A,	/* 分析的矩阵 */
    Int *fset,			/* 0:(A->ncol)-1的子集 */
    size_t fsize,		/* fset的大小 */
    int pack,			/* if TRUE, 填充L的列 */
    /* ---- in/out --- */
    sparse_factor *L,	/* 因式分解，在输出中删除条目 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Lx, *Lz ;
    Int i, j, k, row, parent, p, pend, pdest, ncol, apacked, sorted, nrow, nf,
	use_fset, mark, jj, stype, xtype ;
    Int *Ap, *Ai, *Anz, *Li, *Lp, *Lnz, *Flag, *Head, *Link, *Anext, *Iwork ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    ncol = A->ncol ;
    nrow = A->nrow ;
    stype = A->stype ;
    if (stype > 0)
    {
	/* 上三角对称，不支持 */
	ERROR (SPARSE_INVALID, "symmetric upper not supported ") ;
	return (FALSE) ;
    }
    if (L->is_super)
    {
	/* 不能在超节点或符号分解上操作 */
	ERROR (SPARSE_INVALID, "cannot operate on supernodal L") ;
	return (FALSE) ;
    }
    if (L->n != A->nrow)
    {
	ERROR (SPARSE_INVALID, "A and L dimensions do not match") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    /* s = 2*nrow + (stype ? 0 : ncol) */
    s = SparseCore_mult_size_t (nrow, 2, &ok) ;
    if (stype != 0)
    {
	s = SparseCore_add_size_t (s, ncol, &ok) ;
    }
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (nrow, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;	/* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    Ai = A->i ;
    Ap = A->p ;
    Anz = A->nz ;
    apacked = A->packed ;
    sorted = A->sorted ;

    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;
    Lp = L->p ;
    Lnz = L->nz ;
    xtype = L->xtype ;

    /* 如果L在输入上是单调的，那么它可以在输出上被填充或不填充，
	 * 这取决于填充输入参数。 */

    /* 不能填充一个非单调矩阵 */
    if (!(L->is_monotonic))
    {
	pack = FALSE ;
    }

    pdest = 0 ;

    /* ---------------------------------------------------------------------- */
    /* 获取工作空间 */
    /* ---------------------------------------------------------------------- */

    Flag  = Common->Flag ;	/* 大小为nrow */
    Head  = Common->Head ;	/* 大小为nrow+1 */
    Iwork = Common->Iwork ;
    Link  = Iwork ;		/* 大小为nrow (i/i/l) [ */
    Lnz   = Iwork + nrow ;	/* 如果L没填充，大小为nrow (i/i/l) */
    Anext = Iwork + 2*((size_t) nrow) ;	/* 大小为ncol (i/i/l), unsym */
    for (j = 0 ; j < nrow ; j++)
    {
	Link [j] = EMPTY ;
    }

    /* 在L本身中使用Lnz */
    Lnz = L->nz ;

    /* ---------------------------------------------------------------------- */
    /* 对于非对称情况，将A (:，f)的每一列放入队列 */
    /* ---------------------------------------------------------------------- */

    /* 将基集的每一列放在链接列表中对应的 */
    /* 该列中的最小行索引 */

    if (stype == 0)
    {
	use_fset = (fset != NULL) ;
	if (use_fset)
	{
	    nf = fsize ;
	    /* 这是SparseCore_resymbol中唯一的O(ncol)循环。只需要检查fset。 */
	    for (j = 0 ; j < ncol ; j++)
	    {
		Anext [j] = -2 ;
	    }
	    for (jj = 0 ; jj < nf ; jj++)
	    {
		j = fset [jj] ;
		if (j < 0 || j > ncol || Anext [j] != -2)
		{
		    /* fset中超出范围或重复条目 */
		    ERROR (SPARSE_INVALID, "fset invalid") ;
		    return (FALSE) ;
		}
		/* 将列j标记为已看到 */
		Anext [j] = EMPTY ;
	    }
	}
	else
	{
	    nf = ncol ;
	}
	for (jj = 0 ; jj < nf ; jj++)
	{
	    j = (use_fset) ? (fset [jj]) : jj ;
	    /* 第j列是fset;查找最小的行(如果有的话) */
	    p = Ap [j] ;
	    pend = (apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    if (pend > p)
	    {
		k = Ai [p] ;
		if (!sorted)
		{
		    for ( ; p < pend ; p++)
		    {
			k = MIN (k, Ai [p]) ;
		    }
		}
		/* 将列j放到链接列表k上 */
		Anext [j] = Head [k] ;
		Head [k] = j ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 重新计算符号LDL因子分解 */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < nrow ; k++)
    {

	/* ------------------------------------------------------------------ */
	/* 计算I+F*F'或者I+A的第k列 */
	/* ------------------------------------------------------------------ */

	/* 标记对角项 */
	/* mark = SparseCore_clear_flag (Common) ; */
	SPARSE_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	Flag [k] = mark ;

	if (stype != 0)
	{
	    /* 将A的k列合并为标志(仅限下三角形部分) */
	    p = Ap [k] ;
	    pend = (apacked) ? (Ap [k+1]) : (p + Anz [k]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i > k)
		{
		    Flag [i] = mark ;
		}
	    }
	}
	else
	{
	    /* 对于第一行下标在k行的每一列j */
	    for (j = Head [k] ; j != EMPTY ; j = Anext [j])
	    {
		/* 将A的第j列合并到标记中 */
		p = Ap [j] ;
		pend = (apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    Flag [Ai [p]] = mark ;
		}
	    }
	    /* 清除第k个链表 */
	    Head [k] = EMPTY ;
	}

	/* ------------------------------------------------------------------ */
	/*计算L的第k列的pruned pattern联合的孩子 */
	/* ------------------------------------------------------------------ */

	/* 对于父结点为k的L的每一列j */
	for (j = Link [k] ; j != EMPTY ; j = Link [j])
	{
	    /* 将L的第j列合并到标记中 */
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;
	    p++ ;	    /* 跳过对角项 */
	    for ( ; p < pend ; p++)
	    {
		/* 加到pattern */
		Flag [Li [p]] = mark ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* 修剪L的第k列 */
	/* ------------------------------------------------------------------ */

	p = Lp [k] ;
	pend = p + Lnz [k] ;

	if (pack)
	{
	    /* 向上移动k列 */
	    Lp [k] = pdest ;
	}
	else
	{
	    /* 将列k保留到位，只需降低Lnz [k] */
	    pdest = p ;
	}

	for ( ; p < pend ; p++)
	{
	    row = Li [p] ;
	    if (Flag [row] == mark)
	    {
		/* 保留这一项 */
		Li [pdest] = row ;
		if (xtype == SPARSE_REAL)
		{
		    Lx [pdest] = Lx [p] ;
		}
		pdest++ ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* 为父列准备此列 */
	/* ------------------------------------------------------------------ */

	Lnz [k] = pdest - Lp [k] ;

	/* 父元素是该列对角线后的第一个条目 */
	parent = (Lnz [k] > 1) ? (Li [Lp [k] + 1]) : EMPTY ;

	if (parent != EMPTY)
	{
	    Link [k] = Link [parent] ;
	    Link [parent] = k ;
	}
    }

    /* 使用Iwork完成链接，Lnz(如果需要)和Anext */

    /* ---------------------------------------------------------------------- */
    /* 如果需要，将L转换为打包 */
    /* ---------------------------------------------------------------------- */

    if (pack)
    {
	/* 完成Lp */
	Lp [nrow] = pdest ;
	/* 把L缩小到足够大。它不能失败。 */
	/* 工作空间： none */
	SparseCore_reallocate_factor (Lp [nrow], L, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 清除工作空间 */
    /* ---------------------------------------------------------------------- */

    /* SparseCore_clear_flag (Common) ; */
    SPARSE_CLEAR_FLAG (Common) ;
    return (TRUE) ;
}

#endif
