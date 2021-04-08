/**
 * @file SparseCore_change_factor.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */

/* 更改SparseCore_factor对象的状态
 * numeric/symbolic
 * LL/LDL
 * simplicial/super
 * packed/unpacked,
 * monotonic/non-monotonic
 */
#include "Sparse_internal.h"
#include "SparseCore.h"

static void natural_list (sparse_factor *L) ;

/* ========================================================================== */
/* === 模板 ============================================================= */
/* ========================================================================== */

#define REAL

#include "Sparse_template.h"

/**
 * @brief 更改numeric/symbolic状态
 */
static void TEMPLATE (change_simplicial_numeric)
(
    sparse_factor *L,
    Int to_ll,
    Int to_packed,
    Int *newLi,
    double *newLx,
    double *newLz,
    Int lnz,
    Int grow,
    double grow1,
    Int grow2,
    Int make_ll,
    Int make_monotonic,
    Int make_ldl,
    sparse_common *Common
)
{
    double xlen, dj [1], ljj [1], lj2 [1] ;
    double *Lx, *Lz ;
    Int *Lp, *Li, *Lnz ;
    Int n, j, len, pnew, pold, k, p, pend ;

    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;
    Lnz = L->nz ;

    if (make_ll)
    {
	L->minor = n ;
    }

    if (make_monotonic)
    {

	/* ------------------------------------------------------------------ */
	/* 重新排列列使它们单调 */
	/* ------------------------------------------------------------------ */

	pnew = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    /* 拷贝和压缩列j */
	    len = Lnz [j] ;
	    pold = Lp [j] ;

	    if (make_ll)
	    {

		/* ---------------------------------------------------------- */
		/* 复制并将LDL'转化为LL' */
		/* ---------------------------------------------------------- */

		/* dj = Lx [pold] ; */
		ASSIGN_REAL (dj,0, Lx,pold) ;

		if (IS_LE_ZERO (dj [0]))
		{
		    /* 转换失败了;矩阵不是正定的。
			 * 不要修改该列，以便在需要时通过转换回LDL'恢复LDL'的因数分解。
			 * 继续转换，但标记错误。 */
		    if (L->minor == (size_t) n)
		    {
			ERROR (SPARSE_NOT_POSDEF, "L not positive definite") ;
			L->minor = j ;
		    }
		    for (k = 0 ; k < len ; k++)
		    {
			newLi [pnew + k] = Li [pold + k] ;
			/* newLx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (newLx, newLz, pnew+k, Lx, Lz, pold+k) ;
		    }
		}
		else
		{
		    ljj [0] = sqrt (dj [0]) ;
		    newLi [pnew] = j ;
		    /* newLx [pnew] = ljj ; */
		    ASSIGN_REAL (newLx, pnew, ljj, 0) ;
		    CLEAR_IMAG (newLx, newLz, pnew) ;

		    for (k = 1 ; k < len ; k++)
		    {
			newLi [pnew + k] = Li [pold + k] ;
			/* newLx [pnew + k] = Lx [pold + k] * ljj ; */
			MULT_REAL (newLx, newLz, pnew+k, Lx, Lz, pold+k, ljj,0);
		    }
		}

	    }
	    else if (make_ldl)
	    {

		/* ---------------------------------------------------------- */
		/* 复制并将LL'转换为LDL' */
		/* ---------------------------------------------------------- */

		/* ljj = Lx [pold] ; */
		ASSIGN_REAL (ljj, 0, Lx, pold) ;

		if (ljj [0] <= 0)
		{
		    /* 矩阵不是正定的;按原样复制列 */
		    for (k = 0 ; k < len ; k++)
		    {
			newLi [pnew + k] = Li [pold + k] ;
			/* newLx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (newLx, newLz, pnew+k, Lx, Lz, pold+k) ;
		    }
		}
		else
		{
		    newLi [pnew] = j ;
		    /* newLx [pnew] = ljj*ljj ; */
		    lj2 [0] = ljj [0] * ljj [0] ;
		    ASSIGN_REAL (newLx, pnew, lj2, 0) ;
		    CLEAR_IMAG (newLx, newLz, pnew) ;

		    for (k = 1 ; k < len ; k++)
		    {
			newLi [pnew + k] = Li [pold + k] ;
			/* newLx [pnew + k] = Lx [pold + k] / ljj ; */
			DIV_REAL (newLx, newLz, pnew+k, Lx, Lz, pold+k, ljj,0) ;
		    }
		}

	    }
	    else
	    {

		/* ---------------------------------------------------------- */
		/* 复制并保持LL'或LDL'原样 */
		/* ---------------------------------------------------------- */

		for (k = 0 ; k < len ; k++)
		{
		    newLi [pnew + k] = Li [pold + k] ;
		    /* newLx [pnew + k] = Lx [pold + k] ; */
		    ASSIGN (newLx, newLz, pnew+k, Lx, Lz, pold+k) ;
		}
	    }

	    Lp [j] = pnew ;

	    /* 以双精度计算len以避免整型溢出 */
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (Int) xlen ;
	    }
	    pnew += len ;
	}
	Lp [n] = pnew ;

	/* 释放旧的L->i和L->x，并更换新的 */
	CORE(free) (L->nzmax, sizeof (Int), L->i, Common) ;

#ifdef REAL
	CORE(free) (L->nzmax, sizeof (double), L->x, Common) ;
#else
	CORE(free) (L->nzmax, sizeof (double), L->x, Common) ;
	CORE(free) (L->nzmax, sizeof (double), L->z, Common) ;
#endif

	L->i = newLi ;
	L->x = newLx ;
	L->z = newLz ;
	L->nzmax = lnz ;

	/* 重构链接列表 */
	natural_list (L) ;

    }
    else if (to_packed)
    {

	/* ------------------------------------------------------------------ */
	/* 已经是单调的了，只是把L的列集合起来 */
	/* ------------------------------------------------------------------ */

	pnew = 0 ;

	if (make_ll)
	{

	    /* -------------------------------------------------------------- */
	    /* 压缩并将LDL'转化为LL' */
	    /* -------------------------------------------------------------- */

	    for (j = 0 ; j < n ; j++)
	    {
		/* 压缩 column j */
		pold = Lp [j] ;
		len = Lnz [j] ;

		/* dj = Lx [pold] ; */
		ASSIGN_REAL (dj,0, Lx,pold) ;

		if (IS_LE_ZERO (dj [0]))
		{
		    /* 转换失败了;矩阵不是正定的。
			 * 不要修改该列，以便在需要时通过转换回LDL'恢复LDL'的因数分解。
			 * 继续转换，但标记错误。 */
		    if (L->minor == (size_t) n)
		    {
			ERROR (SPARSE_NOT_POSDEF, "L not positive definite") ;
			L->minor = j ;
		    }
		    for (k = 0 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (Lx, Lz, pnew+k, Lx, Lz, pold+k) ;
		    }
		}
		else
		{
		    ljj [0] = sqrt (dj [0]) ;
		    Li [pnew] = j ;

		    /* Lx [pnew] = ljj ; */
		    ASSIGN_REAL (Lx, pnew, ljj, 0) ;
		    CLEAR_IMAG (Lx, Lz, pnew) ;

		    for (k = 1 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] * ljj ; */
			MULT_REAL (Lx, Lz, pnew+k, Lx, Lz, pold+k, ljj,0) ;
		    }
		}
		Lp [j] = pnew ;
		pnew += len ;
	    }

	}
	else if (make_ldl)
	{

	    /* -------------------------------------------------------------- */
	    /* 压缩并将LL'转化为LDL' */
	    /* -------------------------------------------------------------- */

	    for (j = 0 ; j < n ; j++)
	    {
		/* 压缩第j列 */
		pold = Lp [j] ;
		len = Lnz [j] ;

		/* ljj = Lx [pold] ; */
		ASSIGN_REAL (ljj, 0, Lx, pold) ;
		if (ljj [0] <= 0)
		{
		    /* 矩阵不是正定的;按原样包列 */
		    for (k = 0 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (Lx, Lz, pnew+k, Lx, Lz, pold+k) ;
		    }
		}
		else
		{
		    Li [pnew] = Li [pold] ;

		    /* Lx [pnew] = ljj*ljj ; */
		    lj2 [0] = ljj [0] * ljj [0] ;
		    ASSIGN_REAL (Lx, pnew, lj2, 0) ;
		    CLEAR_IMAG (Lx, Lz, pnew) ;

		    for (k = 1 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] / ljj ; */
			DIV_REAL (Lx, Lz, pnew+k, Lx, Lz, pold+k, ljj,0) ;
		    }
		}
		Lp [j] = pnew ;
		pnew += len ;
	    }

	}
	else
	{

	    /* ---------------------------------------------------------- */
	    /* 压缩并保留LL'或LDL'原样 */
	    /* ---------------------------------------------------------- */

	    for (j = 0 ; j < n ; j++)
	    {
		/* 压缩第j列 */
		pold = Lp [j] ;
		len = Lnz [j] ;
		if (pnew < pold)
		{
		    for (k = 0 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (Lx, Lz, pnew+k, Lx, Lz, pold+k) ;
		    }
		    Lp [j] = pnew ;
		}
		pnew += len ;
	    }
	}

	Lp [n] = pnew ;

    }
    else if (make_ll)
    {

	/* ------------------------------------------------------------------ */
	/* 将LDL'转换为' LL'，但在适当的地方这样做 */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < n ; j++)
	{
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;

	    /* dj = Lx [p] ; */
	    ASSIGN_REAL (dj,0, Lx,p) ;

	    if (IS_LE_ZERO (dj [0]))
	    {
		/* 转换失败了;矩阵不是正定的。
		 * 不要修改该列，以便在需要时通过转换回LDL'恢复LDL'的因数分解。
		 * 继续转换，但标记错误。 */
		if (L->minor == (size_t) n)
		{
		    ERROR (SPARSE_NOT_POSDEF, "L not positive definite") ;
		    L->minor = j ;
		}
	    }
	    else
	    {
		ljj [0] = sqrt (dj [0]) ;
		/* Lx [p] = ljj ; */
		ASSIGN_REAL (Lx,p, ljj,0) ;
		CLEAR_IMAG (Lx, Lz, p) ;

		for (p++ ; p < pend ; p++)
		{
		    /* Lx [p] *= ljj ; */
		    MULT_REAL (Lx,Lz,p, Lx,Lz,p, ljj,0) ;
		}
	    }
	}

    }
    else if (make_ldl)
    {

	/* ------------------------------------------------------------------ */
	/* 将LL'转化为LDL'，在原内存空间进行 */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < n ; j++)
	{
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;

	    /* ljj = Lx [p] ; */
	    ASSIGN_REAL (ljj, 0, Lx, p) ;

	    if (ljj [0] > 0)
	    {
		/* Lx [p] = ljj*ljj ; */
		lj2 [0] = ljj [0] * ljj [0] ;
		ASSIGN_REAL (Lx, p, lj2, 0) ;
		CLEAR_IMAG (Lx, Lz, p) ;

		for (p++ ; p < pend ; p++)
		{
		    /* Lx [p] /= ljj ; */
		    DIV_REAL (Lx,Lz,p, Lx,Lz,p, ljj,0) ;
		}
	    }
	}
    }

    L->is_ll = to_ll ;
}

/**
 * @brief 	超节点L只能是真实的
 * 
 */
static void TEMPLATE (ll_super_to_simplicial_numeric)
(
    sparse_factor *L,
    Int to_packed,
    Int to_ll,
    sparse_common *Common
)
{
    double ljj [1], lj2 [1] ;
    double *Lx ;
    Int *Ls, *Lpi, *Lpx, *Super, *Lp, *Li, *Lnz ;
    Int n, lnz, s, nsuper, p, psi, psx, psend, nsrow, nscol, ii, jj, j, k1, k2,
	q ;

    L->is_ll = to_ll ;

    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lnz = L->nz ;
    lnz = L->nzmax ;

    n = L->n ;
    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;

    p = 0 ;

    for (s = 0 ; s < nsuper ; s++)
    {
	k1 = Super [s] ;
	k2 = Super [s+1] ;
	psi = Lpi [s] ;
	psend = Lpi [s+1] ;
	psx = Lpx [s] ;
	nsrow = psend - psi ;
	nscol = k2 - k1 ;

	for (jj = 0 ; jj < nscol ; jj++)
	{
	    /* L的列j从这里开始 */
	    j = jj + k1 ;

	    if (to_ll)
	    {
		if (to_packed)
		{

		    /* ------------------------------------------------------ */
		    /* 转换到LL'压缩 */
		    /* ------------------------------------------------------ */

		    Lp [j] = p ;
		    for (ii = jj ; ii < nsrow ; ii++)
		    {
			/* 从超节点中获取L(i,j)并存储在第j列中 */
			Li [p] = Ls [psi + ii] ;
			/* Lx [p] = Lx [psx + ii + jj*nsrow] ; */
			q = psx + ii + jj*nsrow ;
			ASSIGN (Lx,-,p, Lx,-,q) ;
			p++ ;
		    }
		    Lnz [j] = p - Lp [j] ;

		}
		else
		{

		    /* ------------------------------------------------------ */
		    /* 转换到LL'不压缩 */
		    /* ------------------------------------------------------ */

		    p = psx + jj + jj*nsrow ;
		    Lp [j] = p ;
		    Li [p] = j ;
		    Lnz [j] = nsrow - jj ;
		    p++ ;
		    for (ii = jj + 1 ; ii < nsrow ; ii++)
		    {
			/* 从超节点中获取L(i,j)并存储在第j列中 */
			Li [psx + ii + jj*nsrow] = Ls [psi + ii] ;
		    }

		}
	    }
	    else
	    {
		if (to_packed)
		{

		    /* ------------------------------------------------------ */
		    /* 转换到LDL'压缩 */
		    /* ------------------------------------------------------ */

		    Lp [j] = p ;
		    /* ljj = Lx [psx + jj + jj*nsrow] ; */
		    ASSIGN_REAL (ljj, 0, Lx, psx + jj + jj*nsrow) ;

		    if (ljj [0] <= 0)
		    {
			/* 矩阵不是正定的;不划分 */
			/* Lx [p] = ljj ; */
			ASSIGN_REAL (Lx, p, ljj, 0) ;
			CLEAR_IMAG (Lx, Lz, p) ;
			ljj [0] = 1 ;
		    }
		    else
		    {
			lj2 [0] = ljj [0] * ljj [0] ;
			/* Lx [p] = ljj*ljj ; */
			ASSIGN_REAL (Lx, p, lj2, 0) ;
			CLEAR_IMAG (Lx, Lz, p) ;
		    }
		    Li [p] = j ;
		    p++ ;
		    for (ii = jj + 1 ; ii < nsrow ; ii++)
		    {
			/* 从超节点中获取L(i,j)并存储在第j列中 */
			Li [p] = Ls [psi + ii] ;

			/* Lx [p] = Lx [psx + ii + jj*nsrow] / ljj ; */
			q = psx + ii + jj*nsrow ;
			DIV_REAL (Lx, Lz, p, Lx, Lz, q, ljj,0) ;
			p++ ;
		    }
		    Lnz [j] = p - Lp [j] ;

		}
		else
		{

		    /* ------------------------------------------------------ */
		    /* 转换到LDL'不压缩 */
		    /* ------------------------------------------------------ */

		    p = psx + jj + jj*nsrow ;
		    Lp [j] = p ;

		    /* ljj = Lx [p] ; */
		    ASSIGN_REAL (ljj,0, Lx,p) ;

		    if (ljj [0] <= 0)
		    {
			/* 矩阵不是正定的;不划分 */
			/* Lx [p] = ljj ; */
			ASSIGN_REAL (Lx, p, ljj, 0) ;
			CLEAR_IMAG (Lx, Lz, p) ;
			ljj [0] = 1 ;
		    }
		    else
		    {
			lj2 [0] = ljj [0] * ljj [0] ;
			/* Lx [p] = ljj*ljj ; */
			ASSIGN_REAL (Lx, p, lj2, 0) ;
			CLEAR_IMAG (Lx, Lz, p) ;
		    }
		    Li [p] = j ;
		    Lnz [j] = nsrow - jj ;
		    p++ ;
		    for (ii = jj + 1 ; ii < nsrow ; ii++)
		    {
			/* 从超节点中获取L(i,j)并存储在第j列中 */
			Li [psx + ii + jj*nsrow] = Ls [psi + ii] ;

			/* Lx [psx + ii + jj*nsrow] /= ljj ; */
			q = psx + ii + jj*nsrow ;
			DIV_REAL (Lx, Lz, q, Lx, Lz, q, ljj,0) ;
		    }
		}
	    }
	}
    }

    if (to_packed)
    {
	Lp [n] = p ;
	/* 减小L->x的大小以匹配L->i。这不能失败。 */
	L->x = CORE(realloc) (lnz, 
		sizeof (double), L->x, &(L->xsize), Common) ;
	Common->status = SPARSE_OK ;
    }
    else
    {
	Lp [n] = Lpx [nsuper] ;
    }
}

#undef PATTERN
#undef REAL

/**
 * @brief 创建一个自然排序的双链接列列表。
 * 
 * @param L 
 */
static void natural_list (sparse_factor *L)
{
    Int head, tail, n, j ;
    Int *Lnext, *Lprev ;
    Lnext = L->next ;
    Lprev = L->prev ;
    n = L->n ;
    head = n+1 ;
    tail = n ;
    Lnext [head] = 0 ;
    Lprev [head] = EMPTY ;
    Lnext [tail] = EMPTY ;
    Lprev [tail] = n-1 ;
    for (j = 0 ; j < n ; j++)
    {
	Lnext [j] = j+1 ;
	Lprev [j] = j-1 ;
    }
    Lprev [0] = head ;
    L->is_monotonic = TRUE ;
}

/* Allocate O(n) arrays for simplicial numeric factorization.  Initializes
 * the link lists only.  Does not allocate the L->i, L->x, or L->z arrays. */

/**
 * @brief 	分配O(n)数组用于简单数值因数分解。
 * 			只初始化链接列表。不分配L->i, L->x, L->z数组。
 * 
 * @param L 
 * @param Common 
 * @return int 
 */
static int allocate_simplicial_numeric
(
    sparse_factor *L,
    sparse_common *Common
)
{
    Int n ;
    Int *Lp, *Lnz, *Lprev, *Lnext ;
    size_t n1, n2 ;

    n = L->n ;

    /* 这不会导致size_t溢出 */
    n1 = ((size_t) n) + 1 ;
    n2 = ((size_t) n) + 2 ;

    Lp = CORE(malloc) (n1, sizeof (Int), Common) ;
    Lnz = CORE(malloc) (n, sizeof (Int), Common) ;
    Lprev = CORE(malloc) (n2, sizeof (Int), Common) ;
    Lnext = CORE(malloc) (n2, sizeof (Int), Common) ;

    if (Common->status < SPARSE_OK)
    {
	CORE(free) (n1, sizeof (Int), Lp,    Common) ;
	CORE(free) (n,   sizeof (Int), Lnz,   Common) ;
	CORE(free) (n2, sizeof (Int), Lprev, Common) ;
	CORE(free) (n2, sizeof (Int), Lnext, Common) ;
	return (FALSE) ;	/* 内存溢出 */
    }

    /* ======================将更改提交到L======================== */

    L->p = Lp ;
    L->nz = Lnz ;
    L->prev = Lprev ;
    L->next = Lnext ;
    /* 按自然顺序初始化列的双链表 */
    natural_list (L) ;
    return (TRUE) ;
}

/**
 * @brief 转化为单纯符号因子，超节点符号因子。不初始化新空间。
 * 
 * @param L 
 * @param Common 
 * @return int 
 */
static int simplicial_symbolic_to_super_symbolic
(
    sparse_factor *L,
    sparse_common *Common
)
{
    Int nsuper, xsize, ssize ;
    Int *Lsuper, *Lpi, *Lpx, *Ls ;
    size_t nsuper1 ;

    xsize  = L->xsize ;
    ssize  = L->ssize ;
    nsuper = L->nsuper ;
    nsuper1 = ((size_t) nsuper) + 1 ;

    /* O(nsuper)大小的数组, nsuper <= n */
    Lsuper = CORE(malloc) (nsuper1, sizeof (Int), Common) ;
    Lpi    = CORE(malloc) (nsuper1, sizeof (Int), Common) ;
    Lpx    = CORE(malloc) (nsuper1, sizeof (Int), Common) ;

    /* O(ssize)大小的数组, ssize <= nnz(L), 并且一般可能会更小 */
    Ls = CORE(malloc) (ssize, sizeof (Int), Common) ;

    if (Common->status < SPARSE_OK)
    {
	CORE(free) (nsuper1, sizeof (Int), Lsuper, Common) ;
	CORE(free) (nsuper1, sizeof (Int), Lpi,    Common) ;
	CORE(free) (nsuper1, sizeof (Int), Lpx,    Common) ;
	CORE(free) (ssize,    sizeof (Int), Ls,     Common) ;
	return (FALSE) ;	/* 内存溢出 */
    }

    /* ========================将更改提交到L======================  */

    L->maxcsize = 0 ;
    L->maxesize = 0 ;

    L->super = Lsuper ;
    L->pi = Lpi ;
    L->px = Lpx ;
    L->s  = Ls ;
    Ls [0] = EMPTY ;	    /* 超节点模式定义 */

    L->is_super = TRUE ;
    L->is_ll = TRUE ;	    /* 不支持超节点LDL */
    L->xtype = SPARSE_PATTERN ;
    L->dtype = DTYPE ;
    L->minor = L->n ;
    return (TRUE) ;
}

/* Convert any factor L to a simplicial symbolic factor, leaving only L->Perm
 * and L->ColCount.  Cannot fail.  Any of the components of L (except Perm and
 * ColCount) may already be free'd.  */

/**
 * @brief 	将任何因子L转换为一个简单符号因子，只留下L->Perm和L->ColCount。
 * 			不能失败。L的任何组成部分(除了Perm和ColCount)可能已经被释放。
 * 
 * @param L 
 * @param to_ll 
 * @param Common 
 */
static void any_to_simplicial_symbolic
(
    sparse_factor *L,
    int to_ll,
    sparse_common *Common
)
{
    Int n, lnz, xs, ss, s, e ;
    size_t n1, n2 ;

    /* ===========================将更改提交到L===================  */

    n = L->n ;
    lnz = L->nzmax ;
    s = L->nsuper + 1 ;
    xs = (L->is_super) ? ((Int) (L->xsize)) : (lnz) ;
    e = 1 ;
    ss = L->ssize ;

    /* 这不会导致size_t溢出 */
    n1 = ((size_t) n) + 1 ;
    n2 = ((size_t) n) + 2 ;

    /* 释放除了符号分析之外的所有参数(Perm，ColCount) */
    L->p     = CORE(free) (n1,  sizeof (Int),      L->p,     Common) ;
    L->i     = CORE(free) (lnz, sizeof (Int),      L->i,     Common) ;
    L->x     = CORE(free) (xs,  e*sizeof (double), L->x,     Common) ;
    L->z     = CORE(free) (lnz, sizeof (double),   L->z,     Common) ;
    L->nz    = CORE(free) (n,   sizeof (Int),      L->nz,    Common) ;
    L->next  = CORE(free) (n2,  sizeof (Int),      L->next,  Common) ;
    L->prev  = CORE(free) (n2,  sizeof (Int),      L->prev,  Common) ;
    L->super = CORE(free) (s,   sizeof (Int),      L->super, Common) ;
    L->pi    = CORE(free) (s,   sizeof (Int),      L->pi,    Common) ;
    L->px    = CORE(free) (s,   sizeof (Int),      L->px,    Common) ;
    L->s     = CORE(free) (ss,  sizeof (Int),      L->s,     Common) ;
    L->nzmax = 0 ;
    L->is_super = FALSE ;
    L->xtype = SPARSE_PATTERN ;
    L->dtype = DTYPE ;
    L->minor = n ;
    L->is_ll = to_ll ;
}

/**
 * @brief 将数值超节点L转换为符号超节点。不能失败。
 * 
 * @param L 
 * @param Common 
 */
static void ll_super_to_super_symbolic
(
    sparse_factor *L,
    sparse_common *Common
)
{

    /* ===========================将更改提交到L===================  */

    /* 释放除超节点数值因子以外的所有因子 */
    L->x = CORE(free) (L->xsize, sizeof (double), L->x, Common) ;
    L->xtype = SPARSE_PATTERN ;
    L->dtype = DTYPE ;
    L->minor = L->n ;
    L->is_ll = TRUE ;	    /* 不支持超节点LDL */
}

/* Convert a simplicial symbolic L to a simplicial numeric L; allocate space
 * for L using L->ColCount from symbolic analysis, and set L to identity.
 *
 * If packed < 0, then this routine is creating a copy of another factor
 * (via SparseCore_copy_factor).  In this case, the space is not initialized. */

/**
 * @brief 	将简单符号L转换为简单数值L;使用符号分析中的L->ColCount为L分配空间，
 * 			并将L设为相等。
 * 
 * 			如果packed<0，那么这个例程创建另一个因子的副本(通过SparseCore_copy_factor)。
 * 			在这种情况下，空间没有初始化。
 * 
 * @param L 
 * @param to_ll 
 * @param packed 
 * @param to_xtype 
 * @param Common 
 */
static void simplicial_symbolic_to_simplicial_numeric
(
    sparse_factor *L,
    int to_ll,
    int packed,
    int to_xtype,
    sparse_common *Common
)
{
    double grow0, grow1, xlen, xlnz ;
    double *Lx, *Lz ;
    Int *Li, *Lp, *Lnz, *ColCount ;
    Int n, grow, grow2, p, j, lnz, len, ok, e ;

    if (!allocate_simplicial_numeric (L, Common))
    {
	return ;	/* 内存溢出 */
    }

    ColCount = L->ColCount ;
    Lnz = L->nz ;
    Lp = L->p ;
    ok = TRUE ;
    n = L->n ;

    if (packed < 0)
    {

	/* ------------------------------------------------------------------ */
	/* 由SparseCore_copy_factor用于分配一个因子对象的副本 */
	/* ------------------------------------------------------------------ */

	lnz = L->nzmax ;
	L->nzmax = 0 ;

    }
    else if (packed)
    {

	/* ------------------------------------------------------------------ */
	/* 压缩LDL'或者LL' */
	/* ------------------------------------------------------------------ */
	lnz = 0 ;
	for (j = 0 ; ok && j < n ; j++)
	{
	    /* 确保 1 < len < n-j */
	    len = ColCount [j] ;
	    len = MAX (1, len) ;
	    len = MIN (len, n-j) ;
	    lnz += len ;
	    ok = (lnz >= 0) ;
	}
	for (j = 0 ; j <= n ; j++)
	{
	    Lp [j] = j ;
	}
	for (j = 0 ; j < n ; j++)
	{
	    Lnz [j] = 1 ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* LDL'不压缩 */
	/* ------------------------------------------------------------------ */

	/* 计算新的lnzmax */
	/* 如果任何参数为NaN，则grow为假 */
	grow0 = Common->grow0 ;
	grow1 = Common->grow1 ;
	grow2 = Common->grow2 ;
	grow0 = IS_NAN (grow0) ? 1 : grow0 ;
	grow1 = IS_NAN (grow1) ? 1 : grow1 ;

	grow = (grow0 >= 1.0) && (grow1 >= 1.0) && (grow2 > 0) ;
	/* 对每个列初始化Lp和Lnz */
	lnz = 0 ;
	for (j = 0 ; ok && j < n ; j++)
	{
	    Lp [j] = lnz ;
	    Lnz [j] = 1 ;

	    /* 确保1 < len < n-j */
	    len = ColCount [j] ;
	    len = MAX (1, len) ;
	    len = MIN (len, n-j) ;

	    /* 以双精度计算len以避免整数溢出 */
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (Int) xlen ;
	    }
	    lnz += len ;
	    ok = (lnz >= 0) ;
	}
	if (ok)
	{
	    Lp [n] = lnz ;
	    if (grow)
	    {
		/* 添加额外的空间 */
		xlnz = (double) lnz ;
		xlnz *= grow0 ;
		xlnz = MIN (xlnz, Size_max) ;
		xlnz = MIN (xlnz, ((double) n * (double) n + (double) n) / 2) ;
		lnz = (Int) xlnz ;
	    }
	}
    }

    lnz = MAX (1, lnz) ;

    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
    }

    /* 指定 L->i, L->x, L->z */
    e = 1 ;
    if (!ok || !CORE(realloc_multiple) (lnz, 1, to_xtype, &(L->i), NULL,
		&(L->x), &(L->z), &(L->nzmax), Common))
    {
	L->p    = CORE(free) (n+1, sizeof (Int),      L->p, Common) ;
	L->nz   = CORE(free) (n,   sizeof (Int),      L->nz, Common) ;
	L->prev = CORE(free) (n+2, sizeof (Int),      L->prev, Common) ;
	L->next = CORE(free) (n+2, sizeof (Int),      L->next, Common) ;
	L->i    = CORE(free) (lnz, sizeof (Int),      L->i, Common) ;
	L->x    = CORE(free) (lnz, e*sizeof (double), L->x, Common) ;
	L->z    = CORE(free) (lnz, sizeof (double),   L->z, Common) ;
	return ;	/* 内存溢出 */
    }

    /* =========================将更改提交到L=====================  */

    /* 将L初始化为单位矩阵 */
    L->xtype = to_xtype ;
    L->dtype = DTYPE ;
    L->minor = n ;

    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;

    if (packed >= 0)
    {
	/* 为LL'或LDL'创建单位对角线 */

	switch (L->xtype)
	{
	    case SPARSE_REAL:
		for (j = 0 ; j < n ; j++)
		{
		    p = Lp [j] ;
		    Li [p] = j ;
		    Lx [p] = 1 ;
		}
		break ;
	}
    }

    L->is_ll = to_ll ;
}

/**
 * @brief 	将LL'改为LDL'， LDL'改为LL'，或者保持原样。
 * 
 * @param L 
 * @param to_ll 
 * @param to_packed 
 * @param to_monotonic 
 * @param Common 
 */
static void change_simplicial_numeric
(
    sparse_factor *L,
    int to_ll,
    int to_packed,
    int to_monotonic,
    sparse_common *Common
)
{
    double grow0, grow1, xlen, xlnz ;
    void *newLi, *newLx, *newLz ;
    double *Lx, *Lz ;
    Int *Lp, *Li, *Lnz ;
    Int make_monotonic, grow2, n, j, lnz, len, grow, ok, make_ll, make_ldl ;
    size_t nzmax0 ;

    make_monotonic = ((to_packed || to_monotonic) && !(L->is_monotonic)) ;
    make_ll  = (to_ll && !(L->is_ll)) ;
    make_ldl = (!to_ll && L->is_ll) ;

    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;
    Lnz = L->nz ;

    grow = FALSE ;
    grow0 = Common->grow0 ;
    grow1 = Common->grow1 ;
    grow2 = Common->grow2 ;
    grow0 = IS_NAN (grow0) ? 1 : grow0 ;
    grow1 = IS_NAN (grow1) ? 1 : grow1 ;
    ok = TRUE ;
    newLi = NULL ;
    newLx = NULL ; 
    newLz = NULL ; 
    lnz = 0 ;

    if (make_monotonic)
    {

	/* ------------------------------------------------------------------ */
	/* 列的顺序，但将重新排序和选择性压缩。 */
	/* ------------------------------------------------------------------ */
	/* 计算新的L->nzmax */
	if (!to_packed)
	{
	    /* 如果任何参数为NaN，则grow为false */
	    grow = (grow0 >= 1.0) && (grow1 >= 1.0) && (grow2 > 0) ;
	}
	for (j = 0 ; ok && j < n ; j++)
	{
	    len = Lnz [j] ;

	    /* 以双精度计算len以避免整数溢出 */
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (Int) xlen ;
	    }

	    lnz += len ;
	    ok = (lnz >= 0) ;
	}

	if (!ok)
	{
	    ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	    return ;
	}

	if (grow)
	{
	    xlnz = (double) lnz ;
	    xlnz *= grow0 ;
	    xlnz = MIN (xlnz, Size_max) ;
	    xlnz = MIN (xlnz, ((double) n * (double) n + (double) n) / 2) ;
	    lnz = (Int) xlnz ;
	}

	lnz = MAX (1, lnz) ;
	nzmax0 = 0 ;

	CORE(realloc_multiple) (lnz, 1, L->xtype, &newLi, NULL,
		&newLx, &newLz, &nzmax0, Common) ;

	if (Common->status < SPARSE_OK)
	{
	    return ;	    /* 内存溢出 */
	}
    }

    /* =========================将更改提交到L=====================  */

    /* ---------------------------------------------------------------------- */
    /* 使用模板例程转换简单的L */
    /* ---------------------------------------------------------------------- */

    switch (L->xtype)
    {

	case SPARSE_REAL:
	    r_change_simplicial_numeric (L, to_ll, to_packed,
		    newLi, newLx, newLz, lnz, grow, grow1, grow2,
		    make_ll, make_monotonic, make_ldl, Common) ;
	    break ;
    }
}

/**
 * @brief 	将一个超节点数值因数分解转换为任何简单数值分解。保留L->xtype不变(实)。
 * 
 * @param L 
 * @param to_packed 
 * @param to_ll 
 * @param Common 
 */
static void ll_super_to_simplicial_numeric
(
    sparse_factor *L,
    int to_packed,
    int to_ll,
    sparse_common *Common
)
{
    Int *Ls, *Lpi, *Lpx, *Super, *Li ;
    Int n, lnz, s, nsuper, psi, psend, nsrow, nscol, k1, k2, erows ;

    n = L->n ;
    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;

    /* Int溢出不能发生，因为超节点L已经存在 */

    if (to_packed)
    {
	/* 计算L中非零的个数，每个超节点都是这种形式
	 *
	 *    l	. . .	    以这个矩阵为例, nscol = 4 (# columns). nsrow = 9.
	 *    l l . .	    "." 在超节点因子中分配，但未使用。它们没有复制到简单因子。
	 *    l l l .	    如果对某些“l”和“e”项分别进行了严格的单形分解或重符号分解，
	 *    l l l l	    则由于数值消去和松弛的超节点合并，它们可能在数值上为零，
	 *    e e e e	    甚至在符号上为零。
	 *    e e e e	    
	 *    e e e e	    
	 *    e e e e	    
	 *    e e e e
	 */
	lnz = 0 ;
	for (s = 0 ; s < nsuper ; s++)
	{
	    k1 = Super [s] ;
	    k2 = Super [s+1] ;
	    psi = Lpi [s] ;
	    psend = Lpi [s+1] ;
	    nsrow = psend - psi ;
	    nscol = k2 - k1 ;
	    erows = nsrow - nscol ;

	    /* 下三角形部分，包括对角线，计算上图中的“l”项。 */
	    lnz += nscol * (nscol+1) / 2 ;

	    /* 矩形部分，对角线下的方块(“e”项) */
	    lnz += nscol * erows ;
	}
    }
    else
    {
	/* Li和Lx的大小相同 */
	lnz = L->xsize ;
    }

    Li = CORE(malloc) (lnz, sizeof (Int), Common) ;
    if (Common->status < SPARSE_OK)
    {
	return ;	/* 内存溢出 */
    }

    if (!allocate_simplicial_numeric (L, Common))
    {
	CORE(free) (lnz, sizeof (Int), Li, Common) ;
	return ;	/* 内存溢出 */
    }

    /* =============================将更改提交到L=================  */

    L->i = Li ;
    L->nzmax = lnz ;

    /* ---------------------------------------------------------------------- */
    /* 使用模板例程转换超节点L */
    /* ---------------------------------------------------------------------- */

    switch (L->xtype)
    {

	case SPARSE_REAL:
	    r_ll_super_to_simplicial_numeric (L, to_packed, to_ll, Common) ;
	    break ;
    }

    /* ---------------------------------------------------------------------- */
    /* 释放L的闲置部分 */
    /* ---------------------------------------------------------------------- */

    L->super = CORE(free) (nsuper+1, sizeof (Int), L->super, Common) ;
    L->pi    = CORE(free) (nsuper+1, sizeof (Int), L->pi, Common) ;
    L->px    = CORE(free) (nsuper+1, sizeof (Int), L->px, Common) ;
    L->s     = CORE(free) (L->ssize, sizeof (Int), L->s, Common) ;

    L->ssize = 0 ;
    L->xsize = 0 ;
    L->nsuper = 0 ;
    L->maxesize = 0 ;
    L->maxcsize = 0 ;

    L->is_super = FALSE ;
}

/**
 * @brief 	通过分配L->x将一个超节点符号分解转化为一个超节点数值分解。L->x的含量未定义。
 * 
 * @param to_xtype 
 * @param L 
 * @param Common 
 * @return int 
 */
static int super_symbolic_to_ll_super
(
    int to_xtype,
    sparse_factor *L,
    sparse_common *Common
)
{
    double *Lx ;
    Int wentry = (to_xtype == SPARSE_REAL) ? 1 : 2 ;
    Lx = CORE(malloc) (L->xsize, wentry * sizeof (double), Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;	/* 内存溢出 */
    }

    /* ==========================将更改提交到L====================  */

    if (L->xsize == 1)
    {
	/* 调用者不会期望访问这个条目，但是一些hnuchol例程可能会。将其设置为零，这样valgrind不会冲突 */
	switch (to_xtype)
	{
	    case SPARSE_REAL:
		Lx [0] = 0 ;
		break ;
	}
    }

    L->x = Lx ;
    L->xtype = to_xtype ;
    L->dtype = DTYPE ;
    L->minor = L->n ;
    return (TRUE) ;
}

/**
 * @brief 	转换一个因子L。有些转换只是分配未初始化的空间，以便以后填充。
 * 			
 * 			如果转换失败，该因子将保留其原始形式，但有一个例外。
 * 			将超节点符号因子转换为单纯数字因子(L=D=I)可以使该因子保持单纯符号形式。
 * 		
 * 			下面列出了为每个转换分配的内存。
 */
int CORE(change_factor)
(
    /* ---- input ---- */
    int to_xtype,		/* 转换为SPARSE_PATTERN, _REAL */
    int to_ll,			/* TRUE: LL', FALSE: LDL' */
    int to_super,		/* TRUE: 超节点, FALSE: 单节点 */
    int to_packed,		/* TRUE: 压缩单节点列, FALSE: 不压缩 */
    int to_monotonic,	/* TRUE: 单节点列排序, FALSE: 不排序 */
    /* ---- in/out --- */
    sparse_factor *L,	/* 待修改因子 */
    /* --------------- */
    sparse_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    if (to_xtype < SPARSE_PATTERN || to_xtype > SPARSE_REAL)
    {
	ERROR (SPARSE_INVALID, "xtype invalid") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* 确保所有参数为TRUE/FALSE */
    to_ll = BOOLEAN (to_ll) ;
    to_super = BOOLEAN (to_super) ;

    /* ---------------------------------------------------------------------- */
    /* 转换 */
    /* ---------------------------------------------------------------------- */

    if (to_xtype == SPARSE_PATTERN)
    {

	/* ------------------------------------------------------------------ */
	/* 转换为符号 */
	/* ------------------------------------------------------------------ */

	if (!to_super)
	{

	    /* -------------------------------------------------------------- */
	    /* 将任何因子转换为简单符号因子 */
	    /* -------------------------------------------------------------- */

	    any_to_simplicial_symbolic (L, to_ll, Common) ;    /* cannot fail */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* 转化为一个超节点符号因子 */
	    /* -------------------------------------------------------------- */

	    if (L->xtype != SPARSE_PATTERN && L->is_super)
	    {
		/* 从超节点数字到超节点符号的转换。这保留了L的符号模式，抛弃了数值 */
		ll_super_to_super_symbolic (L, Common) ;       /* 不能失败 */
	    }
	    else if (L->xtype == SPARSE_PATTERN && !(L->is_super))
	    {
		/* 从简单符号到超节点符号的转换。超节点模式的内容未初始化。不是为终端用户设计的。 */
		simplicial_symbolic_to_super_symbolic (L, Common) ;
	    }
	    else
	    {
		/* 不能从简单的数字到超节点符号的转换 */
		ERROR (SPARSE_INVALID,
			"cannot convert L to supernodal symbolic") ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* 转换为数值 */
	/* ------------------------------------------------------------------ */
	    
	if (to_super)
	{

	    /* -------------------------------------------------------------- */
	    /* 转换到超节点数值因子 */
	    /* -------------------------------------------------------------- */

	    if (L->xtype == SPARSE_PATTERN)
	    {
		if (L->is_super)
		{
		    /* 转换超节点符号到超节点数字。超节点数值的内容未初始化。
			 * 它由SparseCore_super_numeric使用。不是为终端用户设计的。 */
		    super_symbolic_to_ll_super (to_xtype, L, Common) ;
		}
		else
		{
		    /* 转换简单符号到超节点数字。内容没有定义。
			 * 这只由Core/SparseCore_copy_factor使用。不是为终端用户设计的。 */
		    if (!simplicial_symbolic_to_super_symbolic (L, Common))
		    {
			/* 失败了，就转换回简单的符号 */
			any_to_simplicial_symbolic (L, to_ll, Common) ;
		    }
		    else
		    {
			/* 转换到超级符号OK，分配数字部分 */
			super_symbolic_to_ll_super (to_xtype, L, Common) ;
		    }
		}
	    }
	    else
	    {
		/* 如果L已经是超节点数值形式，什么都不做 */
		if (!(L->is_super))
		{
		    ERROR (SPARSE_INVALID,
			"cannot convert simplicial L to supernodal") ;
		}
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* 将任何因子转换为单纯数值 */
	    /* -------------------------------------------------------------- */

	    if (L->xtype == SPARSE_PATTERN && !(L->is_super))
	    {

		/* ---------------------------------------------------------- */
		/* 将单纯符号转换为简单数值(L=I,D=I)*/
		/* ---------------------------------------------------------- */

		simplicial_symbolic_to_simplicial_numeric (L, to_ll, to_packed,
			to_xtype, Common) ;

	    }
	    else if (L->xtype != SPARSE_PATTERN && L->is_super)
	    {

		/* ---------------------------------------------------------- */
		/* 将一个超节点LL'转换为简单数值 */
		/* ---------------------------------------------------------- */

		ll_super_to_simplicial_numeric (L, to_packed, to_ll, Common) ;

	    }
	    else if (L->xtype == SPARSE_PATTERN && L->is_super)
	    {

		/* ---------------------------------------------------------- */
		/* 将超节点符号转换为单纯数值(L=D=I)*/
		/* ---------------------------------------------------------- */

		any_to_simplicial_symbolic (L, to_ll, Common) ;
		/* 如果下列方法失败，则将该因子保留为简单的符号形式 */
		simplicial_symbolic_to_simplicial_numeric (L, to_ll, to_packed,
			to_xtype, Common) ;

	    }
	    else
	    {

		/* ---------------------------------------------------------- */
		/* 改变一个简单的数值因子 */
		/* ---------------------------------------------------------- */

		/* 将LL'改为LDL'， LDL'改为LL'，或者保持原样。
		 * 打包L的列，或者保持原样。确保列是单调的，或者保持原样。 */

		change_simplicial_numeric (L, to_ll, to_packed, to_monotonic,
			Common) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */

    return (Common->status >= SPARSE_OK) ;
}
