/**
 * @file SparseCore_t_transpose.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-19
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_template.h"

/**
 * @brief 	计算F = A', A(:,f)', A(p,f)'
 * 			其中A是非对称的，并且F已经分配好
 * 
 */
static int TEMPLATE (SparseCore_transpose_unsym)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 转置的矩阵 */
    Int *Perm,			/* 大小为nrow，可以为空 */
    Int *fset,			/* 0:(A->ncol)-1的子集 */
    Int nf,				/* 子集的大小 */
    /* ---- output --- */
    sparse_csc *F,	/* F = A', A(:,f)', A(p,f)' */
    /* --------------- */
    sparse_common *Common
)
{
    double *Ax, *Az, *Fx, *Fz ;
    Int *Ap, *Anz, *Ai, *Fp, *Fnz, *Fj, *Wi, *Iwork ;
    Int j, p, pend, nrow, ncol, Apacked, use_fset, fp, Fpacked, jj, permute ;

    /* 检查输入，确保A和F的xtype匹配(如果是模式版本则忽略) */
    if (!XTYPE_OK (A->xtype))
    {
	ERROR (SPARSE_INVALID, "real/complex mismatch") ;
	return (FALSE) ;
    }

    /*
     * 输入
     */
    use_fset = (fset != NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;

    Ap = A->p ;		/* A的列指针，大小为A->ncol+1 */
    Ai = A->i ;		/* A的行索引，大小为nz = Ap [A->ncol] */
    Ax = A->x ;		/* A的实际值，大小为nz */
    Az = A->z ;		/* A的imag values，大小为nz */
    Anz = A->nz ;
    Apacked = A->packed ;

    permute = (Perm != NULL) ;

    Fp = F->p ;		/* F的行指针，大小为A->nrow+1 */
    Fj = F->i ;		/* F的列索引，大小为nz */
    Fx = F->x ;		/* F的实际值，大小为nz */
    Fz = F->z ;		/* F的imag values，大小为nz */
    Fnz = F->nz ;
    Fpacked = F->packed ;

    nf = (use_fset) ? nf : ncol ;

    /* 
     * 获取工作空间 
     */
    Iwork = Common->Iwork ;
    Wi = Iwork ;		/* 大小为nrow (i/l/l) */

    /* 
     * 转置
     */
    for (jj = 0 ; jj < nf ; jj++)
    {
	j = (use_fset) ? (fset [jj]) : jj ;
	p = Ap [j] ;
	pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	for ( ; p < pend ; p++)
	{
	    fp = Wi [Ai [p]]++ ;
	    Fj [fp] = j ;
	    ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
	}
    }

    return (TRUE) ;
}


/**
 * @brief 	计算F=A',A(p,p)'
 * 			其中A是对称的，并且F已经分配好
 * 
 */
static int TEMPLATE (SparseCore_transpose_sym)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 转置的矩阵 */
    Int *Perm,			/* 大小为n，可以为空 */
    /* ---- output --- */
    sparse_csc *F,	/* F = A', A(p,p)' */
    /* --------------- */
    sparse_common *Common
)
{
    double *Ax, *Az, *Fx, *Fz ;
    Int *Ap, *Anz, *Ai, *Fp, *Fj, *Wi, *Pinv, *Iwork ;
    Int p, pend, packed, fp, upper, permute, jold, n, i, j, iold ;

    /* 检查输入，确定A和F的xtype是相同的(模式版本时忽略) */
    if (!XTYPE_OK (A->xtype))
    {
	ERROR (SPARSE_INVALID, "real/complex mismatch") ;
	return (FALSE) ;
    }

    /* 
	 * 输入
	 */
    permute = (Perm != NULL) ;
    n = A->nrow ;
    Ap = A->p ;		/* A的列指针，大小为A->ncol+1 */
    Ai = A->i ;		/* A的行索引，大小为nz = Ap [A->ncol] */
    Ax = A->x ;		/* A的实际值，大小为nz */
    Az = A->z ;		/* A的imag values，大小为nz */
    Anz = A->nz ;
    packed = A->packed ;
    upper = (A->stype > 0) ;

    Fp = F->p ;		/* F的行指针，大小为A->nrow+1 */
    Fj = F->i ;		/* F的列索引，大小为nz */
    Fx = F->x ;		/* F的实际值，大小为nz */
    Fz = F->z ;		/* F的imag values，大小为nz */

	/*
	 * 获取工作空间
	 */
    Iwork = Common->Iwork ;
    Wi = Iwork ;		/* 大小为n (i/l/l) */
    Pinv = Iwork + n ;	/* 如果perm为空，则不使用，大小为n (i/i/l) */

    /* 
    /* 转置
    */
    if (permute)
    {
	if (upper)
	{
	    /* 排列上部分 */
	    for (j = 0 ; j < n ; j++)
	    {
		jold = Perm [j] ;
		p = Ap [jold] ;
		pend = (packed) ? Ap [jold+1] : p + Anz [jold] ;
		for ( ; p < pend ; p++)
		{
		    iold = Ai [p] ;
		    if (iold <= jold)
		    {
			i = Pinv [iold] ;
			if (i < j)
			{
			    fp = Wi [i]++ ;
			    Fj [fp] = j ;
			    ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
			}
			else
			{
			    fp = Wi [j]++ ;
			    Fj [fp] = i ;
			    ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
			}
		    }
		}
	    }
	}
	else
	{
	    /* 排列下部分 */
	    for (j = 0 ; j < n ; j++)
	    {
		jold = Perm [j] ;
		p = Ap [jold] ;
		pend = (packed) ? Ap [jold+1] : p + Anz [jold] ;
		for ( ; p < pend ; p++)
		{
		    iold = Ai [p] ;
		    if (iold >= jold)
		    {
			i = Pinv [iold] ;
			if (i > j)
			{
			    fp = Wi [i]++ ;
			    Fj [fp] = j ;
			    ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
			}
			else
			{
			    fp = Wi [j]++ ;
			    Fj [fp] = i ;
			    ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
			}
		    }
		}
	    }
	}
    }
    else
    {
	if (upper)
	{
	    /* 不排列 上部分 */
	    for (j = 0 ; j < n ; j++)
	    {
		p = Ap [j] ;
		pend = (packed) ? Ap [j+1] : p + Anz [j] ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i <= j)
		    {
			fp = Wi [i]++ ;
			Fj [fp] = j ;
			ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
		    }
		}
	    }
	}
	else
	{
	    /* 不排列 下部分 */
	    for (j = 0 ; j < n ; j++)
	    {
		p = Ap [j] ;
		pend = (packed) ? Ap [j+1] : p + Anz [j] ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i >= j)
		    {
			fp = Wi [i]++ ;
			Fj [fp] = j ;
			ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
		    }
		}
	    }
	}
    }

    return (TRUE) ;
}

#undef PATTERN
#undef REAL
