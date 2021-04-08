/**
 * @file SparseCore_norm.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-20
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_internal.h"
#include "SparseCore.h"

#define REAL
#include "SparseCore_t_sdmult.c"

/**
 * @brief 计算绝对值
 * 
 * @param xtype 
 * @param Ax 
 * @param Az 
 * @param p 
 * @param Common 
 * @return double 
 */
static double abs_value
(
    int xtype,
    double *Ax,
    double *Az,
    Int p,
    sparse_common *Common
)
{
    double s = 0 ;
    switch (xtype)
    {
	case SPARSE_PATTERN:
	    s = 1 ;
	    break ;

	case SPARSE_REAL:
	    s = fabs (Ax [p]) ;
	    break ;

    }
    return (s) ;
}

/**
 * @brief 稠密矩阵的范数
 * 
 */
double SparseCore_norm_dense
(
    /* ---- input ---- */
    dense_array *X,	/* 需要计算范数的稠密矩阵 */
    int norm,			/* 0，1，2对应范数 */
    /* --------------- */
    sparse_common *Common
)
{
    double xnorm, s, x, z ;
    double *Xx, *Xz, *W ;
    Int nrow, ncol, d, i, j, use_workspace, xtype ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    
    Common->status = SPARSE_OK ;
    ncol = X->ncol ;
    if (norm < 0 || norm > 2 || (norm == 2 && ncol > 1))
    {
		return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    nrow = X->nrow ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;
    xtype = X->xtype ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    W = NULL ;
    use_workspace = (norm == 0 && ncol > 4) ;
    if (use_workspace)
    {
		SparseCore_allocate_work (0, 0, nrow, Common) ;
		W = Common->Xwork ;
		if (Common->status < SPARSE_OK)
		{
			/* 无工作空间 */
			use_workspace = FALSE ;
		}
    }


    /* ---------------------------------------------------------------------- */
    /* 计算范数 */
    /* ---------------------------------------------------------------------- */

    xnorm = 0 ;

    if (use_workspace)
    {

	/* ------------------------------------------------------------------ */
	/* 无穷大范数=使用stride-1访问X的最大行的和 */
	/* ------------------------------------------------------------------ */

	/* 这比stride-d快，但是需要O(nrow)的工作空间 */
	for (j = 0 ; j < ncol ; j++)
	{
	    for (i = 0 ; i < nrow ; i++)
	    {
		W [i] += abs_value (xtype, Xx, Xz, i+j*d, Common) ;
	    }
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    s = W [i] ;
	    if ((IS_NAN (s) || s > xnorm) && !IS_NAN (xnorm))
	    {
		xnorm = s ;
	    }
	    W [i] = 0 ;
	}

    }
    else if (norm == 0)
    {

	/* ------------------------------------------------------------------ */
	/* 无穷大范数=使用stride-d访问X的最大行的和 */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < nrow ; i++)
	{
	    s = 0 ;
	    for (j = 0 ; j < ncol ; j++)
	    {
		s += abs_value (xtype, Xx, Xz, i+j*d, Common) ;
	    }
	    if ((IS_NAN (s) || s > xnorm) && !IS_NAN (xnorm))
	    {
		xnorm = s ;
	    }
	}

    }
    else if (norm == 1)
    {

	/* ------------------------------------------------------------------ */
	/* 一范数 最大列和 */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < ncol ; j++)
	{
	    s = 0 ;
	    for (i = 0 ; i < nrow ; i++)
	    {
		s += abs_value (xtype, Xx, Xz, i+j*d, Common) ;
	    }
	    if ((IS_NAN (s) || s > xnorm) && !IS_NAN (xnorm))
	    {
		xnorm = s ;
	    }
	}
    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* 二范数sqrt (sum (X.^2)) */
	/* ------------------------------------------------------------------ */

	switch (xtype)
	{

	    case SPARSE_REAL:
		for (i = 0 ; i < nrow ; i++)
		{
		    x = Xx [i] ;
		    xnorm += x*x ;
		}
		break ; 
        default:
            printf("NOT A REAL MATRIX !\n");
            break;
	}   

	xnorm = sqrt (xnorm) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */

    return (xnorm) ;
}


/**
 * @brief 计算稀疏矩阵的范数
 * 
 */
double SparseCore_norm_sparse
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要计算范数的稀疏矩阵 */
    int norm,			/* norm: 0: inf. norm, 1: 1-norm */
    /* --------------- */
    sparse_common *Common
)
{
    double anorm, s ;
    double *Ax, *Az, *W ;
    Int *Ap, *Ai, *Anz ;
    Int i, j, p, pend, nrow, ncol, packed, xtype ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    
    Common->status = SPARSE_OK ;
    ncol = A->ncol ;
    nrow = A->nrow ;
    if (norm < 0 || norm > 1)
    {
	
	return (EMPTY) ;
    }
    if (A->stype && nrow != ncol)
    {
	
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    packed = A->packed ;
    xtype = A->xtype ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    W = NULL ;
    if (A->stype || norm == 0)
    {
	SparseCore_allocate_work (0, 0, nrow, Common) ;
	W = Common->Xwork ;
	if (Common->status < SPARSE_OK)
	{
	    /* 内存溢出 */
	    return (EMPTY) ;
	}
	
    }

    /* ---------------------------------------------------------------------- */
    /* 计算范数 */
    /* ---------------------------------------------------------------------- */

    anorm = 0 ;

    if (A->stype > 0)
    {

	/* ------------------------------------------------------------------ */
	/* A是对称的，存储上三角形部分 */
	/* ------------------------------------------------------------------ */

	/* 无穷范数 = 1范数 = 最大行（列）和 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		s = abs_value (xtype, Ax, Az, p, Common) ;
		if (i == j)
		{
		    W [i] += s ;
		}
		else if (i < j)
		{
		    W [i] += s ;
		    W [j] += s ;
		}
	    }
	}

    }
    else if (A->stype < 0)
    {

	/* ------------------------------------------------------------------ */
	/* A是对称的，存储下三角形部分 */
	/* ------------------------------------------------------------------ */

	/* 无穷范数 = 1范数 = 最大行（列）和 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		s = abs_value (xtype, Ax, Az, p, Common) ;
		if (i == j)
		{
		    W [i] += s ;
		}
		else if (i > j)
		{
		    W [i] += s ;
		    W [j] += s ;
		}
	    }
	}

    }
    else if (norm == 0)
    {

	/* ------------------------------------------------------------------ */
	/* A是非对称的，计算无穷范数 */
	/* ------------------------------------------------------------------ */

	/* 无穷范数 = 最大行和 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		W [Ai [p]] += abs_value (xtype, Ax, Az, p, Common) ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* A是非对称的，计算1范数 */
	/* ------------------------------------------------------------------ */

	/* 1范数 = 最大列和 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    if (xtype == SPARSE_PATTERN)
	    {
		s = pend - p ;
	    }
	    else
	    {
		s = 0 ;
		for ( ; p < pend ; p++)
		{
		    s += abs_value (xtype, Ax, Az, p, Common) ;
		}
	    }
	    if ((IS_NAN (s) || s > anorm) && !IS_NAN (anorm))
	    {
		anorm = s ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 计算最大行和 */
    /* ---------------------------------------------------------------------- */

    if (A->stype || norm == 0)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    s = W [i] ;
	    if ((IS_NAN (s) || s > anorm) && !IS_NAN (anorm))
	    {
		anorm = s ;
	    }
	    W [i] = 0 ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */

    return (anorm) ;
}


/**
 * @brief 	稀疏矩阵 乘 稠密矩阵 ：
 * 			Y = alpha*(A*X) + beta*Y or Y = alpha*(A'*X) + beta*Y,
 * 			其中A是稀疏矩阵，X、Y是稠密矩阵.
 * 
 */
int SparseCore_sdmult
(
    /* ---- input ---- */
    sparse_csc *A,	/* 进行乘法的稀疏矩阵 */
    int transpose,		/* A是否转置：0代表不转置，1表示转置 */
    double alpha [2],   /* A的标量因子 */
    double beta [2],    /* Y的标量因子 */
    dense_array *X,	/* 进行乘法的稠密矩阵 */
    /* ---- in/out --- */
    dense_array *Y,	/* 得到的结果——稠密矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    double *w ;
    size_t nx, ny ;
    Int e ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    ny = transpose ? A->ncol : A->nrow ;	/* Y的长度 */
    nx = transpose ? A->nrow : A->ncol ;	/* X的长度 */
    if (X->nrow != nx || X->ncol != Y->ncol || Y->nrow != ny)
    {
	/* X，Y的维度错误 */
	return (FALSE) ;
    }
    if (A->xtype != X->xtype || A->xtype != Y->xtype)
    {
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    w = NULL ;
    e = (A->xtype == SPARSE_REAL ? 1:2) ;
    if (A->stype && X->ncol >= 4)
    {
	w = SparseCore_malloc (nx, 4*e*sizeof (double), Common) ;
    }
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* Y = alpha*op(A)*X + beta*Y通过模板程序计算*/
    /* ---------------------------------------------------------------------- */

    switch (A->xtype)
    {

	case SPARSE_REAL:
	    r_SparseCore_sdmult (A, transpose, alpha, beta, X, Y, w) ;
	    break ;
    }

    /* ---------------------------------------------------------------------- */
    /* 释放工作空间 */
    /* ---------------------------------------------------------------------- */

    SparseCore_free (4*nx, e*sizeof (double), w, Common) ;

    return (TRUE) ;
}

