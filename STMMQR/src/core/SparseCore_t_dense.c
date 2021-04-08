/**
 * @file t_hnucore_dense.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-20
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_template.h"

/**
 * @brief 将稀疏矩阵转换成稠密矩阵
 * 
 */
static dense_array *TEMPLATE (SparseCore_sparse_to_dense)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要转换的稀疏矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Ax, *Xx, *Az, *Xz ;
    Int *Ap, *Ai, *Anz ;
    dense_array *X ;
    Int i, j, p, pend, nrow, ncol, packed ;

    /*
     * 得到输入
     */
    nrow = A->nrow ;
    ncol = A->ncol ;
    packed = A->packed ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;

    /*
     * 分配结果
     */
    X = SparseCore_zeros (nrow, ncol, XTYPE2, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 超出内存 */
    }
    Xx = X->x ;
    Xz = X->z ;

    /*
     * 复制到稠密矩阵之中
     */
    if (A->stype < 0)
    {
	/* A是一个存储下三角部分的对称矩阵，但是X的两个部分都在 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i >= j)
		{
		    ASSIGN2 (Xx, Xz, i+j*nrow, Ax, Az, p) ;
		    ASSIGN2_CONJ (Xx, Xz, j+i*nrow, Ax, Az, p) ;
		}
	    }
	}
    }
    else if (A->stype > 0)
    {
	/* A是一个存储上三角部分的对称矩阵，但是X的两个部分都在 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i <= j)
		{
		    ASSIGN2 (Xx, Xz, i+j*nrow, Ax, Az, p) ;
		    ASSIGN2_CONJ (Xx, Xz, j+i*nrow, Ax, Az, p) ;
		}
	    }
	}
    }
    else
    {
	/* A和X的两个部分都被给出 */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		ASSIGN2 (Xx, Xz, i+j*nrow, Ax, Az, p) ;
	    }
	}
    }

    return (X) ;
}


#ifndef PATTERN

/* SPARSE_PATTERN没有xtype的稠密矩阵 */

/**
 * @brief 稠密矩阵转换成稀疏矩阵
 * 
 */
static sparse_csc *TEMPLATE (SparseCore_dense_to_sparse)
(
    /* ---- input ---- */
    dense_array *X,	/* 需要转换的稠密矩阵 */
    int values,		    /* 如果要复制为TRUE，否则为FALSE */
    /* --------------- */
    sparse_common *Common
)
{
    double *Xx, *Cx, *Xz, *Cz ;
    Int *Ci, *Cp ;
    sparse_csc *C ;
    Int i, j, p, d, nrow, ncol, nz ;

    /*
     * 得到输入
     */
    nrow = X->nrow ;
    ncol = X->ncol ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;

    /*
     * 计算结果中的非零元数量
     */
    nz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    if (ENTRY_IS_NONZERO (Xx, Xz, i+j*d))
	    {
		nz++ ;
	    }
	}
    }

    /*
     * 为结果C分配内存空间
     */
    C = SparseCore_allocate_sparse (nrow, ncol, nz, TRUE, TRUE, 0,
	    values ? XTYPE : SPARSE_PATTERN, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 超出内存 */
    }
    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;
    Cz = C->z ;

    /*
     * 将稠密矩阵X拷贝到稀疏矩阵C中
     */
    p = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	Cp [j] = p ;
	for (i = 0 ; i < nrow ; i++)
	{
	    if (ENTRY_IS_NONZERO (Xx, Xz, i+j*d))
	    {
		Ci [p] = i ;
		if (values)
		{
		    ASSIGN (Cx, Cz, p, Xx, Xz, i+j*d) ;
		}
		p++ ;
	    }
	}
    }
    Cp [ncol] = nz ;

    return (C) ;
}

/**
 * @brief Y = X X和Y均已经分配好空间
 * 
 */
static int TEMPLATE (SparseCore_copy_dense2)
(
    /* ---- input ---- */
    dense_array *X,	/* 需要拷贝的稠密矩阵X */
    /* ---- output --- */
    dense_array *Y	/* 稠密矩阵X的复制 */
)
{
    double *Xx, *Xz, *Yx, *Yz ;
    Int i, j, nrow, ncol, dy, dx ;

    /*
     * 得到输入
     */
    Xx = X->x ;
    Xz = X->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    dx = X->d ;
    dy = Y->d ;
    nrow = X->nrow ;
    ncol = X->ncol ;

    /*
     * 拷贝
     */
    CLEAR (Yx, Yz, 0) ;
    for (j = 0 ; j < ncol ; j++)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    ASSIGN (Yx, Yz, i+j*dy, Xx, Xz, i+j*dx) ;
	}
    }
    return (TRUE) ;
}

#endif

#undef PATTERN
#undef REAL
