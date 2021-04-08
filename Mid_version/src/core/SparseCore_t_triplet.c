/**
 * @file SparseCore_t_triplet.c
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
 * @brief 三元组转换为稀疏矩阵
 * 
 */
static size_t TEMPLATE (SparseCore_triplet_to_sparse)
(
    sparse_triplet *T,	/* [in] 	输入的三元组 */
    sparse_csc *R,	/* [in/out]	由三元组转换得到的稀疏矩阵 */
    sparse_common *Common
)
{
    double *Rx, *Rz, *Tx, *Tz ;
    Int *Wj, *Rp, *Ri, *Rnz, *Ti, *Tj  ;
    Int i, j, p, p1, p2, pdest, pj, k, stype, nrow, ncol, nz ;
    size_t anz ;

    /* 
	 * 在输入时使用Wj(i/l/l)作为临时指针包含Rp的副本 
	 */
    Wj = Common->Iwork ;	/* 大小为MAX (nrow,ncol). */

    Rp = R->p ;
    Ri = R->i ;
    Rnz = R->nz ;
    Rx = R->x ;
    Rz = R->z ;

    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    Tz = T->z ;
    nz = T->nnz ;
    nrow = T->nrow ;
    ncol = T->ncol ;
    stype = SIGN (T->stype) ;

	/*
	 * 按行结构存储
	 */
    /* 如果Ti是无序的，这部分将占用运行时间 */
    if (stype > 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < j)
	    {
		/* 将三元组(j,i,x)存入R的第i列 */
		p = Wj [i]++ ;
		Ri [p] = j ;
	    }
	    else
	    {
		/* 将三元组(i,j,x)存入R的第j列 */
		p = Wj [j]++ ;
		Ri [p] = i ;
	    }
	    ASSIGN (Rx, Rz, p, Tx, Tz, k) ;
	}
    }
    else if (stype < 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i > j)
	    {
		/* 将三元组(j,i,x)存入R的第i列 */
		p = Wj [i]++ ;
		Ri [p] = j ;
	    }
	    else
	    {
		/* 将三元组(i,j,x)存入R的第j列 */
		p = Wj [j]++ ;
		Ri [p] = i ;
	    }
	    ASSIGN (Rx, Rz, p, Tx, Tz, k) ;
	}
    }
    else
    {
	for (k = 0 ; k < nz ; k++)
	{
	    /* 将三元组(i,j,x)存入R的第i列 */
	    p = Wj [Ti [k]]++ ;
	    Ri [p] = Tj [k] ;
	    ASSIGN (Rx, Rz, p, Tx, Tz, k) ;
	}
    }

    /* 
	 * 使用大小为ncol的Wj (i/l/l)来跟踪每行中的重复项 
	 */
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = EMPTY ;
    }

    anz = 0 ;
    for (i = 0 ; i < nrow ; i++)
    {
	p1 = Rp [i] ;
	p2 = Rp [i+1] ;
	pdest = p1 ;
	/* 此时Wj[j] < p1对所有的列j都成立，因为Ri/Rx是以行方式存储的 */
	for (p = p1 ; p < p2 ; p++)
	{
	    j = Ri [p] ;
	    pj = Wj [j] ;
	    if (pj >= p1)
	    {
		/* 这个列索引j已经在第i行的位置pj处；
		 * 把重复项加起来 */
		/* Rx [pj] += Rx [p] ; */
		ASSEMBLE (Rx, Rz, pj, Rx, Rz, p) ;
	    }
	    else
	    {
		/* 保留条目并在Wj[j]中跟踪上面的情况 */
		Wj [j] = pdest ;
		if (pdest != p)
		{
		    Ri [pdest] = j ;
		    ASSIGN (Rx, Rz, pdest, Rx, Rz, p) ;
		}
		pdest++ ;
	    }
	}
	Rnz [i] = pdest - p1 ;
	anz += (pdest - p1) ;
    }

    return (anz) ;
}

#undef PATTERN
#undef REAL