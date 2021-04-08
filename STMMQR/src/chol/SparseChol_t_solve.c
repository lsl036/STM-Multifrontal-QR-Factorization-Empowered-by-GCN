/**
 * @file t_SparseChol_solve.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-22
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_template.h"

/* LL': 求解非单位对角的Lx=b */
#define LL
#include "SparseChol_t_lsolve.c"

/* LDL':求解LDx=b */
#define LD
#include "SparseChol_t_lsolve.c"

/* LDL':求解单位对角的Lx=b */
#include "SparseChol_t_lsolve.c"

/* LL':求解非单位对角的L'x=b */
#define LL
#include "SparseChol_t_ltsolve.c"

/* LDL':求解DL'x=b */
#define LD
#include "SparseChol_t_ltsolve.c"

/* LDL':求解单位对角的L'x=b */
#include "SparseChol_t_ltsolve.c"

/**
 * @brief 	解出Dx=b的LDL'分解，其中Y在输入上为b'，在输出上为x'。
 * 			右手边的数量(nrhs)并不受限制，即使是Yseti存在。
 */
static void TEMPLATE (ldl_dsolve)
(
    sparse_factor *L,
    dense_array *Y,		/* n*n，前导维度为nr */
    Int *Yseti, Int ysetlen
)
{
    double d [1] ;
    double *Lx, *Yx, *Yz ;
    Int *Lp ;
    Int n, nrhs, k, p, k1, k2, kk, kkiters ;

    nrhs = Y->nrow ;
    n = L->n ;
    Lp = L->p ;
    Lx = L->x ;
    Yx = Y->x ;
    Yz = Y->z ;
    kkiters = Yseti ? ysetlen : n ;
    for (kk = 0 ; kk < kkiters ; kk++)
    {
        k = Yseti ? Yseti [kk] : kk ;
	k1 = k*nrhs ;
	k2 = (k+1)*nrhs ;
	ASSIGN_REAL (d,0, Lx,Lp[k]) ;
	for (p = k1 ; p < k2 ; p++)
	{
	    DIV_REAL (Yx,Yz,p, Yx,Yz,p, d,0) ;
	}
    }
}

/**
 * @brief 	解一个线性系统，其中Y'包含输入的右边(数组转置)和输出的解。
 * 			没有应用排列;这些肯定已经应用到Y上了。
 * 
 * 			Yseti [0..ysetlen-1]是SparseChol_lsolve_pattern的可选索引列表。
 * 			求解只在与Yseti中的项相对应的L列上执行。如果为NULL则忽略。
 * 			如果存在，大多数函数要求Y'由单个密集列组成。
 */
static void TEMPLATE (simplicial_solver)
(
    int sys,		    	/* 求解系统 */
    sparse_factor *L,	    /* 使用的分解，简单的LL'或LDL' */
    dense_array *Y,	    /* 输入右侧，输出解 */
    Int *Yseti, Int ysetlen
)
{
    if (L->is_ll)
    {
	/* 分解为LL' */
	if (sys == SPARSE_A || sys == SPARSE_LDLt)
	{
	    /* 求解Ax=b或者LL'x=b */
	    TEMPLATE (ll_lsolve_k) (L, Y, Yseti, ysetlen) ;
	    TEMPLATE (ll_ltsolve_k) (L, Y, Yseti, ysetlen) ;
	}
	else if (sys == SPARSE_L || sys == SPARSE_LD)
	{
	    /* 求解Lx=b */
	    TEMPLATE (ll_lsolve_k) (L, Y, Yseti, ysetlen) ;
	}
	else if (sys == SPARSE_Lt || sys == SPARSE_DLt)
	{
	    /* 求解L'x=b */
	    TEMPLATE (ll_ltsolve_k) (L, Y, Yseti, ysetlen) ;
	}
    }
    else
    {
	/* 分解为LDL' */
	if (sys == SPARSE_A || sys == SPARSE_LDLt)
	{
	    /* 求解Ax=b或者LDL'x=b */
	    TEMPLATE (ldl_lsolve_k) (L, Y, Yseti, ysetlen) ;
	    TEMPLATE (ldl_dltsolve_k) (L, Y, Yseti, ysetlen) ;
	}
	else if (sys == SPARSE_LD)
	{
	    /* 求解LDx=b */
	    TEMPLATE (ldl_ldsolve_k) (L, Y, Yseti, ysetlen) ;
	}
	else if (sys == SPARSE_L)
	{
	    /* 求解Lx=b */
	    TEMPLATE (ldl_lsolve_k) (L, Y, Yseti, ysetlen) ;
	}
	else if (sys == SPARSE_Lt)
	{
	    /* 求解L'x=b */
	    TEMPLATE (ldl_ltsolve_k) (L, Y, Yseti, ysetlen) ;
	}
	else if (sys == SPARSE_DLt)
	{
	    /* 求解DL'x=b */
	    TEMPLATE (ldl_dltsolve_k) (L, Y, Yseti, ysetlen) ;
	}
	else if (sys == SPARSE_D)
	{
	    /* 求解Dx=b */
	    TEMPLATE (ldl_dsolve) (L, Y, Yseti, ysetlen) ;
	}
    }
}

#undef PATTERN
#undef REAL