/**
 * @file SparseChol_super_solve.c
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

/* ========================================================================== */
/* === 模板 ============================================================= */
/* ========================================================================== */

#define REAL
#include "SparseChol_t_super_solve.c"

/**
 * @brief 
 * 
 * @return int 成功则为TRUE，发生BLAS溢出则为FLASE
 */
int SparseChol_super_lsolve
(
    /* ---- input ---- */
    sparse_factor *L,	/* 用于正向求解的因子 */
    /* ---- output ---- */
    dense_array *X,	/* 输入为b，输出为Lx=b的解 */
    /* ---- workspace ---- */
    dense_array *E,	/* 大小为nrhs*(L->maxesize)的工作空间 */
    /* --------------- */
    sparse_common *Common
)
{
    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    if (L->xtype != X->xtype)
    {
	ERROR (SPARSE_INVALID, "L and X must have the same xtype") ;
	return (FALSE) ;
    }
    if (L->xtype != E->xtype)
    {
	ERROR (SPARSE_INVALID, "L and E must have the same xtype") ;
	return (FALSE) ;
    }
    if (X->d < X->nrow || L->n != X->nrow)
    {
	ERROR (SPARSE_INVALID, "X and L dimensions must match") ;
	return (FALSE) ;
    }
    if (E->nzmax < X->ncol * (L->maxesize))
    {
	ERROR (SPARSE_INVALID, "workspace E not large enough") ;
	return (FALSE) ;
    }
    if (!(L->is_ll) || !(L->is_super))
    {
	ERROR (SPARSE_INVALID, "L not supernodal") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;
    if (L->n == 0 || X->ncol == 0)
    {
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 使用模板例程求解Lx=b */
    /* ---------------------------------------------------------------------- */

    switch (L->xtype)
    {

	case SPARSE_REAL:
	    r_SparseChol_super_lsolve (L, X, E, Common) ;
	    break ;

    }

    if (CHECK_BLAS_INT && !Common->blas_ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large for the BLAS") ;
    }
    return (Common->blas_ok) ;
}

/**
 * @brief   解L'x=b，其中x和b的大小为n×nrhs。b被解x覆盖，
 *          在输入时b以共主顺序存储，前导维数为d，在输出时x以共主顺序存储。
 *          
 *          工作空间E的内容在输入和输出上都没有定义。
 *          不需要额外工作空间
 * 
 * @return int 成功则为TRUE，发生BLAS溢出则为FLASE
 */
int SparseChol_super_ltsolve
(
    /* ---- input ---- */
    sparse_factor *L,	/* 用于反向求解的因子 */
    /* ---- output ---- */
    dense_array *X,	/* 输入为b,输出为L'x=b的解 */
    /* ---- workspace ---- */
    dense_array *E,	/* 工作空间大小为nrhs*(L->maxesize) */
    /* --------------- */
    sparse_common *Common
)
{
    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    
    if (L->xtype != X->xtype)
    {
	ERROR (SPARSE_INVALID, "L and X must have the same xtype") ;
	return (FALSE) ;
    }
    if (L->xtype != E->xtype)
    {
	ERROR (SPARSE_INVALID, "L and E must have the same xtype") ;
	return (FALSE) ;
    }
    if (X->d < X->nrow || L->n != X->nrow)
    {
	ERROR (SPARSE_INVALID, "X and L dimensions must match") ;
	return (FALSE) ;
    }
    if (E->nzmax < X->ncol * (L->maxesize))
    {
	ERROR (SPARSE_INVALID, "workspace E not large enough") ;
	return (FALSE) ;
    }
    if (!(L->is_ll) || !(L->is_super))
    {
	ERROR (SPARSE_INVALID, "L not supernodal") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;
    
    if (L->n == 0 || X->ncol == 0)
    {
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 使用模板例程求解Lx=b */
    /* ---------------------------------------------------------------------- */

    switch (L->xtype)
    {

	case SPARSE_REAL:
	    r_SparseChol_super_ltsolve (L, X, E, Common) ;
	    break ;

	
    }

    if (CHECK_BLAS_INT && !Common->blas_ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large for the BLAS") ;
    }
    return (Common->blas_ok) ;
}
#endif
#endif