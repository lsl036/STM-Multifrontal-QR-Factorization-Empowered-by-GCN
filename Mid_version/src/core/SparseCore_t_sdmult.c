/**
 * @file SparseCore_t_sdmult.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-19
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_template.h"

#undef ADVANCE
#ifdef REAL
#define ADVANCE(x,z,d) x += d
#endif

/**
 * @brief 稀疏矩阵与稠密矩阵相乘
 * 
 */
static void TEMPLATE (SparseCore_sdmult)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 做乘法的稀疏矩阵 */
    int transpose,		/* A是否转置，0不转置 */
    double alpha [2],   /* A的标量因子 */
    double beta [2],    /* Y的标量因子 */
    dense_array *X,	/* 做乘法的稠密矩阵 */
    /* ---- in/out --- */
    dense_array *Y,	/* 结果稠密矩阵 */
    /* -- workspace -- */
    double *W			/* 大小为4*nx，如果需要，是复数的两倍 */
)
{

    double yx [8], xx [8], ax [2] ;


    double *Ax, *Az, *Xx, *Xz, *Yx, *Yz, *w, *Wz ;
    Int *Ap, *Ai, *Anz ;
    size_t nx, ny, dx, dy ;
    Int packed, nrow, ncol, j, k, p, pend, kcol, i ;

	/*
	 * 得到输入
	 */
    ny = transpose ? A->ncol : A->nrow ;	/* Y的长度 */
    nx = transpose ? A->nrow : A->ncol ;	/* X的长度 */

    nrow = A->nrow ;
    ncol = A->ncol ;

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    Az  = A->z ;
    packed = A->packed ;
    Xx = X->x ;
    Xz = X->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    kcol = X->ncol ;
    dy = Y->d ;
    dx = X->d ;
    w = W ;
    Wz = W + 4*nx ;

	/*
	 * Y = beta * Y
	 */
    if (ENTRY_IS_ZERO (beta, betaz, 0))
    {
	for (k = 0 ; k < kcol ; k++)
	{
	    for (i = 0 ; i < ((Int) ny) ; i++)
	    {
		/* y [i] = 0. ; */
		CLEAR (Yx, Yz, i) ;
	    }
	    /* y += dy ; */
	    ADVANCE (Yx,Yz,dy) ;
	}
    }
    else if (!ENTRY_IS_ONE (beta, betaz, 0))
    {
	for (k = 0 ; k < kcol ; k++)
	{
	    for (i = 0 ; i < ((Int) ny) ; i++)
	    {
		/* y [i] *= beta [0] ; */
		MULT (Yx,Yz,i, Yx,Yz,i, beta,betaz, 0) ;
	    }
	    /* y += dy ; */
	    ADVANCE (Yx,Yz,dy) ;
	}
    }

    if (ENTRY_IS_ZERO (alpha, alphaz, 0))
    {
	return ;
    }

	/*
	 * Y += alpha * op(A) * X ，其中op(A) = A或者A'
	 */
    Yx = Y->x ;
    Yz = Y->z ;

    k = 0 ;

    if (A->stype == 0)
    {

	if (transpose)
	{

		/*
		 * 不对称情况下 Y += alpha * A' * x
		 */
	    if (kcol % 4 == 1)
	    {

		for (j = 0 ; j < ncol ; j++)
		{
		    /* yj = 0. ; */
		    CLEAR (yx, yz, 0) ;
		    p = Ap [j] ;
		    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		    for ( ; p < pend ; p++)
		    {
			/* yj += conj(Ax [p]) * x [Ai [p]] ; */
			i = Ai [p] ;
			ASSIGN_CONJ (ax,az,0, Ax,Az,p) ;
			MULTADD (yx,yz,0, ax,az,0, Xx,Xz,i) ;
		    }
		    /* y [j] += alpha [0] * yj ; */
		    MULTADD (Yx,Yz,j, alpha,alphaz,0, yx,yz,0) ;
		}
		/* y += dy ; */
		/* x += dx ; */
		ADVANCE (Yx,Yz,dy) ;
		ADVANCE (Xx,Xz,dx) ;
		k++ ;

	    }
	    else if (kcol % 4 == 2)
	    {

		for (j = 0 ; j < ncol ; j++)
		{
		    /* yj0 = 0. ; */
		    /* yj1 = 0. ; */
		    CLEAR (yx,yz,0) ;
		    CLEAR (yx,yz,1) ;

		    p = Ap [j] ;
		    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			/* aij = conj (Ax [p]) ; */
			ASSIGN_CONJ (ax,az,0, Ax,Az,p) ;

			/* yj0 += aij * x [i   ] ; */
			/* yj1 += aij * x [i+dx] ; */
			MULTADD (yx,yz,0, ax,az,0, Xx,Xz,i) ;
			MULTADD (yx,yz,1, ax,az,0, Xx,Xz,i+dx) ;
		    }
		    /* y [j   ] += alpha [0] * yj0 ; */
		    /* y [j+dy] += alpha [0] * yj1 ; */
		    MULTADD (Yx,Yz,j,      alpha,alphaz,0, yx,yz,0) ;
		    MULTADD (Yx,Yz,j+dy,   alpha,alphaz,0, yx,yz,1) ;
		}
		/* y += 2*dy ; */
		/* x += 2*dx ; */
		ADVANCE (Yx,Yz,2*dy) ;
		ADVANCE (Xx,Xz,2*dx) ;
		k += 2 ;

	    }
	    else if (kcol % 4 == 3)
	    {

		for (j = 0 ; j < ncol ; j++)
		{
		    /* yj0 = 0. ; */
		    /* yj1 = 0. ; */
		    /* yj2 = 0. ; */
		    CLEAR (yx,yz,0) ;
		    CLEAR (yx,yz,1) ;
		    CLEAR (yx,yz,2) ;

		    p = Ap [j] ;
		    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			/* aij = conj (Ax [p]) ; */
			ASSIGN_CONJ (ax,az,0, Ax,Az,p) ;

			/* yj0 += aij * x [i     ] ; */
			/* yj1 += aij * x [i+  dx] ; */
			/* yj2 += aij * x [i+2*dx] ; */
			MULTADD (yx,yz,0, ax,az,0, Xx,Xz,i) ;
			MULTADD (yx,yz,1, ax,az,0, Xx,Xz,i+dx) ;
			MULTADD (yx,yz,2, ax,az,0, Xx,Xz,i+2*dx) ;
		    }
		    /* y [j     ] += alpha [0] * yj0 ; */
		    /* y [j+  dy] += alpha [0] * yj1 ; */
		    /* y [j+2*dy] += alpha [0] * yj2 ; */
		    MULTADD (Yx,Yz,j,      alpha,alphaz,0, yx,yz,0) ;
		    MULTADD (Yx,Yz,j+dy,   alpha,alphaz,0, yx,yz,1) ;
		    MULTADD (Yx,Yz,j+2*dy, alpha,alphaz,0, yx,yz,2) ;
		}
		/* y += 3*dy ; */
		/* x += 3*dx ; */
		ADVANCE (Yx,Yz,3*dy) ;
		ADVANCE (Xx,Xz,3*dx) ;
		k += 3 ;
	    }

	    for ( ; k < kcol ; k += 4)
	    {
		for (j = 0 ; j < ncol ; j++)
		{
		    /* yj0 = 0. ; */
		    /* yj1 = 0. ; */
		    /* yj2 = 0. ; */
		    /* yj3 = 0. ; */
		    CLEAR (yx,yz,0) ;
		    CLEAR (yx,yz,1) ;
		    CLEAR (yx,yz,2) ;
		    CLEAR (yx,yz,3) ;

		    p = Ap [j] ;
		    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			/* aij = conj(Ax [p]) ; */
			ASSIGN_CONJ (ax,az,0, Ax,Az,p) ;

			/* yj0 += aij * x [i     ] ; */
			/* yj1 += aij * x [i+  dx] ; */
			/* yj2 += aij * x [i+2*dx] ; */
			/* yj3 += aij * x [i+3*dx] ; */
			MULTADD (yx,yz,0, ax,az,0, Xx,Xz,i) ;
			MULTADD (yx,yz,1, ax,az,0, Xx,Xz,i+dx) ;
			MULTADD (yx,yz,2, ax,az,0, Xx,Xz,i+2*dx) ;
			MULTADD (yx,yz,3, ax,az,0, Xx,Xz,i+3*dx) ;

		    }
		    /* y [j     ] += alpha [0] * yj0 ; */
		    /* y [j+  dy] += alpha [0] * yj1 ; */
		    /* y [j+2*dy] += alpha [0] * yj2 ; */
		    /* y [j+3*dy] += alpha [0] * yj3 ; */
		    MULTADD (Yx,Yz,j,      alpha,alphaz,0, yx,yz,0) ;
		    MULTADD (Yx,Yz,j+dy,   alpha,alphaz,0, yx,yz,1) ;
		    MULTADD (Yx,Yz,j+2*dy, alpha,alphaz,0, yx,yz,2) ;
		    MULTADD (Yx,Yz,j+3*dy, alpha,alphaz,0, yx,yz,3) ;
		}
		/* y += 4*dy ; */
		/* x += 4*dx ; */
		ADVANCE (Yx,Yz,4*dy) ;
		ADVANCE (Xx,Xz,4*dx) ;
	    }

	}
	else
	{

		/*
		 * 不对称情况下 Y += alpha * A * x
		 */
	    if (kcol % 4 == 1)
	    {

		for (j = 0 ; j < ncol ; j++)
		{
		    /*  xj = alpha [0] * x [j] ; */
		    MULT (xx,xz,0, alpha,alphaz,0, Xx,Xz,j) ;

		    p = Ap [j] ;
		    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		    for ( ; p < pend ; p++)
		    {
			/* y [Ai [p]] += Ax [p] * xj ; */
			i = Ai [p] ;
			MULTADD (Yx,Yz,i, Ax,Az,p, xx,xz,0) ;
		    }
		}
		/* y += dy ; */
		/* x += dx ; */
		ADVANCE (Yx,Yz,dy) ;
		ADVANCE (Xx,Xz,dx) ;
		k++ ;

	    }
	    else if (kcol % 4 == 2)
	    {

		for (j = 0 ; j < ncol ; j++)
		{
		    /* xj0 = alpha [0] * x [j   ] ; */
		    /* xj1 = alpha [0] * x [j+dx] ; */
		    MULT (xx,xz,0, alpha,alphaz,0, Xx,Xz,j) ;
		    MULT (xx,xz,1, alpha,alphaz,0, Xx,Xz,j+dx) ;

		    p = Ap [j] ;
		    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i   ] += aij * xj0 ; */
			/* y [i+dy] += aij * xj1 ; */
			MULTADD (Yx,Yz,i,    ax,az,0, xx,xz,0) ;
			MULTADD (Yx,Yz,i+dy, ax,az,0, xx,xz,1) ;
		    }
		}
		/* y += 2*dy ; */
		/* x += 2*dx ; */
		ADVANCE (Yx,Yz,2*dy) ;
		ADVANCE (Xx,Xz,2*dx) ;
		k += 2 ;

	    }
	    else if (kcol % 4 == 3)
	    {

		for (j = 0 ; j < ncol ; j++)
		{
		    /* xj0 = alpha [0] * x [j     ] ; */
		    /* xj1 = alpha [0] * x [j+  dx] ; */
		    /* xj2 = alpha [0] * x [j+2*dx] ; */
		    MULT (xx,xz,0, alpha,alphaz,0, Xx,Xz,j) ;
		    MULT (xx,xz,1, alpha,alphaz,0, Xx,Xz,j+dx) ;
		    MULT (xx,xz,2, alpha,alphaz,0, Xx,Xz,j+2*dx) ;

		    p = Ap [j] ;
		    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i     ] += aij * xj0 ; */
			/* y [i+  dy] += aij * xj1 ; */
			/* y [i+2*dy] += aij * xj2 ; */
			MULTADD (Yx,Yz,i,      ax,az,0, xx,xz,0) ;
			MULTADD (Yx,Yz,i+dy,   ax,az,0, xx,xz,1) ;
			MULTADD (Yx,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
		    }
		}
		/* y += 3*dy ; */
		/* x += 3*dx ; */
		ADVANCE (Yx,Yz,3*dy) ;
		ADVANCE (Xx,Xz,3*dx) ;
		k += 3 ;
	    }

	    for ( ; k < kcol ; k += 4)
	    {
		for (j = 0 ; j < ncol ; j++)
		{
		    /* xj0 = alpha [0] * x [j     ] ; */
		    /* xj1 = alpha [0] * x [j+  dx] ; */
		    /* xj2 = alpha [0] * x [j+2*dx] ; */
		    /* xj3 = alpha [0] * x [j+3*dx] ; */
		    MULT (xx,xz,0, alpha,alphaz,0, Xx,Xz,j) ;
		    MULT (xx,xz,1, alpha,alphaz,0, Xx,Xz,j+dx) ;
		    MULT (xx,xz,2, alpha,alphaz,0, Xx,Xz,j+2*dx) ;
		    MULT (xx,xz,3, alpha,alphaz,0, Xx,Xz,j+3*dx) ;

		    p = Ap [j] ;
		    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i     ] += aij * xj0 ; */
			/* y [i+  dy] += aij * xj1 ; */
			/* y [i+2*dy] += aij * xj2 ; */
			/* y [i+3*dy] += aij * xj3 ; */
			MULTADD (Yx,Yz,i,      ax,az,0, xx,xz,0) ;
			MULTADD (Yx,Yz,i+dy,   ax,az,0, xx,xz,1) ;
			MULTADD (Yx,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
			MULTADD (Yx,Yz,i+3*dy, ax,az,0, xx,xz,3) ;
		    }
		}
		/* y += 4*dy ; */
		/* x += 4*dx ; */
		ADVANCE (Yx,Yz,4*dy) ;
		ADVANCE (Xx,Xz,4*dx) ;
	    }
	}

    }
    else
    {

	/*
	 * 对称情况下（上三角或者下三角） Y += alpha * (A or A') * x。
	 * 
	 * 只使用了A的上三角或者下三角部分以及对角线元素。由于x和y都是在
	 * 内层循环中写入的，因此如果直接使用x，可能会造成缓存冲突。所以，
	 * 如果x有四列或者更多列，则每次获取x的四列作为副本。
	 */
	if (kcol % 4 == 1)
	{

	    for (j = 0 ; j < ncol ; j++)
	    {
		/* yj = 0. ; */
		CLEAR (yx,yz,0) ;

		/* xj = alpha [0] * x [j] ; */
		MULT (xx,xz,0, alpha,alphaz,0, Xx,Xz,j) ;

		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i == j)
		    {
			/* y [i] += Ax [p] * xj ; */
			MULTADD (Yx,Yz,i, Ax,Az,p, xx,xz,0) ;
		    }
		    else if ((A->stype > 0 && i < j) || (A->stype < 0 && i > j))
		    {
			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i] += aij * xj ; */
			/* yj    += aij * x [i] ; */
			MULTADD     (Yx,Yz,i, ax,az,0, xx,xz,0) ;
			MULTADDCONJ (yx,yz,0, ax,az,0, Xx,Xz,i) ;


		    }
		}
		/* y [j] += alpha [0] * yj ; */
		MULTADD (Yx,Yz,j, alpha,alphaz,0, yx,yz,0) ;

	    }
	    /* y += dy ; */
	    /* x += dx ; */
	    ADVANCE (Yx,Yz,dy) ;
	    ADVANCE (Xx,Xz,dx) ;
	    k++ ;

	}
	else if (kcol % 4 == 2)
	{

	    for (j = 0 ; j < ncol ; j++)
	    {
		/* yj0 = 0. ; */
		/* yj1 = 0. ; */
		CLEAR (yx,yz,0) ;
		CLEAR (yx,yz,1) ;

		/* xj0 = alpha [0] * x [j   ] ; */
		/* xj1 = alpha [0] * x [j+dx] ; */
		MULT (xx,xz,0, alpha,alphaz,0, Xx,Xz,j) ;
		MULT (xx,xz,1, alpha,alphaz,0, Xx,Xz,j+dx) ;

		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i == j)
		    {
			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i   ] += aij * xj0 ; */
			/* y [i+dy] += aij * xj1 ; */
			MULTADD (Yx,Yz,i,    ax,az,0, xx,xz,0) ;
			MULTADD (Yx,Yz,i+dy, ax,az,0, xx,xz,1) ;

		    }
		    else if ((A->stype > 0 && i < j) || (A->stype < 0 && i > j))
		    {
			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i   ] += aij * xj0 ; */
			/* y [i+dy] += aij * xj1 ; */
			/* yj0 += aij * x [i   ] ; */
			/* yj1 += aij * x [i+dx] ; */
			MULTADD     (Yx,Yz,i,    ax,az,0, xx,xz,0) ;
			MULTADD     (Yx,Yz,i+dy, ax,az,0, xx,xz,1) ;
			MULTADDCONJ (yx,yz,0,    ax,az,0, Xx,Xz,i) ;
			MULTADDCONJ (yx,yz,1,    ax,az,0, Xx,Xz,i+dx) ;

		    }
		}
		/* y [j   ] += alpha [0] * yj0 ; */
		/* y [j+dy] += alpha [0] * yj1 ; */
		MULTADD (Yx,Yz,j,    alpha,alphaz,0, yx,yz,0) ;
		MULTADD (Yx,Yz,j+dy, alpha,alphaz,0, yx,yz,1) ;

	    }
	    /* y += 2*dy ; */
	    /* x += 2*dx ; */
	    ADVANCE (Yx,Yz,2*dy) ;
	    ADVANCE (Xx,Xz,2*dx) ;
	    k += 2 ;

	}
	else if (kcol % 4 == 3)
	{

	    for (j = 0 ; j < ncol ; j++)
	    {
		/* yj0 = 0. ; */
		/* yj1 = 0. ; */
		/* yj2 = 0. ; */
		CLEAR (yx,yz,0) ;
		CLEAR (yx,yz,1) ;
		CLEAR (yx,yz,2) ;

		/* xj0 = alpha [0] * x [j     ] ; */
		/* xj1 = alpha [0] * x [j+  dx] ; */
		/* xj2 = alpha [0] * x [j+2*dx] ; */
		MULT (xx,xz,0, alpha,alphaz,0, Xx,Xz,j) ;
		MULT (xx,xz,1, alpha,alphaz,0, Xx,Xz,j+dx) ;
		MULT (xx,xz,2, alpha,alphaz,0, Xx,Xz,j+2*dx) ;

		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i == j)
		    {

			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i     ] += aij * xj0 ; */
			/* y [i+  dy] += aij * xj1 ; */
			/* y [i+2*dy] += aij * xj2 ; */
			MULTADD (Yx,Yz,i,      ax,az,0, xx,xz,0) ;
			MULTADD (Yx,Yz,i+dy,   ax,az,0, xx,xz,1) ;
			MULTADD (Yx,Yz,i+2*dy, ax,az,0, xx,xz,2) ;

		    }
		    else if ((A->stype > 0 && i < j) || (A->stype < 0 && i > j))
		    {

			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i     ] += aij * xj0 ; */
			/* y [i+  dy] += aij * xj1 ; */
			/* y [i+2*dy] += aij * xj2 ; */
			/* yj0 += aij * x [i     ] ; */
			/* yj1 += aij * x [i+  dx] ; */
			/* yj2 += aij * x [i+2*dx] ; */
			MULTADD     (Yx,Yz,i,      ax,az,0, xx,xz,0) ;
			MULTADD     (Yx,Yz,i+dy,   ax,az,0, xx,xz,1) ;
			MULTADD     (Yx,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
			MULTADDCONJ (yx,yz,0,      ax,az,0, Xx,Xz,i) ;
			MULTADDCONJ (yx,yz,1,      ax,az,0, Xx,Xz,i+dx) ;
			MULTADDCONJ (yx,yz,2,      ax,az,0, Xx,Xz,i+2*dx) ;

		    }
		}
		/* y [j     ] += alpha [0] * yj0 ; */
		/* y [j+  dy] += alpha [0] * yj1 ; */
		/* y [j+2*dy] += alpha [0] * yj2 ; */
		MULTADD (Yx,Yz,j,      alpha,alphaz,0, yx,yz,0) ;
		MULTADD (Yx,Yz,j+dy,   alpha,alphaz,0, yx,yz,1) ;
		MULTADD (Yx,Yz,j+2*dy, alpha,alphaz,0, yx,yz,2) ;

	    }
	    /* y += 3*dy ; */
	    /* x += 3*dx ; */
	    ADVANCE (Yx,Yz,3*dy) ;
	    ADVANCE (Xx,Xz,3*dx) ;

	    k += 3 ;
	}

	/* 将X的四列以行形式复制到W中 */
	for ( ; k < kcol ; k += 4)
	{

	    for (j = 0 ; j < ncol ; j++)
	    {
		/* w [4*j  ] = x [j     ] ; */
		/* w [4*j+1] = x [j+  dx] ; */
		/* w [4*j+2] = x [j+2*dx] ; */
		/* w [4*j+3] = x [j+3*dx] ; */
		ASSIGN (w,Wz,4*j  , Xx,Xz,j     ) ;
		ASSIGN (w,Wz,4*j+1, Xx,Xz,j+dx  ) ;
		ASSIGN (w,Wz,4*j+2, Xx,Xz,j+2*dx) ;
		ASSIGN (w,Wz,4*j+3, Xx,Xz,j+3*dx) ;
	    }

	    for (j = 0 ; j < ncol ; j++)
	    {
		/* yj0 = 0. ; */
		/* yj1 = 0. ; */
		/* yj2 = 0. ; */
		/* yj3 = 0. ; */
		CLEAR (yx,yz,0) ;
		CLEAR (yx,yz,1) ;
		CLEAR (yx,yz,2) ;
		CLEAR (yx,yz,3) ;

		/* xj0 = alpha [0] * w [4*j  ] ; */
		/* xj1 = alpha [0] * w [4*j+1] ; */
		/* xj2 = alpha [0] * w [4*j+2] ; */
		/* xj3 = alpha [0] * w [4*j+3] ; */
		MULT (xx,xz,0, alpha,alphaz,0, w,Wz,4*j) ;
		MULT (xx,xz,1, alpha,alphaz,0, w,Wz,4*j+1) ;
		MULT (xx,xz,2, alpha,alphaz,0, w,Wz,4*j+2) ;
		MULT (xx,xz,3, alpha,alphaz,0, w,Wz,4*j+3) ;

		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i == j)
		    {
			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i     ] += aij * xj0 ; */
			/* y [i+  dy] += aij * xj1 ; */
			/* y [i+2*dy] += aij * xj2 ; */
			/* y [i+3*dy] += aij * xj3 ; */
			MULTADD (Yx,Yz,i     , ax,az,0, xx,xz,0) ;
			MULTADD (Yx,Yz,i+dy  , ax,az,0, xx,xz,1) ;
			MULTADD (Yx,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
			MULTADD (Yx,Yz,i+3*dy, ax,az,0, xx,xz,3) ;

		    }
		    else if ((A->stype > 0 && i < j) || (A->stype < 0 && i > j))
		    {
			/* aij = Ax [p] ; */
			ASSIGN (ax,az,0, Ax,Az,p) ;

			/* y [i     ] += aij * xj0 ; */
			/* y [i+  dy] += aij * xj1 ; */
			/* y [i+2*dy] += aij * xj2 ; */
			/* y [i+3*dy] += aij * xj3 ; */
			/* yj0 += aij * w [4*i  ] ; */
			/* yj1 += aij * w [4*i+1] ; */
			/* yj2 += aij * w [4*i+2] ; */
			/* yj3 += aij * w [4*i+3] ; */
			MULTADD     (Yx,Yz,i,      ax,az,0, xx,xz,0) ;
			MULTADD     (Yx,Yz,i+dy,   ax,az,0, xx,xz,1) ;
			MULTADD     (Yx,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
			MULTADD     (Yx,Yz,i+3*dy, ax,az,0, xx,xz,3) ;
			MULTADDCONJ (yx,yz,0,     ax,az,0, w,Wz,4*i) ;
			MULTADDCONJ (yx,yz,1,     ax,az,0, w,Wz,4*i+1) ;
			MULTADDCONJ (yx,yz,2,     ax,az,0, w,Wz,4*i+2) ;
			MULTADDCONJ (yx,yz,3,     ax,az,0, w,Wz,4*i+3) ;

		    }
		}
		/* y [j     ] += alpha [0] * yj0 ; */
		/* y [j+  dy] += alpha [0] * yj1 ; */
		/* y [j+2*dy] += alpha [0] * yj2 ; */
		/* y [j+3*dy] += alpha [0] * yj3 ; */
		MULTADD (Yx,Yz,j     , alpha,alphaz,0, yx,yz,0) ;
		MULTADD (Yx,Yz,j+dy  , alpha,alphaz,0, yx,yz,1) ;
		MULTADD (Yx,Yz,j+2*dy, alpha,alphaz,0, yx,yz,2) ;
		MULTADD (Yx,Yz,j+3*dy, alpha,alphaz,0, yx,yz,3) ;

	    }
	    /* y += 4*dy ; */
	    /* x += 4*dx ; */
	    ADVANCE (Yx,Yz,4*dy) ;
	    ADVANCE (Xx,Xz,4*dx) ;

	}
    }
}

#undef PATTERN
#undef REAL