/**
 * @file t_SparseChol_super_solve.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-22
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_template.h"

/**
 * @brief 
 * 
 */
static void TEMPLATE (SparseChol_super_lsolve)
(
    /* ---- input ---- */
    sparse_factor *L,	/* 用于正向求解的因子 */
    /* ---- output ---- */
    dense_array *X,	/* 输入为b，输出为Lx=b的解 */
    /* ---- workspace ---- */
    dense_array *E,	/* 大小为nrhs*(L->maxesize)工作空间 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Lx, *Xx, *Ex ;
    double minus_one [2], one [2] ;
    Int *Lpi, *Lpx, *Ls, *Super ;
    Int nsuper, k1, k2, psi, psend, psx, nsrow, nscol, ii, s,
	nsrow2, n, ps2, j, i, d, nrhs ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    nrhs = X->ncol ;
    Ex = E->x ;
    Xx = X->x ;
    n = L->n ;
    d = X->d ;

    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;
    Lx = L->x ;
    minus_one [0] = -1.0 ;
    minus_one [1] = 0 ;
    one [0] = 1.0 ;
    one [1] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* 求解Lx=b */
    /* ---------------------------------------------------------------------- */

    if (nrhs == 1)
    {

		for (s = 0 ; s < nsuper ; s++)
		{
			k1 = Super [s] ;
			k2 = Super [s+1] ;
			psi = Lpi [s] ;
			psend = Lpi [s+1] ;
			psx = Lpx [s] ;
			nsrow = psend - psi ;
			nscol = k2 - k1 ;
			nsrow2 = nsrow - nscol ;
			ps2 = psi + nscol ;
			

			/* L1大小为nscol*nscol, 是具有非单位对角线的下三角矩阵.	     
			 * L2大小为nsrow2*nscol.  L1和L2前导维度为nsrow. 
			 * x1大小为nscol*nsrow, 前导维度为n.
			 * E大小为nsrow2*1, 前导维度为nsrow2.
			 */

			/* 把X集合成E */
			for (ii = 0 ; ii < nsrow2 ; ii++)
			{
			/* Ex [ii] = Xx [Ls [ps2 + ii]] ; */
			ASSIGN (Ex,-,ii, Xx,-,Ls [ps2 + ii]) ;
			}

			/* solve L1*x1 (that is, x1 = L1\x1) */
			BLAS_dtrsv ("L", "N", "N",
			nscol,			    
			Lx + ENTRY_SIZE*psx, nsrow,
			Xx + ENTRY_SIZE*k1, 1) ;    

			/* E = E - L2*x1 */
			BLAS_dgemv ("N",
			nsrow2, nscol,		   
			minus_one,		    
			Lx + ENTRY_SIZE*(psx + nscol),  
			nsrow,
			Xx + ENTRY_SIZE*k1, 1,	    
			one,			  
			Ex, 1) ;		   


			/* 把E传回X中 */
			for (ii = 0 ; ii < nsrow2 ; ii++)
			{
			/* Xx [Ls [ps2 + ii]] = Ex [ii] ; */
			ASSIGN (Xx,-,Ls [ps2 + ii], Ex,-,ii) ;
			}
		}
    }
    else
    {

		for (s = 0 ; s < nsuper ; s++)
		{
			k1 = Super [s] ;
			k2 = Super [s+1] ;
			psi = Lpi [s] ;
			psend = Lpi [s+1] ;
			psx = Lpx [s] ;
			nsrow = psend - psi ;
			nscol = k2 - k1 ;
			nsrow2 = nsrow - nscol ;
			ps2 = psi + nscol ;

			/* E大小为nsrow2*nrhs，前导维为nsrow2。 */

			/* 把X集合成E */
			for (ii = 0 ; ii < nsrow2 ; ii++)
			{
			i = Ls [ps2 + ii] ;
			for (j = 0 ; j < nrhs ; j++)
			{
				/* Ex [ii + j*nsrow2] = Xx [i + j*d] ; */
				ASSIGN (Ex,-,ii+j*nsrow2, Xx,-,i+j*d) ;
			}
			}

			/* 求解L1*x1 */
			BLAS_dtrsm ("L", "L", "N", "N",
			nscol, nrhs,			
			one,				
			Lx + ENTRY_SIZE*psx, nsrow,	
			Xx + ENTRY_SIZE*k1, d) ;	

			/* E = E - L2*x1 */
			if (nsrow2 > 0)
			{
			BLAS_dgemm ("N", "N",
				nsrow2, nrhs, nscol,	    
				minus_one,			    
				Lx + ENTRY_SIZE*(psx + nscol),  
				nsrow,
				Xx + ENTRY_SIZE*k1, d,	    
				one,			    
				Ex, nsrow2) ;		   
			}

			/* 把E传回X中 */
			for (ii = 0 ; ii < nsrow2 ; ii++)
			{
			i = Ls [ps2 + ii] ;
			for (j = 0 ; j < nrhs ; j++)
			{
				/* Xx [i + j*d] = Ex [ii + j*nsrow2] ; */
				ASSIGN (Xx,-,i+j*d, Ex,-,ii+j*nsrow2) ;
			}
			}
	}
    }
}

/**
 * @brief 
 * 
 */
static void TEMPLATE (SparseChol_super_ltsolve)
(
    /* ---- input ---- */
    sparse_factor *L,	/* 用于正向求解的因子 */
    /* ---- output ---- */
    dense_array *X,	/* 输入为b，输出为Lx=b的解 */
    /* ---- workspace ---- */
    dense_array *E,	/* 大小为nrhs*(L->maxesize)工作空间 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Lx, *Xx, *Ex ;
    double minus_one [2], one [2] ;
    Int *Lpi, *Lpx, *Ls, *Super ;
    Int nsuper, k1, k2, psi, psend, psx, nsrow, nscol, ii, s,
	nsrow2, n, ps2, j, i, d, nrhs ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    nrhs = X->ncol ;
    Ex = E->x ;
    Xx = X->x ;
    n = L->n ;
    d = X->d ;

    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;
    Lx = L->x ;
    minus_one [0] = -1.0 ;
    minus_one [1] = 0 ;
    one [0] = 1.0 ;
    one [1] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* 求解L'x=b */
    /* ---------------------------------------------------------------------- */

    if (nrhs == 1)
    {

		for (s = nsuper-1 ; s >= 0 ; s--)
		{
			k1 = Super [s] ;
			k2 = Super [s+1] ;
			psi = Lpi [s] ;
			psend = Lpi [s+1] ;
			psx = Lpx [s] ;
			nsrow = psend - psi ;
			nscol = k2 - k1 ;
			nsrow2 = nsrow - nscol ;
			ps2 = psi + nscol ;
			

			/* L1大小为nscol*nscol, 是具有非单位对角线的下三角矩阵.
			* L2大小为nsrow2*nscol.  L1和L2前导维度为nsrow. 
			* x1大小为nscol*nsrow, 前导维度为n.
			* E大小为nsrow2*1, 前导维度为nsrow2.
			*/

			/*把X集合成E */
			for (ii = 0 ; ii < nsrow2 ; ii++)
			{
			/* Ex [ii] = Xx [Ls [ps2 + ii]] ; */
			ASSIGN (Ex,-,ii, Xx,-,Ls [ps2 + ii]) ;
			}

			/* x1 = x1 - L2'*E */
			BLAS_dgemv ("C",
			nsrow2, nscol,		    
			minus_one,		    
			Lx + ENTRY_SIZE*(psx + nscol),   
			nsrow,
			Ex, 1,			    
			one,			    
			Xx + ENTRY_SIZE*k1, 1) ;    

			/* solve L1'*x1 */
			BLAS_dtrsv ("L", "C", "N",
			nscol,			    
			Lx + ENTRY_SIZE*psx, nsrow,	   
			Xx + ENTRY_SIZE*k1, 1) ;	    

		}
    }
    else
    {

		for (s = nsuper-1 ; s >= 0 ; s--)
		{
			k1 = Super [s] ;
			k2 = Super [s+1] ;
			psi = Lpi [s] ;
			psend = Lpi [s+1] ;
			psx = Lpx [s] ;
			nsrow = psend - psi ;
			nscol = k2 - k1 ;
			nsrow2 = nsrow - nscol ;
			ps2 = psi + nscol ;

			/*把X集合成E */
			for (ii = 0 ; ii < nsrow2 ; ii++)
			{
			i = Ls [ps2 + ii] ;
			for (j = 0 ; j < nrhs ; j++)
			{
				/* Ex [ii + j*nsrow2] = Xx [i + j*d] ; */
				ASSIGN (Ex,-,ii+j*nsrow2, Xx,-,i+j*d) ;
			}
			}

			/* x1 = x1 - L2'*E */
			if (nsrow2 > 0)
			{
			BLAS_dgemm ("C", "N",
				nscol, nrhs, nsrow2,	
				minus_one,		
				Lx + ENTRY_SIZE*(psx + nscol),  
				nsrow,
				Ex, nsrow2,		
				one,			
				Xx + ENTRY_SIZE*k1, d) ;	
			}

			/* solve L1'*x1 */
			BLAS_dtrsm ("L", "L", "C", "N",
			nscol,	nrhs,			
			one,				
			Lx + ENTRY_SIZE*psx, nsrow,	
			Xx + ENTRY_SIZE*k1, d) ;	

		}
    }
}

#undef PATTERN
#undef REAL