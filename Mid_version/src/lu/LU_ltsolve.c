/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_ltsolve.c
 * BRIEF:   LU,L的转置求解方程组
 *****************************************************************************/


#include "SparseLU_internal.h"
#include "SparseLU_function.h"



//   LU_ltsolve or LU_lhsolve
GLOBAL double
#ifdef CONJUGATE_SOLVE
LU_lhsolve			
#else
LU_ltsolve			
#endif
(
    NumericType *Numeric,
    Entry X [ ],		
    Int Pattern [ ]		
)
{
    Entry xk ;
    Entry *xp, *Lval ;
    Int k, deg, *ip, j, row, *Lpos, *Lilen, kstart, kend, *Lip, llen,
	lp, pos, npiv, n1, *Li ;

    

    if (Numeric->n_row != Numeric->n_col) return (0.) ;
    npiv = Numeric->npiv ;
    Lpos = Numeric->Lpos ;
    Lilen = Numeric->Lilen ;
    Lip = Numeric->Lip ;
    kstart = npiv ;
    n1 = Numeric->n1 ;

    
    
    

    for (kend = npiv-1 ; kend >= n1 ; kend = kstart-1)
    {

	
	
	

	
	kstart = kend ;
	while (kstart >= 0 && Lip [kstart] > 0)
	{
	    kstart-- ;
	}

	

	
	
	

	deg = 0 ;
	for (k = kstart ; k <= kend ; k++)
	{

	    
	    
	    

	    
	    pos = Lpos [k] ;
	    if (pos != EMPTY)
	    {
		Pattern [pos] = Pattern [--deg] ;
	    }

	    
	    lp = Lip [k] ;
	    if (k == kstart)
	    {
		lp = -lp ;
	    }
	    ip = (Int *) (Numeric->Memory + lp) ;
	    llen = Lilen [k] ;
	    for (j = 0 ; j < llen ; j++)
	    {
		row = *ip++ ;
		Pattern [deg++] = row ;
	    }

	}
	

	
	
	

	for (k = kend ; k >= kstart ; k--)
	{

	    
	    
	    

	    lp = Lip [k] ;
	    if (k == kstart)
	    {
		lp = -lp ;
	    }
	    llen = Lilen [k] ;
	    xp = (Entry *) (Numeric->Memory + lp + UNITS (Int, llen)) ;
	    xk = X [k] ;
	    for (j = 0 ; j < deg ; j++)
	    {
#ifdef CONJUGATE_SOLVE
		
		MULT_SUB_CONJ (xk, X [Pattern [j]], *xp) ;
#else
		
		MULT_SUB (xk, X [Pattern [j]], *xp) ;
#endif

		xp++ ;
	    }
	    X [k] = xk ;

	    
	    
	    

	    
	    deg -= llen ;

	    
	    pos = Lpos [k] ;
	    if (pos != EMPTY)
	    {
		Pattern [deg++] = Pattern [pos] ;
		Pattern [pos] = k ;
	    }
	}
    }

    
    
    

    for (k = n1 - 1 ; k >= 0 ; k--)
    {
	deg = Lilen [k] ;
	if (deg > 0)
	{
	    xk = X [k] ;
	    lp = Lip [k] ;
	    Li = (Int *) (Numeric->Memory + lp) ;
	    lp += UNITS (Int, deg) ;
	    Lval = (Entry *) (Numeric->Memory + lp) ;
	    for (j = 0 ; j < deg ; j++)
	    {
#ifdef CONJUGATE_SOLVE
		
		MULT_SUB_CONJ (xk, X [Li [j]], Lval [j]) ;
#else
		
		MULT_SUB (xk, X [Li [j]], Lval [j]) ;
#endif
	    }
	    X [k] = xk ;
	}
    }

    return (MULTSUB_FLOPS * ((double) Numeric->lnz)) ;
}