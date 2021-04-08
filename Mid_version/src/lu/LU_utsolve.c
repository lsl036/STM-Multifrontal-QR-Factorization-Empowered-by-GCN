/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_utsolve.c
 * BRIEF:   LU U的转置求解器
 *****************************************************************************/

#include "SparseLU_internal.h"
#include "SparseLU_function.h"




//   LU_utsolve or LU_uhsolve
GLOBAL double
#ifdef CONJUGATE_SOLVE
LU_uhsolve			
#else
LU_utsolve			
#endif
(
    NumericType *Numeric,
    Entry X [ ],		
    Int Pattern [ ]		
)
{
    
    
    

    Entry xk ;
    Entry *xp, *D, *Uval ;
    Int k, deg, j, *ip, col, *Upos, *Uilen, kstart, kend, up,
	*Uip, n, uhead, ulen, pos, npiv, n1, *Ui ;

    
    
    

    if (Numeric->n_row != Numeric->n_col) return (0.) ;
    n = Numeric->n_row ;
    npiv = Numeric->npiv ;
    Upos = Numeric->Upos ;
    Uilen = Numeric->Uilen ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;
    kend = 0 ;
    n1 = Numeric->n1 ;

    
    
    

    for (k = 0 ; k < n1 ; k++)
    {

#ifndef NO_DIVIDE_BY_ZERO
	
#ifdef CONJUGATE_SOLVE
	
	DIV_CONJ (xk, X [k], D [k]) ;
#else
	
	DIV (xk, X [k], D [k]) ;
#endif
#else
	
	if (IS_NONZERO (D [k]))
	{
#ifdef CONJUGATE_SOLVE
	    
	    DIV_CONJ (xk, X [k], D [k]) ;
#else
	    
	    DIV (xk, X [k], D [k]) ;
#endif
	}
#endif

	X [k] = xk ;
	deg = Uilen [k] ;
	if (deg > 0 && IS_NONZERO (xk))
	{
	    up = Uip [k] ;
	    Ui = (Int *) (Numeric->Memory + up) ;
	    up += UNITS (Int, deg) ;
	    Uval = (Entry *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		
#ifdef CONJUGATE_SOLVE
		
		MULT_SUB_CONJ (X [Ui [j]], xk, Uval [j]) ;
#else
		
		MULT_SUB (X [Ui [j]], xk, Uval [j]) ;
#endif
	    }
	}
    }

    
    
    

    for (kstart = n1 ; kstart < npiv ; kstart = kend + 1)
    {

	
	
	

	
	kend = kstart ;
	while (kend < npiv && Uip [kend+1] > 0)
	{
	    kend++ ;
	}

	
	
	

	k = kend+1 ;

	
	
	

	if (k == npiv)
	{
	    deg = Numeric->ulen ;
	    if (deg > 0)
	    {
		
		for (j = 0 ; j < deg ; j++)
		{
		    Pattern [j] = Numeric->Upattern [j] ;
		}
	    }
	}
	else
	{
	    up = -Uip [k] ;
	    deg = Uilen [k] ;
	    ip = (Int *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		col = *ip++ ;
		
		Pattern [j] = col ;
	    }
	}

	
	uhead = n ;

	for (k = kend ; k > kstart ; k--)
	{
	    

	    
	    
	    
	    ulen = Uilen [k] ;
	    
	    for (j = 0 ; j < ulen ; j++)
	    {
		Pattern [--uhead] = Pattern [--deg] ;
	    }
	    

	    pos = Upos [k] ;
	    if (pos != EMPTY)
	    {
		
		
		Pattern [deg++] = Pattern [pos] ;
		Pattern [pos] = k ;
	    }
	}

	

	
	
	


	for (k = kstart ; k <= kend ; k++)
	{

	    
	    
	    

	    
	    pos = Upos [k] ;
	    if (pos != EMPTY)
	    {
		
		
		Pattern [pos] = Pattern [--deg] ;
	    }

	    up = Uip [k] ;
	    ulen = Uilen [k] ;
	    if (k > kstart)
	    {
		
		for (j = 0 ; j < ulen ; j++)
		{
		    Pattern [deg++] = Pattern [uhead++] ;
		}

	    }

	    
	    
	    

#ifndef NO_DIVIDE_BY_ZERO
	    
#ifdef CONJUGATE_SOLVE
	    
	    DIV_CONJ (xk, X [k], D [k]) ;
#else
	    
	    DIV (xk, X [k], D [k]) ;
#endif
#else
	    
	    if (IS_NONZERO (D [k]))
	    {
#ifdef CONJUGATE_SOLVE
		
		DIV_CONJ (xk, X [k], D [k]) ;
#else
		
		DIV (xk, X [k], D [k]) ;
#endif
	    }
#endif

	    X [k] = xk ;
	    if (IS_NONZERO (xk))
	    {
		if (k == kstart)
		{
		    up = -up ;
		    xp = (Entry *) (Numeric->Memory + up + UNITS (Int, ulen)) ;
		}
		else
		{
		    xp = (Entry *) (Numeric->Memory + up) ;
		}
		for (j = 0 ; j < deg ; j++)
		{
#ifdef CONJUGATE_SOLVE
		    
		    MULT_SUB_CONJ (X [Pattern [j]], xk, *xp) ;
#else
		    
		    MULT_SUB (X [Pattern [j]], xk, *xp) ;
#endif
		    xp++ ;
		}
	    }
	}
    }

#ifndef NO_DIVIDE_BY_ZERO
    for (k = npiv ; k < n ; k++)
    {
	
	
	
	DIV (xk, X [k], D [k]) ;
	X [k] = xk ;
    }
#endif

    return (DIV_FLOPS * ((double) n) + MULTSUB_FLOPS * ((double) Numeric->unz));
}
