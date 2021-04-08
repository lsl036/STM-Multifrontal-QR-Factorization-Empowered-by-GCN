/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    SparseLU_solve.c
 * BRIEF:   SparseLU求解器
 *****************************************************************************/





#include "SparseLU_internal.h"
#include "SparseLU_function.h"

GLOBAL Int
#ifdef WSOLVE
SparseLU_wsolve
#else
SparseLU_solve
#endif
(
    Int sys,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    double Xx [ ],
#ifdef COMPLEX
    double Xz [ ],
#endif
    const double Bx [ ],
#ifdef COMPLEX
    const double Bz [ ],
#endif
    void *NumericHandle,
    const double Control [SparseLU_CONTROL],
    double User_Info [SparseLU_INFO]
#ifdef WSOLVE
    , Int Pattern [ ],
    double W [ ]
#endif
)
{
    
    
    

    double Info2 [SparseLU_INFO], stats [2] ;
    double *Info ;
    NumericType *Numeric ;
    Int n, i, irstep, status ;
#ifndef WSOLVE
    Int *Pattern, wsize ;
    double *W ;
#endif

    
    
    

    irstep = GET_CONTROL (SparseLU_IRSTEP, SparseLU_DEFAULT_IRSTEP) ; //internel

    if (User_Info != (double *) NULL)
    {
	
	Info = User_Info ;
	
	for (i = SparseLU_IR_TAKEN ; i <= SparseLU_SOLVE_TIME ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }
    else
    {
	
	Info = Info2 ;
	for (i = 0 ; i < SparseLU_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }

    Info [SparseLU_STATUS] = SparseLU_OK ;
    Info [SparseLU_SOLVE_FLOPS] = 0 ;

    Numeric = (NumericType *) NumericHandle ;
    if (!LU_valid_numeric (Numeric))
    {
	Info [SparseLU_STATUS] = SparseLU_ERROR_invalid_Numeric_object ;
	return (SparseLU_ERROR_invalid_Numeric_object) ;
    }

    Info [SparseLU_NROW] = Numeric->n_row ;
    Info [SparseLU_NCOL] = Numeric->n_col ;

    if (Numeric->n_row != Numeric->n_col)
    {
	
	Info [SparseLU_STATUS] = SparseLU_ERROR_invalid_system ;
	return (SparseLU_ERROR_invalid_system) ;
    }
    n = Numeric->n_row ;
    if (Numeric->nnzpiv < n
	|| SCALAR_IS_ZERO (Numeric->rcond) || SCALAR_IS_NAN (Numeric->rcond)) //lu_version.h
    {
	
	
	irstep = 0 ;
    }

    if (!Xx || !Bx)
    {
	Info [SparseLU_STATUS] = SparseLU_ERROR_argument_missing ;
	return (SparseLU_ERROR_argument_missing) ;
    }

    if (sys >= SparseLU_Pt_L)
    {
	
	irstep = 0 ;
    }

    
    
    

#ifdef WSOLVE

    if (!W || !Pattern)
    {
	Info [SparseLU_STATUS] = SparseLU_ERROR_argument_missing ;
	return (SparseLU_ERROR_argument_missing) ;
    }

#else

#ifdef COMPLEX
    if (irstep > 0)
    {
	wsize = 10*n ;		
    }
    else
    {
	wsize = 4*n ;		
    }
#else
    if (irstep > 0)
    {
	wsize = 5*n ;		
    }
    else
    {
	wsize = n ;		
    }
#endif

    Pattern = (Int *) LU_malloc (n, sizeof (Int)) ;
    W = (double *) LU_malloc (wsize, sizeof (double)) ;
    if (!W || !Pattern)
    {
	
	Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
	(void) LU_free ((void *) W) ;
	(void) LU_free ((void *) Pattern) ;
	return (SparseLU_ERROR_out_of_memory) ;
    }

#endif	

    
    
    

    status = LU_solve (sys, Ap, Ai, Ax, Xx, Bx,
#ifdef COMPLEX
	Az, Xz, Bz,
#endif
	Numeric, irstep, Info, Pattern, W) ;

    
    
    

#ifndef WSOLVE
    (void) LU_free ((void *) W) ;
    (void) LU_free ((void *) Pattern) ;
    
#endif

    
    
    

    Info [SparseLU_STATUS] = status ;
    if (status >= 0)
    {
	Info [SparseLU_SOLVE_WALLTIME] = stats [0] ;
	Info [SparseLU_SOLVE_TIME] = stats [1] ;
    } 

    return (status) ;
}
