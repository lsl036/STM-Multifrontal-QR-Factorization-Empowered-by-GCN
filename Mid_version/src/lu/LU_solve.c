/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_solve.c
 * BRIEF:   LU求解器
 *****************************************************************************/



#include "SparseLU_internal.h"
#include "SparseLU_function.h"

PRIVATE Int do_step
(
    double omega [3],
    Int step,
    const double B2 [ ],
    Entry X [ ],
    const Entry W [ ],
    const double Y [ ],
    const double Z2 [ ],
    Entry S [ ],
    Int n,
    double Info [SparseLU_INFO]
) ;





GLOBAL Int LU_solve
(
    Int sys,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    double Xx [ ],
    const double Bx [ ],
#ifdef COMPLEX
    const double Az [ ],
    double Xz [ ],
    const double Bz [ ],
#endif
    NumericType *Numeric,
    Int irstep,
    double Info [SparseLU_INFO],
    Int Pattern [ ],		
    double SolveWork [ ]	
				
)
{
    
    
    

    Entry axx, wi, xj, zi, xi, aij, bi ;
    double omega [3], d, z2i, yi, flops ;
    Entry *W, *Z, *S, *X ;
    double *Z2, *Y, *B2, *Rs ;
    Int *Rperm, *Cperm, i, n, p, step, j, nz, status, p2, do_scale ;
#ifdef COMPLEX
    Int AXsplit ;
    Int Bsplit ;
#endif
#ifndef NRECIPROCAL
    Int do_recip = Numeric->do_recip ;
#endif

    
    
    

    nz = 0 ;
    omega [0] = 0. ;
    omega [1] = 0. ;
    omega [2] = 0. ;
    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;
    Rs = Numeric->Rs ;		
    do_scale = (Rs != (double *) NULL) ;
    flops = 0 ;
    Info [SparseLU_SOLVE_FLOPS] = 0 ;
    Info [SparseLU_IR_TAKEN] = 0 ;
    Info [SparseLU_IR_ATTEMPTED] = 0 ;

    
    
    n = Numeric->n_row ;
    if (Numeric->nnzpiv < n
	|| SCALAR_IS_ZERO (Numeric->rcond) || SCALAR_IS_NAN (Numeric->rcond))
    {
	
	
	
	status = SparseLU_WARNING_singular_matrix ;
	irstep = 0 ;
    }
    else
    {
	status = SparseLU_OK ;
    }
    irstep = MAX (0, irstep) ;			

    W = (Entry *) SolveWork ;			

    Z = (Entry *) NULL ;	
    S = (Entry *) NULL ;
    Y = (double *) NULL ;
    Z2 = (double *) NULL ;
    B2 = (double *) NULL ;

#ifdef COMPLEX
    if (irstep > 0)
    {
	if (!Ap || !Ai || !Ax)
	{
	    return (SparseLU_ERROR_argument_missing) ;
	}
	
	AXsplit = SPLIT (Az) || SPLIT(Xz);
	Z = (Entry *) (SolveWork + 4*n) ;	
	S = (Entry *) (SolveWork + 6*n) ;	
	Y = (double *) (SolveWork + 8*n) ;	
	B2 = (double *) (SolveWork + 9*n) ;	
	Z2 = (double *) Z ;		
    }
    else
    {
      
      AXsplit = SPLIT(Xz);
    }
    Bsplit = SPLIT (Bz);

    if (AXsplit)
    {
	X = (Entry *) (SolveWork + 2*n) ;	
    }
    else
    {
	X = (Entry *) Xx ;			
    }
#else
    X = (Entry *) Xx ;				
    if (irstep > 0)
    {
	if (!Ap || !Ai || !Ax)
	{
	    return (SparseLU_ERROR_argument_missing) ;
	}
	Z = (Entry *) (SolveWork + n) ;		
	S = (Entry *) (SolveWork + 2*n) ;	
	Y = (double *) (SolveWork + 3*n) ;	
	B2 = (double *) (SolveWork + 4*n) ;	
	Z2 = (double *) Z ;		
    }
#endif

    
    
    

    if (sys == SparseLU_A)
    {

	
	
	

	if (irstep > 0)
	{

	    
	    
	    

	    nz = Ap [n] ;
	    Info [SparseLU_NZ] = nz ;

	    
	    
	    for (i = 0 ; i < n ; i++)
	    {
		Y [i] = 0. ;
	    }
	    flops += (ABS_FLOPS + 1) * nz ;
	    p2 = Ap [n] ;
	    for (p = 0 ; p < p2 ; p++)
	    {
		
	        ASSIGN (aij, Ax, Az, p, AXsplit) ;
		ABS (d, aij) ;
		Y [Ai [p]] += d ;
	    }

	    
	    flops += ABS_FLOPS * n ;
	    for (i = 0 ; i < n ; i++)
	    {
		
		ASSIGN (bi, Bx, Bz, i, Bsplit) ;
		ABS (B2 [i], bi) ;
	    }

	    
	    if (do_scale)
	    {
		
		
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    
		    for (i = 0 ; i < n ; i++)
		    {
			Y [i]  *= Rs [i] ;
			B2 [i] *= Rs [i] ;
		    }
		}
		else
#endif
		{
		    
		    for (i = 0 ; i < n ; i++)
		    {
			Y [i]  /= Rs [i] ;
			B2 [i] /= Rs [i] ;
		    }
		}

		flops += 2 * n ;
	    }

	}

	for (step = 0 ; step <= irstep ; step++)
	{

	    
	    
	    
	    
	    
	    

	    if (step == 0)
	    {
		if (do_scale)
		{
		    
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    ASSIGN (X [i], Bx, Bz, i, Bsplit) ;
			    SCALE (X [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    ASSIGN (X [i], Bx, Bz, i, Bsplit) ;
			    SCALE_DIV (X [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		    for (i = 0 ; i < n ; i++)
		    {
			W [i] = X [Rperm [i]] ;
		    }
		}
		else
		{
		    
		    for (i = 0 ; i < n ; i++)
		    {
			
			ASSIGN (W [i], Bx, Bz, Rperm [i], Bsplit) ;
		    }
		}
	    }
	    else
	    {
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSIGN (Z [i], Bx, Bz, i, Bsplit) ;
		}
		flops += MULTSUB_FLOPS * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    xi = X [i] ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			
			ASSIGN (aij, Ax, Az, p, AXsplit) ;
			MULT_SUB (Z [Ai [p]], aij, xi) ;
		    }
		}
		
		if (do_scale)
		{
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE (Z [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (Z [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}
		for (i = 0 ; i < n ; i++)
		{
		    W [i] = Z [Rperm [i]] ;
		}
	    }

	    flops += LU_lsolve (Numeric, W, Pattern) ; // lu_core.c
	    flops += LU_usolve (Numeric, W, Pattern) ; // lu_core.c

	    if (step == 0)
	    {
		for (i = 0 ; i < n ; i++)
		{
		    X [Cperm [i]] = W [i] ;
		}
	    }
	    else
	    {
		flops += ASSEMBLE_FLOPS * n ;
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSEMBLE (X [Cperm [i]], W [i]) ;
		}
	    }

	    
	    
	    

	    if (irstep > 0)
	    {

		
		
		
		
		

		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSIGN (W [i], Bx, Bz, i, Bsplit) ;
		    Z2 [i] = 0. ;
		}
		flops += (MULT_FLOPS + DECREMENT_FLOPS + ABS_FLOPS + 1) * nz ;
		for (j = 0 ; j < n ; j++)
		{
		    xj = X [j] ;
		    p2 = Ap [j+1] ;
		    for (p = Ap [j] ; p < p2 ; p++)
		    {
			i = Ai [p] ;

			
			ASSIGN (aij, Ax, Az, p, AXsplit) ;
			MULT (axx, aij, xj) ;

			
			DECREMENT (W [i], axx) ;

			
			ABS (d, axx) ;
			Z2 [i] += d ;
		    }
		}

		
		if (do_scale)
		{
		    
		    
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE (W [i], Rs [i]) ;
			    Z2 [i] *= Rs [i] ;
			}
		    }
		    else
#endif
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (W [i], Rs [i]) ;
			    Z2 [i] /= Rs [i] ;
			}
		    }
		    flops += (SCALE_FLOPS + 1) * n ;
		}

		flops += (2*ABS_FLOPS + 5) * n ;
		if (do_step (omega, step, B2, X, W, Y, Z2, S, n, Info))
		{
		    
		    break ;
		}

	    }

	}

    }
    else if (sys == SparseLU_At)
    {

	
	
	

	

	if (irstep > 0)
	{

	    
	    
	    

	    nz = Ap [n] ;
	    Info [SparseLU_NZ] = nz ;

	    
	    

	    if (do_scale)
	    {
		flops += (ABS_FLOPS + 2) * nz ;
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    
		    for (i = 0 ; i < n ; i++)
		    {
			yi = 0. ;
			p2 = Ap [i+1] ;
			for (p = Ap [i] ; p < p2 ; p++)
			{
			    
			    
			    ASSIGN (aij, Ax, Az, p, AXsplit) ;
			    ABS (d, aij) ;
			    yi += (d * Rs [Ai [p]]) ;
			}
			Y [i] = yi ;
		    }
		}
		else
#endif
		{
		    
		    for (i = 0 ; i < n ; i++)
		    {
			yi = 0. ;
			p2 = Ap [i+1] ;
			for (p = Ap [i] ; p < p2 ; p++)
			{
			    
			    
			    ASSIGN (aij, Ax, Az, p, AXsplit) ;
			    ABS (d, aij) ;
			    yi += (d / Rs [Ai [p]]) ;
			}
			Y [i] = yi ;
		    }
		}
	    }
	    else
	    {
		
		flops += (ABS_FLOPS + 1) * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    yi = 0. ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			
			
			ASSIGN (aij, Ax, Az, p, AXsplit) ;
			ABS (d, aij) ;
			yi += d ;
		    }
		    Y [i] = yi ;
		}
	    }

	    
	    for (i = 0 ; i < n ; i++)
	    {
		
		ASSIGN (bi, Bx, Bz, i, Bsplit) ;
		ABS (B2 [i], bi) ;
	    }

	}

	for (step = 0 ; step <= irstep ; step++)
	{

	    
	    
	    
	    
	    
	    

	    if (step == 0)
	    {
		
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSIGN (W [i], Bx, Bz, Cperm [i], Bsplit) ;
		}
	    }
	    else
	    {
		
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSIGN (Z [i], Bx, Bz, i, Bsplit) ;
		}
		flops += MULTSUB_FLOPS * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    zi = Z [i] ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			
			ASSIGN (aij, Ax, Az, p, Bsplit) ;
			MULT_SUB_CONJ (zi, X [Ai [p]], aij) ;
		    }
		    Z [i] = zi ;
		}
		
		for (i = 0 ; i < n ; i++)
		{
		    W [i] = Z [Cperm [i]] ;
		}
	    }

	    flops += LU_uhsolve (Numeric, W, Pattern) ;
	    flops += LU_lhsolve (Numeric, W, Pattern) ;

	    if (step == 0)
	    {

		
		

		
		for (i = 0 ; i < n ; i++)
		{
		    X [Rperm [i]] = W [i] ;
		}
		if (do_scale)
		{
		    
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE (X [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (X [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}

	    }
	    else
	    {

		
		for (i = 0 ; i < n ; i++)
		{
		    Z [Rperm [i]] = W [i] ;
		}
		if (do_scale)
		{
		    
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE (Z [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (Z [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}

		flops += ASSEMBLE_FLOPS * n ;
		
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSEMBLE (X [i], Z [i]) ;	
		}
	    }

	    
	    
	    

	    if (irstep > 0)
	    {

		
		
		
		
		

		flops += (MULT_FLOPS + DECREMENT_FLOPS + ABS_FLOPS + 1) * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSIGN (wi, Bx, Bz, i, Bsplit) ;
		    z2i = 0. ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			
			ASSIGN (aij, Ax, Az, p, AXsplit) ;
			MULT_CONJ (axx, X [Ai [p]], aij) ;

			
			DECREMENT (wi, axx) ;

			
			ABS (d, axx) ;
			z2i += d ;
		    }
		    W [i] = wi ;
		    Z2 [i] = z2i ;
		}

		flops += (2*ABS_FLOPS + 5) * n ;
		if (do_step (omega, step, B2, X, W, Y, Z2, S, n, Info))
		{
		    
		    break ;
		}

	    }

	}

    }
    else if (sys == SparseLU_Aat)
    {

	
	
	

	

	if (irstep > 0)
	{

	    
	    
	    

	    nz = Ap [n] ;
	    Info [SparseLU_NZ] = nz ;

	    
	    

	    if (do_scale)
	    {
		flops += (ABS_FLOPS + 2) * nz ;
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    
		    for (i = 0 ; i < n ; i++)
		    {
			yi = 0. ;
			p2 = Ap [i+1] ;
			for (p = Ap [i] ; p < p2 ; p++)
			{
			    
			    
			    ASSIGN (aij, Ax, Az, p, AXsplit) ;
			    ABS (d, aij) ;
			    yi += (d * Rs [Ai [p]]) ;
			}
			Y [i] = yi ;
		    }
		}
		else
#endif
		{
		    
		    for (i = 0 ; i < n ; i++)
		    {
			yi = 0. ;
			p2 = Ap [i+1] ;
			for (p = Ap [i] ; p < p2 ; p++)
			{
			    
			    
			    ASSIGN (aij, Ax, Az, p, AXsplit) ;
			    ABS (d, aij) ;
			    yi += (d / Rs [Ai [p]]) ;
			}
			Y [i] = yi ;
		    }
		}
	    }
	    else
	    {
		
		flops += (ABS_FLOPS + 1) * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    yi = 0. ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			
			
			ASSIGN (aij, Ax, Az, p, AXsplit) ;
			ABS (d, aij) ;
			yi += d ;
		    }
		    Y [i] = yi ;
		}
	    }

	    
	    for (i = 0 ; i < n ; i++)
	    {
		
		ASSIGN (bi, Bx, Bz, i, Bsplit) ;
		ABS (B2 [i], bi) ;
	    }

	}

	for (step = 0 ; step <= irstep ; step++)
	{

	    
	    
	    
	    
	    
	    

	    if (step == 0)
	    {
		
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSIGN (W [i], Bx, Bz, Cperm [i], Bsplit) ;
		}
	    }
	    else
	    {
		
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSIGN (Z [i], Bx, Bz, i, Bsplit) ;
		}
		flops += MULTSUB_FLOPS * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    zi = Z [i] ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			
			ASSIGN (aij, Ax, Az, p, AXsplit) ;
			MULT_SUB (zi, aij, X [Ai [p]]) ;
		    }
		    Z [i] = zi ;
		}
		
		for (i = 0 ; i < n ; i++)
		{
		    W [i] = Z [Cperm [i]] ;
		}
	    }

	    flops += LU_utsolve (Numeric, W, Pattern) ;
	    flops += LU_ltsolve (Numeric, W, Pattern) ;

	    if (step == 0)
	    {

		
		

		
		for (i = 0 ; i < n ; i++)
		{
		    X [Rperm [i]] = W [i] ;
		}
		if (do_scale)
		{
		    
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE (X [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (X [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}

	    }
	    else
	    {

		
		for (i = 0 ; i < n ; i++)
		{
		    Z [Rperm [i]] = W [i] ;
		}
		if (do_scale)
		{
		    
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE (Z [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (Z [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}

		flops += ASSEMBLE_FLOPS * n ;
		
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSEMBLE (X [i], Z [i]) ;	
		}
	    }

	    
	    
	    

	    if (irstep > 0)
	    {

		
		
		
		
		

		flops += (MULT_FLOPS + DECREMENT_FLOPS + ABS_FLOPS + 1) * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    
		    ASSIGN (wi, Bx, Bz, i, Bsplit) ;
		    z2i = 0. ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			
			ASSIGN (aij, Ax, Az, p, AXsplit) ;
			MULT (axx, aij, X [Ai [p]]) ;

			
			DECREMENT (wi, axx) ;

			
			ABS (d, axx) ;
			z2i += d ;
		    }
		    W [i] = wi ;
		    Z2 [i] = z2i ;
		}

		flops += (2*ABS_FLOPS + 5) * n ;
		if (do_step (omega, step, B2, X, W, Y, Z2, S, n, Info))
		{
		    
		    break ;
		}

	    }

	}

    }
    else if (sys == SparseLU_Pt_L)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (X [i], Bx, Bz, Rperm [i], Bsplit) ;
	}
	flops = LU_lsolve (Numeric, X, Pattern) ;//要隔
	status = SparseLU_OK ;

    }
    else if (sys == SparseLU_L)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (X [i], Bx, Bz, i, Bsplit) ;
	}
	flops = LU_lsolve (Numeric, X, Pattern) ;//要隔
	status = SparseLU_OK ;

    }
    else if (sys == SparseLU_Lt_P)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (W [i], Bx, Bz, i, Bsplit) ;
	}
	flops = LU_lhsolve (Numeric, W, Pattern) ;//要隔
	for (i = 0 ; i < n ; i++)
	{
	    X [Rperm [i]] = W [i] ;
	}
	status = SparseLU_OK ;

    }
    else if (sys == SparseLU_Lat_P)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (W [i], Bx, Bz, i, Bsplit) ;
	}
	flops = LU_ltsolve (Numeric, W, Pattern) ;//要隔
	for (i = 0 ; i < n ; i++)
	{
	    X [Rperm [i]] = W [i] ;
	}
	status = SparseLU_OK ;

    }
    else if (sys == SparseLU_Lt)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (X [i], Bx, Bz, i, Bsplit) ;
	}
	flops = LU_lhsolve (Numeric, X, Pattern) ;//要隔
	status = SparseLU_OK ;

    }
    else if (sys == SparseLU_Lat)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (X [i], Bx, Bz, i, Bsplit) ;
	}
	flops = LU_ltsolve (Numeric, X, Pattern) ;//要隔
	status = SparseLU_OK ;

    }
    else if (sys == SparseLU_U_Qt)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (W [i], Bx, Bz, i, Bsplit) ;
	}
	flops = LU_usolve (Numeric, W, Pattern) ;//要隔
	for (i = 0 ; i < n ; i++)
	{
	    X [Cperm [i]] = W [i] ;
	}

    }
    else if (sys == SparseLU_U)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (X [i], Bx, Bz, i, Bsplit) ;
	}
	flops = LU_usolve (Numeric, X, Pattern) ;//要隔

    }
    else if (sys == SparseLU_Q_Ut)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (X [i], Bx, Bz, Cperm [i], Bsplit) ;
	}
	flops = LU_uhsolve (Numeric, X, Pattern) ;//要隔

    }
    else if (sys == SparseLU_Q_Uat)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (X [i], Bx, Bz, Cperm [i], Bsplit) ;
	}
	flops = LU_utsolve (Numeric, X, Pattern) ;//要隔

    }
    else if (sys == SparseLU_Ut)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	  ASSIGN (X [i], Bx, Bz, i, Bsplit) ;
	}
	flops = LU_uhsolve (Numeric, X, Pattern) ; //要隔

    }
    else if (sys == SparseLU_Uat)
    {

	
	
	

	for (i = 0 ; i < n ; i++)
	{
	    
	    ASSIGN (X [i], Bx, Bz, i, Bsplit) ; //version
	}
	flops = LU_utsolve (Numeric, X, Pattern) ; //要隔

    }
    else
    {
	return (SparseLU_ERROR_invalid_system) ;
    }

#ifdef COMPLEX
    
    if (AXsplit)
    {
	for (i = 0 ; i < n ; i++)
	{
	    Xx [i] = REAL_COMPONENT (X [i]) ;
	    Xz [i] = IMAG_COMPONENT (X [i]) ;
	}
    }
#endif

    
    
    Info [SparseLU_SOLVE_FLOPS] = flops ;
    return (status) ;
}








PRIVATE Int do_step		
(
    double omega [3],
    Int step,			
    const double B2 [ ],	
    Entry X [ ],
    const Entry W [ ],
    const double Y [ ],
    const double Z2 [ ],
    Entry S [ ],
    Int n,
    double Info [SparseLU_INFO]
)
{
    double last_omega [3], tau, nctau, d1, wd1, d2, wd2, xi, yix, wi, xnorm ;
    Int i ;

    
    

    nctau = 1000 * n * DBL_EPSILON ;

    
    

    
    
    

    last_omega [0] = omega [0] ;
    last_omega [1] = omega [1] ;
    last_omega [2] = omega [2] ;

    
    
    

    
    xnorm = 0.0 ;
    for (i = 0 ; i < n ; i++)
    {
	
	ABS (xi, X [i]) ;
	if (SCALAR_IS_NAN (xi))
	{
	    xnorm = xi ;
	    break ;
	}
	
	xnorm = MAX (xnorm, xi) ;
    }

    omega [1] = 0. ;
    omega [2] = 0. ;
    for (i = 0 ; i < n ; i++)
    {
	yix = Y [i] * xnorm ;
	tau = (yix + B2 [i]) * nctau ;
	d1 = Z2 [i] + B2 [i] ;
	
	ABS (wi, W [i]) ;
	if (SCALAR_IS_NAN (d1))
	{
	    omega [1] = d1 ;
	    omega [2] = d1 ;
	    break ;
	}
	if (SCALAR_IS_NAN (tau))
	{
	    omega [1] = tau ;
	    omega [2] = tau ;
	    break ;
	}
	if (d1 > tau)		
	{
	    wd1 = wi / d1 ;
	    omega [1] = MAX (omega [1], wd1) ;
	}
	else if (tau > 0.0)	
	{
	    d2 = Z2 [i] + yix ;
	    wd2 = wi / d2 ;
	    omega [2] = MAX (omega [2], wd2) ;
	}
    }

    omega [0] = omega [1] + omega [2] ;
    Info [SparseLU_OMEGA1] = omega [1] ;
    Info [SparseLU_OMEGA2] = omega [2] ;

    
    
    

    Info [SparseLU_IR_TAKEN] = step ;
    Info [SparseLU_IR_ATTEMPTED] = step ;

    if (SCALAR_IS_NAN (omega [0]))
    {
	
	return (TRUE) ;
    }

    if (omega [0] < DBL_EPSILON)    
    {
	
	return (TRUE) ;
    }

    
    
    

    
    if (step > 0 && omega [0] > last_omega [0] / 2)
    {
	
	if (omega [0] > last_omega [0])
	{
	    
	    
	    for (i = 0 ; i < n ; i++)
	    {
		X [i] = S [i] ;
	    }
	    Info [SparseLU_OMEGA1] = last_omega [1] ;
	    Info [SparseLU_OMEGA2] = last_omega [2] ;
	}
	Info [SparseLU_IR_TAKEN] = step - 1 ;
	
	return (TRUE) ;
    }

    
    
    

    for (i = 0 ; i < n ; i++)
    {
	S [i] = X [i] ;
    }

    
    
    

    return (FALSE) ;
}
