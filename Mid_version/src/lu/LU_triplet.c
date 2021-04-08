/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_triplet.c
 * BRIEF:   LU三元组存储格式
 *****************************************************************************/

#include "SparseLU_internal.h"
#include "SparseLU_function.h"

#ifdef DO_MAP
#ifdef DO_VALUES
GLOBAL Int LU_triplet_map_x
#else
GLOBAL Int LU_triplet_map_nox
#endif
#else
#ifdef DO_VALUES
GLOBAL Int LU_triplet_nomap_x
#else
GLOBAL Int LU_triplet_nomap_nox
#endif
#endif
(
    Int n_row,
    Int n_col,
    Int nz,
    const Int Ti [ ],		
    const Int Tj [ ],		
    Int Ap [ ],			
    Int Ai [ ],			
    Int Rp [ ],			
    Int Rj [ ],			
    Int W [ ],			
    Int RowCount [ ]		
#ifdef DO_VALUES
    , const double Tx [ ]	
    , double Ax [ ]		
    , double Rx [ ]		
#endif
#ifdef DO_MAP
    , Int Map [ ]		
    , Int Map2 [ ]		
#endif
)
{

    
    
    

    Int i, j, k, p, cp, p1, p2, pdest, pj ;
#ifdef DO_MAP
    Int duplicates ;
#endif

    
    
    

    
    for (i = 0 ; i < n_row ; i++)
    {
	W [i] = 0 ;
    }

    for (k = 0 ; k < nz ; k++)
    {
	i = Ti [k] ;
	j = Tj [k] ;
	if (i < 0 || i >= n_row || j < 0 || j >= n_col)
	{
	    return (SparseLU_ERROR_invalid_matrix) ;
	}
	W [i]++ ;
    }

    
    
    

    Rp [0] = 0 ;
    for (i = 0 ; i < n_row ; i++)
    {
	Rp [i+1] = Rp [i] + W [i] ;
	W [i] = Rp [i] ;
    }

    

    
    
    

    for (k = 0 ; k < nz ; k++)
    {
	p = W [Ti [k]]++ ;
#ifdef DO_MAP
	Map [k] = p ;
#endif
	Rj [p] = Tj [k] ;
#ifdef DO_VALUES
	Rx [p] = Tx [k] ;
#endif
    }

    

    
    
    

    

    for (j = 0 ; j < n_col ; j++)
    {
	W [j] = EMPTY ;
    }

#ifdef DO_MAP
    duplicates = FALSE ;
#endif

    for (i = 0 ; i < n_row ; i++)
    {
	p1 = Rp [i] ;
	p2 = Rp [i+1] ;
	pdest = p1 ;
	
	
	for (p = p1 ; p < p2 ; p++)
	{
	    j = Rj [p] ;
	    ASSERT (j >= 0 && j < n_col) ;
	    pj = W [j] ;
	    if (pj >= p1)
	    {
		
		ASSERT (pj < p) ;
		ASSERT (Rj [pj] == j) ;
#ifdef DO_MAP
		Map2 [p] = pj ;
		duplicates = TRUE ;
#endif
#ifdef DO_VALUES
		
		Rx [pj] += Rx [p] ;
#endif
	    }
	    else
	    {
		
		
		W [j] = pdest ;
#ifdef DO_MAP
		Map2 [p] = pdest ;
#endif
		
		if (pdest != p)
		{
		    Rj [pdest] = j ;
#ifdef DO_VALUES
		    Rx [pdest] = Rx [p] ;
#endif
		}
		pdest++ ;
	    }
	}
	RowCount [i] = pdest - p1 ;
    }

    

    
    
    

#ifdef DO_MAP
    if (duplicates)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    Map [k] = Map2 [Map [k]] ;
	}
    }
#endif

    

    
    
    

    
    for (j = 0 ; j < n_col ; j++)
    {
	W [j] = 0 ;
    }

    for (i = 0 ; i < n_row ; i++)
    {
	for (p = Rp [i] ; p < Rp [i] + RowCount [i] ; p++)
	{
	    j = Rj [p] ;
	    ASSERT (j >= 0 && j < n_col) ;
	    W [j]++ ;
	}
    }

    
    
    

    Ap [0] = 0 ;
    for (j = 0 ; j < n_col ; j++)
    {
	Ap [j+1] = Ap [j] + W [j] ;
    }
    

    for (j = 0 ; j < n_col ; j++)
    {
	W [j] = Ap [j] ;
    }

    
    
    

    for (i = 0 ; i < n_row ; i++)
    {
	for (p = Rp [i] ; p < Rp [i] + RowCount [i] ; p++)
	{
	    cp = W [Rj [p]]++ ;
#ifdef DO_MAP
	    Map2 [p] = cp ;
#endif
	    Ai [cp] = i ;
#ifdef DO_VALUES
	    Ax [cp] = Rx [p] ;
#endif
	}
    }

    
    
    

#ifdef DO_MAP
    for (k = 0 ; k < nz ; k++)
    {
	Map [k] = Map2 [Map [k]] ;
    }
#endif

    

    return (SparseLU_OK) ;
}
