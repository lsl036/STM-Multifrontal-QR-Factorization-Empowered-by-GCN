/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_analyze.c
 * BRIEF:   LU分解分析
 *****************************************************************************/

#include "SparseLU_internal.h"
#include "SparseLU_function.h"



GLOBAL Int LU_analyze
(
    Int n_row,		
    Int n_col,
    Int Ai [ ],		

    Int Ap [ ],		

    Int Up [ ],		

    Int fixQ,


    Int W [ ],		
    Int Link [ ],	

   
    Int Front_ncols [ ],	
    Int Front_nrows [ ],	
    Int Front_npivcol [ ],	
    Int Front_parent [ ],	
    Int *nfr_out,

    Int *p_ncompactions		
)
{

	
    
    

    Int j, j3, col, k, row, parent, j2, pdest, p, p2, thickness, npivots, nfr,
	i, *Winv, kk, npiv, jnext, krow, knext, pfirst, jlast, ncompactions,
	*Front_stack, *Front_order, *Front_child, *Front_sibling,
	Wflag, npivcol, fallrows, fallcols, fpiv, frows, fcols, *Front_size ;

    nfr = 0 ;

    
    
    

#pragma ivdep
    for (j = 0 ; j < n_col ; j++)
    {
	Link [j] = EMPTY ;
	W [j] = EMPTY ;
	Up [j] = EMPTY ;

	
	Front_npivcol [j] = 0 ;		
	Front_nrows [j] = 0 ;		
	Front_ncols [j] = 0 ;		
	Front_parent [j] = EMPTY ;	
    }

   
    krow = 0 ;
    pfirst = Ap [0] ;
    jlast = EMPTY ;
    jnext = EMPTY ;
    Wflag = 0 ;

    
    ASSERT (pfirst >= n_col) ;	

    
    pdest = 0 ;
    ncompactions = 0 ;

    
    
    

    for (j = 0 ; j < n_col ; j = jnext)
    {

	

	if (pdest + (n_col-j) > pfirst)
	{
	   

	    pdest = 0 ;
	    ncompactions++ ;
	    for (row = 0 ; row < j ; row++)
	    {
		if (Up [row] != EMPTY)
		{
		    
		    p = Up [row] ;
		    p2 = p + (Front_ncols [row] - Front_npivcol [row]) ;
		    Up [row] = pdest ;
		    for ( ; p < p2 ; p++)
		    {
			Ai [pdest++] = Ai [p] ;
		    }
		}
	    }

	}

	if (pdest + (n_col-j) > pfirst)
	{
	   
	    return (FALSE) ;	
	}

	

	if (jlast != EMPTY && Link [j] == jlast)
	{
	    

	    Up [j] = Up [jlast] ;
	    Up [jlast] = EMPTY ;

	    
	    parent = n_col ;
	    for (p = Up [j] ; p < pdest ; )
	    {
		j3 = Ai [p] ;
		if (j == j3)
		{
		    Ai [p] = Ai [--pdest] ;
		}
		else
		{
		    if (j3 < parent)
		    {
			parent = j3 ;
		    }
		    p++ ;
		}
	    }

	    
	    Link [j] = Link [jlast] ;
	    thickness = (Front_nrows [jlast] - Front_npivcol [jlast]) ;

	}
	else
	{
	    Up [j] = pdest ;
	    parent = n_col ;
	    
	    thickness = 0 ;
	    Wflag = j ;
	}

	
	W [j] = Wflag ;

	

	jnext = n_col ;
	for (knext = krow ; knext < n_row ; knext++)
	{
	    jnext = Ai [Ap [knext]] ;
	    if (jnext != j)
	    {
		break ;
	    }
	}

	
	if (knext == n_row)
	{
	    jnext = n_col ;
	}

	

	for (k = krow ; k < knext ; k++)
	{
	    p = Ap [k] ;
	    p2 = Ap [k+1] ;

	    
	    for ( ; p < p2 ; p++)
	    {
		
		col = Ai [p] ;
		if (W [col] != Wflag)
		{
		    Ai [pdest++] = col ;
		    
		    W [col] = Wflag ;
		    if (col < parent)
		    {
			parent = col ;
		    }
		}
	    }
	    thickness++ ;
	}



	krow = knext ;
	pfirst = Ap [knext] ;

	
	for (k = Link [j] ; k != EMPTY ; k = Link [k])
	{
	   
	    p = Up [k] ;
	    p2 = p + (Front_ncols [k] - Front_npivcol [k]) ;
	    for ( ; p < p2 ; p++)
	    {
		
		col = Ai [p] ;
		if (W [col] != Wflag)
		{
		    Ai [pdest++] = col ;
		    
		    W [col] = Wflag ;
		    if (col < parent)
		    {
			parent = col ;
		    }
		}
	    }

	    
	    Up [k] = EMPTY ;
	    thickness += (Front_nrows [k] - Front_npivcol [k]) ;
	}

	

	for (j2 = j+1 ; j2 < jnext ; j2++)
	{
	    ASSERT (j2 >= 0 && j2 < n_col) ;
	    if (W [j2] != Wflag || Link [j2] != EMPTY)
	    {
		break ;
	    }
	}

	
	jnext = j2 ;
	j2-- ;

	npivots = j2-j+1 ;

	

	if (j2 > j)
	{
	   

	    parent = n_col ;
	    p2 = pdest ;
	    pdest = Up [j] ;
	    for (p = Up [j] ; p < p2 ; p++)
	    {
		col = Ai [p] ;
		if (col > j2)
		{
		   
		    Ai [pdest++] = col ;
		    if (col < parent)
		    {
			parent = col ;
		    }
		}
	    }
	}

	if (parent == n_col)
	{
	   
	    parent = EMPTY ;
	}

	

	npivcol = npivots ;
	fallrows = thickness ;
	fallcols = npivots + pdest - Up [j] ;

	
	fpiv = MIN (npivcol, fallrows) ;

	
	frows = fallrows - fpiv ;
	fcols = fallcols - fpiv ;

	if (frows == 0 || fcols == 0)
	{
	    
	    Up [j] = EMPTY ;
	    parent = EMPTY ;
	}

	Front_npivcol [j] = npivots ;
	Front_nrows [j] = fallrows ;
	Front_ncols [j] = fallcols ;
	Front_parent [j] = parent ;

	
	nfr++ ;


	if (parent != EMPTY)
	{
	    Link [j] = Link [parent] ;
	    Link [parent] = j ;
	}

	ASSERT (jnext > j) ;

	jlast = j ;
    }

   

    *nfr_out = nfr ;

    Front_order = W ;	

    if (fixQ)
    {
	
	k = 0 ;
	
#pragma novector
	for (j = 0 ; j < n_col ; j++)
	{
	    if (Front_npivcol [j] > 0)
	    {
		Front_order [j] = k++ ;
	    }
	    else
	    {
		Front_order [j] = EMPTY ;
	    }
	}
    }
    else
    {

	
	Front_child = Ap ;
	Front_sibling = Link ;

	
	Front_stack = Ai ;
	Front_size = Front_stack + n_col ;

	LU_fsize (n_col, Front_size, Front_nrows, Front_ncols,
	    Front_parent, Front_npivcol) ;

	AMD_postorder (n_col, Front_parent, Front_npivcol, Front_size,
	    Front_order, Front_child, Front_sibling, Front_stack) ;

	
	Winv = Ai ;
	for (k = 0 ; k < nfr ; k++)
	{
	    Winv [k] = EMPTY ;
	}

	
	for (i = 0 ; i < n_col ; i++)
	{
	    k = Front_order [i] ;
	    if (k != EMPTY)
	    {
		Winv [k] = i ;
	    }
	}

	
	kk = 0 ;
	for (k = 0 ; k < nfr ; k++)
	{
	    i = Winv [k] ;
	    for (npiv = 0 ; npiv < Front_npivcol [i] ; npiv++)
	    {
		Up [kk] = i + npiv ;
		kk++ ;
	    }
	}
	ASSERT (kk == n_col) ;

	
    }

  
   

    LU_apply_order (Front_npivcol, Front_order, Ai, n_col, nfr) ;
    LU_apply_order (Front_nrows,   Front_order, Ai, n_col, nfr) ;
    LU_apply_order (Front_ncols,   Front_order, Ai, n_col, nfr) ;
    LU_apply_order (Front_parent,  Front_order, Ai, n_col, nfr) ;


    for (i = 0 ; i < nfr ; i++)
    {
	parent = Front_parent [i] ;
	if (parent != EMPTY)
	{
	    ASSERT (parent >= 0 && parent < n_col) ;
	    ASSERT (Front_order [parent] >= 0 && Front_order [parent] < nfr) ;
	    Front_parent [i] = Front_order [parent] ;
	}
    }



    *p_ncompactions = ncompactions ;
    return (TRUE) ;
}




GLOBAL void LU_apply_order
(
    Int Front [ ],	   
    const Int Order [ ],    
			
    Int Temp [ ],	  
    Int nn,	
    Int nfr		 
)
{
    Int i, k ;
    for (i = 0 ; i < nn ; i++)
    {
	k = Order [i] ;
	ASSERT (k >= EMPTY && k < nfr) ;
	if (k != EMPTY)
	{
	    Temp [k] = Front [i] ;
	}
    }

    for (k = 0 ; k < nfr ; k++)
    {
	Front [k] = Temp [k] ;
    }
}


#include "Sparse.h"

#define SPARSE_start       SparseCore_start
#define SPARSE_transpose   SparseCore_transpose
#define SPARSE_analyze     SparseChol_analyze
#define SPARSE_free_sparse SparseCore_free_sparse
#define SPARSE_free_factor SparseCore_free_factor
#define SPARSE_finish      SparseCore_finish
#define SPARSE_print_common      SparseCore_print_common

int LU_chol
(

    Int nrow,             
    Int ncol,             
    Int symmetric,         
    Int Ap [ ],            
    Int Ai [ ],          

    Int Perm [ ],          

    void *user_params,     
    double user_info [3]   
                               
)
{
    double dmax, flops, c, lnz ;
    sparse_csc Amatrix, *A, *AT, *S ;
    sparse_factor *L ;
    sparse_common cm ;
    Int *P, *ColCount ;
    Int k, ordering_option, print_level, *params ;

    params = (Int *) user_params ;
    ordering_option = params [0] ;
    print_level = params [1] - 1 ;
    params [2] = -1 ;

    if (Ap == NULL || Ai == NULL || Perm == NULL || nrow < 0 || ncol < 0)
    {

        return (FALSE) ;
    }
    if (nrow != ncol)
    {
        symmetric = FALSE ;
    }


    SPARSE_start (&cm) ;
    cm.supernodal = SPARSE_SIMPLICIAL ;
    cm.print = print_level ;


    switch (ordering_option)
    {

        default:
        case SparseLU_ORDERING_AMD:
      
            cm.nmethods = 1 ;
            cm.method [0].ordering = symmetric ? SPARSE_AMD : SPARSE_COLAMD ;
            cm.postorder = TRUE ;
            break ;

        case SparseLU_ORDERING_NONE:
        case SparseLU_ORDERING_GIVEN:
        case SparseLU_ORDERING_USER:
       
            cm.nmethods = 1 ;
            cm.method [0].ordering = SPARSE_NATURAL ;
            cm.postorder = FALSE ;
            break ;
    }


    A = &Amatrix ;
    A->nrow = nrow ;    
    A->ncol = ncol ;
    A->nzmax = Ap [ncol] ;       
    A->packed = TRUE ;             
    if (symmetric)
    {
        A->stype = 1 ;              
    }
    else
    {
        A->stype = 0 ;             
    }
    A->itype = SPARSE_INT ;
    A->xtype = SPARSE_PATTERN ;
    A->dtype = SPARSE_DOUBLE ;
    A->nz = NULL ;
    A->p = Ap ;                 
    A->i = Ai ;                     
    A->x = NULL ;                  
    A->z = NULL ;
    A->sorted = FALSE ;           

    if (symmetric)
    {
      
        AT = NULL ;
        S = A ;
    }
    else
    {
       
        AT = SPARSE_transpose (A, 0, &cm) ;
        S = AT ;
    }


    L = SPARSE_analyze (S, &cm) ;
    SPARSE_free_sparse (&AT, &cm) ;
    if (L == NULL)
    {
        return (FALSE) ;
    }


    switch (L->ordering)
    {

        case SPARSE_AMD:
        case SPARSE_COLAMD:
            params [2] = SparseLU_ORDERING_AMD ;
            break ;

        case SPARSE_GIVEN:
        case SPARSE_NATURAL:
        default:
            params [2] = SparseLU_ORDERING_NONE ;
            break ;
    }

   
    P = L->Perm ;
    ColCount = L->ColCount ;
    dmax = 1 ;
    lnz = 0 ;
    flops = 0 ;
    for (k = 0 ; k < ncol ; k++)
    {
        Perm [k] = P [k] ;
        c = ColCount [k] ;
        if (c > dmax) dmax = c ;
        lnz += c ;
        flops += c*c ;
    }
    user_info [0] = dmax ;
    user_info [1] = lnz ;
    user_info [2] = flops ;

    SPARSE_free_factor (&L, &cm) ;
    if (print_level > 0) 
    {
        SPARSE_print_common ("for SparseLU", &cm) ;
    }
    SPARSE_finish (&cm) ;
    return (TRUE) ;

}



GLOBAL void *LU_free (void *p)
{
    if (p)
    {

        SparseBase_free (p) ;

    #if defined (LU_MALLOC_COUNT)
    
        LU_malloc_count-- ;
    #endif

    }

    return ((void *) NULL) ;
}



GLOBAL void LU_fsize
(
    Int nn,
    Int Fsize [ ],
    Int Fnrows [ ],
    Int Fncols [ ],
    Int Parent [ ],
    Int Npiv [ ]
)
{
    Int j, parent, frsize, r, c ;

    for (j = 0 ; j < nn ; j++)
    {
	Fsize [j] = EMPTY ;
    }

    
    
    
    for (j = 0 ; j < nn ; j++)
    {
	if (Npiv [j] > 0)
	{
	    
	    parent = Parent [j] ;
	    r = Fnrows [j] ;
	    c = Fncols [j] ;
	    frsize = r * c ;
	    
	    if (INT_OVERFLOW (((double) r) * ((double) c)))
	    {
		
		frsize = Int_MAX ;
	    }
	    Fsize [j] = MAX (Fsize [j], frsize) ;
	    if (parent != EMPTY)
	    {
		
		Fsize [parent] = MAX (Fsize [parent], Fsize [j]) ;
	    }
	}
    }
}




GLOBAL Int LU_is_permutation
(
    const Int P [ ],	
    Int W [ ],		
    Int n,
    Int r
)
{
    Int i, k ;

    if (!P)
    {
	
	return (TRUE) ;
    }

    for (i = 0 ; i < n ; i++)
    {
	W [i] = FALSE ;
    }
    for (k = 0 ; k < r ; k++)
    {
	i = P [k] ;
	if (i < 0 || i >= n)
	{
	    return (FALSE) ;
	}
	if (W [i])
	{
	    return (FALSE) ;
	}
	W [i] = TRUE ;
    }
    return (TRUE) ;
}



#if defined (LU_MALLOC_COUNT)



GLOBAL Int LU_malloc_count = 0 ;

#endif

#ifdef LU_TCOV_TEST

GLOBAL int lu_fail, lu_fail_lo, lu_fail_hi ;
GLOBAL int lu_realloc_fail, lu_realloc_lo, lu_realloc_hi ;
#endif

GLOBAL void *LU_malloc
(
    Int n_objects,
    size_t size_of_object
)
{
    size_t size ;
    void *p ;

#ifdef LU_TCOV_TEST
    
    
    lu_fail-- ;
    if (lu_fail <= lu_fail_hi && lu_fail >= lu_fail_lo)
    {
	return ((void *) NULL) ;
    }
#endif

    p = SparseBase_malloc (n_objects, size_of_object) ;

#if defined (LU_MALLOC_COUNT)
    if (p)
    {
	
	
	LU_malloc_count++ ;
    }
#endif

    return (p) ;
}



GLOBAL void *LU_realloc
(
    void *p,
    Int n_objects,
    size_t size_of_object
)
{
    size_t size ;
    void *p2 ;

#ifdef LU_TCOV_TEST
    
    
    lu_realloc_fail-- ;
    if (lu_realloc_fail <= lu_realloc_hi &&
	lu_realloc_fail >= lu_realloc_lo)
    {
	return ((void *) NULL) ;
    }
#endif

    
    n_objects = MAX (1, n_objects) ;

    size = (size_t) n_objects ;
    ASSERT (size_of_object > 1) ;
    if (size > Int_MAX / size_of_object)
    {
	
	return ((void *) NULL) ;
    }
    size *= size_of_object ;

    p2 = SparseBase_config.realloc_func (p, size) ;

#if defined (LU_MALLOC_COUNT)
    
    if (p == (void *) NULL && p2 != (void *) NULL)
    {
	LU_malloc_count++ ;
    }
#endif

    return (p2) ;
}
















PRIVATE void create_row_form
(
    
    Int n_row,		    
    Int n_col,
    const Int Ap [ ],	    
    const Int Ai [ ],	    
    Int Rdeg [ ],	    

    
    Int Rp [ ],		    
    Int Ri [ ],		    

    
    Int W [ ]		    
)
{
    Int row, col, p, p2 ;

    
    Rp [0] = 0 ;
    W [0] = 0 ;
    for (row = 0 ; row < n_row ; row++)
    {
	Rp [row+1] = Rp [row] + Rdeg [row] ;
	W [row] = Rp [row] ;
    }

    
    for (col = 0 ; col < n_col ; col++)
    {
	p2 = Ap [col+1] ;
	for (p = Ap [col] ; p < p2 ; p++)
	{
	    Ri [W [Ai [p]]++] = col ;
	}
    }
}





PRIVATE int order_singletons	
(
    Int k,	    
    Int head,
    Int tail,
    Int Next [ ],
    Int Xdeg [ ], Int Xperm [ ], const Int Xp [ ], const Int Xi [ ],
    Int Ydeg [ ], Int Yperm [ ], const Int Yp [ ], const Int Yi [ ]
)
{
    Int xpivot, x, y, ypivot, p, p2, deg ;

    while (head != EMPTY)
    {
	
	xpivot = head ;
	head = Next [xpivot] ;
	if (head == EMPTY) tail = EMPTY ;

	if (Xdeg [xpivot] != 1)
	{
	    
	    continue ;
	}

	

	ypivot = EMPTY ;
	p2 = Xp [xpivot+1] ;
	for (p = Xp [xpivot] ; p < p2 ; p++)
	{
	    y = Xi [p] ;
	    if (Ydeg [y] >= 0)
	    {
		
		ypivot = y ;
		break ;
	    }
	}

	
	p2 = Yp [ypivot+1] ;
	for (p = Yp [ypivot] ; p < p2 ; p++)
	{
	    x = Yi [p] ;
	    if (Xdeg [x] < 0) continue ;
	    if (x == xpivot) continue ;
	    deg = --(Xdeg [x]) ;
	    if (deg == 1)
	    {
		
		Next [x] = EMPTY ;
		if (head == EMPTY)
		{
		    head = x ;
		}
		else
		{
		    Next [tail] = x ;
		}
		tail = x ;
	    }
	}

	
	Xdeg [xpivot] = FLIP (1) ;
	Ydeg [ypivot] = FLIP (Ydeg [ypivot]) ;

	
	Xperm [k] = xpivot ;
	Yperm [k] = ypivot ;
	k++ ;
    }

    return (k) ;
}





PRIVATE Int find_any_singletons	    
(
    
    Int n_row,
    Int n_col,
    const Int Ap [ ],	    
    const Int Ai [ ],	    

    
    Int Cdeg [ ],	    
    Int Rdeg [ ],	    

    
    Int Cperm [ ],	    
    Int Rperm [ ],	    
    Int *p_n1r,		    
    Int *p_n1c,		    

    
    Int Rp [ ],		    
    Int Ri [ ],		    
    Int W [ ],		    
    Int Next [ ]	    
)
{
    Int n1, col, row, row_form, head, tail, n1r, n1c ;

    
    
    

    n1 = 0 ;
    n1r = 0 ;
    n1c = 0 ;
    row_form = FALSE ;

    head = EMPTY ;
    tail = EMPTY ;
    for (col = n_col-1 ; col >= 0 ; col--)
    {
	if (Cdeg [col] == 1)
	{
	    
	    if (head == EMPTY) tail = col ;
	    Next [col] = head ;
	    head = col ;
	}
    }

    if (head != EMPTY)
    {

	
	
	

	create_row_form (n_row, n_col, Ap, Ai, Rdeg, Rp, Ri, W) ;
	row_form = TRUE ;

	
	
	

	n1 = order_singletons (0, head, tail, Next,
		Cdeg, Cperm, Ap, Ai,
		Rdeg, Rperm, Rp, Ri
		) ;
	n1c = n1 ;
    }

    
    
    

    head = EMPTY ;
    tail = EMPTY ;
    for (row = n_row-1 ; row >= 0 ; row--)
    {
	if (Rdeg [row] == 1)
	{
	    
	    if (head == EMPTY) tail = row ;
	    Next [row] = head ;
	    head = row ;
	}
    }

    if (head != EMPTY)
    {

	
	
	

	if (!row_form)
	{
	    create_row_form (n_row, n_col, Ap, Ai, Rdeg, Rp, Ri, W) ;
	}

	
	
	

	n1 = order_singletons (n1, head, tail, Next,
		Rdeg, Rperm, Rp, Ri,
		Cdeg, Cperm, Ap, Ai
		) ;
	n1r = n1 - n1c ;
    }
    *p_n1r = n1r ;
    *p_n1c = n1c ;
    return (n1) ;
}





PRIVATE Int find_user_singletons	
(
    
    Int n_row,
    Int n_col,
    const Int Ap [ ],	    
    const Int Ai [ ],	    
    const Int Quser [ ],    

    
    Int Cdeg [ ],	    
    Int Rdeg [ ],	    

    
    Int Cperm [ ],	    
    Int Rperm [ ],	    
    Int *p_n1r,		    
    Int *p_n1c,		    

    
    Int Rp [ ],		    
    Int Ri [ ],		    
    Int W [ ]		    
)
{
    Int n1, col, row, p, p2, pivcol, pivrow, found, k, n1r, n1c ;

    n1 = 0 ;
    n1r = 0 ;
    n1c = 0 ;
    *p_n1r = 0 ;
    *p_n1c = 0 ;

    
    pivcol = Quser [0] ;
    found = (Cdeg [pivcol] == 1) ;
    if (!found)
    {
	
	for (p = Ap [pivcol] ; p < Ap [pivcol+1] ; p++)
	{
	    if (Rdeg [Ai [p]] == 1)
	    {
		found = TRUE ;
		break ;
	    }
	}
    }

    if (!found)
    {
	
	return (0) ;
    }

    
    create_row_form (n_row, n_col, Ap, Ai, Rdeg, Rp, Ri, W) ;

    n1 = 0 ;

    for (k = 0 ; k < n_col ; k++)
    {
	pivcol = Quser [k] ;
	pivrow = EMPTY ;

	
	
	

	found = (Cdeg [pivcol] == 1) ;

	if (found)
	{

	    
	    
	    

	    

	    p2 = Ap [pivcol+1] ;
	    for (p = Ap [pivcol] ; p < p2 ; p++)
	    {
		row = Ai [p] ;
		if (Rdeg [row] >= 0)
		{
		    
		    pivrow = row ;
		    break ;
		}
	    }

	    
	    p2 = Rp [pivrow+1] ;
	    for (p = Rp [pivrow] ; p < p2 ; p++)
	    {
		col = Ri [p] ;
		if (Cdeg [col] < 0) continue ;
		Cdeg [col]-- ;
	    }

	    
	    Cdeg [pivcol] = FLIP (1) ;
	    Rdeg [pivrow] = FLIP (Rdeg [pivrow]) ;
	    n1c++ ;

	}
	else
	{

	    
	    
	    

	    p2 = Ap [pivcol+1] ;
	    for (p = Ap [pivcol] ; p < p2 ; p++)
	    {
		pivrow = Ai [p] ;
		if (Rdeg [pivrow] == 1)
		{
		    found = TRUE ;
		    break ;
		}
	    }

	    if (!found)
	    {
		break ;
	    }

	    
	    p2 = Ap [pivcol+1] ;
	    for (p = Ap [pivcol] ; p < p2 ; p++)
	    {
		row = Ai [p] ;
		if (Rdeg [row] < 0) continue ;
		Rdeg [row]-- ;
	    }

	    
	    Cdeg [pivcol] = FLIP (Cdeg [pivcol]) ;
	    Rdeg [pivrow] = FLIP (1) ;
	    n1r++ ;
	}

	
	Cperm [k] = pivcol ;
	Rperm [k] = pivrow ;
	n1++ ;

    }

    *p_n1r = n1r ;
    *p_n1c = n1c ;
    return (n1) ;
}







PRIVATE Int finish_permutation
(
    Int n1,
    Int nx,
    Int Xdeg [ ],
    const Int Xuser [ ],
    Int Xperm [ ],
    Int *p_max_deg
)
{
    Int nempty, x, deg, s, max_deg, k ;
    nempty = 0 ;
    s = n1 ;
    max_deg = 0 ;
    for (k = 0 ; k < nx ; k++)
    {
	x = (Xuser != (Int *) NULL) ? Xuser [k] : k ;
	deg = Xdeg [x] ;
	if (deg == 0)
	{
	    
	    nempty++ ;
	    Xperm [nx - nempty] = x ;
	}
	else if (deg > 0)
	{
	    
	    ASSERT (s < nx - nempty) ;
	    Xperm [s++] = x ;
	    max_deg = MAX (max_deg, deg) ;
	}
	else
	{
	    
	    Xdeg [x] = FLIP (deg) ;
	}
    }
    ASSERT (s == nx - nempty) ;
    *p_max_deg = max_deg ;
    return (nempty) ;
}





GLOBAL Int LU_singletons
(

    
    Int n_row,
    Int n_col,
    const Int Ap [ ],	    
    const Int Ai [ ],	    
    const Int Quser [ ],    
    Int strategy,	    
    Int do_singletons,      

    
    Int Cdeg [ ],	
    Int Cperm [ ],	
    Int Rdeg [ ],	
    Int Rperm [ ],	
    Int InvRperm [ ],	
    Int *p_n1,		
    Int *p_n1c,		
    Int *p_n1r,		
    Int *p_nempty_col,	
    Int *p_nempty_row,	
    Int *p_is_sym,	
    Int *p_max_rdeg,	

    
    Int Rp [ ],		
    Int Ri [ ],		
    Int W [ ],		
    Int Next [ ]	
)
{
    Int n1, s, col, row, p, p1, p2, cdeg, last_row, is_sym, k,
	nempty_row, nempty_col, max_cdeg, max_rdeg, n1c, n1r ;

    
    
    

    
    
    

    if (Ap [0] != 0 || Ap [n_col] < 0)
    {
	return (SparseLU_ERROR_invalid_matrix) ;
    }
    for (row = 0 ; row < n_row ; row++)
    {
	Rdeg [row] = 0 ;
    }
    for (col = 0 ; col < n_col ; col++)
    {
	p1 = Ap [col] ;
	p2 = Ap [col+1] ;
	cdeg = p2 - p1 ;
	if (cdeg < 0)
	{
	    return (SparseLU_ERROR_invalid_matrix) ;
	}
	last_row = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    row = Ai [p] ;
	    if (row <= last_row || row >= n_row)
	    {
		return (SparseLU_ERROR_invalid_matrix) ;
	    }
	    Rdeg [row]++ ;
	    last_row = row ;
	}
	Cdeg [col] = cdeg ;
    }

    
    
    

    if (!do_singletons)
    {
        
        n1 = 0 ;
        n1r = 0 ;
        n1c = 0 ;
    }
    else if (Quser != (Int *) NULL)
    {
	
	if (strategy == SparseLU_STRATEGY_UNSYMMETRIC)
	{
	    
	    n1 = find_user_singletons (n_row, n_col, Ap, Ai, Quser,
		    Cdeg, Rdeg, Cperm, Rperm, &n1r, &n1c, Rp, Ri, W) ;
	}
	else
	{
	    
	    n1 = 0 ;
	    n1r = 0 ;
	    n1c = 0 ;
	}
    }
    else
    {
	
	n1 = find_any_singletons (n_row, n_col, Ap, Ai,
		Cdeg, Rdeg, Cperm, Rperm, &n1r, &n1c, Rp, Ri, W, Next) ;
    }

    
    
    

    nempty_col = finish_permutation (n1, n_col, Cdeg, Quser, Cperm, &max_cdeg) ;

    
    
    

    if (Quser != (Int *) NULL && strategy == SparseLU_STRATEGY_SYMMETRIC)
    {
	
	ASSERT (n_row == n_col) ;
	nempty_row = finish_permutation (n1, n_row, Rdeg, Quser, Rperm,
	    &max_rdeg) ;
    }
    else
    {
	
	nempty_row = finish_permutation (n1, n_row, Rdeg, (Int *) NULL, Rperm,
	    &max_rdeg) ;
    }

    
    
    

    for (k = 0 ; k < n_row ; k++)
    {
	ASSERT (Rperm [k] >= 0 && Rperm [k] < n_row) ;
	InvRperm [Rperm [k]] = k ;
    }

    
    
    

    

    if (n_row == n_col && nempty_row == nempty_col)
    {
	
	is_sym = TRUE ;
	for (s = n1 ;  is_sym &&
            
            s < n_col - nempty_col ; s++)
	{
	    if (Cperm [s] != Rperm [s])
	    {
		is_sym = FALSE ;
		
	    }
	}
    }
    else
    {
	is_sym = FALSE ;
    }

    *p_n1 = n1 ;
    *p_n1r = n1r ;
    *p_n1c = n1c ;
    *p_is_sym = is_sym ;
    *p_nempty_col = nempty_col ;
    *p_nempty_row = nempty_row ;
    *p_max_rdeg = max_rdeg ;
    return (SparseLU_OK) ;
}
























#define ALIVE	(0)
#define DEAD	(-1)


#define DEAD_PRINCIPAL		(-1)
#define DEAD_NON_PRINCIPAL	(-2)


#define ROW_IS_DEAD(r)			ROW_IS_MARKED_DEAD (Row[r].shared2.mark)
#define ROW_IS_MARKED_DEAD(row_mark)	(row_mark < ALIVE)
#define ROW_IS_ALIVE(r)			(Row [r].shared2.mark >= ALIVE)
#define COL_IS_DEAD(c)			(Col [c].start < ALIVE)
#define COL_IS_ALIVE(c)			(Col [c].start >= ALIVE)
#define COL_IS_DEAD_PRINCIPAL(c)	(Col [c].start == DEAD_PRINCIPAL)
#define KILL_ROW(r)			{ Row [r].shared2.mark = DEAD ; }
#define KILL_PRINCIPAL_COL(c)		{ Col [c].start = DEAD_PRINCIPAL ; }
#define KILL_NON_PRINCIPAL_COL(c)	{ Col [c].start = DEAD_NON_PRINCIPAL ; }









PRIVATE Int init_rows_cols
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int p []
    
) ;

PRIVATE void init_scoring
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    double knobs [COLAMD_KNOBS],
    Int *p_n_row2,
    Int *p_n_col2,
    Int *p_max_deg
    
    
    , Int *p_ndense_row		
    , Int *p_nempty_row		
    , Int *p_nnewlyempty_row	
    , Int *p_ndense_col		
    , Int *p_nempty_col		
    , Int *p_nnewlyempty_col	
) ;

PRIVATE Int find_ordering
(
    Int n_row,
    Int n_col,
    Int Alen,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int n_col2,
    Int max_deg,
    Int pfree
    
    
    , Int Front_npivcol [ ]
    , Int Front_nrows [ ]
    , Int Front_ncols [ ]
    , Int Front_parent [ ]
    , Int Front_cols [ ]
    , Int *p_nfr
    , Int aggressive
    , Int InFront [ ]
    
) ;





PRIVATE void detect_super_cols
(
    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int row_start,
    Int row_length
) ;

PRIVATE Int garbage_collection
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int *pfree
) ;

PRIVATE Int clear_mark
(
    Int n_row,
    Colamd_Row Row []
) ;
















GLOBAL void LU_colamd_set_defaults
(
    

    double knobs [COLAMD_KNOBS]		
)
{
    

    Int i ;

#if 0
    if (!knobs)
    {
	return ;			
    }
#endif
    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
	knobs [i] = 0 ;
    }
    knobs [COLAMD_DENSE_ROW] = 0.2 ;	
    knobs [COLAMD_DENSE_COL] = 0.2 ;	
    knobs [COLAMD_AGGRESSIVE] = TRUE ;	
}









GLOBAL Int LU_colamd		
(
    

    Int n_row,			
    Int n_col,			
    Int Alen,			
    Int A [],			
    Int p [],			
    double knobs [COLAMD_KNOBS],
    Int stats [COLAMD_STATS]	

    
    
    , Int Front_npivcol [ ]	
    , Int Front_nrows [ ]	
    , Int Front_ncols [ ]	
    , Int Front_parent [ ]	
    , Int Front_cols [ ]	
    , Int *p_nfr		
    , Int InFront [ ]		
    
)
{
    

    Int row ;			
    Int i ;			
    Int nnz ;			
    Int Row_size ;		
    Int Col_size ;		
#if 0
    Int need ;			
#endif
    Colamd_Row *Row ;		
    Colamd_Col *Col ;		
    Int n_col2 ;		
    Int n_row2 ;		
    Int ngarbage ;		
    Int max_deg ;		
    Int aggressive ;		
#if 0
    double default_knobs [COLAMD_KNOBS] ;	
#endif

    
    
    Int ndense_row, nempty_row, parent, ndense_col,
	nempty_col, k, col, nfr, *Front_child, *Front_sibling, *Front_stack,
	*Front_order, *Front_size ;
    Int nnewlyempty_col, nnewlyempty_row ;
    

    

#if 0
    if (!stats)
    {
	return (FALSE) ;	
    }
#endif

    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

#if 0
    if (!A)		
    {
	
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
	return (FALSE) ;
    }

    if (!p)		
    {
	
	stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
	return (FALSE) ;
    }

    if (n_row < 0)	
    {
	
	stats [COLAMD_STATUS] = COLAMD_ERROR_nrow_negative ;
	stats [COLAMD_INFO1] = n_row ;
	return (FALSE) ;
    }

    if (n_col < 0)	
    {
	
	stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
	stats [COLAMD_INFO1] = n_col ;
	return (FALSE) ;
    }
#endif

    nnz = p [n_col] ;

#if 0
    if (nnz < 0)	
    {
	
	stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
	stats [COLAMD_INFO1] = nnz ;
	return (FALSE) ;
    }

    if (p [0] != 0)	
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero	;
	stats [COLAMD_INFO1] = p [0] ;
	return (FALSE) ;
    }
#endif

    

#if 0
    if (!knobs)
    {
	
	LU_colamd_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }
#endif

    ASSERT (knobs != (double *) NULL) ;

    
    
    aggressive = (knobs [COLAMD_AGGRESSIVE] != 0) ;
    

    

    Col_size = LU_COLAMD_C (n_col) ;
    Row_size = LU_COLAMD_R (n_row) ;

#if 0
    need = MAX (2*nnz, 4*n_col) + n_col + Col_size + Row_size ;
    if (need > Alen)
    {
	
	
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_too_small ;
	stats [COLAMD_INFO1] = need ;
	stats [COLAMD_INFO2] = Alen ;
	return (FALSE) ;
    }
#endif

    Alen -= Col_size + Row_size ;
    Col = (Colamd_Col *) &A [Alen] ;
    Row = (Colamd_Row *) &A [Alen + Col_size] ;

    

    

    i = init_rows_cols (n_row, n_col, Row, Col, A, p) ;

#if 0
    if (!i)
    {
	
	return (FALSE) ;
    }
#endif

    

    for (col = 0 ; col < n_col ; col++)
    {
	Front_npivcol [col] = 0 ;
	Front_nrows [col] = 0 ;
	Front_ncols [col] = 0 ;
	Front_parent [col] = EMPTY ;
	Front_cols [col] = EMPTY ;
    }

    

    init_scoring (n_row, n_col, Row, Col, A, p, knobs,
	&n_row2, &n_col2, &max_deg
	
	
	, &ndense_row, &nempty_row, &nnewlyempty_row
	, &ndense_col, &nempty_col, &nnewlyempty_col
	
	) ;
    ASSERT (n_row2 == n_row - nempty_row - nnewlyempty_row - ndense_row) ;
    ASSERT (n_col2 == n_col - nempty_col - nnewlyempty_col - ndense_col) ;

    

    ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
	n_col2, max_deg, 2*nnz
	
	
	, Front_npivcol, Front_nrows, Front_ncols, Front_parent, Front_cols
	, &nfr, aggressive, InFront
	
	) ;

    
    

    
    
    Front_child   = A ;
    Front_sibling = Front_child + nfr ;
    Front_stack   = Front_sibling + nfr ;
    Front_order   = Front_stack + nfr ;
    Front_size    = Front_order + nfr ;

    LU_fsize (nfr, Front_size, Front_nrows, Front_ncols,
	    Front_parent, Front_npivcol) ;

    AMD_postorder (nfr, Front_parent, Front_npivcol, Front_size,
	Front_order, Front_child, Front_sibling, Front_stack) ;

    

    
    LU_apply_order (Front_npivcol, Front_order, A, nfr, nfr) ;
    LU_apply_order (Front_nrows,   Front_order, A, nfr, nfr) ;
    LU_apply_order (Front_ncols,   Front_order, A, nfr, nfr) ;
    LU_apply_order (Front_parent,  Front_order, A, nfr, nfr) ;
    LU_apply_order (Front_cols,    Front_order, A, nfr, nfr) ;

    
    for (i = 0 ; i < nfr ; i++)
    {
	parent = Front_parent [i] ;
	if (parent != EMPTY)
	{
	    Front_parent [i] = Front_order [parent] ;
	}
    }

    
    for (row = 0 ; row < n_row ; row++)
    {
	i = InFront [row] ;
	ASSERT (i >= EMPTY && i < nfr) ;
	if (i != EMPTY)
	{
	    InFront [row] = Front_order [i] ;
	}
    }

    

    

    
    for (i = 0 ; i < n_col ; i++)
    {
	A [i] = EMPTY ;
    }
    k = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	for (col = Front_cols [i] ; col != EMPTY ; col = Col [col].nextcol)
	{
	    p [k] = col ;
	    A [col] = k ;
	    k++ ;
	}
    }

    

    if (n_col2 < n_col)
    {
	for (col = 0 ; col < n_col ; col++)
	{
	    if (A [col] == EMPTY)
	    {
		k = Col [col].shared2.order ;
		p [k] = col ;
		A [col] = k ;
	    }
	}
    }

    

    

    
    
    stats [COLAMD_DENSE_ROW] = ndense_row ;
    stats [COLAMD_EMPTY_ROW] = nempty_row ;
    stats [COLAMD_NEWLY_EMPTY_ROW] = nnewlyempty_row ;
    stats [COLAMD_DENSE_COL] = ndense_col ;
    stats [COLAMD_EMPTY_COL] = nempty_col ;
    stats [COLAMD_NEWLY_EMPTY_COL] = nnewlyempty_col ;
    
    stats [COLAMD_DEFRAG_COUNT] = ngarbage ;
    *p_nfr = nfr ;
    return (TRUE) ;
}









PRIVATE Int init_rows_cols	
(
    

    Int n_row,			
    Int n_col,			
    Colamd_Row Row [],		
    Colamd_Col Col [],		
    Int A [],			
    Int p []			

)
{
    

    Int col ;			
    Int row ;			
    Int *cp ;			
    Int *cp_end ;		

    

    for (col = 0 ; col < n_col ; col++)
    {
	Col [col].start = p [col] ;
	Col [col].length = p [col+1] - p [col] ;

#if 0
	if (Col [col].length < 0)
	{
	    
	    stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
	    stats [COLAMD_INFO1] = col ;
	    stats [COLAMD_INFO2] = Col [col].length ;
	    return (FALSE) ;
	}
#endif

	
	Col [col].shared1.thickness = 1 ;
	Col [col].shared2.score = 0 ;
	Col [col].shared3.prev = EMPTY ;
	Col [col].shared4.degree_next = EMPTY ;

	
	
	Col [col].nextcol = EMPTY ;
	Col [col].lastcol = col ;
	
    }

    

    

    
    
    
    

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].length = 0 ;
	
	
	
	
	
	
	Row [row].thickness = 1 ;
	Row [row].front = EMPTY ;
	
    }

    for (col = 0 ; col < n_col ; col++)
    {

	cp = &A [p [col]] ;
	cp_end = &A [p [col+1]] ;

	while (cp < cp_end)
	{
	    row = *cp++ ;

#if 0
	    
	    if (row < 0 || row >= n_row)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		
		
		
		
		return (FALSE) ;
	    }
#endif

	    ASSERT (row >= 0 && row < n_row) ;

#if 0
	    
	    
	    if (row <= last_row)
	    {
		
		
		stats [COLAMD_STATUS] = COLAMD_ERROR_jumbled_matrix ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;

		return (FALSE) ;
	    }
	    
#endif

	    ASSERT (row > last_row) ;

	    
	    
	    Row [row].length++ ;
	    
	}
    }

    

    
    
    Row [0].start = p [n_col] ;
    Row [0].shared1.p = Row [0].start ;
    
    
    
    
    for (row = 1 ; row < n_row ; row++)
    {
	Row [row].start = Row [row-1].start + Row [row-1].length ;
	Row [row].shared1.p = Row [row].start ;
	
	
	
	
    }

    

    
    
    

	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		A [(Row [*cp++].shared1.p)++] = col ;
	    }
	}

    

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].shared2.mark = 0 ;
	Row [row].shared1.degree = Row [row].length ;
    }

    
    
    

    return (TRUE) ;
}








PRIVATE void init_scoring
(
    

    Int n_row,			
    Int n_col,			
    Colamd_Row Row [],		
    Colamd_Col Col [],		
    Int A [],			
    Int head [],		
    double knobs [COLAMD_KNOBS],
    Int *p_n_row2,		
    Int *p_n_col2,		
    Int *p_max_deg		
    
    
    , Int *p_ndense_row		
    , Int *p_nempty_row		
    , Int *p_nnewlyempty_row	
    , Int *p_ndense_col		
    , Int *p_nempty_col		
    , Int *p_nnewlyempty_col	
    
)
{
    

    Int c ;			
    Int r, row ;		
    Int *cp ;			
    Int deg ;			
    Int *cp_end ;		
    Int *new_cp ;		
    Int col_length ;		
    Int score ;			
    Int n_col2 ;		
    Int n_row2 ;		
    Int dense_row_count ;	
    Int dense_col_count ;	
    Int min_score ;		
    Int max_deg ;		
    Int next_col ;		

    
    
    Int ndense_row ;		
    Int nempty_row ;		
    Int nnewlyempty_row ;	
    Int ndense_col ;		
    Int nempty_col ;		
    Int nnewlyempty_col ;	
    Int ne ;
    

    

    
    
    
    
    dense_row_count =
	SparseLU_DENSE_DEGREE_THRESHOLD (knobs [COLAMD_DENSE_ROW], n_col) ;
    dense_col_count =
	SparseLU_DENSE_DEGREE_THRESHOLD (knobs [COLAMD_DENSE_COL], n_row) ;
    
    dense_row_count = MAX (0, MIN (dense_row_count, n_col)) ;
    dense_col_count = MAX (0, MIN (dense_col_count, n_row)) ;
    

    max_deg = 0 ;
    n_col2 = n_col ;
    n_row2 = n_row ;

    
    
    ndense_col = 0 ;
    nempty_col = 0 ;
    nnewlyempty_col = 0 ;
    ndense_row = 0 ;
    nempty_row = 0 ;
    nnewlyempty_row = 0 ;
    

    

    

#if 0
    
    
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	deg = Col [c].length ;
	if (deg == 0)
	{
	    
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	    
	    
	    nempty_col++ ;
	    
	}
    }
#endif

    

#if 0
    for (r = 0 ; r < n_row ; r++)
    {
	deg = Row [r].shared1.degree ;
	if (deg == 0)
	{
	    
	    nempty_row++ ;
	}
    }
#endif

    

    
    for (c = n_col-1 ; c >= 0 ; c--)
    {

	
#if 0
	
	
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
#endif
	

	deg = Col [c].length ;
	if (deg > dense_col_count)
	{
	    
	    Col [c].shared2.order = --n_col2 ;
	    
	    
	    ndense_col++ ;
	    
	    
	    cp = &A [Col [c].start] ;
	    cp_end = cp + Col [c].length ;
	    while (cp < cp_end)
	    {
		Row [*cp++].shared1.degree-- ;
	    }
	    KILL_PRINCIPAL_COL (c) ;
	}
    }

    

    

    ne = 0 ;
    for (r = 0 ; r < n_row ; r++)
    {
	deg = Row [r].shared1.degree ;
	
	
	if (deg > dense_row_count)
	{
	    
	    
	    ndense_row++ ;
	}
	if (deg == 0)
	{
	    
	    ne++ ;
	}
	
	if (deg > dense_row_count || deg == 0)
	{
	    
	    KILL_ROW (r) ;
	    
	    
	    Row [r].thickness = 0 ;
	    
	    --n_row2 ;
	}
	else
	{
	    
	    max_deg = MAX (max_deg, deg) ;
	}
    }
    nnewlyempty_row = ne - nempty_row ;

    

    
    
    
    

    
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	score = 0 ;
	cp = &A [Col [c].start] ;
	new_cp = cp ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    
	    row = *cp++ ;
	    
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }
	    
	    *new_cp++ = row ;
	    
	    score += Row [row].shared1.degree - 1 ;
	    
	    score = MIN (score, n_col) ;
	}
	
	col_length = (Int) (new_cp - &A [Col [c].start]) ;
	if (col_length == 0)
	{
	    
	    
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	    
	    
	    nnewlyempty_col++ ;
	    
	}
	else
	{
	    
	    Col [c].length = col_length ;
	    Col [c].shared2.score = score ;
	}
    }

    
    
    
    

    

    
    for (c = 0 ; c <= n_col ; c++)
    {
	head [c] = EMPTY ;
    }
    min_score = n_col ;
    
    
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	
	if (COL_IS_ALIVE (c))
	{

	    

	    score = Col [c].shared2.score ;

	    
	    next_col = head [score] ;
	    Col [c].shared3.prev = EMPTY ;
	    Col [c].shared4.degree_next = next_col ;

	    
	    
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = c ;
	    }
	    head [score] = c ;

	    
	    min_score = MIN (min_score, score) ;

	}
    }

    

    *p_n_col2 = n_col2 ;
    *p_n_row2 = n_row2 ;
    *p_max_deg = max_deg ;

    
    
    *p_ndense_row = ndense_row ;
    *p_nempty_row = nempty_row ;	
    *p_nnewlyempty_row = nnewlyempty_row ;
    *p_ndense_col = ndense_col ;
    *p_nempty_col = nempty_col ;	
    *p_nnewlyempty_col = nnewlyempty_col ;
    
}








PRIVATE Int find_ordering	
(
    

    Int n_row,			
    Int n_col,			
    Int Alen,			
    Colamd_Row Row [],		
    Colamd_Col Col [],		
    Int A [],			
    Int head [],		
    Int n_col2,			
    Int max_deg,		
    Int pfree			
    
    
    , Int Front_npivcol [ ]
    , Int Front_nrows [ ]
    , Int Front_ncols [ ]
    , Int Front_parent [ ]
    , Int Front_cols [ ]
    , Int *p_nfr		
    , Int aggressive
    , Int InFront [ ]
    
)
{
    

    Int k ;			
    Int pivot_col ;		
    Int *cp ;			
    Int *rp ;			
    Int pivot_row ;		
    Int *new_cp ;		
    Int *new_rp ;		
    Int pivot_row_start ;	
    Int pivot_row_degree ;	
    Int pivot_row_length ;	
    Int pivot_col_score ;	
    Int needed_memory ;		
    Int *cp_end ;		
    Int *rp_end ;		
    Int row ;			
    Int col ;			
    Int max_score ;		
    Int cur_score ;		
    unsigned Int hash ;		
    Int head_column ;		
    Int first_col ;		
    Int tag_mark ;		
    Int row_mark ;		
    Int set_difference ;	
    Int min_score ;		
    Int col_thickness ;		
    Int max_mark ;		
    Int pivot_col_thickness ;	
    Int prev_col ;		
    Int next_col ;		
    Int ngarbage ;		

    
    
    Int pivot_row_thickness ;	
    Int nfr = 0 ;		
    Int child ;
    

    

    max_mark = MAX_MARK (n_col) ;	
    tag_mark = clear_mark (n_row, Row) ;
    min_score = 0 ;
    ngarbage = 0 ;

    for (row = 0 ; row < n_row ; row++)
    {
	InFront [row] = EMPTY ;
    }

    

    for (k = 0 ; k < n_col2 ; )
    {

	

	

	
	while (head [min_score] == EMPTY && min_score < n_col)
	{
	    min_score++ ;
	}
	pivot_col = head [min_score] ;
	next_col = Col [pivot_col].shared4.degree_next ;
	head [min_score] = next_col ;
	if (next_col != EMPTY)
	{
	    Col [next_col].shared3.prev = EMPTY ;
	}

	
	pivot_col_score = Col [pivot_col].shared2.score ;

	
	Col [pivot_col].shared2.order = k ;

	
	pivot_col_thickness = Col [pivot_col].shared1.thickness ;
	
	
	k += pivot_col_thickness ;
	

	

	needed_memory = MIN (pivot_col_score, n_col - k) ;
	if (pfree + needed_memory >= Alen)
	{
	    pfree = garbage_collection (n_row, n_col, Row, Col, A, &A [pfree]) ;
	    ngarbage++ ;
	    
	    ASSERT (pfree + needed_memory < Alen) ;
	    
	    tag_mark = clear_mark (n_row, Row) ;
	}

	

	
	pivot_row_start = pfree ;

	
	pivot_row_degree = 0 ;

	
	
	pivot_row_thickness = 0 ;
	

	
	
	Col [pivot_col].shared1.thickness = -pivot_col_thickness ;

	
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    
	    row = *cp++ ;
	    
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }

	    
	    
	    
	    
	    pivot_row_thickness += Row [row].thickness ;
	    

	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		
		col = *rp++ ;
		
		col_thickness = Col [col].shared1.thickness ;
		if (col_thickness > 0 && COL_IS_ALIVE (col))
		{
		    
		    Col [col].shared1.thickness = -col_thickness ;
		    
		    A [pfree++] = col ;
		    pivot_row_degree += col_thickness ;
		    
		    
		    
		}
		
		
		
	    }
	}

	
	
	
	
	

	
	Col [pivot_col].shared1.thickness = pivot_col_thickness ;	
	max_deg = MAX (max_deg, pivot_row_degree) ;

	

	
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    
	    row = *cp++ ;

	    
	    if (ROW_IS_ALIVE (row))
	    {
		if (Row [row].front != EMPTY)
		{
		    
		    
		    child = Row [row].front ;
		    Front_parent [child] = nfr ;
		}
		else
		{
		    
		    InFront [row] = nfr ;
		}
	    }

	    KILL_ROW (row) ;

	    
	    
	    Row [row].thickness = 0 ;
	    
	}

	

	pivot_row_length = pfree - pivot_row_start ;
	if (pivot_row_length > 0)
	{
	    
	    pivot_row = A [Col [pivot_col].start] ;
	}
	else
	{
	    
	    pivot_row = EMPTY ;
	}

	

	
	
	
	
	
	

	
	
	
	
	
	
	
	
	

	

	

	
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;

	    
	    col_thickness = -Col [col].shared1.thickness ;
	    Col [col].shared1.thickness = col_thickness ;

	    

	    cur_score = Col [col].shared2.score ;
	    prev_col = Col [col].shared3.prev ;
	    next_col = Col [col].shared4.degree_next ;
	    if (prev_col == EMPTY)
	    {
		head [cur_score] = next_col ;
	    }
	    else
	    {
		Col [prev_col].shared4.degree_next = next_col ;
	    }
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = prev_col ;
	    }

	    

	    cp = &A [Col [col].start] ;
	    cp_end = cp + Col [col].length ;
	    while (cp < cp_end)
	    {
		
		row = *cp++ ;
		row_mark = Row [row].shared2.mark ;
		
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    continue ;
		}
		set_difference = row_mark - tag_mark ;
		
		if (set_difference < 0)
		{
		    set_difference = Row [row].shared1.degree ;
		}
		
		set_difference -= col_thickness ;

		
		if (set_difference == 0 && aggressive)
		{
		    
		    if (Row [row].front != EMPTY)
		    {
			
			child = Row [row].front ;
			Front_parent [child] = nfr ;
		    }
		    else
		    {
			
			InFront [row] = nfr ;
		    }

		    KILL_ROW (row) ;

		    
		    
		    pivot_row_thickness += Row [row].thickness ;
		    Row [row].thickness = 0 ;

		}
		else
		{
		    
		    Row [row].shared2.mark = set_difference + tag_mark ;
		}
	    }
	}

	

	
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    
	    col = *rp++ ;
	    hash = 0 ;
	    cur_score = 0 ;
	    cp = &A [Col [col].start] ;
	    
	    new_cp = cp ;
	    cp_end = cp + Col [col].length ;

	    while (cp < cp_end)
	    {
		
		row = *cp++ ;
		row_mark = Row [row].shared2.mark ;
		
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    
		    
		    
		    continue ;
		}
		
		
		
		
		*new_cp++ = row ;
		
		hash += row ;
		
		cur_score += row_mark - tag_mark ;
		
		cur_score = MIN (cur_score, n_col) ;
	    }

	    
	    Col [col].length = (Int) (new_cp - &A [Col [col].start]) ;

	    

	    if (Col [col].length == 0)
	    {
		
		KILL_PRINCIPAL_COL (col) ;
		pivot_row_degree -= Col [col].shared1.thickness ;
		
		Col [col].shared2.order = k ;
		
		k += Col [col].shared1.thickness ;

		
		
		pivot_col_thickness += Col [col].shared1.thickness ;

		
		Col [Col [col].lastcol].nextcol = Front_cols [nfr] ;
		Front_cols [nfr] = col ;
		

	    }
	    else
	    {
		

		
		Col [col].shared2.score = cur_score ;

		
		
		hash %= n_col + 1 ;

		head_column = head [hash] ;
		if (head_column > EMPTY)
		{
		    
		    
		    first_col = Col [head_column].shared3.headhash ;
		    Col [head_column].shared3.headhash = col ;
		}
		else
		{
		    
		    first_col = - (head_column + 2) ;
		    head [hash] = - (col + 2) ;
		}
		Col [col].shared4.hash_next = first_col ;

		
		Col [col].shared3.hash = (Int) hash ;
	    }
	}

	

	

	detect_super_cols (

		Col, A, head, pivot_row_start, pivot_row_length) ;

	

	KILL_PRINCIPAL_COL (pivot_col) ;

	
	
	
	Col [Col [pivot_col].lastcol].nextcol = Front_cols [nfr] ;
	Front_cols [nfr] = pivot_col ;
	

	

	tag_mark += (max_deg + 1) ;
	if (tag_mark >= max_mark)
	{
	    tag_mark = clear_mark (n_row, Row) ;
	}

	

	
	rp = &A [pivot_row_start] ;
	
	new_rp = rp ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    
	    if (COL_IS_DEAD (col))
	    {
		continue ;
	    }
	    *new_rp++ = col ;
	    
	    A [Col [col].start + (Col [col].length++)] = pivot_row ;

	    
	    
	    
	    cur_score = Col [col].shared2.score + pivot_row_degree ;

	    
	    
	    
	    max_score = n_col - k - Col [col].shared1.thickness ;

	    
	    cur_score -= Col [col].shared1.thickness ;

	    
	    cur_score = MIN (cur_score, max_score) ;

	    
	    Col [col].shared2.score = cur_score ;

	    

	    next_col = head [cur_score] ;
	    Col [col].shared4.degree_next = next_col ;
	    Col [col].shared3.prev = EMPTY ;
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = col ;
	    }
	    head [cur_score] = col ;

	    
	    min_score = MIN (min_score, cur_score) ;

	}

	
	
	
	

	
	Front_npivcol [nfr] = pivot_col_thickness ;

	
	Front_nrows [nfr] = pivot_row_thickness ;

	
	Front_ncols [nfr] = pivot_col_thickness + pivot_row_degree ;

	Front_parent [nfr] = EMPTY ;

	pivot_row_thickness -= pivot_col_thickness ;
	pivot_row_thickness = MAX (0, pivot_row_thickness) ;
	

	

	if (pivot_row_degree > 0
	
	
	&& pivot_row_thickness > 0
	
	)
	{
	    
	    
	    Row [pivot_row].start  = pivot_row_start ;
	    Row [pivot_row].length = (Int) (new_rp - &A[pivot_row_start]) ;
	    ASSERT (Row [pivot_row].length > 0) ;
	    Row [pivot_row].shared1.degree = pivot_row_degree ;
	    Row [pivot_row].shared2.mark = 0 ;
	    
	    
	    Row [pivot_row].thickness = pivot_row_thickness ;
	    Row [pivot_row].front = nfr ;
	    
	    
	}

	
	
	nfr++ ; 
	

    }

    

    
    
    *p_nfr = nfr ;
    

    return (ngarbage) ;
}








PRIVATE void detect_super_cols
(
    

    Colamd_Col Col [],		
    Int A [],			
    Int head [],		
    Int row_start,		
    Int row_length		
)
{
    

    Int hash ;			
    Int *rp ;			
    Int c ;			
    Int super_c ;		
    Int *cp1 ;			
    Int *cp2 ;			
    Int length ;		
    Int prev_c ;		
    Int i ;			
    Int *rp_end ;		
    Int col ;			
    Int head_column ;		
    Int first_col ;		

    

    rp = &A [row_start] ;
    rp_end = rp + row_length ;
    while (rp < rp_end)
    {
	col = *rp++ ;
	if (COL_IS_DEAD (col))
	{
	    continue ;
	}

	
	hash = Col [col].shared3.hash ;
	ASSERT (hash <= n_col) ;

	

	head_column = head [hash] ;
	if (head_column > EMPTY)
	{
	    first_col = Col [head_column].shared3.headhash ;
	}
	else
	{
	    first_col = - (head_column + 2) ;
	}

	

	for (super_c = first_col ; super_c != EMPTY ;
	    super_c = Col [super_c].shared4.hash_next)
	{
	    ASSERT (COL_IS_ALIVE (super_c)) ;
	    ASSERT (Col [super_c].shared3.hash == hash) ;
	    length = Col [super_c].length ;

	    
	    prev_c = super_c ;

	    

	    for (c = Col [super_c].shared4.hash_next ;
		c != EMPTY ; c = Col [c].shared4.hash_next)
	    {
		ASSERT (c != super_c) ;
		ASSERT (COL_IS_ALIVE (c)) ;
		ASSERT (Col [c].shared3.hash == hash) ;

		
		if (Col [c].length != length ||
		    Col [c].shared2.score != Col [super_c].shared2.score)
		{
		    prev_c = c ;
		    continue ;
		}

		
		cp1 = &A [Col [super_c].start] ;
		cp2 = &A [Col [c].start] ;

		for (i = 0 ; i < length ; i++)
		{
		    
		    ASSERT (ROW_IS_ALIVE (*cp1))  ;
		    ASSERT (ROW_IS_ALIVE (*cp2))  ;
		    
		    
		    if (*cp1++ != *cp2++)
		    {
			break ;
		    }
		}

		
		if (i != length)
		{
		    prev_c = c ;
		    continue ;
		}

		

		ASSERT (Col [c].shared2.score == Col [super_c].shared2.score) ;

		Col [super_c].shared1.thickness += Col [c].shared1.thickness ;
		Col [c].shared1.parent = super_c ;
		KILL_NON_PRINCIPAL_COL (c) ;

		Col [c].shared2.order = EMPTY ;
		
		Col [prev_c].shared4.hash_next = Col [c].shared4.hash_next ;

		
		
		
		ASSERT (Col [super_c].lastcol >= 0) ;
		ASSERT (Col [super_c].lastcol < n_col) ;
		Col [Col [super_c].lastcol].nextcol = c ;
		Col [super_c].lastcol = Col [c].lastcol ;
		

	    }
	}

	

	if (head_column > EMPTY)
	{
	    
	    Col [head_column].shared3.headhash = EMPTY ;
	}
	else
	{
	    
	    head [hash] = EMPTY ;
	}
    }
}








PRIVATE Int garbage_collection  
(
    

    Int n_row,			
    Int n_col,			
    Colamd_Row Row [],		
    Colamd_Col Col [],		
    Int A [],			
    Int *pfree			
)
{
    

    Int *psrc ;			
    Int *pdest ;		
    Int j ;			
    Int r ;			
    Int c ;			
    Int length ;		

    

    pdest = &A[0] ;
    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    psrc = &A [Col [c].start] ;

	    
	    ASSERT (pdest <= psrc) ;
	    Col [c].start = (Int) (pdest - &A [0]) ;
	    length = Col [c].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		r = *psrc++ ;
		if (ROW_IS_ALIVE (r))
		{
		    *pdest++ = r ;
		}
	    }
	    Col [c].length = (Int) (pdest - &A [Col [c].start]) ;
	}
    }

    

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    if (Row [r].length == 0)
	    {
		
		
		KILL_ROW (r) ;
	    }
	    else
	    {
		
		psrc = &A [Row [r].start] ;
		Row [r].shared2.first_column = *psrc ;
		
		*psrc = ONES_COMPLEMENT (r) ;
	    }
	}
    }

    

    psrc = pdest ;
    while (psrc < pfree)
    {
	
	if (*psrc++ < 0)
	{
	    psrc-- ;
	    
	    r = ONES_COMPLEMENT (*psrc) ;
	    
	    *psrc = Row [r].shared2.first_column ;
	    
	    Row [r].start = (Int) (pdest - &A [0]) ;
	    length = Row [r].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		c = *psrc++ ;
		if (COL_IS_ALIVE (c))
		{
		    *pdest++ = c ;
		}
	    }
	    Row [r].length = (Int) (pdest - &A [Row [r].start]) ;

	}
    }

    

    return ((Int) (pdest - &A [0])) ;
}








PRIVATE Int clear_mark	
(
    

    Int n_row,		
    Colamd_Row Row []	
)
{
    

    Int r ;

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    Row [r].shared2.mark = 0 ;
	}
    }

    
    return (1) ;
    

}

