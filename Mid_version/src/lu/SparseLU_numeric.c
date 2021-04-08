/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_numeric.c
 * BRIEF:   LU数值分解
 *****************************************************************************/



#include "SparseLU_internal.h"
#include "SparseLU_function.h"

PRIVATE Int work_alloc
(
    WorkType *Work,
    SymbolicType *Symbolic
) ;

PRIVATE void free_work
(
    WorkType *Work
) ;

PRIVATE Int numeric_alloc
(
    NumericType **NumericHandle,
    SymbolicType *Symbolic,
    double alloc_init,
    Int scale
) ;

PRIVATE void error
(
    NumericType **Numeric,
    WorkType *Work
) ;







GLOBAL Int SparseLU_numeric
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    void *SymbolicHandle,
    void **NumericHandle,
    const double Control [SparseLU_CONTROL],
    double User_Info [SparseLU_INFO]
)
{

    
    
    

    double Info2 [SparseLU_INFO], alloc_init, relpt, relpt2, droptol,
	front_alloc_init, stats [2] ;
    double *Info ;
    WorkType WorkSpace, *Work ;
    NumericType *Numeric ;
    SymbolicType *Symbolic ;
    Int n_row, n_col, n_inner, newsize, i, status, *inew, npiv, ulen, scale ;
    Unit *mnew ;

    
    
    

    

    relpt = GET_CONTROL (SparseLU_PIVOT_TOLERANCE,
	SparseLU_DEFAULT_PIVOT_TOLERANCE) ;
    relpt2 = GET_CONTROL (SparseLU_SYM_PIVOT_TOLERANCE,
	SparseLU_DEFAULT_SYM_PIVOT_TOLERANCE) ;
    alloc_init = GET_CONTROL (SparseLU_ALLOC_INIT, SparseLU_DEFAULT_ALLOC_INIT) ;
    front_alloc_init = GET_CONTROL (SparseLU_FRONT_ALLOC_INIT,
	SparseLU_DEFAULT_FRONT_ALLOC_INIT) ;
    scale = GET_CONTROL (SparseLU_SCALE, SparseLU_DEFAULT_SCALE) ;
    droptol = GET_CONTROL (SparseLU_DROPTOL, SparseLU_DEFAULT_DROPTOL) ;

    relpt   = MAX (0.0, MIN (relpt,  1.0)) ;
    relpt2  = MAX (0.0, MIN (relpt2, 1.0)) ;
    droptol = MAX (0.0, droptol) ;
    front_alloc_init = MIN (1.0, front_alloc_init) ;

    if (scale != SparseLU_SCALE_NONE && scale != SparseLU_SCALE_MAX)
    {
	    scale = SparseLU_DEFAULT_SCALE ;
    }

    if (User_Info != (double *) NULL)
    {
        
        Info = User_Info ;
        
        for (i = SparseLU_NUMERIC_SIZE ; i <= SparseLU_MAX_FRONT_NCOLS ; i++)
        {
            Info [i] = EMPTY ;
        }
        for (i = SparseLU_NUMERIC_DEFRAG ; i < SparseLU_IR_TAKEN ; i++)
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

    Symbolic = (SymbolicType *) SymbolicHandle ;
    Numeric = (NumericType *) NULL ;
    if (!LU_valid_symbolic (Symbolic))
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_invalid_Symbolic_object ;
        return (SparseLU_ERROR_invalid_Symbolic_object) ;
    }

    
    if ( alloc_init >= 0
        && Symbolic->amd_lunz > 0)
    {
        alloc_init = (Symbolic->nz + Symbolic->amd_lunz) / Symbolic->lunz_bound;
        alloc_init = MIN (1.0, alloc_init) ;
        alloc_init *= LU_REALLOC_INCREASE ;
    }

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n_inner = MIN (n_row, n_col) ;

    
    if (INT_OVERFLOW (Symbolic->dnum_mem_init_usage * sizeof (Unit)))
    {
        
        
        Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
        return (SparseLU_ERROR_out_of_memory) ;
    }

    Info [SparseLU_STATUS] = SparseLU_OK ;
    Info [SparseLU_NROW] = n_row ;
    Info [SparseLU_NCOL] = n_col ;
    Info [SparseLU_SIZE_OF_UNIT] = (double) (sizeof (Unit)) ;

    if (!Ap || !Ai || !Ax || !NumericHandle)
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_argument_missing ;
        return (SparseLU_ERROR_argument_missing) ;
    }

    Info [SparseLU_NZ] = Ap [n_col] ;
    *NumericHandle = (void *) NULL ;

    
    
    

    

    Work = &WorkSpace ;
    Work->n_row = n_row ;
    Work->n_col = n_col ;
    Work->nfr = Symbolic->nfr ;
    Work->nb = Symbolic->nb ;
    Work->n1 = Symbolic->n1 ;

    if (!work_alloc (Work, Symbolic))
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
        error (&Numeric, Work) ;
        return (SparseLU_ERROR_out_of_memory) ;
    }

    
    
    

    

    


    if (!numeric_alloc (&Numeric, Symbolic, alloc_init, scale))
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
        error (&Numeric, Work) ;
        return (SparseLU_ERROR_out_of_memory) ;
    }

    
    Numeric->relpt = relpt ;
    Numeric->relpt2 = relpt2 ;
    Numeric->droptol = droptol ;
    Numeric->alloc_init = alloc_init ;
    Numeric->front_alloc_init = front_alloc_init ;
    Numeric->scale = scale ;

    
    
    

    

    

    status = LU_kernel (Ap, Ai, Ax, Numeric, Work, Symbolic) ;

    Info [SparseLU_STATUS] = status ;
    if (status < SparseLU_OK)
    {
        
        error (&Numeric, Work) ;
        return (status) ;
    }

    Info [SparseLU_FORCED_UPDATES] = Work->nforced ;
    Info [SparseLU_VARIABLE_INIT] = Numeric->init_usage ;
    if (Symbolic->prefer_diagonal)
    {
	    Info [SparseLU_NOFF_DIAG] = Work->noff_diagonal ;
    }

    npiv = Numeric->npiv ;	
    ulen = Numeric->ulen ;	

    
    
    

    

    free_work (Work) ;

    
    
    

    

    if (npiv < n_row)
    {
        
        inew = (Int *) LU_realloc (Numeric->Lpos, npiv+1, sizeof (Int)) ;
        if (inew)
        {
            Numeric->Lpos = inew ;
        }
        inew = (Int *) LU_realloc (Numeric->Uilen, npiv+1, sizeof (Int)) ;
        if (inew)
        {
            Numeric->Uilen = inew ;
        }
        inew = (Int *) LU_realloc (Numeric->Uip, npiv+1, sizeof (Int)) ;
        if (inew)
        {
            Numeric->Uip = inew ;
        }
    }

    if (npiv < n_col)
    {
        
        inew = (Int *) LU_realloc (Numeric->Upos, npiv+1, sizeof (Int)) ;
        if (inew)
        {
            Numeric->Upos = inew ;
        }
        inew = (Int *) LU_realloc (Numeric->Lilen, npiv+1, sizeof (Int)) ;
        if (inew)
        {
            Numeric->Lilen = inew ;
        }
        inew = (Int *) LU_realloc (Numeric->Lip, npiv+1, sizeof (Int)) ;
        if (inew)
        {
            Numeric->Lip = inew ;
        }
    }

    
    
    

    

    if (ulen > 0 && ulen < n_col)
    {
        inew = (Int *) LU_realloc (Numeric->Upattern, ulen+1, sizeof (Int)) ;
        if (inew)
        {
            Numeric->Upattern = inew ;
        }
    }

    
    
    

    

    newsize = Numeric->ihead ;
    if (newsize < Numeric->size)
    {
        mnew = (Unit *) LU_realloc (Numeric->Memory, newsize, sizeof (Unit)) ;
        if (mnew)
        {
            
            Numeric->Memory = mnew ;
            Numeric->size = newsize ;
        }
    }
    Numeric->ihead = Numeric->size ;
    Numeric->itail = Numeric->ihead ;
    Numeric->tail_usage = 0 ;
    Numeric->ibig = EMPTY ;
    

    
    
    

    LU_set_stats (
	Info,
	Symbolic,
	(double) Numeric->max_usage,	
	(double) Numeric->size,		
	Numeric->flops,			
	(double) Numeric->lnz + n_inner,		
	(double) Numeric->unz + Numeric->nnzpiv,	
	(double) Numeric->maxfrsize,	
	(double) ulen,			
	(double) npiv,			
	(double) Numeric->maxnrows,	
	(double) Numeric->maxncols,	
	scale != SparseLU_SCALE_NONE,
	Symbolic->prefer_diagonal,
	ACTUAL) ;

    Info [SparseLU_ALLOC_INIT_USED] = Numeric->alloc_init ;
    Info [SparseLU_NUMERIC_DEFRAG] = Numeric->ngarbage ;
    Info [SparseLU_NUMERIC_REALLOC] = Numeric->nrealloc ;
    Info [SparseLU_NUMERIC_COSTLY_REALLOC] = Numeric->ncostly ;
    Info [SparseLU_COMPRESSED_PATTERN] = Numeric->isize ;
    Info [SparseLU_LU_ENTRIES] = Numeric->nLentries + Numeric->nUentries +
	    Numeric->npiv ;
    Info [SparseLU_UDIAG_NZ] = Numeric->nnzpiv ;
    Info [SparseLU_RSMIN] = Numeric->rsmin ;
    Info [SparseLU_RSMAX] = Numeric->rsmax ;
    Info [SparseLU_WAS_SCALED] = Numeric->scale ;

    
    Info [SparseLU_ALL_LNZ] = Numeric->all_lnz + n_inner ;
    Info [SparseLU_ALL_UNZ] = Numeric->all_unz + Numeric->nnzpiv ;
    Info [SparseLU_NZDROPPED] = (Numeric->all_lnz - Numeric->lnz)
	+ (Numeric->all_unz - Numeric->unz) ;

    
    if (SCALAR_IS_ZERO (Numeric->min_udiag)
     || SCALAR_IS_ZERO (Numeric->max_udiag)
     ||	SCALAR_IS_NAN (Numeric->min_udiag)
     ||	SCALAR_IS_NAN (Numeric->max_udiag))
    {
        
        Numeric->rcond = 0.0 ;
    }
    else
    {
        
        
        Numeric->rcond = Numeric->min_udiag / Numeric->max_udiag ;
    }
    Info [SparseLU_UMIN]  = Numeric->min_udiag ;
    Info [SparseLU_UMAX]  = Numeric->max_udiag ;
    Info [SparseLU_RCOND] = Numeric->rcond ;

    if (Numeric->nnzpiv < n_inner
    || SCALAR_IS_ZERO (Numeric->rcond) || SCALAR_IS_NAN (Numeric->rcond))
    {
        
        status = SparseLU_WARNING_singular_matrix ;
        Info [SparseLU_STATUS] = status ;
    }

    Numeric->valid = NUMERIC_VALID ;
    *NumericHandle = (void *) Numeric ;

    
    return (status) ;

}








PRIVATE Int numeric_alloc
(
    NumericType **NumericHandle,
    SymbolicType *Symbolic,
    double alloc_init,
    Int scale
)
{
    double nsize, bsize ;
    Int n_row, n_col, n_inner, min_usage, trying ;
    NumericType *Numeric ;

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n_inner = MIN (n_row, n_col) ;
    *NumericHandle = (NumericType *) NULL ;

    
    Numeric = (NumericType *) LU_malloc (1, sizeof (NumericType)) ;

    if (!Numeric)
    {
	return (FALSE) ;	
    }
    Numeric->valid = 0 ;
    *NumericHandle = Numeric ;

    
    Numeric->D = (Entry *) LU_malloc (n_inner+1, sizeof (Entry)) ;
    Numeric->Rperm = (Int *) LU_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Cperm = (Int *) LU_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Lpos = (Int *) LU_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Lilen = (Int *) LU_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Lip = (Int *) LU_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Upos = (Int *) LU_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Uilen = (Int *) LU_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Uip = (Int *) LU_malloc (n_row+1, sizeof (Int)) ;

    
    if (scale != SparseLU_SCALE_NONE)
    {
	Numeric->Rs = (double *) LU_malloc (n_row, sizeof (double)) ;
    }
    else
    {
	Numeric->Rs = (double *) NULL ;
    }

    Numeric->Memory = (Unit *) NULL ;

    
    Numeric->Upattern = (Int *) NULL ;	

    if (!Numeric->D || !Numeric->Rperm || !Numeric->Cperm || !Numeric->Upos ||
	!Numeric->Lpos || !Numeric->Lilen || !Numeric->Uilen || !Numeric->Lip ||
	!Numeric->Uip || (scale != SparseLU_SCALE_NONE && !Numeric->Rs))
    {
	return (FALSE) ;	
    }

    
    
    

    if (alloc_init < 0)
    {
	
	nsize = -alloc_init ;
    }
    else
    {
	
	nsize = (alloc_init * Symbolic->num_mem_usage_est) + 1 ;
    }
    min_usage = Symbolic->num_mem_init_usage ;

    
    nsize = MAX (min_usage, nsize) ;

    
    
    bsize = ((double) Int_MAX) / sizeof (Unit) - 1 ;
    nsize = MIN (nsize, bsize) ;

    Numeric->size = (Int) nsize ;

    
    
    trying = TRUE ;
    while (trying)
    {
	Numeric->Memory = (Unit *) LU_malloc (Numeric->size, sizeof (Unit)) ;
	if (Numeric->Memory)
	{
	    return (TRUE) ;
	}
	
	
	trying = Numeric->size > min_usage ;
	Numeric->size = (Int)
	    (LU_REALLOC_REDUCTION * ((double) Numeric->size)) ;
	Numeric->size = MAX (min_usage, Numeric->size) ;
    }

    return (FALSE) ;	
}








PRIVATE Int work_alloc
(
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    Int n_row, n_col, nn, maxnrows, maxncols, nfr, ok, maxnrc, n1 ;

    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nn = MAX (n_row, n_col) ;
    nfr = Work->nfr ;
    n1 = Symbolic->n1 ;

    maxnrows = Symbolic->maxnrows + Symbolic->nb ;
    maxnrows = MIN (n_row, maxnrows) ;
    maxncols = Symbolic->maxncols + Symbolic->nb ;
    maxncols = MIN (n_col, maxncols) ;
    maxnrc = MAX (maxnrows, maxncols) ;

    
    
    Work->Wx = (Entry *) LU_malloc (maxnrows + 1, sizeof (Entry)) ;
    Work->Wy = (Entry *) LU_malloc (maxnrows + 1, sizeof (Entry)) ;
    Work->Frpos    = (Int *) LU_malloc (n_row + 1, sizeof (Int)) ;
    Work->Lpattern = (Int *) LU_malloc (n_row + 1, sizeof (Int)) ;
    Work->Fcpos = (Int *) LU_malloc (n_col + 1, sizeof (Int)) ;
    Work->Wp = (Int *) LU_malloc (nn + 1, sizeof (Int)) ;
    Work->Wrp = (Int *) LU_malloc (MAX (n_col,maxnrows) + 1, sizeof (Int)) ;
    Work->Frows = (Int *) LU_malloc (maxnrows + 1, sizeof (Int)) ;
    Work->Wm    = (Int *) LU_malloc (maxnrows + 1, sizeof (Int)) ;
    Work->Fcols = (Int *) LU_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Wio   = (Int *) LU_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Woi   = (Int *) LU_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Woo = (Int *) LU_malloc (maxnrc + 1, sizeof (Int));
    Work->elen = (n_col - n1) + (n_row - n1) + MIN (n_col-n1, n_row-n1) + 1 ;
    Work->E = (Int *) LU_malloc (Work->elen, sizeof (Int)) ;
    Work->Front_new1strow = (Int *) LU_malloc (nfr + 1, sizeof (Int)) ;

    ok = (Work->Frpos && Work->Fcpos && Work->Lpattern
	&& Work->Wp && Work->Wrp && Work->Frows && Work->Fcols
	&& Work->Wio && Work->Woi && Work->Woo && Work->Wm
	&& Work->E && Work->Front_new1strow && Work->Wx && Work->Wy) ;

    
    if (Symbolic->prefer_diagonal)
    {
	Work->Diagonal_map  = (Int *) LU_malloc (nn, sizeof (Int)) ;
	Work->Diagonal_imap = (Int *) LU_malloc (nn, sizeof (Int)) ;
	ok = ok && Work->Diagonal_map && Work->Diagonal_imap ;
    }
    else
    {
	
	Work->Diagonal_map = (Int *) NULL ;
	Work->Diagonal_imap = (Int *) NULL ;
    }

    
    Work->Upattern = (Int *) LU_malloc (n_col + 1, sizeof (Int)) ;
    ok = ok && Work->Upattern ;

    
    Work->Flublock = (Entry *) NULL ;
    Work->Flblock  = (Entry *) NULL ;
    Work->Fublock  = (Entry *) NULL ;
    Work->Fcblock  = (Entry *) NULL ;

    return (ok) ;
}






PRIVATE void free_work
(
    WorkType *Work
)
{
    if (Work)
    {
	
	Work->Wx = (Entry *) LU_free ((void *) Work->Wx) ;
	Work->Wy = (Entry *) LU_free ((void *) Work->Wy) ;
	Work->Frpos = (Int *) LU_free ((void *) Work->Frpos) ;
	Work->Fcpos = (Int *) LU_free ((void *) Work->Fcpos) ;
	Work->Lpattern = (Int *) LU_free ((void *) Work->Lpattern) ;
	Work->Upattern = (Int *) LU_free ((void *) Work->Upattern) ;
	Work->Wp = (Int *) LU_free ((void *) Work->Wp) ;
	Work->Wrp = (Int *) LU_free ((void *) Work->Wrp) ;
	Work->Frows = (Int *) LU_free ((void *) Work->Frows) ;
	Work->Fcols = (Int *) LU_free ((void *) Work->Fcols) ;
	Work->Wio = (Int *) LU_free ((void *) Work->Wio) ;
	Work->Woi = (Int *) LU_free ((void *) Work->Woi) ;
	Work->Woo = (Int *) LU_free ((void *) Work->Woo) ;
	Work->Wm = (Int *) LU_free ((void *) Work->Wm) ;
	Work->E = (Int *) LU_free ((void *) Work->E) ;
	Work->Front_new1strow =
	    (Int *) LU_free ((void *) Work->Front_new1strow) ;

	
	Work->Diagonal_map = (Int *) LU_free ((void *) Work->Diagonal_map) ;
	Work->Diagonal_imap = (Int *) LU_free ((void *) Work->Diagonal_imap) ;
    }
}








PRIVATE void error
(
    NumericType **Numeric,
    WorkType *Work
)
{
    free_work (Work) ;
    SparseLU_free_numeric ((void **) Numeric) ;
    ASSERT (LU_malloc_count == init_count) ;
}
