/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    SparseLU_kernel.c
 * BRIEF:   稀疏LU核函数
 *****************************************************************************/

#include "SparseLU_internal.h"
#include "SparseLU_function.h"


#define DO(action) { if (! (action)) { return (SparseLU_ERROR_out_of_memory) ; }}

GLOBAL Int LU_kernel
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{

    
    
    

    Int j, f1, f2, chain, nchains, *Chain_start, status, fixQ, evaporate,
	*Front_npivcol, jmax, nb, drop ;

    
    
    

    if (!LU_kernel_init (Ap, Ai, Ax, Numeric, Work, Symbolic))
    {
		
		
		
		
		
		return (SparseLU_ERROR_different_pattern) ;
    }

    
    
    

    nchains = Symbolic->nchains ;
    Chain_start = Symbolic->Chain_start ;
    Front_npivcol = Symbolic->Front_npivcol ;
    nb = Symbolic->nb ;
    fixQ = Symbolic->fixQ ;
    // drop = Numeric->droptol > 0.0 ;
	drop = 0;

    
    
    
    for (chain = 0 ; chain < nchains ; ++chain)
    {
		f1 = Chain_start [chain] ;
		f2 = Chain_start [chain+1] - 1 ;

		
		
		
		DO (LU_start_front (chain, Numeric, Work, Symbolic)) ;

		
		
		
		for (Work->frontid = f1 ; Work->frontid <= f2 ; Work->frontid++)
		{

			
			
			

			Work->ncand = Front_npivcol [Work->frontid] ;
			Work->lo = Work->nextcand ;
			Work->hi = Work->nextcand + Work->ncand - 1 ;
			jmax = MIN (MAX_CANDIDATES, Work->ncand) ;
			if (fixQ)
			{
				
				jmax = 1 ;
			}
			for (j = 0 ; j < jmax ; j++)
			{
				Work->Candidates [j] = Work->nextcand++ ;
			}
			Work->nCandidates = jmax ;

			
			
			

			while (Work->ncand > 0)
			{
				
				
				
				status = LU_local_search (Numeric, Work, Symbolic) ;
				if (status == SparseLU_ERROR_different_pattern)
				{
					
					
					return (SparseLU_ERROR_different_pattern) ;
				}
				if (status == SparseLU_WARNING_singular_matrix)
				{
					
					continue ;
				}

				
				
				

				if (Work->do_update)
				{
					LU_blas3_update (Work) ;
					if (drop)
					{
						DO (LU_store_lu_drop (Numeric, Work)) ;
					}
					else
					{
						DO (LU_store_lu (Numeric, Work)) ;
					}
				}

				
				
				

				if (Work->do_extend)
				{
					
					DO (LU_extend_front (Numeric, Work)) ;
				}
				else
				{
					
					DO (LU_create_element (Numeric, Work, Symbolic)) ;
					DO (LU_init_front (Numeric, Work)) ;
				}

				
				
				

				if (fixQ)
				{
					LU_assemble_fixq (Numeric, Work) ;
				}
				else
				{
					LU_assemble (Numeric, Work) ;
				}

				
				
				

				LU_scale_column (Numeric, Work) ;

				
				
				

				evaporate = Work->fnrows == 0 || Work->fncols == 0 ;
				if (Work->fnpiv >= nb || evaporate)
				{
					LU_blas3_update (Work) ;
					if (drop)
					{
						DO (LU_store_lu_drop (Numeric, Work)) ;
					}
					else
					{
						DO (LU_store_lu (Numeric, Work)) ;
					}

				}

				Work->pivrow_in_front = FALSE ;
				Work->pivcol_in_front = FALSE ;

				
				
				

				if (evaporate)
				{
					
					(void) LU_create_element (Numeric, Work, Symbolic) ;
					Work->fnrows = 0 ;
					Work->fncols = 0 ;
				}
			}
		}

		

		LU_blas3_update (Work) ;
		if (drop)
		{
			DO (LU_store_lu_drop (Numeric, Work)) ;
		}
		else
		{
			DO (LU_store_lu (Numeric, Work)) ;
		}
		Work->fnrows_new = Work->fnrows ;
		Work->fncols_new = Work->fncols ;
		DO (LU_create_element (Numeric, Work, Symbolic)) ;

		
		
		

		Work->fnrows = 0 ;
		Work->fncols = 0 ;
    }

    
    
    

    LU_kernel_wrapup (Numeric, Symbolic, Work) ;

    
    return (SparseLU_OK) ;
}









PRIVATE Int packsp	
(
    Int pnew,		
    Int *p_p,		
    Int *p_len,		
    Int drop,		
    double droptol,	
    Unit *Memory	
)
{
    Entry x ;
    double s ;
    Entry *Bx, *Bx2 ;
    Int p, i, len, len_new, *Bi, *Bi2 ;

    
    p = *p_p ;
    len = *p_len ;
    Bi = (Int   *) (Memory + p) ; p += UNITS (Int,   len) ;
    Bx = (Entry *) (Memory + p) ; p += UNITS (Entry, len) ;

    

    
    len_new = 0 ;
    for (p = 0 ; p < len ; p++)
    {
	i = Bi [p] ;
	x = Bx [p] ;
	
	if (IS_ZERO (x)) continue ;
	if (drop)
	{
	    APPROX_ABS (s, x) ;
	    if (s <= droptol) continue ;
	}
	
	if (len_new != p)
	{
	    Bi [len_new] = i ;
	    Bx [len_new] = x ;
	}
	len_new++ ;
    }

    

    
    *p_p = pnew ;
    *p_len = len_new ;
    Bi2 = (Int   *) (Memory + pnew) ; pnew += UNITS (Int,   len_new) ;
    Bx2 = (Entry *) (Memory + pnew) ; pnew += UNITS (Entry, len_new) ;

    
    for (p = 0 ; p < len_new ; p++)
    {
	Bi2 [p] = Bi [p] ;
    }
    for (p = 0 ; p < len_new ; p++)
    {
	Bx2 [p] = Bx [p] ;
    }

    
    return (pnew) ;
}






GLOBAL Int LU_kernel_init
(
    const Int Ap [ ],		
    const Int Ai [ ],
    const double Ax [ ],
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    
    
    

    Entry x, pivot_value ;
    double unused = 0, rsmin, rsmax, rs, droptol ;
    Entry *D, *C, *Lval, **Rpx ;
    double *Rs ;
    Int row, k, oldcol, size, e, p1, p2, p, nz, *Rows, *Cols, *E, i, *Upos,
	*Lpos, n_row, n_col, *Wp, *Cperm_init, *Frpos, *Fcpos, *Row_degree, nn,
	*Row_tlen, *Col_degree, *Col_tlen, oldrow, newrow, ilast, *Wrp,
	*Rperm_init, col, n_inner, prefer_diagonal, *Diagonal_map, nempty,
	*Diagonal_imap, fixQ, rdeg, cdeg, nempty_col, *Esize, esize, pnew,
	*Lip, *Uip, *Lilen, *Uilen, llen, pa, *Cdeg, *Rdeg, n1, clen, do_scale,
	lnz, unz, lip, uip, k1, *Rperm, *Cperm, pivcol, *Li, lilen, drop,
	**Rpi, nempty_row, dense_row_threshold, empty_elements, rpi, rpx ;
    Element *ep ;
    Unit *Memory ;
#ifndef NRECIPROCAL
    Int do_recip = FALSE ;
#endif

    
    
    

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    nempty_col = Symbolic->nempty_col ;
    nempty_row = Symbolic->nempty_row ;
    nempty = MIN (nempty_row, nempty_col) ;
    Cperm_init = Symbolic->Cperm_init ;
    Rperm_init = Symbolic->Rperm_init ;
    Cdeg = Symbolic->Cdeg ;
    Rdeg = Symbolic->Rdeg ;
    n1 = Symbolic->n1 ;
    dense_row_threshold = Symbolic->dense_row_threshold ;
    Work->nforced = 0 ;
    Work->ndiscard = 0 ;
    Work->noff_diagonal = 0 ;

    nz = Ap [n_col] ;
    if (nz < 0 || Ap [0] != 0 || nz != Symbolic->nz)
    {
	return (FALSE) ;	
    }

    prefer_diagonal = Symbolic->prefer_diagonal ;
    Diagonal_map = Work->Diagonal_map ;
    Diagonal_imap = Work->Diagonal_imap ;

    
    
    

    LU_mem_init_memoryspace (Numeric) ;

    
    
    

    
    Work->fnpiv = 0 ;
    Work->fncols = 0 ;
    Work->fnrows = 0 ;
    Work->fncols_max = 0 ;
    Work->fnrows_max = 0 ;
    Work->fnzeros = 0 ;
    Work->fcurr_size = 0 ;
    Work->fnr_curr = 0 ;
    Work->fnc_curr = 0 ;

    Work->nz = nz ;
    Work->prior_element = EMPTY ;
    Work->ulen = 0 ;
    Work->llen = 0 ;
    Work->npiv = n1 ;
    Work->frontid = 0 ;
    Work->nextcand = n1 ;

    Memory = Numeric->Memory ;
    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;
    Row_degree = Numeric->Rperm ;
    Col_degree = Numeric->Cperm ;
    
    Row_tlen   = Numeric->Uilen ;
    
    Col_tlen   = Numeric->Lilen ;

    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    Lilen = Numeric->Lilen ;
    Uilen = Numeric->Uilen ;

    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;
    Wp = Work->Wp ;
    Wrp = Work->Wrp ;

    D = Numeric->D ;
    Upos = Numeric->Upos ;
    Lpos = Numeric->Lpos ;
    for (k = 0 ; k < n_inner ; k++)
    {
	CLEAR (D [k]) ;
    }

    Rs = Numeric->Rs ;

    for (row = 0 ; row <= n_row ; row++)
    {
	Lpos [row] = EMPTY ;
	
	
	Row_tlen [row] = 0 ;
	
    }

    for (col = 0 ; col <= n_col ; col++)
    {
	Upos [col] = EMPTY ;
	
	
	Col_tlen [col] = 0 ;
	Fcpos [col] = EMPTY ;
	Wrp [col] = 0 ;
    }
    Work->Wrpflag = 1 ;

    
    for (i = 0 ; i <= nn ; i++)
    {
	Wp [i] = EMPTY ;
    }
    

    
    

    

    
    Work->cdeg0 = 1 ;
    Work->rdeg0 = 1 ;

    fixQ = Symbolic->fixQ ;

    E = Work->E ;

    Numeric->n_row = n_row ;
    Numeric->n_col = n_col ;
    Numeric->npiv = 0 ;
    Numeric->nnzpiv = 0 ;
    Numeric->min_udiag = 0.0 ;
    Numeric->max_udiag = 0.0 ;
    Numeric->rcond = 0.0 ;
    Numeric->isize = 0 ;
    Numeric->nLentries = 0 ;
    Numeric->nUentries = 0 ;
    Numeric->lnz = 0 ;
    Numeric->unz = 0 ;
    Numeric->all_lnz = 0 ;
    Numeric->all_unz = 0 ;
    Numeric->maxfrsize = 0 ;
    Numeric->maxnrows = 0 ;
    Numeric->maxncols = 0 ;
    Numeric->flops = 0. ;
    Numeric->n1 = n1 ;
    droptol = Numeric->droptol ;
    drop = (droptol > 0) ;

    
    
    

    

    do_scale = (Numeric->scale != SparseLU_SCALE_NONE) ;

    if (do_scale)
    {
	int do_max = Numeric->scale == SparseLU_SCALE_MAX ;
	for (row = 0 ; row < n_row ; row++)
	{
	    Rs [row] = 0.0 ;
	}
	for (col = 0 ; col < n_col ; col++)
	{
	    ilast = EMPTY ;
	    p1 = Ap [col] ;
	    p2 = Ap [col+1] ;
	    if (p1 > p2)
	    {
		
		return (FALSE) ;
	    }
	    for (p = p1 ; p < p2 ; p++)
	    {
		Entry aij ;
		double value ;
		row = Ai [p] ;
		if (row <= ilast || row >= n_row)
		{
		    
		    return (FALSE) ;
		}
		ASSIGN (aij, Ax, Az, p, split) ;
		APPROX_ABS (value, aij) ;
		rs = Rs [row] ;
		if (!SCALAR_IS_NAN (rs))
		{
		    if (SCALAR_IS_NAN (value))
		    {
			
			Rs [row] = value ;
		    }
		    else if (do_max)
		    {
			Rs [row] = MAX (rs, value) ;
		    }
		    else
		    {
			Rs [row] += value ;
		    }
		}
		ilast = row ;
	    }
	}
	for (row = 0 ; row < n_row ; row++)
	{
	    rs = Rs [row] ;
	    if (SCALAR_IS_ZERO (rs) || SCALAR_IS_NAN (rs))
	    {
		
		Rs [row] = 1.0 ;
	    }
	}
	rsmin = Rs [0] ;
	rsmax = Rs [0] ;
	for (row = 0 ; row < n_row ; row++)
	{
	    rsmin = MIN (rsmin, Rs [row]) ;
	    rsmax = MAX (rsmax, Rs [row]) ;
	}
#ifndef NRECIPROCAL
	
	do_recip = (rsmin >= RECIPROCAL_TOLERANCE) ;
	if (do_recip)
	{
	    
	    for (row = 0 ; row < n_row ; row++)
	    {
		Rs [row] = 1.0 / Rs [row] ;
	    }
	}
#endif
    }
    else
    {
	
	rsmin = -1 ;
	rsmax = -1 ;
#ifndef NRECIPROCAL
	do_recip = FALSE ;
#endif
	
	if (AMD_valid (n_row, n_col, Ap, Ai) != AMD_OK)
	{
	    
	    return (FALSE) ;
	}
    }

    Numeric->rsmin = rsmin ;
    Numeric->rsmax = rsmax ;
#ifndef NRECIPROCAL
    Numeric->do_recip = do_recip ;
#else
    Numeric->do_recip = FALSE ;
#endif

    
    
    

    for (newrow = 0 ; newrow < n_row ; newrow++)
    {
	oldrow = Rperm_init [newrow] ;
	Frpos [oldrow] = newrow ;
    }

    
    
    

    if (prefer_diagonal)
    {

	for (k = 0 ; k < nn ; k++)
	{
	    newrow = Symbolic->Diagonal_map [k] ;
	    Diagonal_map [k] = newrow ;
	    Diagonal_imap [newrow] = k ;
	}
    }

    
    
    

    rpi = LU_mem_alloc_tail_block (Numeric, UNITS (Int *, n_row+1)) ;
    rpx = LU_mem_alloc_tail_block (Numeric, UNITS (Entry *, n_row+1)) ;
    if (!rpi || !rpx)
    {
	
	
	return (FALSE) ;	
    }
    Rpi = (Int   **) (Memory + rpx) ;
    Rpx = (Entry **) (Memory + rpi) ;

    
    
    

    for (k = 0 ; k < n1 ; k++)
    {
	lnz = Cdeg [k] - 1 ;
	unz = Rdeg [k] - 1 ;

	size = UNITS (Int, lnz) + UNITS (Entry, lnz)
	     + UNITS (Int, unz) + UNITS (Entry, unz) ;
	p = LU_mem_alloc_head_block (Numeric, size) ;
	if (!p)
	{
	    
	    return (FALSE) ;	
	}

	Numeric->all_lnz += lnz ;
	Numeric->all_unz += unz ;

	
	lip = p ;
	p += UNITS (Int, lnz) ;
	p += UNITS (Entry, lnz) ;

	
	uip = p ;
	Rpi [k] = (Int *) (Memory + p) ;
	p += UNITS (Int, unz) ;
	Rpx [k] = (Entry *) (Memory + p) ;
	

	
	Lip [k] = lip ;
	Lilen [k] = lnz ;

	
	Uip [k] = uip ;
	Uilen [k] = unz ;

	Wp [k] = unz ;

	
	k1 = ONES_COMPLEMENT (k) ;
	Rperm [k] = k1 ;			
	Cperm [k] = k1 ;			
    }

    
    
    

    e = 0 ;
    E [e] = 0 ;
    Work->Flublock = (Entry *) NULL ;
    Work->Flblock  = (Entry *) NULL ;
    Work->Fublock  = (Entry *) NULL ;
    Work->Fcblock  = (Entry *) NULL ;

    
    
    

    Esize = Symbolic->Esize ;
    empty_elements = FALSE  ;
    for (k = n1 ; k < n_col - nempty_col ; k++)
    {
	e = k - n1 + 1 ;
	ASSERT (e < Work->elen) ;
	esize = Esize ? Esize [k-n1] : Cdeg [k] ;
	if (esize > 0)
	{
	    
	    E [e] = LU_mem_alloc_element (Numeric, esize, 1, &Rows, &Cols, &C,
		&size, &ep) ;
	    if (E [e] <= 0)
	    {
		
		return (FALSE) ;	
	    }
	    Cols [0] = k ;
	}
	else
	{
	    
	    E [e] = 0 ;
	    empty_elements = TRUE  ;
	}
    }

    
    
    

    if (Esize)
    {
	for (k = n1 ; k < n_row - nempty_row ; k++)
	{
	    rdeg = Rdeg [k] ;
	    if (rdeg > dense_row_threshold)
	    {
		
		e++ ;
		E [e] = LU_mem_alloc_element (Numeric, 1, rdeg, &Rows, &Cols,
		    &C, &size, &ep) ;
		if (E [e] <= 0)
		{
		    
		    return (FALSE) ;	
		}
		Rows [0] = k ;
		Rpi [k] = Cols ;
		Rpx [k] = C ;
		Wp [k] = rdeg ;
	    }
	}
    }

    
    Work->nel = e ;

    
    
    

    for (k = 0 ; k < n1 ; k++)
    {
	pivcol = Cperm_init [k] ;
	p2 = Ap [pivcol+1] ;

	
	p = Lip [k] ;
	Li = (Int *) (Memory + p) ;
	lilen = Lilen [k] ;
	p += UNITS (Int, lilen) ;
	Lval = (Entry *) (Memory + p) ;

	llen = 0 ;

	for (pa = Ap [pivcol] ; pa < p2 ; pa++)
	{
	    oldrow = Ai [pa] ;
	    newrow = Frpos [oldrow] ;
	    ASSIGN (x, Ax, Az, pa, split) ;

	    
	    if (do_scale)
	    {
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    SCALE (x, Rs [oldrow]) ;
		}
		else
#endif
		{
		    SCALE_DIV (x, Rs [oldrow]) ;
		}
	    }

	    if (newrow == k)
	    {
		
		D [k] = x ;
	    }
	    else if (newrow < k)
	    {
		
		if (--(Wp [newrow]) < 0)
		{
		    
		    return (FALSE) ;	
		}
		*(Rpi [newrow]++) = k ;
		*(Rpx [newrow]++) = x ;
	    }
	    else
	    {
		
		if (llen >= lilen)
		{
		    return (FALSE) ;	
		}
		Li   [llen] = newrow ;
		Lval [llen] = x ;
		llen++ ;
	    }
	}

	if (llen != lilen)
	{
	    
	    return (FALSE) ;	
	}

	
	if (llen > 0)
	{
	    pivot_value = D [k] ;
	    LU_scale (llen, pivot_value, Lval) ;
	}

    }

    
    
    

    
    for (k = n1 ; k < n_col ; k++)
    {
	
	oldcol = Cperm_init [k] ;

	ASSERT (oldcol >= 0 && oldcol < n_col) ;

	p2 = Ap [oldcol+1] ;

	cdeg = Cdeg [k] ;
	ASSERT (cdeg >= 0) ;

	
	Col_degree [k] = fixQ ? 0 : cdeg ;

	
	e = k - n1 + 1 ;
	if (k < n_col - nempty_col)
	{
	    esize = Esize ? Esize [k-n1] : cdeg ;
	    if (E [e])
	    {
		Int ncols, nrows ;
		Unit *pp ;
		pp = Memory + E [e] ;
		GET_ELEMENT (ep, pp, Cols, Rows, ncols, nrows, C) ;
		ASSERT (ncols == 1) ;
		ASSERT (nrows == esize) ;
		ASSERT (Cols [0] == k) ;
	    }
	}
	else
	{
	    ASSERT (cdeg == 0) ;
	    esize = 0 ;
	}

	clen = 0 ;

	for (pa = Ap [oldcol] ; pa < p2 ; pa++)
	{
	    oldrow = Ai [pa] ;
	    newrow = Frpos [oldrow] ;
	    ASSIGN (x, Ax, Az, pa, split) ;

	    
	    if (do_scale)
	    {
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    
		    SCALE (x, Rs [oldrow]) ;
		}
		else
#endif
		{
		    
		    SCALE_DIV (x, Rs [oldrow]) ;
		}
	    }

	    rdeg = Rdeg [newrow] ;
	    if (newrow < n1 || rdeg > dense_row_threshold)
	    {
		
		if (--(Wp [newrow]) < 0)
		{
		    return (FALSE) ;	
		}
		*(Rpi [newrow]++) = k ;
		*(Rpx [newrow]++) = x ;
	    }
	    else
	    {
		
		if (clen >= esize)
		{
		    return (FALSE) ;	
		}
		Rows [clen] = newrow ;
		C    [clen] = x ;
		clen++ ;
	    }
	}

	if (clen != esize)
	{
	    
	    return (FALSE) ;	
	}
    }

    
    
    

    LU_mem_free_tail_block (Numeric, rpi) ;
    LU_mem_free_tail_block (Numeric, rpx) ;

    
    
    

    if (n1 > 0)
    {
	pnew = Lip [0] ;
	for (k = 0 ; k < n1 ; k++)
	{
	    pnew = packsp (pnew, &Lip [k], &Lilen [k], drop, droptol, Memory) ;
	    Numeric->lnz += Lilen [k] ;
	    pnew = packsp (pnew, &Uip [k], &Uilen [k], drop, droptol, Memory) ;
	    Numeric->unz += Uilen [k] ;
	}
	
	Numeric->ihead = pnew ;
    }

    
    
    

    for (k = 0 ; k < n1 ; k++)
    {
	if (Wp [k] != 0)
	{
	    
	    return (FALSE) ;	
	}
    }

    for (k = n1 ; k < n_row ; k++)
    {
	rdeg = Rdeg [k] ;
	Row_degree [k] = rdeg ;
	if (rdeg > dense_row_threshold && Wp [k] != 0)
	{
	    
	    return (FALSE) ;	
	}
    }

    Col_degree [n_col] = 0 ;

    
    
    

    if (empty_elements)
    {
	Int e2 = 0 ;
	for (e = 1 ; e <= Work->nel ; e++)
	{
	    if (E [e])
	    {
		e2++ ;
		E [e2] = E [e] ;
	    }
	}
	Work->nel = e2 ;
    }

    for (e = Work->nel + 1 ; e < Work->elen ; e++)
    {
	E [e] = 0 ;
    }

    
    for (row = 0 ; row <= n_row ; row++)
    {
	Frpos [row] = EMPTY ;
    }

    
    for (i = 0 ; i <= nn ; i++)
    {
	Wp [i] = EMPTY ;
    }

    
    
    

    

    (void) LU_tuple_lengths (Numeric, Work, &unused) ;
    if (!LU_build_tuples (Numeric, Work))
    {
	
	
	
	return (FALSE) ;	
    }

    Numeric->init_usage = Numeric->max_usage ;

    
    
    

    for (i = 0 ; i <= Symbolic->nfr ; i++)
    {
	Work->Front_new1strow [i] = Symbolic->Front_1strow [i] ;
    }

    return (TRUE) ;
}



GLOBAL void LU_kernel_wrapup
(
    NumericType *Numeric,
    SymbolicType *Symbolic,
    WorkType *Work
)
{

    
    
    

    Entry pivot_value ;
    double d ;
    Entry *D ;
    Int i, k, col, row, llen, ulen, *ip, *Rperm, *Cperm, *Lilen, npiv, lp,
	*Uilen, *Lip, *Uip, *Cperm_init, up, pivrow, pivcol, *Lpos, *Upos, *Wr,
	*Wc, *Wp, *Frpos, *Fcpos, *Row_degree, *Col_degree, *Rperm_init,
	n_row, n_col, n_inner, zero_pivot, nan_pivot, n1 ;

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n_inner = MIN (n_row, n_col) ;
    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;
    Lilen = Numeric->Lilen ;
    Uilen = Numeric->Uilen ;
    Upos = Numeric->Upos ;
    Lpos = Numeric->Lpos ;
    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;

    npiv = Work->npiv ;
    Numeric->npiv = npiv ;
    Numeric->ulen = Work->ulen ;

    

    
    
    

    for (k = 0 ; k < npiv ; k++)
    {
	pivot_value = D [k] ;
	ABS (d, pivot_value) ;
	zero_pivot = SCALAR_IS_ZERO (d) ;
	nan_pivot = SCALAR_IS_NAN (d) ;

	if (!zero_pivot)
	{
	    
	    Numeric->nnzpiv++ ;
	}

	if (k == 0)
	{
	    Numeric->min_udiag = d ;
	    Numeric->max_udiag = d ;
	}
	else
	{
	    

	    if (SCALAR_IS_NONZERO (Numeric->min_udiag))
	    {
		if (zero_pivot || nan_pivot)
		{
		    Numeric->min_udiag = d ;
		}
		else if (!SCALAR_IS_NAN (Numeric->min_udiag))
		{
		    
		    Numeric->min_udiag = MIN (Numeric->min_udiag, d) ;
		}
	    }

	    

	    if (nan_pivot)
	    {
		Numeric->max_udiag = d ;
	    }
	    else if (!SCALAR_IS_NAN (Numeric->max_udiag))
	    {
		
		Numeric->max_udiag = MAX (Numeric->max_udiag, d) ;
	    }
	}
    }

    
    
    

    Col_degree = Cperm ;	
    Row_degree = Rperm ;	

    if (npiv < n_row)
    {
	
	k = npiv ;
	for (row = 0 ; row < n_row ; row++)
	{
	    if (NON_PIVOTAL_ROW (row))
	    {
		Rperm [row] = ONES_COMPLEMENT (k) ;
		Lpos [row] = EMPTY ;
		Uip [row] = EMPTY ;
		Uilen [row] = 0 ;
		k++ ;
	    }
	}
    }

    if (npiv < n_col)
    {
	
	k = npiv ;
	for (col = 0 ; col < n_col ; col++)
	{
	    if (NON_PIVOTAL_COL (col))
	    {
		Cperm [col] = ONES_COMPLEMENT (k) ;
		Upos [col] = EMPTY ;
		Lip [col] = EMPTY ;
		Lilen [col] = 0 ;
		k++ ;
	    }
	}
    }

    if (npiv < n_inner)
    {
	
	for (k = npiv ; k < n_inner ; k++)
	{
	    CLEAR (D [k]) ;
	}
    }

    
    if (Numeric->ulen > 0)
    {
	Numeric->Upattern = Work->Upattern ;
	Work->Upattern = (Int *) NULL ;
    }
    if (Numeric->nnzpiv < n_inner && !SCALAR_IS_NAN (Numeric->min_udiag))
    {
	
	Numeric->min_udiag = 0.0 ;
    }

    
    
    

    Frpos = Work->Frpos ;	
    Fcpos = Work->Fcpos ;	
    Wp = Work->Wp ;		
    
    Wr = Work->Lpattern ;	
    Wc = Work->Wrp ;		

    
    
    

    

    for (pivrow = 0 ; pivrow < n_row ; pivrow++)
    {
	k = Rperm [pivrow] ;
	k = ONES_COMPLEMENT (k) ;
	Wp [k] = pivrow ;
	Frpos [pivrow] = k ;
    }
    for (k = 0 ; k < n_row ; k++)
    {
	Rperm [k] = Wp [k] ;
    }

    
    
    

    

    for (pivcol = 0 ; pivcol < n_col ; pivcol++)
    {
	k = Cperm [pivcol] ;
	k = ONES_COMPLEMENT (k) ;
	Wp [k] = pivcol ;
	
	Fcpos [pivcol] = k ;
    }
    for (k = 0 ; k < n_col ; k++)
    {
	Cperm [k] = Wp [k] ;
    }

    
    
    

    for (k = 0 ; k < npiv ; k++)
    {
	pivrow = Rperm [k] ;
	Wr [k] = Uilen [pivrow] ;
	Wp [k] = Uip [pivrow] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	Uilen [k] = Wr [k] ;
	Uip [k] = Wp [k] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	pivrow = Rperm [k] ;
	Wp [k] = Lpos [pivrow] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	Lpos [k] = Wp [k] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	pivcol = Cperm [k] ;
	Wc [k] = Lilen [pivcol] ;
	Wp [k] = Lip [pivcol] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	Lilen [k] = Wc [k] ;
	Lip [k] = Wp [k] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	pivcol = Cperm [k] ;
	Wp [k] = Upos [pivcol] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	Upos [k] = Wp [k] ;
    }

    
    
    

    Upos [npiv] = EMPTY ;
    Lpos [npiv] = EMPTY ;
    Uip [npiv] = EMPTY ;
    Lip [npiv] = EMPTY ;
    Uilen [npiv] = 0 ;
    Lilen [npiv] = 0 ;

    
    
    

    n1 = Symbolic->n1 ;

    for (k = 0 ; k < n1 ; k++)
    {
	
	ulen = Uilen [k] ;
	if (ulen > 0)
	{
	    up = Uip [k] ;
	    ip = (Int *) (Numeric->Memory + up) ;
	    for (i = 0 ; i < ulen ; i++)
	    {
		col = *ip ;
		*ip++ = Fcpos [col] ;
	    }
	}
    }

    for (k = n1 ; k < npiv ; k++)
    {
	up = Uip [k] ;
	if (up < 0)
	{
	    
	    ulen = Uilen [k] ;
	    if (ulen > 0)
	    {
		up = -up ;
		ip = (Int *) (Numeric->Memory + up) ;
		for (i = 0 ; i < ulen ; i++)
		{
		    col = *ip ;
		    *ip++ = Fcpos [col] ;
		}
	    }
	}
    }

    ulen = Numeric->ulen ;
    if (ulen > 0)
    {
	
	for (i = 0 ; i < ulen ; i++)
	{
	    col = Numeric->Upattern [i] ;
	    Numeric->Upattern [i] = Fcpos [col] ;
	}
    }

    

    
    
    

    for (k = 0 ; k < n1 ; k++)
    {
	llen = Lilen [k] ;
	if (llen > 0)
	{
	    lp = Lip [k] ;
	    ip = (Int *) (Numeric->Memory + lp) ;
	    for (i = 0 ; i < llen ; i++)
	    {
		row = *ip ;
		*ip++ = Frpos [row] ;
	    }
	}
    }

    for (k = n1 ; k < npiv ; k++)
    {
	llen = Lilen [k] ;
	if (llen > 0)
	{
	    lp = Lip [k] ;
	    if (lp < 0)
	    {
		
		lp = -lp ;
	    }
	    ip = (Int *) (Numeric->Memory + lp) ;
	    for (i = 0 ; i < llen ; i++)
	    {
		row = *ip ;
		*ip++ = Frpos [row] ;
	    }
	}
    }

    

    
    
    

    Cperm_init = Symbolic->Cperm_init ;
    Rperm_init = Symbolic->Rperm_init ;

    for (k = 0 ; k < n_row ; k++)
    {
	Rperm [k] = Rperm_init [Rperm [k]] ;
    }

    for (k = 0 ; k < n_col ; k++)
    {
	Cperm [k] = Cperm_init [Cperm [k]] ;
    }

    
    
}






GLOBAL Int LU_mem_alloc_element
(
    NumericType *Numeric,
    Int nrows,
    Int ncols,
    Int **Rows,
    Int **Cols,
    Entry **C,
    Int *size,
    Element **epout
)
{

    Element *ep ;
    Unit *p ;
    Int i ;

    ASSERT (Numeric != (NumericType *) NULL) ;
    ASSERT (Numeric->Memory != (Unit *) NULL) ;

    *size = GET_ELEMENT_SIZE (nrows, ncols) ;
    if (INT_OVERFLOW (DGET_ELEMENT_SIZE (nrows, ncols) + 1))
    {
	
	return (0) ;	
    }

    i = LU_mem_alloc_tail_block (Numeric, *size) ;
    (*size)++ ;
    if (!i)
    {
	return (0) ;	
    }
    p = Numeric->Memory + i ;

    ep = (Element *) p ;

    
    p += UNITS (Element, 1) ;		
    *Cols = (Int *) p ;			
    *Rows = *Cols + ncols ;		
    p += UNITS (Int, ncols + nrows) ;
    *C = (Entry *) p ;			

    ep->nrows = nrows ;		
    ep->ncols = ncols ;
    ep->nrowsleft = nrows ;
    ep->ncolsleft = ncols ;
    ep->cdeg = 0 ;
    ep->rdeg = 0 ;
    ep->next = EMPTY ;

    *epout = ep ;

    
    return (i) ;
}




GLOBAL Int LU_mem_alloc_head_block
(
    NumericType *Numeric,
    Int nunits
)
{
    Int p, usage ;

    if (nunits > (Numeric->itail - Numeric->ihead))
    {
	return (0) ;
    }

    
    p = Numeric->ihead ;
    Numeric->ihead += nunits ;

    usage = Numeric->ihead + Numeric->tail_usage ;
    Numeric->max_usage = MAX (Numeric->max_usage, usage) ;
    return (p) ;
}





GLOBAL Int LU_mem_alloc_tail_block
(
    NumericType *Numeric,
    Int nunits
)
{
    Int bigsize, usage ;
    Unit *p, *pnext, *pbig ;

    bigsize = 0 ;
    pbig = (Unit *) NULL ;

    if (Numeric->ibig != EMPTY)
    {
	pbig = Numeric->Memory + Numeric->ibig ;
	bigsize = -pbig->header.size ;
    }

    if (pbig && bigsize >= nunits)
    {

	
	p = pbig ;
	pnext = p + 1 + bigsize ;
	
	
	
	bigsize -= nunits + 1 ;

	if (bigsize < 4)
	{
	    
	    
	    p->header.size = -p->header.size ;
	    
	    Numeric->ibig = EMPTY ;

	}
	else
	{

	    
	    p->header.size = nunits ;
	    
	    Numeric->ibig += nunits + 1 ;
	    pbig = Numeric->Memory + Numeric->ibig ;
	    pbig->header.size = -bigsize ;
	    pbig->header.prevsize = nunits ;
	    pnext->header.prevsize = bigsize ;
	}

    }
    else
    {

	
	pnext = Numeric->Memory + Numeric->itail ;
	if ((nunits + 1) > (Numeric->itail - Numeric->ihead))
	{
	    return (0) ;
	}
	Numeric->itail -= (nunits + 1) ;
	p = Numeric->Memory + Numeric->itail ;
	p->header.size = nunits ;
	p->header.prevsize = 0 ;
	pnext->header.prevsize = nunits ;
    }

    Numeric->tail_usage += p->header.size + 1 ;
    usage = Numeric->ihead + Numeric->tail_usage ;
    Numeric->max_usage = MAX (Numeric->max_usage, usage) ;

    
    
    return ((p - Numeric->Memory) + 1) ;
}



GLOBAL void LU_mem_free_tail_block
(
    NumericType *Numeric,
    Int i
)
{
    Unit *pprev, *pnext, *p, *pbig ;
    Int sprev ;

    if (i == EMPTY || i == 0) return ;	

    
    
    

    p = Numeric->Memory + i ;

    p-- ;	

    Numeric->tail_usage -= p->header.size + 1 ;

    
    
    

    pnext = p + 1 + p->header.size ;

    if (pnext->header.size < 0)
    {
	
	p->header.size += (-(pnext->header.size)) + 1 ;
    }

    
    
    

    if (p > Numeric->Memory + Numeric->itail)
    {
	pprev = p - 1 - p->header.prevsize ;
	sprev = pprev->header.size ;
	if (sprev < 0)
	{
	    
	    pprev->header.size = p->header.size + (-sprev) + 1 ;
	    p = pprev ;
	    
	}
    }

    
    
    

    pnext = p + 1 + p->header.size ;

    if (p == Numeric->Memory + Numeric->itail)
    {
	
	Numeric->itail = pnext - Numeric->Memory ;
	pnext->header.prevsize = 0 ;
	if (Numeric->ibig != EMPTY && Numeric->ibig <= Numeric->itail)
	{
	    
	    Numeric->ibig = EMPTY ;
	}
    }
    else
    {
	
	if (Numeric->ibig == EMPTY)
	{
	    Numeric->ibig = p - Numeric->Memory ;
	}
	else
	{
	    pbig = Numeric->Memory + Numeric->ibig ;
	    if (-(pbig->header.size) < p->header.size)
	    {
		Numeric->ibig = p - Numeric->Memory ;
	    }
	}
	
	pnext->header.prevsize = p->header.size ;
	p->header.size = -(p->header.size) ;
    }

}



GLOBAL void LU_mem_init_memoryspace
(
    NumericType *Numeric
)
{
    Unit *p ;

    Numeric->ngarbage = 0 ;
    Numeric->nrealloc = 0 ;
    Numeric->ncostly = 0 ;
    Numeric->ibig = EMPTY ;
    Numeric->ihead = 0 ;
    Numeric->itail = Numeric->size ;

    
    Numeric->itail -= 2 ;
    p = Numeric->Memory + Numeric->itail ;
    Numeric->tail_usage = 2 ;
    p->header.prevsize = 0 ;
    p->header.size = 1 ;

    
    
    Numeric->ihead++ ;

    
    Numeric->max_usage = 3 ;
    Numeric->init_usage = Numeric->max_usage ;

    
    

}



GLOBAL void LU_scale
(
    Int n,
    Entry pivot,
    Entry X [ ]
)
{
    Entry x ;
    double s ;
    Int i ;

    
    
    

    APPROX_ABS (s, pivot) ;

    if (s < RECIPROCAL_TOLERANCE || IS_NAN (pivot))
    {
	
	
	

	

	for (i = 0 ; i < n ; i++)
	{
	    
	    x = X [i] ;

#ifndef NO_DIVIDE_BY_ZERO
	    if (IS_NONZERO (x))
	    {
		DIV (X [i], x, pivot) ;
	    }
#else
	    
	    if (IS_NONZERO (x) && IS_NONZERO (pivot))
	    {
		DIV (X [i], x, pivot) ;
	    }
#endif

	}

    }
    else
    {

	
	
	

	

	for (i = 0 ; i < n ; i++)
	{
	    
	    x = X [i] ;
	    DIV (X [i], x, pivot) ;
	}
    }
}









PRIVATE void shift_pivot_row (Entry *Fd, Entry *Fs, Entry *Fe, Int len, Int d)
{
    Int j ;
    for (j = 0 ; j < len ; j++)
    {
		Fd [j]   = Fs [j*d] ;
		Fs [j*d] = Fe [j*d] ;
    }
}





GLOBAL void LU_scale_column
(
    NumericType *Numeric,
    WorkType *Work
)
{
    
    
    

    Entry pivot_value ;
    Entry *Fcol, *Flublock, *Flblock, *Fublock, *Fcblock ;
    Int k, k1, fnr_curr, fnrows, fncols, *Frpos, *Fcpos, pivrow, pivcol,
	*Frows, *Fcols, fnc_curr, fnpiv, *Row_tuples, nb,
	*Col_tuples, *Rperm, *Cperm, fspos, col2, row2 ;

    
    
    

    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    fnpiv = Work->fnpiv ;

    

    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;

    

    Flublock = Work->Flublock ;
    Flblock  = Work->Flblock ;
    Fublock  = Work->Fublock ;
    Fcblock  = Work->Fcblock ;

    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;
    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;
    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    pivrow = Work->pivrow ;
    pivcol = Work->pivcol ;

    ASSERT (pivrow >= 0 && pivrow < Work->n_row) ;
    ASSERT (pivcol >= 0 && pivcol < Work->n_col) ;

    Row_tuples = Numeric->Uip ;
    Col_tuples = Numeric->Lip ;
    nb = Work->nb ;

    
    
    

    
    
    

    

    fspos = Fcpos [pivcol] ;

    
    fncols = --(Work->fncols) ;

    if (fspos != fncols * fnr_curr)
    {

	Int fs = fspos / fnr_curr ;

	
	
	

	
	{
	    
	    
	    Int i ;
	    Entry *Fs, *Fe ;
	    Fs = Fcblock + fspos ;
	    Fe = Fcblock + fncols * fnr_curr ;
#pragma ivdep
	    for (i = 0 ; i < fnrows ; i++)
	    {
		Fs [i] = Fe [i] ;
	    }
	}

	
	{
	    
	    
	    Int i ;
	    Entry *Fs, *Fe ;
	    Fs = Fublock + fs ;
	    Fe = Fublock + fncols ;
#pragma ivdep
	    for (i = 0 ; i < fnpiv ; i++)
	    {
		Fs [i * fnc_curr] = Fe [i * fnc_curr] ;
	    }
	}

	
	col2 = Fcols [fncols] ;
	Fcols [fs] = col2 ;
	Fcpos [col2] = fspos ;
    }

    
    Fcpos [pivcol] = EMPTY ;

    
    
    

    fspos = Frpos [pivrow] ;

    
    fnrows = --(Work->fnrows) ;

    if (fspos == fnrows)
    {

	
	
	

	
	{
	    Int j ;
	    Entry *Fd, *Fs ;
	    Fd = Fublock + fnpiv * fnc_curr ;
	    Fs = Fcblock + fspos ;
#pragma ivdep
	    for (j = 0 ; j < fncols ; j++)
	    {
		Fd [j] = Fs [j * fnr_curr] ;
	    }
	}

	
	if (Work->pivrow_in_front)
	{
	    Int j ;
	    Entry *Fd, *Fs ;
	    Fd = Flublock + fnpiv ;
	    Fs = Flblock  + fspos ;
#pragma ivdep
	    for (j = 0 ; j <= fnpiv ; j++)
	    {
		Fd [j * nb] = Fs [j * fnr_curr] ;
	    }
	}
	else
	{
	    Int j ;
	    Entry *Fd, *Fs ;
	    Fd = Flublock + fnpiv ;
	    Fs = Flblock  + fspos ;
#pragma ivdep
	    for (j = 0 ; j < fnpiv ; j++)
	    {
		CLEAR (Fd [j * nb]) ;
	    }
	    Fd [fnpiv * nb] = Fs [fnpiv * fnr_curr] ;
	}
    }
    else
    {

	
	
	
	

	

	
	{
	    
	    
	    
	    Entry *Fd, *Fs, *Fe ;
	    Fd = Fublock + fnpiv * fnc_curr ;
	    Fs = Fcblock + fspos ;
	    Fe = Fcblock + fnrows ;
	    shift_pivot_row (Fd, Fs, Fe, fncols, fnr_curr) ;
	}

	
	if (Work->pivrow_in_front)
	{
	    
	    
	    
	    Int j ;
	    Entry *Fd, *Fs, *Fe ;
	    Fd = Flublock + fnpiv ;
	    Fs = Flblock  + fspos ;
	    Fe = Flblock  + fnrows ;
#pragma ivdep
	    for (j = 0 ; j <= fnpiv ; j++)
	    {
		Fd [j * nb]       = Fs [j * fnr_curr] ;
		Fs [j * fnr_curr] = Fe [j * fnr_curr] ;
	    }
	}
	else
	{
	    Int j ;
	    Entry *Fd, *Fs, *Fe ;
	    Fd = Flublock + fnpiv ;
	    Fs = Flblock  + fspos ;
	    Fe = Flblock  + fnrows ;
#pragma ivdep
	    for (j = 0 ; j < fnpiv ; j++)
	    {
		ASSERT (IS_ZERO (Fs [j * fnr_curr])) ;
		CLEAR (Fd [j * nb]) ;
		Fs [j * fnr_curr] = Fe [j * fnr_curr] ;
	    }
	    Fd [fnpiv * nb]       = Fs [fnpiv * fnr_curr] ;
	    Fs [fnpiv * fnr_curr] = Fe [fnpiv * fnr_curr] ;
	}

	
	row2 = Frows [fnrows] ;
	Frows [fspos] = row2 ;
	Frpos [row2] = fspos ;

    }
    
    Frpos [pivrow] = EMPTY ;

    
    
    
    
    

    
    
    

    
    

    
    
    

    
    Fcol = Flblock + fnpiv * fnr_curr ;
    
    pivot_value = Flublock [fnpiv + fnpiv * nb] ;

    
    k = Work->npiv + fnpiv ;

    LU_scale (fnrows, pivot_value, Fcol) ;

    
    
    

    LU_mem_free_tail_block (Numeric, Row_tuples [pivrow]) ;
    LU_mem_free_tail_block (Numeric, Col_tuples [pivcol]) ;

    Row_tuples [pivrow] = 0 ;
    Col_tuples [pivcol] = 0 ;
    
    k1 = ONES_COMPLEMENT (k) ;
    Rperm [pivrow] = k1 ;			
    Cperm [pivcol] = k1 ;			

    
    
    

    
    Work->Pivrow [fnpiv] = pivrow ;
    Work->Pivcol [fnpiv] = pivcol ;

    
    
    

    

    Work->fnpiv++ ;

}




GLOBAL Int LU_valid_symbolic
(
    SymbolicType *Symbolic
)
{
    
    
    

    if (!Symbolic)
    {
		return (FALSE) ;
    }

    if (Symbolic->valid != SYMBOLIC_VALID)
    {
	
	return (FALSE) ;
    }

    if (!Symbolic->Cperm_init || !Symbolic->Rperm_init ||
	!Symbolic->Front_npivcol || !Symbolic->Front_1strow ||
	!Symbolic->Front_leftmostdesc ||
	!Symbolic->Front_parent || !Symbolic->Chain_start ||
	!Symbolic->Chain_maxrows || !Symbolic->Chain_maxcols ||
	Symbolic->n_row <= 0 || Symbolic->n_col <= 0)
    {
	return (FALSE) ;
    }

    return (TRUE) ;
}





GLOBAL Int LU_valid_numeric
(
    NumericType *Numeric
)
{
    
    
    

    if (!Numeric)
    {
	return (FALSE) ;
    }

    if (Numeric->valid != NUMERIC_VALID)
    {
	
	return (FALSE) ;
    }

    if (Numeric->n_row <= 0 || Numeric->n_col <= 0 || !Numeric->D ||
	!Numeric->Rperm || !Numeric->Cperm ||
	!Numeric->Lpos || !Numeric->Upos ||
	!Numeric->Lilen || !Numeric->Uilen || !Numeric->Lip || !Numeric->Uip ||
	!Numeric->Memory || (Numeric->ulen > 0 && !Numeric->Upattern))
    {
	return (FALSE) ;
    }

    return (TRUE) ;
}



GLOBAL Int LU_transpose
(
    Int n_row,			
    Int n_col,
    const Int Ap [ ],		
    const Int Ai [ ],		
    const double Ax [ ],	

    const Int P [ ],	
			
			

    const Int Q [ ],	
			
			
    Int nq,		

			
    Int Rp [ ],		
    Int Ri [ ],		
    double Rx [ ],	

    Int W [ ],		

    Int check		
)
{

    
    
    

    Int i, j, k, p, bp, newj, do_values ;

    
    
    

    if (check)
    {
	
	
	if (!Ai || !Ap || !Ri || !Rp || !W)
	{
	    return (SparseLU_ERROR_argument_missing) ;
	}
	if (n_row <= 0 || n_col <= 0)		
	{
	    return (SparseLU_ERROR_n_nonpositive) ;
	}
	if (!LU_is_permutation (P, W, n_row, n_row) ||
	    !LU_is_permutation (Q, W, nq, nq))
	{
	    return (SparseLU_ERROR_invalid_permutation) ;
	}
	if (AMD_valid (n_row, n_col, Ap, Ai) != AMD_OK)
	{
	    return (SparseLU_ERROR_invalid_matrix) ;
	}
    }

    
    
    

    

    for (i = 0 ; i < n_row ; i++)
    {
	W [i] = 0 ;
	Rp [i] = 0 ;
    }

    if (Q != (Int *) NULL)
    {
	for (newj = 0 ; newj < nq ; newj++)
	{
	    j = Q [newj] ;
	    ASSERT (j >= 0 && j < n_col) ;
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		ASSERT (i >= 0 && i < n_row) ;
		W [i]++ ;
	    }
	}
    }
    else
    {
	for (j = 0 ; j < n_col ; j++)
	{
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		ASSERT (i >= 0 && i < n_row) ;
		W [i]++ ;
	    }
	}
    }

    
    
    

    if (P != (Int *) NULL)
    {
	Rp [0] = 0 ;
	for (k = 0 ; k < n_row ; k++)
	{
	    i = P [k] ;
	    ASSERT (i >= 0 && i < n_row) ;
	    Rp [k+1] = Rp [k] + W [i] ;
	}
	for (k = 0 ; k < n_row ; k++)
	{
	    i = P [k] ;
	    ASSERT (i >= 0 && i < n_row) ;
	    W [i] = Rp [k] ;
	}
    }
    else
    {
	Rp [0] = 0 ;
	for (i = 0 ; i < n_row ; i++)
	{
	    Rp [i+1] = Rp [i] + W [i] ;
	}
	for (i = 0 ; i < n_row ; i++)
	{
	    W [i] = Rp [i] ;
	}
    }
    ASSERT (Rp [n_row] <= Ap [n_col]) ;

    

    
    
    

    do_values = Ax && Rx ;
    {
	if (Q != (Int *) NULL)
	{
	    if (do_values)
	    {
		{
		    
		    for (newj = 0 ; newj < nq ; newj++)
		    {
			j = Q [newj] ;
			ASSERT (j >= 0 && j < n_col) ;
			for (p = Ap [j] ; p < Ap [j+1] ; p++)
			{
			    bp = W [Ai [p]]++ ;
			    Ri [bp] = newj ;
			    Rx [bp] = Ax [p] ;
			}
		    }
		}
	    }
	    else
	    {
		
		for (newj = 0 ; newj < nq ; newj++)
		{
		    j = Q [newj] ;
		    ASSERT (j >= 0 && j < n_col) ;
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			Ri [W [Ai [p]]++] = newj ;
		    }
		}
	    }
	}
	else
	{
	    if (do_values)
	    {
		{
		    
		    for (j = 0 ; j < n_col ; j++)
		    {
			for (p = Ap [j] ; p < Ap [j+1] ; p++)
			{
			    bp = W [Ai [p]]++ ;
			    Ri [bp] = j ;
			    Rx [bp] = Ax [p] ;
			}
		    }
		}
	    }
	    else
	    {
		
		for (j = 0 ; j < n_col ; j++)
		{
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			Ri [W [Ai [p]]++] = j ;
		    }
		}
	    }
	}
    }

    return (SparseLU_OK) ;
}



GLOBAL void LU_set_stats
(
    double Info [ ],
    SymbolicType *Symbolic,
    double max_usage,		
    double num_mem_size,	
    double flops,		
    double lnz,			
    double unz,			
    double maxfrsize,		
    double ulen,		
    double npiv,		
    double maxnrows,		
    double maxncols,		
    Int scale,			
    Int prefer_diagonal,	
    Int what			
)
{

    double sym_size, work_usage, nn, n_row, n_col, n_inner, num_On_size1,
	num_On_size2, num_usage, sym_maxncols, sym_maxnrows, elen, n1 ;

    n_col = Symbolic->n_col ;
    n_row = Symbolic->n_row ;
    n1 = Symbolic->n1 ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    sym_maxncols = MIN (Symbolic->maxncols + Symbolic->nb, n_col) ;
    sym_maxnrows = MIN (Symbolic->maxnrows + Symbolic->nb, n_row) ;
    elen = (n_col - n1) + (n_row - n1) + MIN (n_col - n1, n_row - n1) + 1 ;

    
    sym_size = LU_symbolic_usage (Symbolic->n_row, Symbolic->n_col,
	Symbolic->nchains, Symbolic->nfr, Symbolic->esize, prefer_diagonal) ;

    
    
    num_On_size1 =
	DUNITS (NumericType, 1)		
	+ DUNITS (Entry, n_inner+1)	
	+ 4 * DUNITS (Int, n_row+1)	
	+ 4 * DUNITS (Int, n_col+1)	
	+ (scale ? DUNITS (Entry, n_row) : 0) ;   

    
    
    num_On_size2 =
	DUNITS (NumericType, 1)		
	+ DUNITS (Entry, n_inner+1)	
	+ DUNITS (Int, n_row+1)		
	+ DUNITS (Int, n_col+1)		
	+ 6 * DUNITS (Int, npiv+1)	
	+ (scale ? DUNITS (Entry, n_row) : 0) ;	    

    
    Info [SparseLU_VARIABLE_PEAK + what] = max_usage ;

    
    Info [SparseLU_VARIABLE_FINAL + what] = num_mem_size ;

    
    Info [SparseLU_NUMERIC_SIZE + what] =
	num_On_size2
	+ num_mem_size		
	+ DUNITS (Int, ulen+1) ;

    
    Info [SparseLU_MAX_FRONT_SIZE + what] = maxfrsize ;
    Info [SparseLU_MAX_FRONT_NROWS + what] = maxnrows ;
    Info [SparseLU_MAX_FRONT_NCOLS + what] = maxncols ;

    
    work_usage =
	
	2 * DUNITS (Entry, sym_maxnrows + 1)	
	+ 2 * DUNITS (Int, n_row+1)		
	+ 2 * DUNITS (Int, n_col+1)		
	+ DUNITS (Int, nn + 1)			
	+ DUNITS (Int, MAX (n_col, sym_maxnrows) + 1)	
	+ 2 * DUNITS (Int, sym_maxnrows + 1)	
	+ 3 * DUNITS (Int, sym_maxncols + 1)	
	+ DUNITS (Int, MAX (sym_maxnrows, sym_maxncols) + 1)	
	+ DUNITS (Int, elen)			
	+ DUNITS (Int, Symbolic->nfr + 1)	
	+ ((n_row == n_col) ? (2 * DUNITS (Int, nn)) : 0) ;  

    
    num_usage =
	sym_size	
	+ num_On_size1	
	+ work_usage	
	+ max_usage ;	

    
    Info [SparseLU_PEAK_MEMORY + what] =
	MAX (Symbolic->peak_sym_usage, num_usage) ;

    Info [SparseLU_FLOPS + what] = flops ;
    Info [SparseLU_LNZ + what] = lnz ;
    Info [SparseLU_UNZ + what] = unz ;
}







GLOBAL Int LU_tuple_lengths	    
(
    NumericType *Numeric,
    WorkType *Work,
    double *p_dusage		    
)
{
    
    
    

    double dusage ;
    Int e, nrows, ncols, nel, i, *Rows, *Cols, row, col, n_row, n_col, *E,
	*Row_degree, *Row_tlen, *Col_degree, *Col_tlen, usage, n1 ;
    Element *ep ;
    Unit *p ;

    
    
    

    E = Work->E ;
    Row_degree = Numeric->Rperm ;   
    Col_degree = Numeric->Cperm ;   
    Row_tlen   = Numeric->Uilen ;
    Col_tlen   = Numeric->Lilen ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    n1 = Work->n1 ;
    nel = Work->nel ;

    

    
    
    

    for (e = 1 ; e <= nel ; e++)	
    {
	if (E [e])
	{
	    p = Numeric->Memory + E [e] ;
	    GET_ELEMENT_PATTERN (ep, p, Cols, Rows, ncols) ;
	    nrows = ep->nrows ;
	    for (i = 0 ; i < nrows ; i++)
	    {
		row = Rows [i] ;
		ASSERT (row == EMPTY || (row >= n1 && row < n_row)) ;
		if (row >= n1)
		{
		    ASSERT (NON_PIVOTAL_ROW (row)) ;
		    Row_tlen [row] ++ ;
		}
	    }
	    for (i = 0 ; i < ncols ; i++)
	    {
		col = Cols [i] ;
		ASSERT (col == EMPTY || (col >= n1 && col < n_col)) ;
		if (col >= n1)
		{
		    ASSERT (NON_PIVOTAL_COL (col)) ;
		    Col_tlen [col] ++ ;
		}
	    }
	}
    }

    
    

    
    
    
    usage = 0 ;
    dusage = 0 ;

    for (col = n1 ; col < n_col ; col++)
    {
	if (NON_PIVOTAL_COL (col))
	{
	    usage  += 1 +  UNITS (Tuple, TUPLES (Col_tlen [col])) ;
	    dusage += 1 + DUNITS (Tuple, TUPLES (Col_tlen [col])) ;
	}
    }

    for (row = n1 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row))
	{
	    usage  += 1 +  UNITS (Tuple, TUPLES (Row_tlen [row])) ;
	    dusage += 1 + DUNITS (Tuple, TUPLES (Row_tlen [row])) ;
	}
    }
    *p_dusage = dusage ;
    return (usage) ;
}



GLOBAL double LU_symbolic_usage
(
    Int n_row,
    Int n_col,
    Int nchains,
    Int nfr,
    Int esize,	    
    Int prefer_diagonal
)
{
    double units ;

    units =
	DUNITS (SymbolicType, 1)	
	+ 2 * DUNITS (Int, n_col+1)	
	+ 2 * DUNITS (Int, n_row+1)	
	+ 3 * DUNITS (Int, nchains+1)	
	+ 4 * DUNITS (Int, nfr+1) ;	

    
    units += DUNITS (Int, esize) ;	

    
    if (prefer_diagonal)
    {
	units += DUNITS (Int, n_col+1) ;    
    }

    return (units) ;
}


//   LU_lsolve
GLOBAL double LU_lsolve
(
    NumericType *Numeric,
    Entry X [ ],		
    Int Pattern [ ]		
)
{
    Entry xk ;
    Entry *xp, *Lval ;
    Int k, deg, *ip, j, row, *Lpos, *Lilen, *Lip, llen, lp, newLchain,
	pos, npiv, n1, *Li ;

    

    if (Numeric->n_row != Numeric->n_col) return (0.) ;
    npiv = Numeric->npiv ;
    Lpos = Numeric->Lpos ;
    Lilen = Numeric->Lilen ;
    Lip = Numeric->Lip ;
    n1 = Numeric->n1 ;

    
    
    

    for (k = 0 ; k < n1 ; k++)
    {
	xk = X [k] ;
	deg = Lilen [k] ;
	if (deg > 0 && IS_NONZERO (xk))
	{
	    lp = Lip [k] ;
	    Li = (Int *) (Numeric->Memory + lp) ;
	    lp += UNITS (Int, deg) ;
	    Lval = (Entry *) (Numeric->Memory + lp) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		
		MULT_SUB (X [Li [j]], xk, Lval [j]) ;
	    }
	}
    }

    
    
    

    deg = 0 ;

    for (k = n1 ; k < npiv ; k++)
    {

	
	
	

	lp = Lip [k] ;
	newLchain = (lp < 0) ;
	if (newLchain)
	{
	    lp = -lp ;
	    deg = 0 ;
	}

	
	pos = Lpos [k] ;
	if (pos != EMPTY)
	{

	    Pattern [pos] = Pattern [--deg] ;
	}

	
	ip = (Int *) (Numeric->Memory + lp) ;
	llen = Lilen [k] ;
	for (j = 0 ; j < llen ; j++)
	{
	    row = *ip++ ;
	    Pattern [deg++] = row ;
	}

	
	
	

	xk = X [k] ;
	if (IS_NONZERO (xk))
	{
	    xp = (Entry *) (Numeric->Memory + lp + UNITS (Int, llen)) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		
		MULT_SUB (X [Pattern [j]], xk, *xp) ;
		xp++ ;
	    }
	}
    }

    return (MULTSUB_FLOPS * ((double) Numeric->lnz)) ;
}

//   LU_usolve
GLOBAL double LU_usolve
(
    NumericType *Numeric,
    Entry X [ ],		
    Int Pattern [ ]		
)
{
    
    
    

    Entry xk ;
    Entry *xp, *D, *Uval ;
    Int k, deg, j, *ip, col, *Upos, *Uilen, pos,
	*Uip, n, ulen, up, newUchain, npiv, n1, *Ui ;

    
    
    

    if (Numeric->n_row != Numeric->n_col) return (0.) ;
    n = Numeric->n_row ;
    npiv = Numeric->npiv ;
    Upos = Numeric->Upos ;
    Uilen = Numeric->Uilen ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;
    n1 = Numeric->n1 ;

    
    
    

#ifndef NO_DIVIDE_BY_ZERO
    
    for (k = n-1 ; k >= npiv ; k--)
    {
	
	xk = X [k] ;
	
	DIV (X [k], xk, D [k]) ;
    }
#else
    
#endif

    deg = Numeric->ulen ;
    if (deg > 0)
    {
	
	for (j = 0 ; j < deg ; j++)
	{
	    Pattern [j] = Numeric->Upattern [j] ;
	}
    }

    
    
    

    for (k = npiv-1 ; k >= n1 ; k--)
    {

	
	
	

	up = Uip [k] ;
	ulen = Uilen [k] ;
	newUchain = (up < 0) ;
	if (newUchain)
	{
	    up = -up ;
	    xp = (Entry *) (Numeric->Memory + up + UNITS (Int, ulen)) ;
	}
	else
	{
	    xp = (Entry *) (Numeric->Memory + up) ;
	}

	xk = X [k] ;
	for (j = 0 ; j < deg ; j++)
	{
	    
	    MULT_SUB (xk, X [Pattern [j]], *xp) ;
	    xp++ ;
	}

#ifndef NO_DIVIDE_BY_ZERO
	
	
	DIV (X [k], xk, D [k]) ;
#else
	
	if (IS_NONZERO (D [k]))
	{
	    
	    DIV (X [k], xk, D [k]) ;
	}
#endif

	
	
	

	if (k == n1) break ;

	if (newUchain)
	{
	    
	    deg = ulen ;
	    ip = (Int *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		col = *ip++ ;
		Pattern [j] = col ;
	    }
	}
	else
	{
	    deg -= ulen ;
	    pos = Upos [k] ;
	    if (pos != EMPTY)
	    {
		
		Pattern [deg++] = Pattern [pos] ;
		Pattern [pos] = k ;
	    }
	}
    }

    
    
    

    for (k = n1 - 1 ; k >= 0 ; k--)
    {
	deg = Uilen [k] ;
	xk = X [k] ;
	if (deg > 0)
	{
	    up = Uip [k] ;
	    Ui = (Int *) (Numeric->Memory + up) ;
	    up += UNITS (Int, deg) ;
	    Uval = (Entry *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		

		MULT_SUB (xk, X [Ui [j]], Uval [j]) ;
	    }
	}

#ifndef NO_DIVIDE_BY_ZERO
	
	
	DIV (X [k], xk, D [k]) ;
#else
	
	if (IS_NONZERO (D [k]))
	{
	    
	    DIV (X [k], xk, D [k]) ;
	}
#endif
    }

    return (DIV_FLOPS * ((double) n) + MULTSUB_FLOPS * ((double) Numeric->unz));
}


GLOBAL void LU_blas3_update (WorkType *Work)
{
    
    
    
    Entry *L, *U, *C, *LU ;
    Int i, j, s, k, m, n, d, nb, dc ;
    Int blas_ok = TRUE ;

    k = Work->fnpiv ;
    if (k == 0)
    {
		
		return ;
    }

    m = Work->fnrows ;
    n = Work->fncols ;
    d = Work->fnr_curr ;
    dc = Work->fnc_curr ;
    nb = Work->nb ;
    C = Work->Fcblock ;	    
    L =	Work->Flblock ;	    
    U = Work->Fublock ;	    
    LU = Work->Flublock ;   

    if (k == 1)
    {
		BLAS_GER (m, n, L, U, C, d) ;
    }
    else
    {
		
		BLAS_TRSM_RIGHT (n, k, LU, nb, U, dc) ;

		
		
		BLAS_GEMM (m, n, k, L, U, dc, C, d) ;
	}
}




