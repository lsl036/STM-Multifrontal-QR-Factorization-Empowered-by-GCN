
/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_base.c
 * BRIEF:   LU基本
 *****************************************************************************/

#include "SparseLU_internal.h"
#include "SparseLU_function.h"



GLOBAL Int LU_build_tuples
(
    NumericType *Numeric,
    WorkType *Work
)
{
    
    
    

    Int e, nrows, ncols, nel, *Rows, *Cols, row, col, n_row, n_col, *E,
	*Row_tuples, *Row_degree, *Row_tlen,
	*Col_tuples, *Col_degree, *Col_tlen, n1 ;
    Element *ep ;
    Unit *p ;
    Tuple tuple, *tp ;

    
    
    

    E = Work->E ;
    Col_degree = Numeric->Cperm ;	
    Row_degree = Numeric->Rperm ;	
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nel = Work->nel ;
    n1 = Work->n1 ;

    
    
    

    
    
    

    for (row = n1 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row))
	{
	    Row_tuples [row] = LU_mem_alloc_tail_block (Numeric,
		UNITS (Tuple, TUPLES (Row_tlen [row]))) ;
	    if (!Row_tuples [row])
	    {
		
		return (FALSE) ;	
	    }
	    Row_tlen [row] = 0 ;
	}
    }

    
    
    for (col = n_col-1 ; col >= n1 ; col--)
    {
	if (NON_PIVOTAL_COL (col))
	{
	    Col_tuples [col] = LU_mem_alloc_tail_block (Numeric,
		UNITS (Tuple, TUPLES (Col_tlen [col]))) ;
	    if (!Col_tuples [col])
	    {
		
		return (FALSE) ;	
	    }
	    Col_tlen [col] = 0 ;
	}
    }

    
    
    

    
    for (e = 1 ; e <= nel ; e++)
    {
	p = Numeric->Memory + E [e] ;
	GET_ELEMENT_PATTERN (ep, p, Cols, Rows, ncols) ;
	nrows = ep->nrows ;
	tuple.e = e ;
	for (tuple.f = 0 ; tuple.f < ncols ; tuple.f++)
	{
	    col = Cols [tuple.f] ;
	    tp = ((Tuple *) (Numeric->Memory + Col_tuples [col]))
		+ Col_tlen [col]++ ;
	    *tp = tuple ;
	}
	for (tuple.f = 0 ; tuple.f < nrows ; tuple.f++)
	{
	    row = Rows [tuple.f] ;
	    tp = ((Tuple *) (Numeric->Memory + Row_tuples [row]))
		+ Row_tlen [row]++ ;
	    *tp = tuple ;
	}
    }

    
    
    
    return (TRUE) ;
}







PRIVATE void copy_column (Int len, Entry *X, Entry *Y)
{
    Int i ;
#pragma ivdep
    for (i = 0 ; i < len ; i++)
    {
	Y [i] = X [i] ;
    }
}





GLOBAL Int LU_create_element
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    
    
    

    Int j, col, row, *Fcols, *Frows, fnrows, fncols, *Cols, len, needunits, t1,
	t2, size, e, i, *E, *Fcpos, *Frpos, *Rows, eloc, fnr_curr, f,
	got_memory, *Row_tuples, *Row_degree, *Row_tlen, *Col_tuples, max_mark,
	*Col_degree, *Col_tlen, nn, n_row, n_col, r2, c2, do_Fcpos ;
    Entry *C, *Fcol ;
    Element *ep ;
    Unit *p, *Memory ;
    Tuple *tp, *tp1, *tp2, tuple, *tpend ;

    
    
    

    Row_degree = Numeric->Rperm ;
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Col_degree = Numeric->Cperm ;
    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nn = MAX (n_row, n_col) ;
    Fcols = Work->Fcols ;
    Frows = Work->Frows ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Memory = Numeric->Memory ;
    fncols = Work->fncols ;
    fnrows = Work->fnrows ;

    tp = (Tuple *) NULL ;
    tp1 = (Tuple *) NULL ;
    tp2 = (Tuple *) NULL ;

    
    
    

    if (!Symbolic->fixQ)
    {
		
	#pragma ivdep
		for (j = 0 ; j < fncols ; j++)
		{
			
			Col_degree [Fcols [j]] += fnrows ;
		}
    }

    
    
    

#pragma ivdep
    for (i = 0 ; i < fnrows ; i++)
    {
		
		Row_degree [Frows [i]] += fncols ;
    }

    
    
    

    E = Work->E ;
    max_mark = MAX_MARK (nn) ;

    if (!Work->pivcol_in_front)
    {
	
	Work->cdeg0 += (nn + 1) ;
	if (Work->cdeg0 >= max_mark)
	{
	    
	    Work->cdeg0 = 1 ;
	#pragma ivdep
	    for (e = 1 ; e <= Work->nel ; e++)
	    {
			if (E [e])
			{
				ep = (Element *) (Memory + E [e]) ;
				ep->cdeg = 0 ;
			}
	    }
	}
    }

    if (!Work->pivrow_in_front)
    {
		
		Work->rdeg0 += (nn + 1) ;
		if (Work->rdeg0 >= max_mark)
		{
			
			Work->rdeg0 = 1 ;
		#pragma ivdep
			for (e = 1 ; e <= Work->nel ; e++)
			{
				if (E [e])
				{
					ep = (Element *) (Memory + E [e]) ;
					ep->rdeg = 0 ;
				}
			}
		}
    }

    
    
    

    if (!Work->pivrow_in_front)
    {
	#pragma ivdep
		for (j = 0 ; j < fncols ; j++)
		{
			Fcpos [Fcols [j]] = EMPTY ;
		}
    }

    if (!Work->pivcol_in_front)
    {
	#pragma ivdep
		for (i = 0 ; i < fnrows ; i++)
		{
			Frpos [Frows [i]] = EMPTY ;
		}
    }

    if (fncols <= 0 || fnrows <= 0)
    {
		
		Work->prior_element = EMPTY ;
		return (TRUE) ;
    }

    
    
    

    needunits = 0 ;
    got_memory = FALSE ;
    eloc = LU_mem_alloc_element (Numeric, fnrows, fncols, &Rows, &Cols, &C,
	&needunits, &ep) ;

    
    if (Work->do_grow)
    {
	
	r2 = fnrows ;
	c2 = fncols ;
	do_Fcpos = FALSE ;
    }
    else
    {
	
	r2 = MAX (fnrows, Work->fnrows_new + 1) ;
	c2 = MAX (fncols, Work->fncols_new + 1) ;
	
	do_Fcpos = Work->pivrow_in_front ;
    }

    if (!eloc)
    {
	
	
	
	if (!LU_get_memory (Numeric, Work, needunits, r2, c2, do_Fcpos))
	{
	    
	    return (FALSE) ;	
	}
	got_memory = TRUE ;
	Memory = Numeric->Memory ;
	eloc = LU_mem_alloc_element (Numeric, fnrows, fncols, &Rows, &Cols, &C,
	    &needunits, &ep) ;
	if (!eloc)
	{
	    
	    return (FALSE) ;	
	}
    }

    e = ++(Work->nel) ;	
    Work->prior_element = e ;

    E [e] = eloc ;

    if (Work->pivcol_in_front)
    {
	
	ep->cdeg = Work->cdeg0 ;
    }

    if (Work->pivrow_in_front)
    {
	
	ep->rdeg = Work->rdeg0 ;
    }

    
    
    

#pragma ivdep
    for (i = 0 ; i < fnrows ; i++)
    {
	Rows [i] = Frows [i] ;
    }
#pragma ivdep
    for (i = 0 ; i < fncols ; i++)
    {
	Cols [i] = Fcols [i] ;
    }
    Fcol = Work->Fcblock ;
    fnr_curr = Work->fnr_curr ;
    for (j = 0 ; j < fncols ; j++)
    {
		copy_column (fnrows, Fcol, C) ;
		Fcol += fnr_curr ;
		C += fnrows ;
    }

    
    
    

    tuple.e = e ;

    if (got_memory)
    {

	
	
	

	
	for (tuple.f = 0 ; tuple.f < fncols ; tuple.f++)
	{
	    col = Fcols [tuple.f] ;
	    tp = ((Tuple *) (Memory + Col_tuples [col])) + Col_tlen [col]++ ;
	    *tp = tuple ;
	}

	

	
	for (tuple.f = 0 ; tuple.f < fnrows ; tuple.f++)
	{
	    row = Frows [tuple.f] ;
	    tp = ((Tuple *) (Memory + Row_tuples [row])) + Row_tlen [row]++ ;
	    *tp = tuple ;
	}

    }
    else
    {

	
	
	

	

	for (tuple.f = 0 ; tuple.f < fncols ; tuple.f++)
	{
	    col = Fcols [tuple.f] ;
	    t1 = Col_tuples [col] ;

	    size = 0 ;
	    len = 0 ;

	    if (t1)
	    {
		p = Memory + t1 ;
		tp = (Tuple *) p ;
		size = GET_BLOCK_SIZE (p) ;
		len = Col_tlen [col] ;
		tp2 = tp + len ;
	    }

	    needunits = UNITS (Tuple, len + 1) ;
	    if (needunits > size && t1)
	    {
		
		tp1 = tp ;
		tp2 = tp ;
		tpend = tp + len ;
		for ( ; tp < tpend ; tp++)
		{
		    e = tp->e ;
		    if (!E [e]) continue ;   
		    f = tp->f ;
		    p = Memory + E [e] ;
		    ep = (Element *) p ;
		    p += UNITS (Element, 1) ;
		    Cols = (Int *) p ;
		    ;
		    if (Cols [f] == EMPTY) continue ;	
		    *tp2++ = *tp ;	
		}
		len = tp2 - tp1 ;
		Col_tlen [col] = len ;
		needunits = UNITS (Tuple, len + 1) ;
	    }

	    if (needunits > size)
	    {
		
		needunits = MIN (2*needunits, (Int) UNITS (Tuple, nn)) ;
		t2 = LU_mem_alloc_tail_block (Numeric, needunits) ;
		if (!t2)
		{
		    
		    
		    
		    
		    return (LU_get_memory (Numeric, Work, 0, r2, c2,do_Fcpos));
		}
		Col_tuples [col] = t2 ;
		tp2 = (Tuple *) (Memory + t2) ;
		if (t1)
		{
		    for (i = 0 ; i < len ; i++)
		    {
			*tp2++ = *tp1++ ;
		    }
		    LU_mem_free_tail_block (Numeric, t1) ;
		}
	    }

	    
	    Col_tlen [col]++ ;
	    *tp2 = tuple ;
	}

	
	
	

	for (tuple.f = 0 ; tuple.f < fnrows ; tuple.f++)
	{
	    row = Frows [tuple.f] ;
	    t1 = Row_tuples [row] ;

	    size = 0 ;
	    len = 0 ;
	    if (t1)
	    {
		p = Memory + t1 ;
		tp = (Tuple *) p ;
		size = GET_BLOCK_SIZE (p) ;
		len = Row_tlen [row] ;
		tp2 = tp + len ;
	    }

	    needunits = UNITS (Tuple, len + 1) ;

	    if (needunits > size && t1)
	    {
		
		tp1 = tp ;
		tp2 = tp ;
		tpend = tp + len ;
		for ( ; tp < tpend ; tp++)
		{
		    e = tp->e ;
		    if (!E [e])
		    {
			continue ;	
		    }
		    f = tp->f ;
		    p = Memory + E [e] ;
		    ep = (Element *) p ;
		    p += UNITS (Element, 1) ;
		    Cols = (Int *) p ;
		    Rows = Cols + (ep->ncols) ;
		    if (Rows [f] == EMPTY) continue ;	
		    *tp2++ = *tp ;	
		}
		len = tp2 - tp1 ;
		Row_tlen [row] = len ;
		needunits = UNITS (Tuple, len + 1) ;
	    }

	    if (needunits > size)
	    {
		

		needunits = MIN (2*needunits, (Int) UNITS (Tuple, nn)) ;
		t2 = LU_mem_alloc_tail_block (Numeric, needunits) ;
		if (!t2)
		{
		    
		    
		    
		    
		    return (LU_get_memory (Numeric, Work, 0, r2, c2,do_Fcpos));
		}
		Row_tuples [row] = t2 ;
		tp2 = (Tuple *) (Memory + t2) ;
		if (t1)
		{
		    for (i = 0 ; i < len ; i++)
		    {
			*tp2++ = *tp1++ ;
		    }
		    LU_mem_free_tail_block (Numeric, t1) ;
		}
	    }

	    
	    Row_tlen [row]++ ;
	    *tp2 = tuple ;
	}

    }

    

    return (TRUE) ;
}







PRIVATE void zero_front (
    Entry *Flblock, Entry *Fublock, Entry *Fcblock,
    Int fnrows, Int fncols, Int fnr_curr, Int fnc_curr,
    Int fnpiv, Int fnrows_extended, Int fncols_extended)
{
    Int j, i ;
    Entry *F, *Fj, *Fi ;

    Fj = Fcblock + fnrows ;
    for (j = 0 ; j < fncols ; j++)
    {
		
		F = Fj ;
		Fj += fnr_curr ;
	#pragma ivdep
		for (i = fnrows ; i < fnrows_extended ; i++)
		{
			
			CLEAR_AND_INCREMENT (F) ;
		}
    }

    Fj -= fnrows ;
    for (j = fncols ; j < fncols_extended ; j++)
    {
	
	F = Fj ;
	Fj += fnr_curr ;
#pragma ivdep
	for (i = 0 ; i < fnrows_extended ; i++)
	{
	    
	    CLEAR_AND_INCREMENT (F) ;
	}
    }

    Fj = Flblock + fnrows ;
    for (j = 0 ; j < fnpiv ; j++)
    {
	
	F = Fj ;
	Fj += fnr_curr ;
#pragma ivdep
	for (i = fnrows ; i < fnrows_extended ; i++)
	{
	    
	    CLEAR_AND_INCREMENT (F) ;
	}
    }

    Fi = Fublock + fncols ;
    for (i = 0 ; i < fnpiv ; i++)
    {
		
		F = Fi ;
		Fi += fnc_curr ;
	#pragma ivdep
		for (j = fncols ; j < fncols_extended ; j++)
		{
			
			CLEAR_AND_INCREMENT (F) ;
		}
    }

}





GLOBAL Int LU_extend_front
(
    NumericType *Numeric,
    WorkType *Work
)
{
    
    
    

    Int j, i, *Frows, row, col, *Wrow, fnr2, fnc2, *Frpos, *Fcpos, *Fcols,
	fnrows_extended, rrdeg, ccdeg, fncols_extended, fnr_curr, fnc_curr,
	fnrows, fncols, pos, fnpiv, *Wm ;
    Entry *Wx, *Wy, *Fu, *Fl ;

    
    
    

    fnpiv = Work->fnpiv ;

    if (Work->do_grow)
    {
	fnr2 = LU_FRONTAL_GROWTH * Work->fnrows_new + 2 ;
	fnc2 = LU_FRONTAL_GROWTH * Work->fncols_new + 2 ;
	if (!LU_grow_front (Numeric, fnr2, fnc2, Work, 1))
	{
	    return (FALSE) ;
	}
    }

    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;

    
    
    

    Frows = Work->Frows ;
    Frpos = Work->Frpos ;
    Fcols = Work->Fcols ;
    Fcpos = Work->Fcpos ;
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    rrdeg = Work->rrdeg ;
    ccdeg = Work->ccdeg ;

    
    
    Work->fscan_col = fncols ;
    Work->NewCols = Fcols ;

    
    
    Work->fscan_row = fnrows ;
    Work->NewRows = Frows ;

    
    
    

    fnrows_extended = fnrows ;
    fncols_extended = fncols ;
    Fl = Work->Flblock + fnpiv * fnr_curr ;

    if (Work->pivcol_in_front)
    {
	
	fnrows_extended += ccdeg ;
	Wy = Work->Wy ;

	for (i = 0 ; i < fnrows_extended ; i++)
	{
	    Fl [i] = Wy [i] ;
	}

    }
    else
    {
	
	Entry *F ;
	Fu = Work->Flublock + fnpiv * Work->nb ;
	Wm = Work->Wm ;
	Wx = Work->Wx ;
	F = Fu ;
	for (i = 0 ; i < fnpiv ; i++)
	{
	    CLEAR_AND_INCREMENT (F) ;
	}
	F = Fl ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    CLEAR_AND_INCREMENT (F) ;
	}
	for (i = 0 ; i < ccdeg ; i++)
	{
	    row = Wm [i] ;
	    pos = Frpos [row] ;
	    if (pos < 0)
	    {
		pos = fnrows_extended++ ;
		Frows [pos] = row ;
		Frpos [row] = pos ;
	    }
	    Fl [pos] = Wx [i] ;
	}
    }

    ASSERT (fnrows_extended <= fnr_curr) ;

    
    
    

    if (Work->pivrow_in_front)
    {
	if (Work->pivcol_in_front)
	{
	    for (j = fncols ; j < rrdeg ; j++)
	    {
		Fcpos [Fcols [j]] = j * fnr_curr ;
	    }
	}
	else
	{
	    
	    Wrow = Work->Wrow ;
	    if (Wrow == Fcols)
	    {
		
		for (j = fncols ; j < rrdeg ; j++)
		{
		    col = Wrow [j] ;
		    
		    Fcpos [col] = j * fnr_curr ;
		}
	    }
	    else
	    {
		for (j = fncols ; j < rrdeg ; j++)
		{
		    col = Wrow [j] ;
		    Fcols [j] = col ;
		    Fcpos [col] = j * fnr_curr ;
		}
	    }
	}
	fncols_extended = rrdeg ;
    }
    else
    {
	Wrow = Work->Wrow ;
	for (j = 0 ; j < rrdeg ; j++)
	{
	    col = Wrow [j] ;
	    if (Fcpos [col] < 0)
	    {
		Fcols [fncols_extended] = col ;
		Fcpos [col] = fncols_extended * fnr_curr ;
		fncols_extended++ ;
	    }
	}
    }

    
    
    

    
    
    

    zero_front (Work->Flblock, Work->Fublock, Work->Fcblock,
	fnrows, fncols, fnr_curr, fnc_curr,
	fnpiv, fnrows_extended, fncols_extended) ;

    
    
    

    Work->fnrows = fnrows_extended ;
    Work->fncols = fncols_extended ;

    return (TRUE) ;

}







GLOBAL void LU_garbage_collection
(
    NumericType *Numeric,
    WorkType *Work,
    Int drnew,	    
    Int dcnew,
    Int do_Fcpos
)
{
    
    
    

    Int size, e, n_row, n_col, nrows, ncols, nrowsleft, ncolsleft, prevsize,
	csize, size2, i2, j2, i, j, cdeg, rdeg, *E, row, col,
	*Rows, *Cols, *Rows2, *Cols2, nel, e2, *Row_tuples, *Col_tuples,
	*Row_degree, *Col_degree ;
    Entry *C, *C1, *C3, *C2 ;
    Unit *psrc, *pdest, *p, *pnext ;
    Element *epsrc, *epdest ;

    
    
    

    Col_degree = Numeric->Cperm ;	
    Row_degree = Numeric->Rperm ;	
    Row_tuples = Numeric->Uip ;
    Col_tuples = Numeric->Lip ;
    E = Work->E ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;

    
    
    

    Numeric->ngarbage++ ;

    
    
    

    
    

    for (row = 0 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row) && Row_tuples [row])
	{
	    p = Numeric->Memory + Row_tuples [row] - 1 ;
	    p->header.size = -p->header.size ;
	    Row_tuples [row] = 0 ;
	}
    }

    for (col = 0 ; col < n_col ; col++)
    {
	if (NON_PIVOTAL_COL (col) && Col_tuples [col])
	{
	    p = Numeric->Memory + Col_tuples [col] - 1 ;
	    p->header.size = -p->header.size ;
	    Col_tuples [col] = 0 ;
	}
    }

    
    
    

    nel = Work->nel ;
    ASSERT (nel < Work->elen) ;

    e2 = 0 ;

    for (e = 0 ; e <= nel ; e++) 
    {
	if (E [e])
	{
	    psrc = Numeric->Memory + E [e] ;
	    psrc-- ;		
	    if (e > 0)
	    {
		e2++ ;	
	    }
	    psrc->header.size = e2  ;	
	    E [e] = 0 ;
	    if (e == Work->prior_element)
	    {
		Work->prior_element = e2 ;
	    }
	}
    }

    
    Work->nel = e2 ;
    nel = Work->nel ;

    
    
    

    
    psrc = Numeric->Memory + Numeric->size - 2 ;
    pdest = psrc ;
    prevsize = psrc->header.prevsize ;

    while (prevsize > 0)
    {

	
	
	
	
	

	size = prevsize ;
	psrc -= (size + 1) ;
	e = psrc->header.size ;
	prevsize = psrc->header.prevsize ;
	

	
	

	if (e == 0)
	{
	    
	    
	    

	    Entry *F1, *F2, *Fsrc, *Fdst ;
	    Int c, r, k, dr, dc, gap, gap1, gap2, nb ;

	    
	    F1 = (Entry *) (psrc + 1) ;

	    
	    k = Work->fnpiv ;
	    dr = Work->fnr_curr ;
	    dc = Work->fnc_curr ;
	    r = Work->fnrows ;
	    c = Work->fncols ;
	    nb = Work->nb ;

	    ASSERT ((dr >= 0 && (dr % 2) == 1) || dr == 0) ;
	    ASSERT (drnew >= 0) ;
	    if (drnew % 2 == 0)
	    {
		
		drnew++ ;
	    }
	    drnew = MIN (dr, drnew) ;
	    ASSERT ((drnew >= 0 && (drnew % 2) == 1) || drnew == 0) ;

	    pnext = pdest ;

	    

	    

	    
	    Fsrc = Work->Flblock ;
	    Fdst = Work->Flblock ;
	    ASSERT (Fdst == F1 + nb*nb) ;
	    gap1 = dr - r ;
	    gap2 = drnew - r ;
	    ASSERT (gap1 >= 0) ;
	    for (j = 0 ; j < k ; j++)
	    {
		for (i = 0 ; i < r ; i++)
		{
		    *Fdst++ = *Fsrc++ ;
		}
		Fsrc += gap1 ;
		Fdst += gap2 ;
	    }
	    ASSERT (Fdst == F1 + nb*nb + drnew*k) ;
	    Fdst += drnew * (nb - k) ;

	    
	    Fsrc = Work->Fublock ;
	    ASSERT (Fdst == F1 + nb*nb + drnew*nb) ;
	    gap1 = dc - c ;
	    gap2 = dcnew - c ;
	    for (i = 0 ; i < k ; i++)
	    {
		for (j = 0 ; j < c ; j++)
		{
		    *Fdst++ = *Fsrc++ ;
		}
		Fsrc += gap1 ;
		Fdst += gap2 ;
	    }
	    ASSERT (Fdst == F1 + nb*nb + drnew*nb + dcnew*k) ;
	    Fdst += dcnew * (nb - k) ;

	    
	    Fsrc = Work->Fcblock ;
	    ASSERT (Fdst == F1 + nb*nb + drnew*nb + nb*dcnew) ;
	    gap1 = dr - r ;
	    gap2 = drnew - r ;
	    for (j = 0 ; j < c ; j++)
	    {
		for (i = 0 ; i < r ; i++)
		{
		    *Fdst++ = *Fsrc++ ;
		}
		Fsrc += gap1 ;
		Fdst += gap2 ;
	    }
	    ASSERT (Fdst == F1 + nb*nb + drnew*nb + nb*dcnew + drnew*c) ;

	    
	    if (do_Fcpos)
	    {
		Int *Fcols, *Fcpos ;
		Fcols = Work->Fcols ;
		Fcpos = Work->Fcpos ;
		for (j = 0 ; j < c ; j++)
		{
		    col = Fcols [j] ;
		    ASSERT (col >= 0 && col < Work->n_col) ;
		    ASSERT (Fcpos [col] == j * dr) ;
		    Fcpos [col] = j * drnew ;
		}
	    }

	    
	    Work->fnr_curr = drnew ;
	    Work->fnc_curr = dcnew ;
	    Work->fcurr_size = (drnew + nb) * (dcnew + nb) ;
	    size = UNITS (Entry, Work->fcurr_size) ;

	    
	    size = MAX (1, size) ;

	    
	    pnext->header.prevsize = size ;
	    pdest -= (size + 1) ;
	    F2 = (Entry *) (pdest + 1) ;

	    ASSERT ((unsigned Int) psrc + 1 + size <= (unsigned Int) pnext) ;
	    ASSERT (psrc <= pdest) ;
	    ASSERT (F1 <= F2) ;

	    
	    Fsrc = F1 + nb*nb + drnew*nb + nb*dcnew + drnew*c ;
	    Fdst = F2 + nb*nb + drnew*nb + nb*dcnew + drnew*c ;
	    gap = drnew - r ;
	    for (j = c-1 ; j >= 0 ; j--)
	    {
		Fsrc -= gap ;
		Fdst -= gap ;
		
		for (i = r-1 ; i >= 0 ; i--)
		{
		    *--Fdst = *--Fsrc ;
		}
	    }
	    ASSERT (Fsrc == F1 + nb*nb + drnew*nb + nb*dcnew) ;
	    ASSERT (Fdst == F2 + nb*nb + drnew*nb + nb*dcnew) ;

	    
	    Fsrc -= dcnew * (nb - k) ;
	    Fdst -= dcnew * (nb - k) ;
	    ASSERT (Fsrc == F1 + nb*nb + drnew*nb + dcnew*k) ;
	    ASSERT (Fdst == F2 + nb*nb + drnew*nb + dcnew*k) ;
	    gap = dcnew - c ;
	    for (i = k-1 ; i >= 0 ; i--)
	    {
		Fsrc -= gap ;
		Fdst -= gap ;
		for (j = c-1 ; j >= 0 ; j--)
		{
		    *--Fdst = *--Fsrc ;
		}
	    }
	    ASSERT (Fsrc == F1 + nb*nb + drnew*nb) ;
	    ASSERT (Fdst == F2 + nb*nb + drnew*nb) ;

	    
	    Fsrc -= drnew * (nb - k) ;
	    Fdst -= drnew * (nb - k) ;
	    ASSERT (Fsrc == F1 + nb*nb + drnew*k) ;
	    ASSERT (Fdst == F2 + nb*nb + drnew*k) ;
	    gap = drnew - r ;
	    for (j = k-1 ; j >= 0 ; j--)
	    {
		Fsrc -= gap ;
		Fdst -= gap ;
		for (i = r-1 ; i >= 0 ; i--)
		{
		    *--Fdst = *--Fsrc ;
		}
	    }
	    ASSERT (Fsrc == F1 + nb*nb) ;
	    ASSERT (Fdst == F2 + nb*nb) ;

	    
	    Fsrc -= nb * (nb - k) ;
	    Fdst -= nb * (nb - k) ;
	    ASSERT (Fsrc == F1 + nb*k) ;
	    ASSERT (Fdst == F2 + nb*k) ;
	    gap = nb - k ;
	    for (j = k-1 ; j >= 0 ; j--)
	    {
		Fsrc -= gap ;
		Fdst -= gap ;
		for (i = k-1 ; i >= 0 ; i--)
		{
		    *--Fdst = *--Fsrc ;
		}
	    }
	    ASSERT (Fsrc == F1) ;
	    ASSERT (Fdst == F2) ;

	    E [0] = (pdest + 1) - Numeric->Memory ;

	    Work->Flublock = (Entry *) (Numeric->Memory + E [0]) ;
	    ASSERT (Work->Flublock == F2) ;
	    Work->Flblock  = Work->Flublock + nb * nb ;
	    Work->Fublock  = Work->Flblock  + drnew * nb ;
	    Work->Fcblock  = Work->Fublock  + nb * dcnew ;

	    pdest->header.prevsize = 0 ;
	    pdest->header.size = size ;

	}
	else if (e > 0)
	{

	    
	    
	    

	    
	    
	    

	    p = psrc + 1 ;
	    GET_ELEMENT (epsrc, p, Cols, Rows, ncols, nrows, C) ;
	    nrowsleft = epsrc->nrowsleft ;
	    ncolsleft = epsrc->ncolsleft ;
	    cdeg = epsrc->cdeg ;
	    rdeg = epsrc->rdeg ;

	    
	    
	    

	    csize = nrowsleft * ncolsleft ;
	    size2 = UNITS (Element, 1)
		  + UNITS (Int, nrowsleft + ncolsleft)
		  + UNITS (Entry, csize) ;

	    pnext = pdest ;
	    pnext->header.prevsize = size2 ;
	    pdest -= (size2 + 1) ;

	    p = pdest + 1 ;
	    epdest = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Cols2 = (Int *) p ;
	    Rows2 = Cols2 + ncolsleft ;
	    p += UNITS (Int, nrowsleft + ncolsleft) ;
	    C2 = (Entry *) p ;

	    
	    
	    

	    

	    if (nrowsleft < nrows || ncolsleft < ncols)
	    {

		
		
		

		C1 = C ;
		C3 = C ;
		for (j = 0 ; j < ncols ; j++)
		{
		    if (Cols [j] >= 0)
		    {
			for (i = 0 ; i < nrows ; i++)
			{
			    if (Rows [i] >= 0)
			    {
				*C3++ = C1 [i] ;
			    }
			}
		    }
		    C1 += nrows ;
		}
	    }

	    
	    C += csize ;
	    C2 += csize ;
	    for (i = 0 ; i < csize ; i++)
	    {
		*--C2 = *--C ;
	    }

	    
	    
	    

	    i2 = nrowsleft ;
	    for (i = nrows - 1 ; i >= 0 ; i--)
	    {
		if (Rows [i] >= 0)
		{
		    Rows2 [--i2] = Rows [i] ;
		}
	    }

	    j2 = ncolsleft ;
	    for (j = ncols - 1 ; j >= 0 ; j--)
	    {
		if (Cols [j] >= 0)
		{
		    Cols2 [--j2] = Cols [j] ;
		}
	    }

	    
	    
	    

	    
	    E [e] = (pdest + 1) - Numeric->Memory ;
	    epdest = (Element *) (pdest + 1) ;

	    epdest->next = EMPTY ;	
	    epdest->ncols = ncolsleft ;
	    epdest->nrows = nrowsleft ;
	    epdest->ncolsleft = ncolsleft ;
	    epdest->nrowsleft = nrowsleft ;
	    epdest->rdeg = rdeg ;
	    epdest->cdeg = cdeg ;

	    pdest->header.prevsize = 0 ;
	    pdest->header.size = size2 ;

	}

    }

    
    
    

    Numeric->itail = pdest - Numeric->Memory ;
    pdest->header.prevsize = 0 ;
    Numeric->ibig = EMPTY ;
    Numeric->tail_usage = Numeric->size - Numeric->itail ;

    
    
    

    for (e = nel+1 ; e < Work->elen ; e++)
    {
	E [e] = 0 ;
    }

}



GLOBAL Int LU_get_memory
(
    NumericType *Numeric,
    WorkType *Work,
    Int needunits,
    Int r2,		
    Int c2,
    Int do_Fcpos
)
{
    double nsize, bsize, tsize ;
    Int i, minsize, newsize, newmem, costly, row, col, *Row_tlen, *Col_tlen,
	n_row, n_col, *Row_degree, *Col_degree ;
    Unit *mnew, *p ;

    
    
    

    n_row = Work->n_row ;
    n_col = Work->n_col ;
    Row_degree = Numeric->Rperm ;	
    Col_degree = Numeric->Cperm ;	
    Row_tlen   = Numeric->Uilen ;
    Col_tlen   = Numeric->Lilen ;

    
    
    

    for (row = 0 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row))
	{
	    Row_tlen [row] = 0 ;
	}
    }
    for (col = 0 ; col < n_col ; col++)
    {
	if (NON_PIVOTAL_COL (col))
	{
	    Col_tlen [col] = 0 ;
	}
    }

    
    
    

    nsize = (double) needunits + 2 ;
    needunits += LU_tuple_lengths (Numeric, Work, &tsize) ;
    nsize += tsize ;
    needunits += 2 ;	

    
    
    

    
    
    

    minsize = Numeric->size + needunits ;
    nsize += (double) Numeric->size ;

    bsize = ((double) Int_MAX) / sizeof (Unit) - 1 ;

    newsize = (Int) (LU_REALLOC_INCREASE * ((double) minsize)) ;
    nsize *= LU_REALLOC_INCREASE ;
    nsize += 1 ;

    if (newsize < 0 || nsize > bsize)
    {
	
	newsize = (Int) bsize ;	
    }
    else
    {
	newsize = MAX (newsize, minsize) ;
    }
    newsize = MAX (newsize, Numeric->size) ;

    
    
    Numeric->ibig = EMPTY ;

    
    
    

    mnew = (Unit *) NULL ;
    while (!mnew)
    {
	mnew = (Unit *) LU_realloc (Numeric->Memory, newsize, sizeof (Unit)) ;
	if (!mnew)
	{
	    if (newsize == minsize)	
	    {
		
		
		
		mnew = Numeric->Memory ;	
		newsize = Numeric->size ;
	    }
	    else
	    {
		
		newsize = (Int) (LU_REALLOC_REDUCTION * ((double) newsize)) ;
		newsize = MAX (minsize, newsize) ;
	    }
	}
    }

    
    costly = (mnew != Numeric->Memory) ;

    
    
    

    Numeric->Memory = mnew ;
    if (Work->E [0])
    {
	Int nb, dr, dc ;
	nb = Work->nb ;
	dr = Work->fnr_curr ;
	dc = Work->fnc_curr ;
	Work->Flublock = (Entry *) (Numeric->Memory + Work->E [0]) ;
	Work->Flblock  = Work->Flublock + nb * nb ;
	Work->Fublock  = Work->Flblock  + dr * nb ;
	Work->Fcblock  = Work->Fublock  + nb * dc ;
    }
    newmem = newsize - Numeric->size ;

    if (newmem >= 2)
    {
	

	
	p = Numeric->Memory + Numeric->size - 2 ;

	
	p->header.size = newmem - 1 ;
	i = Numeric->size - 1 ;
	p += newmem ;

	
	p->header.prevsize = newmem - 1 ;
	p->header.size = 1 ;

	Numeric->size = newsize ;

	
	LU_mem_free_tail_block (Numeric, i) ;

	Numeric->nrealloc++ ;

	if (costly)
	{
	    Numeric->ncostly++ ;
	}

    }
    
    
    

    LU_garbage_collection (Numeric, Work, r2, c2, do_Fcpos) ;

    
    
    

    return (LU_build_tuples (Numeric, Work)) ;
}







PRIVATE void zero_init_front (Int m, Int n, Entry *Fcblock, Int d)
{
    Int i, j ;
    Entry *F, *Fj = Fcblock ;
    for (j = 0 ; j < m ; j++)
    {
	F = Fj ;
	Fj += d ;
	for (i = 0 ; i < n ; i++)
	{
	    
	    CLEAR (*F) ;
	    F++ ;
	}
    }
}





GLOBAL Int LU_init_front
(
    NumericType *Numeric,
    WorkType *Work
)
{
    
    
    

    Int i, j, fnr_curr, row, col, *Frows, *Fcols,
	*Fcpos, *Frpos, fncols, fnrows, *Wrow, fnr2, fnc2, rrdeg, ccdeg, *Wm,
	fnrows_extended ;
    Entry *Fcblock, *Fl, *Wy, *Wx ;

    
    
    

    if (Work->do_grow)
    {
		fnr2 = LU_FRONTAL_GROWTH * Work->fnrows_new + 2 ;
		fnc2 = LU_FRONTAL_GROWTH * Work->fncols_new + 2 ;
		if (!LU_grow_front (Numeric, fnr2, fnc2, Work,
			Work->pivrow_in_front ? 2 : 0))
		{
			
			return (FALSE) ;
		}
    }
    fnr_curr = Work->fnr_curr ;

    
    
    

    

    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;

    Work->fnzeros = 0 ;

    ccdeg = Work->ccdeg ;
    rrdeg = Work->rrdeg ;

    fnrows = Work->fnrows ;
    fncols = Work->fncols ;

    
    

    
    
    

    Fl = Work->Flblock ;

    if (Work->pivcol_in_front)
    {
		
		
		Work->fscan_row = fnrows ;	
		Work->NewRows = Work->Wrp ;
		Wy = Work->Wy ;
		for (i = 0 ; i < fnrows ; i++)
		{
			Fl [i] = Wy [i] ;
		}
		fnrows_extended = fnrows + ccdeg ;
		for (i = fnrows ; i < fnrows_extended ; i++)
		{
			Fl [i] = Wy [i] ;
			
			row = Frows [i] ;
			Work->NewRows [i] = FLIP (row) ;
		}
		fnrows = fnrows_extended ;
    }
    else
    {
		
		Work->fscan_row = 0 ;			
		Work->NewRows = Frows ;
		Wm = Work->Wm ;
		Wx = Work->Wx ;
		for (i = 0 ; i < ccdeg ; i++)
		{
			Fl [i] = Wx [i] ;
			row = Wm [i] ;
			Frows [i] = row ;
			Frpos [row] = i ;
		}
		fnrows = ccdeg ;
    }

    Work->fnrows = fnrows ;

    
    
    

    Wrow = Work->Wrow ;
    if (Work->pivrow_in_front)
    {
		
		Work->fscan_col = fncols ;	
		Work->NewCols = Work->Wp ;
		
		ASSERT (IMPLIES (Work->pivcol_in_front, Wrow == Fcols)) ;
		if (Wrow == Fcols)
		{
			for (j = fncols ; j < rrdeg ; j++)
			{
			col = Wrow [j] ;
			
			
			Work->NewCols [j] = FLIP (col) ;
			Fcpos [col] = j * fnr_curr ;
			}
		}
		else
		{
			for (j = fncols ; j < rrdeg ; j++)
			{
				col = Wrow [j] ;
				Fcols [j] = col ;
				
				Work->NewCols [j] = FLIP (col) ;
				Fcpos [col] = j * fnr_curr ;
			}
		}
    }
    else
    {
		
		Work->fscan_col = 0 ;			
		Work->NewCols = Fcols ;
		for (j = 0 ; j < rrdeg ; j++)
		{
			col = Wrow [j] ;
			Fcols [j] = col ;
			Fcpos [col] = j * fnr_curr ;
		}
    }

    fncols = rrdeg ;
    Work->fncols = fncols ;

    
    
    

    Fcblock = Work->Fcblock ;

    zero_init_front (fncols, fnrows, Fcblock, fnr_curr) ;

    
    
    

    

    return (TRUE) ;

}



GLOBAL Int LU_grow_front
(
    NumericType *Numeric,
    Int fnr2,		
    Int fnc2,
    WorkType *Work,
    Int do_what		
)
{
    
    
    

    double s ;
    Entry *Fcold, *Fcnew ;
    Int j, i, col, *Fcpos, *Fcols, fnrows_max, fncols_max, fnr_curr, nb,
	fnrows_new, fncols_new, fnr_min, fnc_min, minsize,
	newsize, fnrows, fncols, *E, eloc ;

    
    
    

    Fcols = Work->Fcols ;
    Fcpos = Work->Fcpos ;
    E = Work->E ;

    
    
    

    
    nb = Work->nb ;
    fnrows_max = Work->fnrows_max + nb ;
    fncols_max = Work->fncols_max + nb ;

    

    
    fnrows_new = Work->fnrows_new + 1 ;
    fncols_new = Work->fncols_new + 1 ;
    if (fnrows_new % 2 == 0) fnrows_new++ ;
    fnrows_new += nb ;
    fncols_new += nb ;
    fnr_min = MIN (fnrows_new, fnrows_max) ;
    fnc_min = MIN (fncols_new, fncols_max) ;
    minsize = fnr_min * fnc_min ;
    if (INT_OVERFLOW ((double) fnr_min * (double) fnc_min * sizeof (Entry)))
    {
	
	return (FALSE) ;
    }

    
    fnr2 += nb ;
    fnc2 += nb ;
    if (fnr2 % 2 == 0) fnr2++ ;
    fnr2 = MAX (fnr2, fnr_min) ;
    fnc2 = MAX (fnc2, fnc_min) ;
    fnr2 = MIN (fnr2, fnrows_max) ;
    fnc2 = MIN (fnc2, fncols_max) ;

    s = ((double) fnr2) * ((double) fnc2) ;
    if (INT_OVERFLOW (s * sizeof (Entry)))
    {
	
	
	
	double a = 0.9 * sqrt ((Int_MAX / sizeof (Entry)) / s) ;
	fnr2 = MAX (fnr_min, a * fnr2) ;
	fnc2 = MAX (fnc_min, a * fnc2) ;
	
	newsize = fnr2 * fnc2 ;
	if (fnr2 % 2 == 0) fnr2++ ;
	fnc2 = newsize / fnr2 ;
    }

    fnr2 = MAX (fnr2, fnr_min) ;
    fnc2 = MAX (fnc2, fnc_min) ;
    newsize = fnr2 * fnc2 ;

    
    
    

    if (E [0] && do_what != 1)
    {
	
	LU_mem_free_tail_block (Numeric, E [0]) ;
	E [0] = 0 ;
	Work->Flublock = (Entry *) NULL ;
	Work->Flblock  = (Entry *) NULL ;
	Work->Fublock  = (Entry *) NULL ;
	Work->Fcblock  = (Entry *) NULL ;
    }

    
    
    

    eloc = LU_mem_alloc_tail_block (Numeric, UNITS (Entry, newsize)) ;

    if (!eloc)
    {
	
	if (!LU_get_memory (Numeric, Work, 1 + UNITS (Entry, newsize),
	    Work->fnrows, Work->fncols, FALSE))
	{
	    
	    return (FALSE) ;	
	}
	eloc = LU_mem_alloc_tail_block (Numeric, UNITS (Entry, newsize)) ;
    }

    
    while ((fnr2 != fnr_min || fnc2 != fnc_min) && !eloc)
    {
	fnr2 = MIN (fnr2 - 2, fnr2 * LU_REALLOC_REDUCTION) ;
	fnc2 = MIN (fnc2 - 2, fnc2 * LU_REALLOC_REDUCTION) ;
	fnr2 = MAX (fnr_min, fnr2) ;
	fnc2 = MAX (fnc_min, fnc2) ;
	if (fnr2 % 2 == 0) fnr2++ ;
	newsize = fnr2 * fnc2 ;
	eloc = LU_mem_alloc_tail_block (Numeric, UNITS (Entry, newsize)) ;
    }

    
    if (!eloc)
    {
	fnr2 = fnr_min ;
	fnc2 = fnc_min ;
	newsize = minsize ;
	eloc = LU_mem_alloc_tail_block (Numeric, UNITS (Entry, newsize)) ;
    }

    if (!eloc)
    {
	
	return (FALSE) ;
    }

    
    
    

    
    fnr_curr = Work->fnr_curr ;	    
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    Fcold = Work->Fcblock ;

    
    fnr2 -= nb ;
    fnc2 -= nb ;

    
    Work->Flublock = (Entry *) (Numeric->Memory + eloc) ;
    Work->Flblock  = Work->Flublock + nb * nb ;
    Work->Fublock  = Work->Flblock  + nb * fnr2 ;
    Work->Fcblock  = Work->Fublock  + nb * fnc2 ;
    Fcnew = Work->Fcblock ;

    if (E [0])
    {
	
	for (j = 0 ; j < fncols ; j++)
	{
	    col = Fcols [j] ;
	    for (i = 0 ; i < fnrows ; i++)
	    {
		Fcnew [i] = Fcold [i] ;
	    }
	    Fcnew += fnr2 ;
	    Fcold += fnr_curr ;
	    Fcpos [col] = j * fnr2 ;
	}
    }
    else if (do_what == 2)
    {
	
	for (j = 0 ; j < fncols ; j++)
	{
	    col = Fcols [j] ;
	    Fcpos [col] = j * fnr2 ;
	}
    }

    
    LU_mem_free_tail_block (Numeric, E [0]) ;

    
    
    

    E [0] = eloc ;
    Work->fnr_curr = fnr2 ;	    
    Work->fnc_curr = fnc2 ;
    Work->fcurr_size = newsize ;    
    Work->do_grow = FALSE ;	    
    return (TRUE) ;
}



GLOBAL Int LU_start_front    
(
    Int chain,
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    double maxbytes ;
    Int fnrows_max, fncols_max, fnr2, fnc2, fsize, fcurr_size, maxfrsize,
	overflow, nb, f, cdeg ;

    nb = Symbolic->nb ;
    fnrows_max = Symbolic->Chain_maxrows [chain] ;
    fncols_max = Symbolic->Chain_maxcols [chain] ;

    Work->fnrows_max = fnrows_max ;
    Work->fncols_max = fncols_max ;
    Work->any_skip = FALSE ;

    maxbytes = sizeof (Entry) *
	(double) (fnrows_max + nb) * (double) (fncols_max + nb) ;
    fcurr_size = Work->fcurr_size ;

    if (Symbolic->prefer_diagonal)
    {
		
		Int col, tpi, e, *E, *Col_tuples, *Col_tlen, *Cols ;
		Tuple *tp, *tpend ;
		Unit *Memory, *p ;
		Element *ep ;
		E = Work->E ;
		Memory = Numeric->Memory ;
		Col_tuples = Numeric->Lip ;
		Col_tlen = Numeric->Lilen ;
		col = Work->nextcand ;
		tpi = Col_tuples [col] ;
		tp = (Tuple *) Memory + tpi ;
		tpend = tp + Col_tlen [col] ;
		cdeg = 0 ;
		for ( ; tp < tpend ; tp++)
		{
			e = tp->e ;
			if (!E [e]) continue ;
			f = tp->f ;
			p = Memory + E [e] ;
			ep = (Element *) p ;
			p += UNITS (Element, 1) ;
			Cols = (Int *) p ;
			if (Cols [f] == EMPTY) continue ;
			cdeg += ep->nrowsleft ;
		}

		

		
		if (Symbolic->amd_dmax > 0)
		{
			cdeg = MIN (cdeg, Symbolic->amd_dmax) ;
		}

		
		cdeg += 2 ;

		
		cdeg = MIN (cdeg, fnrows_max) ;

    }
    else
    {
		
		cdeg = 0 ;
    }

    

    

    

    
    overflow = INT_OVERFLOW (maxbytes) ;
    if (overflow)
    {
		
		maxfrsize = Int_MAX / sizeof (Entry) ;
    }
    else
    {
		maxfrsize = (fnrows_max + nb) * (fncols_max + nb) ;
    }

    if (Numeric->front_alloc_init < 0)
    {
		
		fsize = -Numeric->front_alloc_init ;
		fsize = MAX (1, fsize) ;
    }
    else
    {
		if (INT_OVERFLOW (Numeric->front_alloc_init * maxbytes))
		{
			
			fsize = Int_MAX / sizeof (Entry) ;
		}
		else
		{
			fsize = Numeric->front_alloc_init * maxfrsize ;
		}

		if (cdeg > 0)
		{
			
			Int fsize2 ;

			
			cdeg += nb ;

			if (INT_OVERFLOW (((double) cdeg * (double) cdeg) * sizeof (Entry)))
			{
				
				fsize2 = Int_MAX / sizeof (Entry) ;
			}
			else
			{
				fsize2 = MAX (cdeg * cdeg, fcurr_size) ;
			}
			fsize = MIN (fsize, fsize2) ;
		}
    }

    fsize = MAX (fsize, 2*nb*nb) ;

    
    ASSERT (!INT_OVERFLOW ((double) fsize * sizeof (Entry))) ;

    Work->fnrows_new = 0 ;
    Work->fncols_new = 0 ;

    

    if (fsize >= maxfrsize && !overflow)
    {
		
		fnr2 = fnrows_max + nb ;
		fnc2 = fncols_max + nb ;
		fsize = maxfrsize ;
    }
    else
    {
		
		if (fnrows_max <= fncols_max)
		{
			fnr2 = (Int) sqrt ((double) fsize) ;
			
			fnr2 = MAX (fnr2, 1) ;
			if (fnr2 % 2 == 0) fnr2++ ;
			fnr2 = MIN (fnr2, fnrows_max + nb) ;
			fnc2 = fsize / fnr2 ;
		}
		else
		{
			fnc2 = (Int) sqrt ((double) fsize) ;
			fnc2 = MIN (fnc2, fncols_max + nb) ;
			fnr2 = fsize / fnc2 ;
			
			fnr2 = MAX (fnr2, 1) ;
			if (fnr2 % 2 == 0)
			{
				fnr2++ ;
				fnc2 = fsize / fnr2 ;
			}
		}
    }
    fnr2 = MIN (fnr2, fnrows_max + nb) ;
    fnc2 = MIN (fnc2, fncols_max + nb) ;

    fnr2 -= nb ;
    fnc2 -= nb ;

    if (fsize > fcurr_size)
    {
		Work->do_grow = TRUE ;
		if (!LU_grow_front (Numeric, fnr2, fnc2, Work, -1))
		{
			
			return (FALSE) ;
		}
    }
    else
    {
		
		Work->fnr_curr = fnr2 ;
		Work->fnc_curr = fnc2 ;
		Work->Flblock  = Work->Flublock + nb * nb ;
		Work->Fublock  = Work->Flblock  + nb * fnr2 ;
		Work->Fcblock  = Work->Fublock  + nb * fnc2 ;
    }

    return (TRUE) ;
}








#define RELAX1 0.25
#define SYM_RELAX1 0.0
#define RELAX2 0.1
#define RELAX3 0.125







PRIVATE void remove_candidate (Int jj, WorkType *Work, SymbolicType *Symbolic)
{

    if (Symbolic->fixQ)
    {
	
	if (Work->ncand > 1)
	{
	    Work->Candidates [0] = Work->nextcand++ ;
	}
	else
	{
	    Work->nCandidates = 0 ;
	}
    }
    else
    {
	
	if (Work->ncand > MAX_CANDIDATES)
	{
	    Work->Candidates [jj] = Work->nextcand++ ;
	}
	else
	{
	    ASSERT (Work->nCandidates == Work->ncand) ;
	    Work->Candidates [jj] = Work->Candidates [Work->ncand - 1] ;
	    Work->Candidates [Work->ncand - 1] = EMPTY ;
	    Work->nCandidates-- ;
	}
    }
    Work->ncand-- ;
}





GLOBAL Int LU_local_search
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    
    
    

    double relax1 ;
    Entry *Flblock, *Fublock, *Fs, *Fcblock, *C, *Wx, *Wy, *Fu, *Flublock,
	*Flu ;
    Int pos, nrows, *Cols, *Rows, e, f, status, max_cdeg, fnzeros, nb, j, col,
	i, row, cdeg_in, rdeg [2][2], fnpiv, nothing [2], new_LUsize,
	pivrow [2][2], pivcol [2], *Wp, *Fcpos, *Frpos, new_fnzeros, cdeg_out,
	*Wm, *Wio, *Woi, *Woo, *Frows, *Fcols, fnrows, fncols, *E, deg, nr_in,
	nc, thiscost, bestcost, nr_out, do_update, extra_cols, extra_rows,
	extra_zeros, relaxed_front, do_extend, fnr_curr, fnc_curr, tpi,
	*Col_tuples, *Col_degree, *Col_tlen, jj, jcand [2], freebie [2],
	did_rowmerge, fnrows_new [2][2], fncols_new [2][2], search_pivcol_out,
	*Diagonal_map, *Diagonal_imap, row2, col2 ;
    Unit *Memory, *p ;
    Tuple *tp, *tpend, *tp1, *tp2 ;
    Element *ep ;

    
    
    

    Memory = Numeric->Memory ;
    E = Work->E ;
    Col_degree = Numeric->Cperm ;

    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;

    Wx = Work->Wx ;
    Wy = Work->Wy ;
    Wp = Work->Wp ;
    Wm = Work->Wm ;
    Woi = Work->Woi ;
    Wio = Work->Wio ;
    Woo = Work->Woo ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    nb = Work->nb ;
    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;
    fnpiv = Work->fnpiv ;
    nothing [0] = EMPTY ;
    nothing [1] = EMPTY ;
    relax1 = (Symbolic->prefer_diagonal) ? SYM_RELAX1 : RELAX1 ;
    fnzeros = Work->fnzeros ;
    new_fnzeros = fnzeros ;
    jj = EMPTY ;

    Fcblock = Work->Fcblock ;	    
    Flblock = Work->Flblock ;	    
    Fublock = Work->Fublock ;	    
    Flublock = Work->Flublock ;	    

    
    max_cdeg = Work->fnrows_max ;
    ASSERT (Work->fnrows_max <= Symbolic->maxnrows) ;
    ASSERT (Work->fncols_max <= Symbolic->maxncols) ;

    if (fnrows == 0 && fncols == 0)
    {
	
	Work->firstsuper = Work->ksuper ;
    }

    
    
    

    

    pivcol [IN] = EMPTY ;
    pivcol [OUT] = EMPTY ;

    cdeg_in = Int_MAX ;
    cdeg_out = Int_MAX ;

    pivrow [IN][IN] = EMPTY ;
    pivrow [IN][OUT] = EMPTY ;
    pivrow [OUT][IN] = EMPTY ;
    pivrow [OUT][OUT] = EMPTY ;

    rdeg [IN][IN] = Int_MAX ;
    rdeg [IN][OUT] = Int_MAX ;
    rdeg [OUT][IN] = Int_MAX ;
    rdeg [OUT][OUT] = Int_MAX ;

    freebie [IN] = FALSE ;
    freebie [OUT] = FALSE ;

    Work->pivot_case = EMPTY ;
    bestcost = EMPTY ;

    nr_out = EMPTY ;
    nr_in = EMPTY ;

    jcand [IN] = EMPTY ;
    jcand [OUT] = EMPTY ;

    fnrows_new [IN][IN] = EMPTY ;
    fnrows_new [IN][OUT] = EMPTY ;
    fnrows_new [OUT][IN] = EMPTY ;
    fnrows_new [OUT][OUT] = EMPTY ;

    fncols_new [IN][IN] = EMPTY ;
    fncols_new [IN][OUT] = EMPTY ;
    fncols_new [OUT][IN] = EMPTY ;
    fncols_new [OUT][OUT] = EMPTY ;

    
    
    
    

    
    
    

    
    
    

    col = Work->Candidates [0] ;

    
    deg = Symbolic->fixQ ? EMPTY : Col_degree [col] ;

    if (Fcpos [col] >= 0)
    {
	
	pivcol [IN] = col ;
	cdeg_in = deg ;		
	jcand [IN] = 0 ;
    }
    else
    {
	
	pivcol [OUT] = col ;
	cdeg_out = deg ;	
	jcand [OUT] = 0 ;
    }

    
    for (j = 1 ; j < Work->nCandidates ; j++)
    {
	col = Work->Candidates [j] ;

	deg = Col_degree [col] ;
	if (Fcpos [col] >= 0)
	{
	    if (deg < cdeg_in || (deg == cdeg_in && col < pivcol [IN]))
	    {
		
		pivcol [IN] = col ;
		cdeg_in = deg ;
		jcand [IN] = j ;
	    }
	}
	else
	{
	    if (deg < cdeg_out || (deg == cdeg_out && col < pivcol [OUT]))
	    {
		
		pivcol [OUT] = col ;
		cdeg_out = deg ;
		jcand [OUT] = j ;
	    }
	}
    }

    cdeg_in = EMPTY ;
    cdeg_out = EMPTY ;

    
    
    

    if (pivcol [IN] != EMPTY)
    {

	
	Fs  = Fcblock  +  Fcpos [pivcol [IN]] ;
	Fu  = Fublock  + (Fcpos [pivcol [IN]] / fnr_curr) ;
	Flu = Flublock + fnpiv * nb ;

	
	
	

	
	for (i = 0 ; i < fnpiv ; i++)
	{
	    Flu [i] = Fu [i*fnc_curr] ;
	}

	
	
	

	

	if (fnpiv > 1)
	{
	    
	    BLAS_TRSV (fnpiv, Flublock, Flu, nb) ;
	}

	
	
	

	for (i = 0 ; i < fnrows ; i++)
	{
	    Wy [i] = Fs [i] ;
	}

	
	
	

	

	
	BLAS_GEMV (fnrows, fnpiv, Flblock, Flu, Wy, fnr_curr) ;

	

	
	
	

	cdeg_in = fnrows ;

	

	ASSERT (pivcol [IN] >= 0 && pivcol [IN] < n_col) ;
	ASSERT (NON_PIVOTAL_COL (pivcol [IN])) ;

	tpi = Col_tuples [pivcol [IN]] ;
	if (tpi)
	{
	    tp = (Tuple *) (Memory + tpi) ;
	    tp1 = tp ;
	    tp2 = tp ;
	    tpend = tp + Col_tlen [pivcol [IN]] ;
	    for ( ; tp < tpend ; tp++)
	    {
		e = tp->e ;
		ASSERT (e > 0 && e <= Work->nel) ;
		if (!E [e]) continue ;	
		f = tp->f ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		p += UNITS (Element, 1) ;
		Cols = (Int *) p ;
		if (Cols [f] == EMPTY) continue ; 
		ASSERT (pivcol [IN] == Cols [f]) ;

		Rows = Cols + ep->ncols ;
		nrows = ep->nrows ;
		p += UNITS (Int, ep->ncols + nrows) ;
		C = ((Entry *) p) + f * nrows ;

		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0) 
		    {
			pos = Frpos [row] ;
			if (pos < 0)
			{
			    
			    if (cdeg_in >= max_cdeg)
			    {
				
				return (SparseLU_ERROR_different_pattern) ;
			    }
			    Frpos [row] = cdeg_in ;
			    Frows [cdeg_in] = row ;
			    Wy [cdeg_in++] = C [i] ;
			}
			else
			{
			    
			    
			    ASSERT (pos < max_cdeg) ;
			    ASSEMBLE (Wy [pos], C [i]) ;
			}
		    }
		}
		*tp2++ = *tp ;	
	    }
	    Col_tlen [pivcol [IN]] = tp2 - tp1 ;
	}

	

	
	
	

	nr_in = cdeg_in - fnrows ;

	
	
	ASSERT (cdeg_in > 0) ;

	
	

	
	
	

	
	
	status = LU_row_search (Numeric, Work, Symbolic,
	    fnrows, cdeg_in, Frows, Frpos,   
	    pivrow [IN], rdeg [IN], Fcols, Wio, nothing, Wy,
	    pivcol [IN], freebie) ;

	
	
	

	if (status == SparseLU_ERROR_different_pattern)
	{
	    
	    return (SparseLU_ERROR_different_pattern) ;
	}

	
	
	

	
	
	

	ASSERT (status != SparseLU_WARNING_singular_matrix) ;

	
	
	

	if (pivrow [IN][IN] != EMPTY)
	{
	    

	    

	    ASSERT (fnrows >= 0 && fncols >= 0) ;

	    ASSERT (cdeg_in > 0) ;
	    nc = rdeg [IN][IN] - fncols ;

	    thiscost =
		
		(nr_in * (fncols - 1)) +
		
		(nc * (cdeg_in - 1)) ;

	    

	    
	    fnrows_new [IN][IN] = (fnrows-1) + nr_in ;
	    fncols_new [IN][IN] = (fncols-1) + nc ;
	    

	    do_extend = TRUE ;

	    
	    
	    fnzeros = Work->fnzeros + fnpiv * (nr_in + nc) ;

	    new_LUsize = (fnpiv+1) * (fnrows + nr_in + fncols + nc) ;

	    

	    
	    do_update = fnpiv > 0 &&
		(((double) fnzeros) / ((double) new_LUsize)) > RELAX2 ;

	    
	    Work->pivot_case = IN_IN ;
	    bestcost = thiscost ;

	    
	    Work->do_extend = do_extend ;
	    Work->do_update = do_update ;
	    new_fnzeros = fnzeros ;

	}

	
	
	

	if (pivrow [IN][OUT] != EMPTY)
	{
	    

	    
	    

	    
	    nc = 0 ;
	    for (i = 0 ; i < rdeg [IN][OUT] ; i++)
	    {
		col = Wio [i] ;
		if (Fcpos [col] < 0) nc++ ;
	    }
	    thiscost =
		
		(nc * fnrows) +
		
		((nr_in-1) * (rdeg [IN][OUT]-1)) ;

	    

	    extra_cols = ((fncols-1) + nc ) - (rdeg [IN][OUT] - 1) ;

	    extra_zeros = (nr_in-1) * extra_cols ;	

	    
	    fnrows_new [IN][OUT] = fnrows + (nr_in-1) ;
	    fncols_new [IN][OUT] = (fncols-1) + nc ;
	    relaxed_front = fnrows_new [IN][OUT] * fncols_new [IN][OUT] ;

	    
	    
	    
	    

	    
	    if (extra_zeros == 0)
	    {
		do_extend = TRUE ;
	    }
	    else
	    {
		do_extend = ((double) extra_zeros) <
		   (relax1 * (double) relaxed_front) ;
	    }

	    if (do_extend)
	    {
		
		thiscost += extra_zeros ;

		
		
		fnzeros = Work->fnzeros + fnpiv * (nr_in + nc) ;

		new_LUsize = (fnpiv+1) * (fnrows + nr_in + fncols + nc) ;


		
		do_update = fnpiv > 0 &&
		    (((double) fnzeros) / ((double) new_LUsize)) > RELAX3 ;
	    }
	    else
	    {
		
		do_update = fnpiv > 0 ;
		fnzeros = 0 ;

		
		fnrows_new [IN][OUT] = cdeg_in - 1 ;
		fncols_new [IN][OUT] = rdeg [IN][OUT] - 1 ;

	    }

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		
		Work->pivot_case = IN_OUT ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}
    }

    
    
    

    search_pivcol_out = (bestcost != 0 && pivcol [OUT] != EMPTY) ;
    if (Symbolic->prefer_diagonal)
    {
	search_pivcol_out = search_pivcol_out && (pivrow [IN][IN] == EMPTY) ;
    }

    if (search_pivcol_out)
    {

	


	
	
	

	
	
	
	

	cdeg_out = 0 ;

	tpi = Col_tuples [pivcol [OUT]] ;
	if (tpi)
	{
	    tp = (Tuple *) (Memory + tpi) ;
	    tp1 = tp ;
	    tp2 = tp ;
	    tpend = tp + Col_tlen [pivcol [OUT]] ;
	    for ( ; tp < tpend ; tp++)
	    {
		e = tp->e ;
		ASSERT (e > 0 && e <= Work->nel) ;
		if (!E [e]) continue ;	
		f = tp->f ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		p += UNITS (Element, 1) ;
		Cols = (Int *) p ;
		if (Cols [f] == EMPTY) continue ; 
		ASSERT (pivcol [OUT] == Cols [f]) ;

		Rows = Cols + ep->ncols ;
		nrows = ep->nrows ;
		p += UNITS (Int, ep->ncols + nrows) ;
		C = ((Entry *) p) + f * nrows ;

		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0) 
		    {
			pos = Wp [row] ;
			if (pos < 0)
			{
			    
			    if (cdeg_out >= max_cdeg)
			    {
				
				return (SparseLU_ERROR_different_pattern) ;
			    }
			    Wp [row] = cdeg_out ;
			    Wm [cdeg_out] = row ;
			    Wx [cdeg_out++] = C [i] ;
			}
			else
			{
			    
			    
			    ASSEMBLE (Wx [pos], C [i]) ;
			}
		    }
		}
		*tp2++ = *tp ;	
	    }
	    Col_tlen [pivcol [OUT]] = tp2 - tp1 ;
	}

	

	
	
	

	
	status = LU_row_search (Numeric, Work, Symbolic,
	    0, cdeg_out, Wm, Wp, 
	    pivrow [OUT], rdeg [OUT], Woi, Woo, pivrow [IN], Wx,
	    pivcol [OUT], freebie) ;

	
	
	

	if (status == SparseLU_ERROR_different_pattern)
	{
	    
	    return (SparseLU_ERROR_different_pattern) ;
	}

	
	
	

	for (i = 0 ; i < cdeg_out ; i++)
	{
	    Wp [Wm [i]] = EMPTY ;	
	}

	
	
	

	if (status == SparseLU_WARNING_singular_matrix)
	{
	    

	    
	    remove_candidate (jcand [OUT], Work, Symbolic) ;
	    Work->ndiscard++ ;

	    
	    Col_tlen [pivcol [OUT]] = 0 ;
	    LU_mem_free_tail_block (Numeric, Col_tuples [pivcol [OUT]]) ;
	    Col_tuples [pivcol [OUT]] = 0 ;

	    
	    return (SparseLU_WARNING_singular_matrix) ;
	}

	

	if (freebie [IN])
	{
	    
	    Woi = Fcols ;
	    rdeg [OUT][IN] = rdeg [IN][IN] ;
	}

	if (freebie [OUT])
	{
	    
	    Woo = Wio ;
	    rdeg [OUT][OUT] = rdeg [IN][OUT] ;
	}

	
	
	

	if (pivrow [OUT][IN] != EMPTY)
	{
	    

	    ASSERT (fnrows > 0 && fncols >= 0) ;

	    did_rowmerge = (cdeg_out == 0) ;
	    if (did_rowmerge)
	    {
		
		
		Wm [0] = pivrow [OUT][IN] ;
		CLEAR (Wx [0]) ;
		cdeg_out = 1 ;
	    }

	    nc = rdeg [OUT][IN] - fncols ;

	    
	    nr_out = 0 ;
	    for (i = 0 ; i < cdeg_out ; i++)
	    {
		row = Wm [i] ;
		if (Frpos [row] < 0 || Frpos [row] >= fnrows) nr_out++ ;
	    }

	    thiscost =
		
		(nr_out * fncols) +
		
		((nc-1) * (cdeg_out-1)) ;

	    

	    extra_rows = ((fnrows-1) + nr_out) - (cdeg_out - 1) ;
	    extra_zeros = (nc-1) * extra_rows ;	

	    
	    fnrows_new [OUT][IN] = (fnrows-1) + nr_out ;
	    fncols_new [OUT][IN] = fncols + (nc-1) ;
	    relaxed_front = fnrows_new [OUT][IN] * fncols_new [OUT][IN] ;

	    
	    
	    
	    
	    if (did_rowmerge)
	    {
		do_extend = FALSE ;
	    }
	    else
	    {
		
		if (extra_zeros == 0)
		{
		    do_extend = TRUE ;
		}
		else
		{
		    do_extend = ((double) extra_zeros) <
		       (relax1 * (double) relaxed_front) ;
		}
	    }

	    if (do_extend)
	    {
		
		thiscost += extra_zeros ;

		
		
		fnzeros = Work->fnzeros + fnpiv * (nr_out + nc) ;

		new_LUsize = (fnpiv+1) * (fnrows + nr_out + fncols + nc) ;

		
		do_update = fnpiv > 0 &&
		    (((double) fnzeros) / ((double) new_LUsize)) > RELAX3 ;
	    }
	    else
	    {
		
		do_update = fnpiv > 0 ;
		fnzeros = 0 ;
		
		fnrows_new [OUT][IN] = cdeg_out - 1 ;
		fncols_new [OUT][IN] = rdeg [OUT][IN] - 1 ;
	    }

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		
		Work->pivot_case = OUT_IN ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}

	
	
	

	if (pivrow [OUT][OUT] != EMPTY)
	{
	    

	    did_rowmerge = (cdeg_out == 0) ;
	    if (did_rowmerge)
	    {
		
		
		Wm [0] = pivrow [OUT][OUT] ;
		CLEAR (Wx [0]) ;
		cdeg_out = 1 ;
		nr_out = 1 ;
	    }

	    if (fnrows == 0 && fncols == 0)
	    {
		
		nc = rdeg [OUT][OUT] ;
		extra_cols = 0 ;
		nr_out = cdeg_out ;
		extra_rows = 0 ;
		extra_zeros = 0 ;

		thiscost = (nc-1) * (cdeg_out-1) ; 

		
		fnrows_new [OUT][OUT] = nr_out-1 ;
		fncols_new [OUT][OUT] = nc-1 ;
		relaxed_front = fnrows_new [OUT][OUT] * fncols_new [OUT][OUT] ;
	    }
	    else
	    {

		
		if (nr_out == EMPTY)
		{
		    nr_out = 0 ;
		    for (i = 0 ; i < cdeg_out ; i++)
		    {
			row = Wm [i] ;
			if (Frpos [row] < 0 || Frpos [row] >= fnrows) nr_out++ ;
		    }
		}

		
		nc = 0 ;
		for (i = 0 ; i < rdeg [OUT][OUT] ; i++)
		{
		    col = Woo [i] ;
		    if (Fcpos [col] < 0) nc++ ;
		}

		extra_cols = (fncols + (nc-1)) - (rdeg [OUT][OUT] - 1) ;
		extra_rows = (fnrows + (nr_out-1)) - (cdeg_out - 1) ;
		extra_zeros = ((nc-1) * extra_rows) + ((nr_out-1) * extra_cols);

		thiscost =
		    
		    ((nc-1) * (cdeg_out-1)) +
		    
		    ((nr_out-1) * (fncols - extra_cols)) ;

		
		fnrows_new [OUT][OUT] = fnrows + (nr_out-1) ;
		fncols_new [OUT][OUT] = fncols + (nc-1) ;
		relaxed_front = fnrows_new [OUT][OUT] * fncols_new [OUT][OUT] ;

	    }

	    
	    
	    
	    
	    if (did_rowmerge)
	    {
		do_extend = FALSE ;
	    }
	    else
	    {
		
		if (extra_zeros == 0)
		{
		    do_extend = TRUE ;
		}
		else
		{
		    do_extend = ((double) extra_zeros) <
		       (relax1 * (double) relaxed_front) ;
		}
	    }

	    if (do_extend)
	    {
		
		thiscost += extra_zeros ;

		
		
		fnzeros = Work->fnzeros + fnpiv * (nr_out + nc) ;

		new_LUsize = (fnpiv+1) * (fnrows + nr_out + fncols + nc) ;

		
		do_update = fnpiv > 0 &&
		    (((double) fnzeros) / ((double) new_LUsize)) > RELAX3 ;
	    }
	    else
	    {
		
		do_update = fnpiv > 0 ;
		fnzeros = 0 ;

		
		fnrows_new [OUT][OUT] = cdeg_out - 1 ;
		fncols_new [OUT][OUT] = rdeg [OUT][OUT] - 1 ;
	    }

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		
		Work->pivot_case = OUT_OUT ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}
    }

    
    

    
    
    

    Work->fnzeros = new_fnzeros ;

    
    
    

    switch (Work->pivot_case)
    {
	case IN_IN:
	    Work->pivcol_in_front = TRUE ;
	    Work->pivrow_in_front = TRUE ;
	    Work->pivcol = pivcol [IN] ;
	    Work->pivrow = pivrow [IN][IN] ;
	    Work->ccdeg = nr_in ;
	    Work->Wrow = Fcols ;
	    Work->rrdeg = rdeg [IN][IN] ;
	    jj = jcand [IN] ;
	    Work->fnrows_new = fnrows_new [IN][IN] ;
	    Work->fncols_new = fncols_new [IN][IN] ;
	    break ;

	case IN_OUT:
	    Work->pivcol_in_front = TRUE ;
	    Work->pivrow_in_front = FALSE ;
	    Work->pivcol = pivcol [IN] ;
	    Work->pivrow = pivrow [IN][OUT] ;
	    Work->ccdeg = nr_in ;
	    Work->Wrow = Wio ;
	    Work->rrdeg = rdeg [IN][OUT] ;
	    jj = jcand [IN] ;
	    Work->fnrows_new = fnrows_new [IN][OUT] ;
	    Work->fncols_new = fncols_new [IN][OUT] ;
	    break ;

	case OUT_IN:
	    Work->pivcol_in_front = FALSE ;
	    Work->pivrow_in_front = TRUE ;
	    Work->pivcol = pivcol [OUT] ;
	    Work->pivrow = pivrow [OUT][IN] ;
	    Work->ccdeg = cdeg_out ;
	    
	    Work->Wrow = Woi ;
	    Work->rrdeg = rdeg [OUT][IN] ;
	    
	    jj = jcand [OUT] ;
	    Work->fnrows_new = fnrows_new [OUT][IN] ;
	    Work->fncols_new = fncols_new [OUT][IN] ;
	    break ;

	case OUT_OUT:
	    Work->pivcol_in_front = FALSE ;
	    Work->pivrow_in_front = FALSE ;
	    Work->pivcol = pivcol [OUT] ;
	    Work->pivrow = pivrow [OUT][OUT] ;
	    Work->ccdeg = cdeg_out ;
	    
	    Work->Wrow = Woo ;
	    Work->rrdeg = rdeg [OUT][OUT] ;
	    jj = jcand [OUT] ;
	    Work->fnrows_new = fnrows_new [OUT][OUT] ;
	    Work->fncols_new = fncols_new [OUT][OUT] ;
	    break ;

    }

    if (!Work->pivcol_in_front && pivcol [IN] != EMPTY)
    {
	
	for (i = fnrows ; i < cdeg_in ; i++)
	{
	    Frpos [Frows [i]] = EMPTY;
	}
    }

    
    
    

    
    
    
    remove_candidate (jj, Work, Symbolic) ;

    
    
    

    Work->do_grow = (Work->fnrows_new + 1 > fnr_curr
		  || Work->fncols_new + 1 > fnc_curr) ;
    if (Work->do_grow)
    {
	
	if (!Work->do_update && fnpiv > 0)
	{
	    
	    Work->nforced++ ;
	    Work->do_update = TRUE ;
	}
    }

    
    
    

    

    
    
    

    

    
    
    

    if (Work->do_extend)
    {
	Work->do_scan2row = (fncols > 0) ;
	Work->do_scan2col = (fnrows > 0) ;
    }
    else
    {
	Work->do_scan2row = (fncols > 0) && Work->pivrow_in_front ;
	Work->do_scan2col = (fnrows > 0) && Work->pivcol_in_front ;
    }

    

    
    
    

    if (Symbolic->prefer_diagonal)
    {
	Diagonal_map = Work->Diagonal_map ;
	Diagonal_imap = Work->Diagonal_imap ;

	row2 = Diagonal_map  [Work->pivcol] ;
	col2 = Diagonal_imap [Work->pivrow] ;

	if (row2 < 0)
	{
	    
	    Work->noff_diagonal++ ;
	    row2 = UNFLIP (row2) ;
	}

	if (row2 != Work->pivrow)
	{
	    

	    Diagonal_map  [Work->pivcol] = FLIP (Work->pivrow) ;
	    Diagonal_imap [Work->pivrow] = Work->pivcol ;

	    Diagonal_map  [col2] = FLIP (row2) ;
	    Diagonal_imap [row2] = col2 ;

	}
    }

    return (SparseLU_OK) ;
}



GLOBAL Int LU_row_search
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic,
    Int cdeg0,			
    Int cdeg1,			
    const Int Pattern [ ],	
    const Int Pos [ ],		
    Int pivrow [2],		
    Int rdeg [2],		
    Int W_i [ ],		
				
    Int W_o [ ],		
				
    Int prior_pivrow [2],	
    const Entry Wxy [ ],	

    Int pivcol,			
    Int freebie [ ]
)
{

    
    
    

    double maxval, toler, toler2, value, pivot [2] ;
    Int i, row, deg, col, *Frpos, fnrows, *E, j, ncols, *Cols, *Rows,
	e, f, Wrpflag, *Fcpos, fncols, tpi, max_rdeg, nans_in_col, was_offdiag,
	diag_row, prefer_diagonal, *Wrp, found, *Diagonal_map ;
    Tuple *tp, *tpend, *tp1, *tp2 ;
    Unit *Memory, *p ;
    Element *ep ;
    Int *Row_tuples, *Row_degree, *Row_tlen ;

    pivot [IN] = 0. ;
    pivot [OUT] = 0. ;

    
    
    

    Row_degree = Numeric->Rperm ;
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Wrp = Work->Wrp ;
    Frpos = Work->Frpos ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    fnrows = Work->fnrows ;

    prefer_diagonal = Symbolic->prefer_diagonal ;
    Diagonal_map = Work->Diagonal_map ;

    if (Diagonal_map)
    {
	diag_row = Diagonal_map [pivcol] ;
	was_offdiag = diag_row < 0 ;
	if (was_offdiag)
	{
	    
	    diag_row = FLIP (diag_row) ;
	}
	ASSERT (diag_row >= 0 && diag_row < Numeric->n_row) ;
    }
    else
    {
	diag_row = EMPTY ;	
	was_offdiag = EMPTY ;	
    }

    
    max_rdeg = Work->fncols_max ;

    
    
    

    maxval = 0.0 ;
    nans_in_col = FALSE ;

    for (i = 0 ; i < cdeg1 ; i++)
    {
	APPROX_ABS (value, Wxy [i]) ;
	if (SCALAR_IS_NAN (value))
	{
	    nans_in_col = TRUE ;
	    maxval = value ;
	    break ;
	}
	
	maxval = MAX (maxval, value) ;
    }

    

    toler = Numeric->relpt * maxval ;
    toler2 = Numeric->relpt2 * maxval ;
    toler2 = was_offdiag ? toler : toler2 ;

    if (SCALAR_IS_NAN (toler) || SCALAR_IS_NAN (toler2))
    {
	nans_in_col = TRUE ;
    }

    if (!nans_in_col)
    {

	
	found = FALSE ;

	if (prefer_diagonal)
	{
	    i = Pos [diag_row] ;
	    if (i >= 0)
	    {
		double a ;

		APPROX_ABS (a, Wxy [i]) ;

		if (SCALAR_IS_NONZERO (a) && a >= toler2)
		{
		    
		    found = TRUE ;
		    if (Frpos [diag_row] >= 0 && Frpos [diag_row] < fnrows)
		    {
			pivrow [IN] = diag_row ;
			pivrow [OUT] = EMPTY ;
		    }
		    else
		    {
			pivrow [IN] = EMPTY ;
			pivrow [OUT] = diag_row ;
		    }
		}
	    }
	}

	
	if (!found)
	{
	    if (cdeg0 > 0)
	    {

		
		for (i = 0 ; i < cdeg0 ; i++)
		{
		    double a ;
		    APPROX_ABS (a, Wxy [i]) ;
		    if (SCALAR_IS_NONZERO (a) && a >= toler)
		    {
			row = Pattern [i] ;
			deg = Row_degree [row] ;
			
			if (deg < rdeg [IN]
			    
			       || (deg == rdeg [IN] && a > pivot [IN])
			    
			    
			   )
			{
			    
			    pivrow [IN] = row ;
			    rdeg [IN] = deg ;
			    pivot [IN] = a ;
			}
		    }
		}
		for ( ; i < cdeg1 ; i++)
		{
		    double a ;
		    APPROX_ABS (a, Wxy [i]) ;
		    if (SCALAR_IS_NONZERO (a) && a >= toler)
		    {
			row = Pattern [i] ;
			deg = Row_degree [row] ;
			
			if (deg < rdeg [OUT]
			    
			       || (deg == rdeg [OUT] && a > pivot [OUT])
			    
			    
			   )
			{
			    
			    pivrow [OUT] = row ;
			    rdeg [OUT] = deg ;
			    pivot [OUT] = a ;
			}
		    }
		}

	    }
	    else
	    {

		
		for (i = 0 ; i < cdeg1 ; i++)
		{
		    double a ;
		    APPROX_ABS (a, Wxy [i]) ;
		    if (SCALAR_IS_NONZERO (a) && a >= toler)
		    {
			row = Pattern [i] ;
			deg = Row_degree [row] ;
			if (Frpos [row] >= 0 && Frpos [row] < fnrows)
			{
			    
			    if (deg < rdeg [IN]
			    
			       || (deg == rdeg [IN] && a > pivot [IN])
			    
			    
			       )
			    {
				
				pivrow [IN] = row ;
				rdeg [IN] = deg ;
				pivot [IN] = a ;
			    }
			}
			else
			{
			    
			    if (deg < rdeg [OUT]
			    
			       || (deg == rdeg[OUT] && a > pivot [OUT])
			    
			    
			       )
			    {
				
				pivrow [OUT] = row ;
				rdeg [OUT] = deg ;
				pivot [OUT] = a ;
			    }
			}
		    }
		}
	    }
	}
    }

    
    
    

    
    

    if (cdeg1 > 0 && pivrow [IN] == EMPTY && pivrow [OUT] == EMPTY)
    {
	

	
	
	row = Pattern [0] ;
	deg = Row_degree [row] ;
	if (Frpos [row] >= 0 && Frpos [row] < fnrows)
	{
	    
	    pivrow [IN] = row ;
	    rdeg [IN] = deg ;
	}
	else
	{
	    
	    pivrow [OUT] = row ;
	    rdeg [OUT] = deg ;
	}

	
	
	
	
	

	for (i = 1 ; i < cdeg1 ; i++)
	{
	    row = Pattern [i] ;
	    deg = Row_degree [row] ;

	    if (Frpos [row] >= 0 && Frpos [row] < fnrows)
	    {
		
		if (deg < rdeg [IN] || (deg == rdeg [IN] && row == diag_row))
		{
		    
		    pivrow [IN] = row ;
		    rdeg [IN] = deg ;
		}
	    }
	    else
	    {
		
		if (deg < rdeg [OUT] || (deg == rdeg [OUT] && row == diag_row))
		{
		    
		    pivrow [OUT] = row ;
		    rdeg [OUT] = deg ;
		}
	    }
	}
    }

    

    

    
    
    

    if (cdeg1  == 0)
    {
	if (fnrows > 0)
	{
	    
	    pivrow [IN] = Work->Frows [0] ;
	}
	else
	{

	    

	    Int *Front_leftmostdesc, *Front_1strow, *Front_new1strow, row1,
		row2, fleftmost, nfr, n_row, frontid ;

	    Front_leftmostdesc = Symbolic->Front_leftmostdesc ;
	    Front_1strow = Symbolic->Front_1strow ;
	    Front_new1strow = Work->Front_new1strow ;
	    nfr = Symbolic->nfr ;
	    n_row = Numeric->n_row ;
	    frontid = Work->frontid ;

	    
	    
	    
	    
	    
	    
	    

	    fleftmost = Front_leftmostdesc [frontid] ;
	    row1 = Front_new1strow [fleftmost] ;
	    row2 = Front_1strow [frontid+1] - 1 ;

	    
	    for (row = row1 ; row <= row2 ; row++)
	    {
		if (NON_PIVOTAL_ROW (row))
		{
		    
		    pivrow [OUT] = row ;
		    break ;
		}
	    }
	    Front_new1strow [fleftmost] = row ;

	    if (pivrow [OUT] == EMPTY)
	    {
		
		row1 = Front_new1strow [nfr] ;
		row2 = n_row-1 ;

		
		for (row = row1 ; row <= row2 ; row++)
		{
		    if (NON_PIVOTAL_ROW (row))
		    {
			
			pivrow [OUT] = row ;
			break ;
		    }
		}
		Front_new1strow [nfr] = row ;
	    }

	    if (pivrow [OUT] == EMPTY)
	    {
		
		
		return (SparseLU_WARNING_singular_matrix) ;
	    }
	}
    }

    
    
    

    if (pivrow [IN] != EMPTY)
    {

	
	freebie [IN] = (pivrow [IN] == prior_pivrow [IN]) && (cdeg1  > 0) ;
	ASSERT (cdeg1  >= 0) ;

	if (!freebie [IN])
	{
	    

	    Fcpos = Work->Fcpos ;
	    fncols = Work->fncols ;

	    Wrpflag = Work->Wrpflag ;

	    
	    
	    

	    rdeg [IN] = fncols ;

	    
	    if (Fcpos [pivcol] < 0)
	    {
		if (rdeg [IN] >= max_rdeg)
		{
		    
		    return (SparseLU_ERROR_different_pattern) ;
		}
		Wrp [pivcol] = Wrpflag ;
		W_i [rdeg [IN]++] = pivcol ;
	    }

	    tpi = Row_tuples [pivrow [IN]] ;
	    if (tpi)
	    {
		tp = (Tuple *) (Memory + tpi) ;
		tp1 = tp ;
		tp2 = tp ;
		tpend = tp + Row_tlen [pivrow [IN]] ;
		for ( ; tp < tpend ; tp++)
		{
		    e = tp->e ;
		    if (!E [e])
		    {
			continue ;	
		    }
		    f = tp->f ;
		    p = Memory + E [e] ;
		    ep = (Element *) p ;
		    p += UNITS (Element, 1) ;
		    Cols = (Int *) p ;
		    ncols = ep->ncols ;
		    Rows = Cols + ncols ;
		    if (Rows [f] == EMPTY)
		    {
			continue ;	
		    }

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if ((col >= 0) && (Wrp [col] != Wrpflag)
			    && Fcpos [col] <0)
			{
			    if (rdeg [IN] >= max_rdeg)
			    {
				
				return (SparseLU_ERROR_different_pattern) ;
			    }
			    Wrp [col] = Wrpflag ;
			    W_i [rdeg [IN]++] = col ;
			}
		    }

		    *tp2++ = *tp ;	
		}
		Row_tlen [pivrow [IN]] = tp2 - tp1 ;
	    }

	    

	    
	    Work->Wrpflag++ ;
	    
	}
    }

    
    
    

    
    

    if (pivrow [OUT] != EMPTY)
    {
	freebie [OUT] = (pivrow [OUT] == prior_pivrow [OUT]) && cdeg1  > 0 ;
	ASSERT (cdeg1  >= 0) ;

	if (!freebie [OUT])
	{

	    Wrpflag = Work->Wrpflag ;

	    
	    
	    

	    rdeg [OUT] = 0 ;

	    
	    if (rdeg [OUT] >= max_rdeg)
	    {
		
		return (SparseLU_ERROR_different_pattern) ;
	    }
	    Wrp [pivcol] = Wrpflag ;
	    W_o [rdeg [OUT]++] = pivcol ;

	    tpi = Row_tuples [pivrow [OUT]] ;
	    if (tpi)
	    {
		tp = (Tuple *) (Memory + tpi) ;
		tp1 = tp ;
		tp2 = tp ;
		tpend = tp + Row_tlen [pivrow [OUT]] ;
		for ( ; tp < tpend ; tp++)
		{
		    e = tp->e ;
		    if (!E [e])
		    {
			continue ;	
		    }
		    f = tp->f ;
		    p = Memory + E [e] ;
		    ep = (Element *) p ;
		    p += UNITS (Element, 1) ;
		    Cols = (Int *) p ;
		    ncols = ep->ncols ;
		    Rows = Cols + ncols ;
		    if (Rows [f] == EMPTY)
		    {
			continue ;	
		    }

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if ((col >= 0) && (Wrp [col] != Wrpflag))
			{
			    if (rdeg [OUT] >= max_rdeg)
			    {
				
				return (SparseLU_ERROR_different_pattern) ;
			    }
			    Wrp [col] = Wrpflag ;
			    W_o [rdeg [OUT]++] = col ;
			}
		    }
		    *tp2++ = *tp ;	
		}
		Row_tlen [pivrow [OUT]] = tp2 - tp1 ;
	    }

	    

	    
	    Work->Wrpflag++ ;
	    

	}

    }

    return (SparseLU_OK) ;
}





