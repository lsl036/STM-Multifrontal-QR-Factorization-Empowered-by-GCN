/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_store_lu.c
 * BRIEF:   LU存储
 *****************************************************************************/



#include "SparseLU_internal.h"
#include "SparseLU_function.h"



#ifdef DROP
GLOBAL Int LU_store_lu_drop
#else
GLOBAL Int LU_store_lu
#endif
(
    NumericType *Numeric,
    WorkType *Work
)
{
    
    
    

    Entry pivot_value ;
#ifdef DROP
    double droptol ;
#endif
    Entry *D, *Lval, *Uval, *Fl1, *Fl2, *Fu1, *Fu2,
	*Flublock, *Flblock, *Fublock ;
    Int i, k, fnr_curr, fnrows, fncols, row, col, pivrow, pivcol, *Frows,
	*Fcols, *Lpattern, *Upattern, *Lpos, *Upos, llen, ulen, fnc_curr, fnpiv,
	uilen, lnz, unz, nb, *Lilen,
	*Uilen, *Lip, *Uip, *Li, *Ui, pivcol_position, newLchain, newUchain,
	pivrow_position, p, size, lip, uip, lnzi, lnzx, unzx, lnz2i, lnz2x,
	unz2i, unz2x, zero_pivot, *Pivrow, *Pivcol, kk,
	Lnz [MAXNB] ;

#ifdef DROP
    Int all_lnz, all_unz ;
    droptol = Numeric->droptol ;
#endif

    
    
    

    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    fnpiv = Work->fnpiv ;

    Lpos = Numeric->Lpos ;
    Upos = Numeric->Upos ;
    Lilen = Numeric->Lilen ;
    Uilen = Numeric->Uilen ;

    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;

    Flublock = Work->Flublock ;
    Flblock  = Work->Flblock ;
    Fublock  = Work->Fublock ;

    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;
    Frows = Work->Frows ;
    Fcols = Work->Fcols ;

    Lpattern = Work->Lpattern ;
    llen = Work->llen ;
    Upattern = Work->Upattern ;
    ulen = Work->ulen ;

    nb = Work->nb ;

    Pivrow = Work->Pivrow ;
    Pivcol = Work->Pivcol ;

    
    
    

    for (kk = 0 ; kk < fnpiv ; kk++)
    {

	
	
	

	k = Work->npiv + kk ;

	
	
	

	pivrow = Pivrow [kk] ;
	pivcol = Pivcol [kk] ;

	
	
	

	
	
	pivrow_position = Lpos [pivrow] ;
	if (pivrow_position != EMPTY)
	{
	    
	    
	    row = Lpattern [--llen] ;
	    
	    Lpattern [pivrow_position] = row ;
	    Lpos [row] = pivrow_position ;
	    Lpos [pivrow] = EMPTY ;
	}

	
	
	

	
	Fl1 = Flublock + kk * nb ;

	
	Fl2 = Flblock + kk * fnr_curr ;

	
	pivot_value = Fl1 [kk] ;

	D [k] = pivot_value ;
	zero_pivot = IS_ZERO (pivot_value) ;

	
	
	

	lnz = 0 ;
	lnz2i = 0 ;
	lnz2x = llen ;

#ifdef DROP
	    all_lnz = 0 ;

	    for (i = kk + 1 ; i < fnpiv ; i++)
	    {
		Entry x ;
		double s ;
		x = Fl1 [i] ;
		if (IS_ZERO (x)) continue ;
		all_lnz++ ;
		APPROX_ABS (s, x) ;
		if (s <= droptol) continue ;
		lnz++ ;
		if (Lpos [Pivrow [i]] == EMPTY) lnz2i++ ;
	    }

	    for (i = 0 ; i < fnrows ; i++)
	    {
		Entry x ;
		double s ;
		x = Fl2 [i] ;
		if (IS_ZERO (x)) continue ;
		all_lnz++ ;
		APPROX_ABS (s, x) ;
		if (s <= droptol) continue ;
		lnz++ ;
		if (Lpos [Frows [i]] == EMPTY) lnz2i++ ;
	    }

#else

	    for (i = kk + 1 ; i < fnpiv ; i++)
	    {
		if (IS_ZERO (Fl1 [i])) continue ;
		lnz++ ;
		if (Lpos [Pivrow [i]] == EMPTY) lnz2i++ ;
	    }

	    for (i = 0 ; i < fnrows ; i++)
	    {
		if (IS_ZERO (Fl2 [i])) continue ;
		lnz++ ;
		if (Lpos [Frows [i]] == EMPTY) lnz2i++ ;
	    }

#endif

	lnz2x += lnz2i ;

	
	if (llen == 0 || zero_pivot)
	{
	    
	    
	    newLchain = TRUE ;
	}
	else
	{
	    newLchain =
		    
		    UNITS (Entry, lnz) + UNITS (Int, lnz)
		<=
		    
		    UNITS (Entry, lnz2x) + UNITS (Int, lnz2i) ;
	}

	if (newLchain)
	{
	    
	    pivrow_position = EMPTY ;

	    
	    for (i = 0 ; i < llen ; i++)
	    {
		row = Lpattern [i] ;
		Lpos [row] = EMPTY ;
	    }
	    llen = 0 ;

	    lnzi = lnz ;
	    lnzx = lnz ;
	}
	else
	{
	    
	    lnzi = lnz2i ;
	    lnzx = lnz2x ;
	}

	
	
	

	size = UNITS (Int, lnzi) + UNITS (Entry, lnzx) ;

	p = LU_mem_alloc_head_block (Numeric, size) ;
	if (!p)
	{
	    Int r2, c2 ;
	    
	    
	    if (Work->do_grow)
	    {
		
		r2 = fnrows ;
		c2 = fncols ;
	    }
	    else
	    {
		
		r2 = MAX (fnrows, Work->fnrows_new + 1) ;
		c2 = MAX (fncols, Work->fncols_new + 1) ;
	    }
	    if (!LU_get_memory (Numeric, Work, size, r2, c2, TRUE))
	    {
		return (FALSE) ;	
	    }
	    p = LU_mem_alloc_head_block (Numeric, size) ;
	    if (!p)
	    {
		return (FALSE) ;	
	    }
	    
	    fnc_curr = Work->fnc_curr ;
	    fnr_curr = Work->fnr_curr ;
	    Flublock = Work->Flublock ;
	    Flblock  = Work->Flblock ;
	    Fublock  = Work->Fublock ;
	    Fl1 = Flublock + kk * nb ;
	    Fl2 = Flblock  + kk * fnr_curr ;
	}

	
	
	

	lip = p ;

	Li = (Int *) (Numeric->Memory + p) ;
	p += UNITS (Int, lnzi) ;
	Lval = (Entry *) (Numeric->Memory + p) ;
	p += UNITS (Entry, lnzx) ;

	for (i = 0 ; i < lnzx ; i++)
	{
	    CLEAR (Lval [i]) ;
	}

	

	if (newLchain)
	{
	    
	    lip = -lip ;

	    ASSERT (llen == 0) ;

#ifdef DROP

		for (i = kk + 1 ; i < fnpiv ; i++)
		{
		    Entry x ;
		    double s ;
		    Int row2, pos ;
		    x = Fl1 [i] ;
		    APPROX_ABS (s, x) ;
		    if (s <= droptol) continue ;
		    row2 = Pivrow [i] ;
		    pos = llen++ ;
		    Lpattern [pos] = row2 ;
		    Lpos [row2] = pos ;
		    Li [pos] = row2 ;
		    Lval [pos] = x ;
		}

		for (i = 0 ; i < fnrows ; i++)
		{
		    Entry x ;
		    double s ;
		    Int row2, pos ;
		    x = Fl2 [i] ;
		    APPROX_ABS (s, x) ;
		    if (s <= droptol) continue ;
		    row2 = Frows [i] ;
		    pos = llen++ ;
		    Lpattern [pos] = row2 ;
		    Lpos [row2] = pos ;
		    Li [pos] = row2 ;
		    Lval [pos] = x ;
		}

#else

		for (i = kk + 1 ; i < fnpiv ; i++)
		{
		    Entry x ;
		    Int row2, pos ;
		    x = Fl1 [i] ;
		    if (IS_ZERO (x)) continue ;
		    row2 = Pivrow [i] ;
		    pos = llen++ ;
		    Lpattern [pos] = row2 ;
		    Lpos [row2] = pos ;
		    Li [pos] = row2 ;
		    Lval [pos] = x ;
		}

		for (i = 0 ; i < fnrows ; i++)
		{
		    Entry x ;
		    Int row2, pos ;
		    x = Fl2 [i] ;
		    if (IS_ZERO (x)) continue ;
		    row2 = Frows [i] ;
		    pos = llen++ ;
		    Lpattern [pos] = row2 ;
		    Lpos [row2] = pos ;
		    Li [pos] = row2 ;
		    Lval [pos] = x ;
		}

#endif

	}
	else
	{
	    ASSERT (llen > 0) ;

#ifdef DROP

		for (i = kk + 1 ; i < fnpiv ; i++)
		{
		    Entry x ;
		    double s ;
		    Int row2, pos ;
		    x = Fl1 [i] ;
		    APPROX_ABS (s, x) ;
		    if (s <= droptol) continue ;
		    row2 = Pivrow [i] ;
		    pos = Lpos [row2] ;
		    if (pos == EMPTY)
		    {
			pos = llen++ ;
			Lpattern [pos] = row2 ;
			Lpos [row2] = pos ;
			*Li++ = row2 ;
		    }
		    Lval [pos] = x ;
		}

		for (i = 0 ; i < fnrows ; i++)
		{
		    Entry x ;
		    double s ;
		    Int row2, pos ;
		    x = Fl2 [i] ;
		    APPROX_ABS (s, x) ;
		    if (s <= droptol) continue ;
		    row2 = Frows [i] ;
		    pos = Lpos [row2] ;
		    if (pos == EMPTY)
		    {
			pos = llen++ ;
			Lpattern [pos] = row2 ;
			Lpos [row2] = pos ;
			*Li++ = row2 ;
		    }
		    Lval [pos] = x ;
		}

#else

		for (i = kk + 1 ; i < fnpiv ; i++)
		{
		    Entry x ;
		    Int row2, pos ;
		    x = Fl1 [i] ;
		    if (IS_ZERO (x)) continue ;
		    row2 = Pivrow [i] ;
		    pos = Lpos [row2] ;
		    if (pos == EMPTY)
		    {
			pos = llen++ ;
			Lpattern [pos] = row2 ;
			Lpos [row2] = pos ;
			*Li++ = row2 ;
		    }
		    Lval [pos] = x ;
		}

		for (i = 0 ; i < fnrows ; i++)
		{
		    Entry x ;
		    Int row2, pos ;
		    x = Fl2 [i] ;
		    if (IS_ZERO (x)) continue ;
		    row2 = Frows [i] ;
		    pos = Lpos [row2] ;
		    if (pos == EMPTY)
		    {
			pos = llen++ ;
			Lpattern [pos] = row2 ;
			Lpos [row2] = pos ;
			*Li++ = row2 ;
		    }
		    Lval [pos] = x ;
		}

#endif

	}

#ifdef DROP
	    Numeric->lnz += lnz ;
	    Numeric->all_lnz += all_lnz ;
	    Lnz [kk] = all_lnz ;

#else

	    Numeric->lnz += lnz ;
	    Numeric->all_lnz += lnz ;
	    Lnz [kk] = lnz ;
#endif

	Numeric->nLentries += lnzx ;
	Work->llen = llen ;
	Numeric->isize += lnzi ;

	
	
	
	

	Lpos [pivrow] = pivrow_position ;	
	Lip [pivcol] = lip ;			
	Lilen [pivcol] = lnzi ;			

    }

    
    
    

    for (kk = 0 ; kk < fnpiv ; kk++)
    {

	
	
	

	k = Work->npiv + kk ;

	
	
	

	pivrow = Pivrow [kk] ;
	pivcol = Pivcol [kk] ;
	
	
	

	zero_pivot = IS_ZERO (D [k]) ;

	
	
	

	
	Fu1 = Flublock + kk ;

	
	Fu2 = Fublock + kk * fnc_curr ;

	unz = 0 ;
	unz2i = 0 ;
	unz2x = ulen ;

	
	pivcol_position = Upos [pivcol] ;
	if (pivcol_position != EMPTY)
	{
	    unz2x-- ;
	}

#ifdef DROP
	    all_unz = 0 ;

	    for (i = kk + 1 ; i < fnpiv ; i++)
	    {
		Entry x ;
		double s ;
		x = Fu1 [i*nb] ;
		if (IS_ZERO (x)) continue ;
		all_unz++ ;
		APPROX_ABS (s, x) ;
		if (s <= droptol) continue ;
		unz++ ;
		if (Upos [Pivcol [i]] == EMPTY) unz2i++ ;
	    }

	    for (i = 0 ; i < fncols ; i++)
	    {
		Entry x ;
		double s ;
		x = Fu2 [i] ;
		if (IS_ZERO (x)) continue ;
		all_unz++ ;
		APPROX_ABS (s, x) ;
		if (s <= droptol) continue ;
		unz++ ;
		if (Upos [Fcols [i]] == EMPTY) unz2i++ ;
	    }

#else

	    for (i = kk + 1 ; i < fnpiv ; i++)
	    {
		if (IS_ZERO (Fu1 [i*nb])) continue ;
		unz++ ;
		if (Upos [Pivcol [i]] == EMPTY) unz2i++ ;
	    }

	    for (i = 0 ; i < fncols ; i++)
	    {
		if (IS_ZERO (Fu2 [i])) continue ;
		unz++ ;
		if (Upos [Fcols [i]] == EMPTY) unz2i++ ;
	    }

#endif

	unz2x += unz2i ;

	ASSERT (IMPLIES (k == 0, ulen == 0)) ;

	
	if (ulen == 0 || zero_pivot)
	{
	    
	    
	    
	    
	    newUchain = TRUE ;
	}
	else
	{
	    newUchain =
		    
		    UNITS (Entry, unz) + UNITS (Int, unz)
		<=
		    
		    UNITS (Entry, unz2x) + UNITS (Int, unz2i) ;

	    
	    
	}

	
	
	

	size = 0 ;
	if (newUchain)
	{
	    
	    size += UNITS (Int, ulen) ;
	    unzx = unz ;
	}
	else
	{
	    unzx = unz2x ;
	}
	size += UNITS (Entry, unzx) ;

	p = LU_mem_alloc_head_block (Numeric, size) ;
	if (!p)
	{
	    Int r2, c2 ;
	    
	    
	    if (Work->do_grow)
	    {
		
		r2 = fnrows ;
		c2 = fncols ;
	    }
	    else
	    {
		
		r2 = MAX (fnrows, Work->fnrows_new + 1) ;
		c2 = MAX (fncols, Work->fncols_new + 1) ;
	    }
	    if (!LU_get_memory (Numeric, Work, size, r2, c2, TRUE))
	    {
		
		return (FALSE) ;	
	    }
	    p = LU_mem_alloc_head_block (Numeric, size) ;
	    if (!p)
	    {
		
		return (FALSE) ;	
	    }
	    
	    fnc_curr = Work->fnc_curr ;
	    fnr_curr = Work->fnr_curr ;
	    Flublock = Work->Flublock ;
	    Flblock  = Work->Flblock ;
	    Fublock  = Work->Fublock ;
	    Fu1 = Flublock + kk ;
	    Fu2 = Fublock  + kk * fnc_curr ;
	}

	
	
	

	uip = p ;

	if (newUchain)
	{
	    
	    uip = -uip ;

	    pivcol_position = EMPTY ;

	    
	    
	    
	    uilen = ulen ;
	    Ui = (Int *) (Numeric->Memory + p) ;
	    Numeric->isize += ulen ;
	    p += UNITS (Int, ulen) ;
	    for (i = 0 ; i < ulen ; i++)
	    {
		col = Upattern [i] ;
		Upos [col] = EMPTY ;
		Ui [i] = col ;
	    }

	    ulen = 0 ;

	}
	else
	{
	    

	    
	    

	    if (pivcol_position != EMPTY)
	    {
		
		
		col = Upattern [--ulen] ;
		Upattern [pivcol_position] = col ;
		Upos [col] = pivcol_position ;
		Upos [pivcol] = EMPTY ;
	    }

	    
	    
	    
	    uilen = unz2i ;

	}

	Uval = (Entry *) (Numeric->Memory + p) ;
	

	for (i = 0 ; i < unzx ; i++)
	{
	    CLEAR (Uval [i]) ;
	}

	if (newUchain)
	{
	    ASSERT (ulen == 0) ;

#ifdef DROP

		for (i = kk + 1 ; i < fnpiv ; i++)
		{
		    Entry x ;
		    double s ;
		    Int col2, pos ;
		    x = Fu1 [i*nb] ;
		    APPROX_ABS (s, x) ;
		    if (s <= droptol) continue ;
		    col2 = Pivcol [i] ;
		    pos = ulen++ ;
		    Upattern [pos] = col2 ;
		    Upos [col2] = pos ;
		    Uval [pos] = x ;
		}

		for (i = 0 ; i < fncols ; i++)
		{
		    Entry x ;
		    double s ;
		    Int col2, pos ;
		    x = Fu2 [i] ;
		    APPROX_ABS (s, x) ;
		    if (s <= droptol) continue ;
		    col2 = Fcols [i] ;
		    pos = ulen++ ;
		    Upattern [pos] = col2 ;
		    Upos [col2] = pos ;
		    Uval [pos] = x ;
		}

#else

		for (i = kk + 1 ; i < fnpiv ; i++)
		{
		    Entry x ;
		    Int col2, pos ;
		    x = Fu1 [i*nb] ;
		    if (IS_ZERO (x)) continue ;
		    col2 = Pivcol [i] ;
		    pos = ulen++ ;
		    Upattern [pos] = col2 ;
		    Upos [col2] = pos ;
		    Uval [pos] = x ;
		}

		for (i = 0 ; i < fncols ; i++)
		{
		    Entry x ;
		    Int col2, pos ;
		    x = Fu2 [i] ;
		    if (IS_ZERO (x)) continue ;
		    col2 = Fcols [i] ;
		    pos = ulen++ ;
		    Upattern [pos] = col2 ;
		    Upos [col2] = pos ;
		    Uval [pos] = x ;
		}

#endif

	}
	else
	{

	    ASSERT (ulen > 0) ;

	    

#ifdef DROP

		for (i = kk + 1 ; i < fnpiv ; i++)
		{
		    Entry x ;
		    double s ;
		    Int col2, pos ;
		    x = Fu1 [i*nb] ;
		    APPROX_ABS (s, x) ;
		    if (s <= droptol) continue ;
		    col2 = Pivcol [i] ;
		    pos = Upos [col2] ;
		    if (pos == EMPTY)
		    {
			pos = ulen++ ;
			Upattern [pos] = col2 ;
			Upos [col2] = pos ;
		    }
		    Uval [pos] = x ;
		}

		for (i = 0 ; i < fncols ; i++)
		{
		    Entry x ;
		    double s ;
		    Int col2, pos ;
		    x = Fu2 [i] ;
		    APPROX_ABS (s, x) ;
		    if (s <= droptol) continue ;
		    col2 = Fcols [i] ;
		    pos = Upos [col2] ;
		    if (pos == EMPTY)
		    {
			pos = ulen++ ;
			Upattern [pos] = col2 ;
			Upos [col2] = pos ;
		    }
		    Uval [pos] = x ;
		}

#else

		for (i = kk + 1 ; i < fnpiv ; i++)
		{
		    Entry x ;
		    Int col2, pos ;
		    x = Fu1 [i*nb] ;
		    if (IS_ZERO (x)) continue ;
		    col2 = Pivcol [i] ;
		    pos = Upos [col2] ;
		    if (pos == EMPTY)
		    {
			pos = ulen++ ;
			Upattern [pos] = col2 ;
			Upos [col2] = pos ;
		    }
		    Uval [pos] = x ;
		}

		for (i = 0 ; i < fncols ; i++)
		{
		    Entry x ;
		    Int col2, pos ;
		    x = Fu2 [i] ;
		    if (IS_ZERO (x)) continue ;
		    col2 = Fcols [i] ;
		    pos = Upos [col2] ;
		    if (pos == EMPTY)
		    {
			pos = ulen++ ;
			Upattern [pos] = col2 ;
			Upos [col2] = pos ;
		    }
		    Uval [pos] = x ;
		}

#endif

	}

#ifdef DROP

	    Numeric->unz += unz ;
	    Numeric->all_unz += all_unz ;
	    
	    Numeric->flops += DIV_FLOPS * Lnz [kk]	
		+ MULTSUB_FLOPS * (Lnz [kk] * all_unz) ;    

#else

	    Numeric->unz += unz ;
	    Numeric->all_unz += unz ;
	    
	    Numeric->flops += DIV_FLOPS * Lnz [kk]	
		+ MULTSUB_FLOPS * (Lnz [kk] * unz) ;    
#endif

	Numeric->nUentries += unzx ;
	Work->ulen = ulen ;

	
	
	

	Upos [pivcol] = pivcol_position ;	
	Uip [pivrow] = uip ;			
	Uilen [pivrow] = uilen ;		

    }

    
    
    

    Work->npiv += fnpiv ;
    Work->fnpiv = 0 ;
    Work->fnzeros = 0 ;
    return (TRUE) ;
}
