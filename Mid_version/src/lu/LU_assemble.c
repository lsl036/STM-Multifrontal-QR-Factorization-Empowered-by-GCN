/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    LU_assemble.c
 * BRIEF:   LU组装
 *****************************************************************************/



#include "SparseLU_internal.h"
#include "SparseLU_function.h"







PRIVATE void row_assemble
(
    Int row,
    NumericType *Numeric,
    WorkType *Work
)
{

    Entry *S, *Fcblock, *Frow ;
    Int tpi, e, *E, *Fcpos, *Frpos, *Row_degree, *Row_tuples, *Row_tlen, rdeg0,
	f, nrows, ncols, *Rows, *Cols, col, ncolsleft, j ;
    Tuple *tp, *tp1, *tp2, *tpend ;
    Unit *Memory, *p ;
    Element *ep ;

#ifndef FIXQ
    Int *Col_degree ;
    Col_degree = Numeric->Cperm ;
#endif

    Row_tuples = Numeric->Uip ;
    tpi = Row_tuples [row] ;
    if (!tpi) return ;

    Memory = Numeric->Memory ;
    E = Work->E ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Row_degree = Numeric->Rperm ;
    Row_tlen   = Numeric->Uilen ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    rdeg0 = Work->rdeg0 ;
    Fcblock = Work->Fcblock ;

    ASSERT (NON_PIVOTAL_ROW (row)) ;

    tp = (Tuple *) (Memory + tpi) ;
    tp1 = tp ;
    tp2 = tp ;
    tpend = tp + Row_tlen [row] ;
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
	Rows = Cols + ep->ncols ;
	if (Rows [f] == EMPTY) continue ;   
	ASSERT (row == Rows [f] && row >= 0 && row < Work->n_row) ;

	if (ep->rdeg == rdeg0)
	{
	    
	    
	    

	    
	    Rows [f] = EMPTY ;

	    nrows = ep->nrows ;
	    ncols = ep->ncols ;

	    p += UNITS (Int, ncols + nrows) ;
	    S = ((Entry *) p) + f ;

	    ncolsleft = ep->ncolsleft ;

	    Frow = Fcblock + Frpos [row] ;

	    Row_degree [row] -= ncolsleft ;

	    if (ncols == ncolsleft)
	    {
		
		
		

#pragma ivdep
		for (j = 0 ; j < ncols ; j++)
		{
		    col = Cols [j] ;
		    ASSERT (col >= 0 && col < Work->n_col) ;
#ifndef FIXQ
		    Col_degree [col] -- ;
#endif
		    
		    ASSEMBLE (Frow [Fcpos [col]], *S) ;
		    S += nrows ;
		}

	    }
	    else
	    {
		
		
		

#pragma ivdep
		for (j = 0 ; j < ncols ; j++)
		{
		    col = Cols [j] ;
		    if (col >= 0)
		    {
			ASSERT (col < Work->n_col) ;
#ifndef FIXQ
			Col_degree [col] -- ;
#endif
			
			ASSEMBLE (Frow [Fcpos [col]], *S) ;
		    }
		    S += nrows ;
		}

	    }
	    ep->nrowsleft-- ;
	    ASSERT (ep->nrowsleft > 0) ;
	}
	else
	{
	    *tp2++ = *tp ;	
	}
    }
    Row_tlen [row] = tp2 - tp1 ;
}





PRIVATE void col_assemble
(
    Int col,
    NumericType *Numeric,
    WorkType *Work
)
{

    Entry *S, *Fcblock, *Fcol ;
    Int tpi, e, *E, *Fcpos, *Frpos, *Row_degree, *Col_tuples, *Col_tlen, cdeg0,
	f, nrows, ncols, *Rows, *Cols, row, nrowsleft, i ;
    Tuple *tp, *tp1, *tp2, *tpend ;
    Unit *Memory, *p ;
    Element *ep ;

#if !defined (FIXQ)
    Int *Col_degree ;
    Col_degree = Numeric->Cperm ;
#endif

    Col_tuples = Numeric->Lip ;
    tpi = Col_tuples [col] ;
    if (!tpi) return ;

    Memory = Numeric->Memory ;
    E = Work->E ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Row_degree = Numeric->Rperm ;
    Col_tlen   = Numeric->Lilen ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    cdeg0 = Work->cdeg0 ;
    Fcblock = Work->Fcblock ;

    tp = (Tuple *) (Memory + tpi) ;
    tp1 = tp ;
    tp2 = tp ;
    tpend = tp + Col_tlen [col] ;
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
	ASSERT (col == Cols [f] && col >= 0 && col < Work->n_col) ;

	if (ep->cdeg == cdeg0)
	{
	    
	    
	    

	    
	    Cols [f] = EMPTY ;

	    nrows = ep->nrows ;
	    ncols = ep->ncols ;
	    Rows = Cols + ncols ;
	    p += UNITS (Int, ncols + nrows) ;
	    S = ((Entry *) p) + f * nrows ;

	    nrowsleft = ep->nrowsleft ;

	    Fcol = Fcblock + Fcpos [col] ;
#ifndef FIXQ
	    Col_degree [col] -= nrowsleft ;
#endif
	    if (nrows == nrowsleft)
	    {
		
		
		

#pragma ivdep
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    ASSERT (row >= 0 && row < Work->n_row) ;
		    Row_degree [row]-- ;
		    
		    ASSEMBLE (Fcol [Frpos [row]], S [i]) ;
		}
	    }
	    else
	    {
		
		
		

#pragma ivdep
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0)
		    {
			ASSERT (row < Work->n_row) ;
			Row_degree [row]-- ;
			
			ASSEMBLE (Fcol [Frpos [row]], S [i]) ;
		    }
		}
	    }
	    ep->ncolsleft-- ;
	    ASSERT (ep->ncolsleft > 0) ;
	}
	else
	{
	    *tp2++ = *tp ;	
	}
    }
    Col_tlen [col] = tp2 - tp1 ;

}






#ifndef FIXQ
GLOBAL void LU_assemble
#else
GLOBAL void LU_assemble_fixq
#endif
(
    NumericType *Numeric,
    WorkType *Work
)
{
    
    
    

    Int e, i, row, col, i2, nrows, ncols, f, tpi, extcdeg, extrdeg, rdeg0,
	cdeg0, son_list, next, nrows_to_assemble,
	ncols_to_assemble, ngetrows, j, j2,
	nrowsleft,	
	ncolsleft,	
	prior_Lson, prior_Uson, *E, *Cols, *Rows, *Wm, *Woo,
	*Row_tuples, *Row_degree, *Row_tlen,
	*Col_tuples, *Col_tlen ;
    Unit *Memory, *p ;
    Element *ep ;
    Tuple *tp, *tp1, *tp2, *tpend ;
    Entry
	*S,		
	*Fcblock,	
	*Fcol ;		
    Int *Frpos,
	*Fcpos,
	fnrows,		
	fncols ;	

#if !defined (FIXQ)
    Int *Col_degree ;
#endif

#if !defined (FIXQ)
    Col_degree = Numeric->Cperm ;   
#endif

    
    
    

    fncols = Work->fncols ;
    fnrows = Work->fnrows ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Row_degree = Numeric->Rperm ;
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    Wm = Work->Wm ;
    Woo = Work->Woo ;
    rdeg0 = Work->rdeg0 ;
    cdeg0 = Work->cdeg0 ;

    Fcblock = Work->Fcblock ;

    
    
    

    ASSERT (fnrows == Work->fnrows_new + 1) ;
    ASSERT (fncols == Work->fncols_new + 1) ;

    Numeric->maxnrows = MAX (Numeric->maxnrows, fnrows) ;
    Numeric->maxncols = MAX (Numeric->maxncols, fncols) ;

    
    Numeric->maxfrsize = MAX (Numeric->maxfrsize, fnrows * fncols) ;

    
    
    

    
    

    son_list = 0 ;	

    
    
    

    if (!Work->do_extend)
    {
	prior_Uson = ( Work->pivcol_in_front && !Work->pivrow_in_front) ;
	prior_Lson = (!Work->pivcol_in_front &&  Work->pivrow_in_front) ;
	if (prior_Uson || prior_Lson)
	{
	    e = Work->prior_element ;
	    if (e != EMPTY)
	    {
		ASSERT (E [e]) ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		ep->next = son_list ;
		son_list = e ;
		ASSERT (E [e]) ;
	    }
	}
    }
    Work->prior_element = EMPTY ;

    
    
    
    

    for (i2 = Work->fscan_row ; i2 < fnrows ; i2++)
    {
	
	row = Work->NewRows [i2] ;
	if (row < 0) row = FLIP (row) ;

	tpi = Row_tuples [row] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tp1 = tp ;
	tp2 = tp ;
	tpend = tp + Row_tlen [row] ;
	for ( ; tp < tpend ; tp++)
	{
	    e = tp->e ;
	    if (!E [e]) continue ;	
	    f = tp->f ;
	    p = Memory + E [e] ;
	    ep = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Rows = ((Int *) p) + ep->ncols ;
	    if (Rows [f] == EMPTY) continue ;	

	    if (ep->cdeg < cdeg0)
	    {
		
		ep->cdeg = ep->nrowsleft + cdeg0 ;
	    }

	    ep->cdeg-- ;	

	    
	    if (ep->cdeg == cdeg0 && ep->next == EMPTY)
	    {
		
		ep->next = son_list ;
		son_list = e ;
	    }

	    *tp2++ = *tp ;	
	}
	Row_tlen [row] = tp2 - tp1 ;
    }

    
    
    
    

    for (j2 = Work->fscan_col ; j2 < fncols ; j2++)
    {
	
	col = Work->NewCols [j2] ;
	if (col < 0) col = FLIP (col) ;

	tpi = Col_tuples [col] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tp1 = tp ;
	tp2 = tp ;
	tpend = tp + Col_tlen [col] ;
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
	    ASSERT (col == Cols [f]) ;

	    if (ep->rdeg < rdeg0)
	    {
		
		ep->rdeg = ep->ncolsleft + rdeg0 ;
	    }

	    ep->rdeg-- ;	

	    if (ep->rdeg == rdeg0 && ep->next == EMPTY)
	    {
		
		ep->next = son_list ;
		son_list = e ;
	    }

	    *tp2++ = *tp ;	
	}
	Col_tlen [col] = tp2 - tp1 ;
    }

    
    
    

    next = EMPTY ;

    for (e = son_list ; e > 0 ; e = next)
    {
	p = Memory + E [e] ;
	GET_ELEMENT (ep, p, Cols, Rows, ncols, nrows, S) ;
	nrowsleft = ep->nrowsleft ;
	ncolsleft = ep->ncolsleft ;
	next = ep->next ;
	ep->next = EMPTY ;

	extrdeg = (ep->rdeg < rdeg0) ? ncolsleft : (ep->rdeg - rdeg0) ;
	extcdeg = (ep->cdeg < cdeg0) ? nrowsleft : (ep->cdeg - cdeg0) ;
	ncols_to_assemble = ncolsleft - extrdeg ;
	nrows_to_assemble = nrowsleft - extcdeg ;

	if (extrdeg == 0 && extcdeg == 0)
	{

	    
	    
	    

	    if (nrows == nrowsleft)
	    {
		
		
		

		
		
#pragma ivdep
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    Row_degree [row] -= ncolsleft ;
		    Wm [i] = Frpos [row] ;
		}

		if (ncols == ncolsleft)
		{
		    
		    
		    

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
#ifndef FIXQ
			Col_degree [col] -= nrowsleft ;
#endif
			Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			for (i = 0 ; i < nrows ; i++)
			{
			    
			    ASSEMBLE (Fcol [Wm [i]], S [i]) ;
			}
			S += nrows ;
		    }


		}
		else
		{
		    
		    
		    

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if (col >= 0)
			{
#ifndef FIXQ
			    Col_degree [col] -= nrowsleft ;
#endif
			    Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			    for (i = 0 ; i < nrows ; i++)
			    {
				
				ASSEMBLE (Fcol [Wm [i]], S [i]) ;
			    }
			}
			S += nrows ;
		    }

		}
		
	    }
	    else
	    {
		
		
		

		
		
		ngetrows = 0 ;
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0)
		    {
			Row_degree [row] -= ncolsleft ;
			Woo [ngetrows] = i ;
			Wm [ngetrows++] = Frpos [row] ;
		    }
		}
		ASSERT (ngetrows == nrowsleft) ;

		if (ncols == ncolsleft)
		{
		    
		    
		    

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
#ifndef FIXQ
			Col_degree [col] -= nrowsleft ;
#endif
			Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			for (i = 0 ; i < nrowsleft ; i++)
			{
			    
			    ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
			}
			S += nrows ;
		    }

		}
		else
		{
		    
		    
		    

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if (col >= 0)
			{
#ifndef FIXQ
			    Col_degree [col] -= nrowsleft ;
#endif
			    Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			    for (i = 0 ; i < nrowsleft ; i++)
			    {
				
				ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
			    }
			}
			S += nrows ;
		    }

		}
		
	    }

	    
	    LU_mem_free_tail_block (Numeric, E [e]) ;
	    E [e] = 0 ;

	}
	else if (extcdeg == 0)
	{

	    
	    
	    

	    if (ncols_to_assemble > 0)
	    {

		Int skip = FALSE ;
		if (ncols_to_assemble * 16 < ncols && nrows == 1)
		{
		    
		    ASSERT (nrowsleft == 1) ;
		    ASSERT (Rows [0] >= 0 && Rows [0] < Work->n_row) ;
		    skip = TRUE ;
		    Work->any_skip = TRUE ;
		}

		if (!skip)
		{

		    if (nrows == nrowsleft)
		    {
			
			
			

			
			
#pragma ivdep
			for (i = 0 ; i < nrows ; i++)
			{
			    row = Rows [i] ;
			    ASSERT (row >= 0 && row < n_row) ;
			    Row_degree [row] -= ncols_to_assemble ;
			    Wm [i] = Frpos [row] ;
			}

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    if ((col >= 0) && (Fcpos [col] >= 0))
			    {
#ifndef FIXQ
				Col_degree [col] -= nrowsleft ;
#endif
				Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
				for (i = 0 ; i < nrows ; i++)
				{
				    
				    ASSEMBLE (Fcol [Wm [i]], S [i]) ;
				}
				
				Cols [j] = EMPTY ;
			    }
			    S += nrows ;
			}


			
		    }
		    else
		    {
			
			
			

			
			
			ngetrows = 0 ;
			for (i = 0 ; i < nrows ; i++)
			{
			    row = Rows [i] ;
			    if (row >= 0)
			    {
				Row_degree [row] -= ncols_to_assemble ;
				ASSERT (row < n_row && Frpos [row] >= 0) ;
				Woo [ngetrows] = i ;
				Wm [ngetrows++] = Frpos [row] ;
			    }
			}
			ASSERT (ngetrows == nrowsleft) ;

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    if ((col >= 0) && (Fcpos [col] >= 0))
			    {
#ifndef FIXQ
				Col_degree [col] -= nrowsleft ;
#endif
				Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
				for (i = 0 ; i < nrowsleft ; i++)
				{
				    
				    ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
				}
				
				Cols [j] = EMPTY ;
			    }
			    S += nrows ;
			}

			
		    }
		    ep->ncolsleft = extrdeg ;
		}
	    }

	}
	else
	{

	    
	    
	    

	    if (nrows_to_assemble > 0)
	    {

		Int skip = FALSE ;
		if (nrows_to_assemble * 16 < nrows && ncols == 1)
		{
		    
		    ASSERT (ncolsleft == 1) ;
		    ASSERT (Cols [0] >= 0 && Cols [0] < Work->n_col) ;
		    Work->any_skip = TRUE ;
		    skip = TRUE ;
		}

		if (!skip)
		{

		    
		    
		    ngetrows = 0 ;
		    for (i = 0 ; i < nrows ; i++)
		    {
			row = Rows [i] ;
			if ((row >= 0) && (Frpos [row] >= 0))
			{
			    ASSERT (row < n_row) ;
			    Row_degree [row] -= ncolsleft ;
			    Woo [ngetrows] = i ;
			    Wm [ngetrows++] = Frpos [row] ;
			    
			    Rows [i] = EMPTY ;
			}
		    }
		    ASSERT (nrowsleft - ngetrows == extcdeg) ;
		    ASSERT (ngetrows == nrows_to_assemble) ;

		    if (ncols == ncolsleft)
		    {
			
			
			

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    ASSERT (col >= 0 && col < n_col) ;
#ifndef FIXQ
			    Col_degree [col] -= nrows_to_assemble ;
#endif
			    Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			    for (i = 0 ; i < nrows_to_assemble ; i++)
			    {
				
				ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
			    }
			    S += nrows ;
			}


		    }
		    else
		    {
			
			
			

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    ASSERT (col < n_col) ;
			    if (col >= 0)
			    {
#ifndef FIXQ
				Col_degree [col] -= nrows_to_assemble ;
#endif
				Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
				for (i = 0 ; i < nrows_to_assemble ; i++)
				{
				    
				    ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
				}
			    }
			    S += nrows ;
			}

		    }

		    

		    ep->nrowsleft = extcdeg ;
		}
	    }
	}
    }

    
    

    

    
    
    

    
    
    

    
    if (Work->any_skip)
    {
	row_assemble (Work->pivrow, Numeric, Work) ;
    }

    if (Work->do_scan2row)
    {
	for (i2 = Work->fscan_row ; i2 < fnrows ; i2++)
	{
	    
	    row = Work->NewRows [i2] ;
	    if (row < 0) row = FLIP (row) ;
	    ASSERT (row >= 0 && row < n_row) ;
	    if (!(row == Work->pivrow && Work->any_skip))
	    {
		
		row_assemble (row, Numeric, Work) ;
	    }
	}
    }

    
    
    

    
    if (Work->any_skip)
    {
	col_assemble (Work->pivcol, Numeric, Work) ;
    }

    if (Work->do_scan2col)
    {

	for (j2 = Work->fscan_col ; j2 < fncols ; j2++)
	{
	    
	    col = Work->NewCols [j2] ;
	    if (col < 0) col = FLIP (col) ;
	    if (!(col == Work->pivcol && Work->any_skip))
	    {
		
		col_assemble (col, Numeric, Work) ;
	    }
	}
    }
}

