/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    SparseLU_base.c
 * BRIEF:   稀疏LU基础
 *****************************************************************************/

#include "SparseLU_internal.h"
#include "SparseLU_function.h"


GLOBAL void SparseLU_defaults
(
    double Control [SparseLU_CONTROL]
)
{
    Int i ;

    if (!Control)
    {
	
	return ;
    }

    for (i = 0 ; i < SparseLU_CONTROL ; i++)
    {
	Control [i] = 0 ;
    }

    
    
    

    
    Control [SparseLU_PRL] = SparseLU_DEFAULT_PRL ;

    Control [SparseLU_DENSE_ROW] = SparseLU_DEFAULT_DENSE_ROW ;
    Control [SparseLU_DENSE_COL] = SparseLU_DEFAULT_DENSE_COL ;
    Control [SparseLU_AMD_DENSE] = SparseLU_DEFAULT_AMD_DENSE ;
    Control [SparseLU_STRATEGY] = SparseLU_DEFAULT_STRATEGY ;
    Control [SparseLU_AGGRESSIVE] = SparseLU_DEFAULT_AGGRESSIVE ;
    Control [SparseLU_SINGLETONS] = SparseLU_DEFAULT_SINGLETONS ;
    Control [SparseLU_ORDERING] = SparseLU_DEFAULT_ORDERING ;
    Control [SparseLU_PIVOT_TOLERANCE] = SparseLU_DEFAULT_PIVOT_TOLERANCE ;
    Control [SparseLU_SYM_PIVOT_TOLERANCE] = SparseLU_DEFAULT_SYM_PIVOT_TOLERANCE;
    Control [SparseLU_BLOCK_SIZE] = SparseLU_DEFAULT_BLOCK_SIZE ;
    Control [SparseLU_ALLOC_INIT] = SparseLU_DEFAULT_ALLOC_INIT ;
    Control [SparseLU_FRONT_ALLOC_INIT] = SparseLU_DEFAULT_FRONT_ALLOC_INIT ;
    Control [SparseLU_SCALE] = SparseLU_DEFAULT_SCALE ;

    
    Control [SparseLU_IRSTEP] = SparseLU_DEFAULT_IRSTEP ;

    
    
    

    
    Control [SparseLU_COMPILED_WITH_BLAS] = 1 ;
}



GLOBAL Int SparseLU_col_to_triplet
(
    Int n_col,
    const Int Ap [ ],
    Int Tj [ ]
)
{

    
    
    

    Int nz, j, p, p1, p2, length ;

    
    
    

    if (!Ap || !Tj)
    {
	return (SparseLU_ERROR_argument_missing) ;
    }
    if (n_col <= 0)
    {
	return (SparseLU_ERROR_n_nonpositive) ;
    }
    if (Ap [0] != 0)
    {
	return (SparseLU_ERROR_invalid_matrix) ;
    }
    nz = Ap [n_col] ;
    if (nz < 0)
    {
	return (SparseLU_ERROR_invalid_matrix) ;
    }

    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	length = p2 - p1 ;
	if (length < 0 || p2 > nz)
	{
	    return (SparseLU_ERROR_invalid_matrix) ;
	}
	for (p = p1 ; p < p2 ; p++)
	{
	    Tj [p] = j ;
	}
    }

    return (SparseLU_OK) ;
}


GLOBAL void SparseLU_free_numeric
(
    void **NumericHandle
)
{

    NumericType *Numeric ;
    if (!NumericHandle)
    {
	return ;
    }
    Numeric = *((NumericType **) NumericHandle) ;
    if (!Numeric)
    {
	return ;
    }

    
    (void) LU_free ((void *) Numeric->D) ;
    (void) LU_free ((void *) Numeric->Rperm) ;
    (void) LU_free ((void *) Numeric->Cperm) ;
    (void) LU_free ((void *) Numeric->Lpos) ;
    (void) LU_free ((void *) Numeric->Lilen) ;
    (void) LU_free ((void *) Numeric->Lip) ;
    (void) LU_free ((void *) Numeric->Upos) ;
    (void) LU_free ((void *) Numeric->Uilen) ;
    (void) LU_free ((void *) Numeric->Uip) ;

    
    (void) LU_free ((void *) Numeric->Rs) ;

    
    (void) LU_free ((void *) Numeric->Upattern) ;

    
    (void) LU_free ((void *) Numeric->Memory) ;
    (void) LU_free ((void *) Numeric) ;

    *NumericHandle = (void *) NULL ;
}



GLOBAL void SparseLU_free_symbolic
(
    void **SymbolicHandle
)
{

    SymbolicType *Symbolic ;
    if (!SymbolicHandle)
    {
	return ;
    }
    Symbolic = *((SymbolicType **) SymbolicHandle) ;
    if (!Symbolic)
    {
	return ;
    }

    (void) LU_free ((void *) Symbolic->Cperm_init) ;
    (void) LU_free ((void *) Symbolic->Rperm_init) ;
    (void) LU_free ((void *) Symbolic->Front_npivcol) ;
    (void) LU_free ((void *) Symbolic->Front_parent) ;
    (void) LU_free ((void *) Symbolic->Front_1strow) ;
    (void) LU_free ((void *) Symbolic->Front_leftmostdesc) ;
    (void) LU_free ((void *) Symbolic->Chain_start) ;
    (void) LU_free ((void *) Symbolic->Chain_maxrows) ;
    (void) LU_free ((void *) Symbolic->Chain_maxcols) ;
    (void) LU_free ((void *) Symbolic->Cdeg) ;
    (void) LU_free ((void *) Symbolic->Rdeg) ;

    
    (void) LU_free ((void *) Symbolic->Esize) ;

    
    (void) LU_free ((void *) Symbolic->Diagonal_map) ;

    (void) LU_free ((void *) Symbolic) ;
    *SymbolicHandle = (void *) NULL ;
}









PRIVATE Int rescale_determinant
(
    Entry *d_mantissa,
    double *d_exponent
)
{
    double d_abs ;

    ABS (d_abs, *d_mantissa) ;

    if (SCALAR_IS_ZERO (d_abs))
    {
	
	*d_exponent = 0 ;
	return (FALSE) ;
    }

    if (SCALAR_IS_NAN (d_abs))
    {
	
	return (FALSE) ;
    }

    while (d_abs < 1.)
    {
	SCALE (*d_mantissa, 10.0) ;
	*d_exponent = *d_exponent - 1.0 ;
	ABS (d_abs, *d_mantissa) ;
    }

    while (d_abs >= 10.)
    {
	SCALE (*d_mantissa, 0.1) ;
	*d_exponent = *d_exponent + 1.0 ;
	ABS (d_abs, *d_mantissa) ;
    }

    return (TRUE) ;
}





GLOBAL Int SparseLU_get_determinant
(
    double *Mx,
    double *Ex,
    void *NumericHandle,
    double User_Info [SparseLU_INFO]
)
{
    
    
    

    Entry d_mantissa, d_tmp ;
    double d_exponent, Info2 [SparseLU_INFO], one [2] = {1.0, 0.0}, d_sign ;
    Entry *D ;
    double *Info, *Rs ;
    NumericType *Numeric ;
    Int i, n, itmp, npiv, *Wi, *Rperm, *Cperm, do_scale ;

#ifndef NRECIPROCAL
    Int do_recip ;
#endif

    
    
    

    if (User_Info != (double *) NULL)
    {
	
	Info = User_Info ;
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

    Numeric = (NumericType *) NumericHandle ;
    if (!LU_valid_numeric (Numeric))
    {
	Info [SparseLU_STATUS] = SparseLU_ERROR_invalid_Numeric_object ;
	return (SparseLU_ERROR_invalid_Numeric_object) ;
    }

    if (Numeric->n_row != Numeric->n_col)
    {
	
	Info [SparseLU_STATUS] = SparseLU_ERROR_invalid_system ;
	return (SparseLU_ERROR_invalid_system) ;
    }

    if (Mx == (double *) NULL)
    {
	Info [SparseLU_STATUS] = SparseLU_ERROR_argument_missing ;
	return (SparseLU_ERROR_argument_missing) ;
    }

    n = Numeric->n_row ;

    
    
    

    Wi = (Int *) LU_malloc (n, sizeof (Int)) ;

    if (!Wi)
    {
	Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
	return (SparseLU_ERROR_out_of_memory) ;
    }

    
    
    

    Rs = Numeric->Rs ;		
    do_scale = (Rs != (double *) NULL) ;

#ifndef NRECIPROCAL
    do_recip = Numeric->do_recip ;
#endif

    d_mantissa = ((Entry *) one) [0] ;
    d_exponent = 0.0 ;
    D = Numeric->D ;

    
    for (i = 0 ; i < n ; i++)
    {
	MULT (d_tmp, d_mantissa, D [i]) ;
	d_mantissa = d_tmp ;

	if (!rescale_determinant (&d_mantissa, &d_exponent))
	{
	    
	    Info [SparseLU_STATUS] = SparseLU_WARNING_singular_matrix ;
	    
	    do_scale = FALSE ;
	    break ;
	}
    }

    
    if (do_scale)
    {
	for (i = 0 ; i < n ; i++)
	{
#ifndef NRECIPROCAL
	    if (do_recip)
	    {
		
		SCALE_DIV (d_mantissa, Rs [i]) ;
	    }
	    else
#endif
	    {
		
		SCALE (d_mantissa, Rs [i]) ;
	    }
	    if (!rescale_determinant (&d_mantissa, &d_exponent))
	    {
		
		Info [SparseLU_STATUS] = SparseLU_WARNING_singular_matrix ;
		break ;
	    }
	}
    }

    
    
    

    npiv = 0 ;
    Rperm = Numeric->Rperm ;

    for (i = 0 ; i < n ; i++)
    {
	Wi [i] = Rperm [i] ;
    }

    for (i = 0 ; i < n ; i++)
    {
	while (Wi [i] != i)
	{
	    itmp = Wi [Wi [i]] ;
	    Wi [Wi [i]] = Wi [i] ;
	    Wi [i] = itmp ;
	    npiv++ ;
	}
    }

    Cperm = Numeric->Cperm ;

    for (i = 0 ; i < n ; i++)
    {
	Wi [i] = Cperm [i] ;
    }

    for (i = 0 ; i < n ; i++)
    {
	while (Wi [i] != i)
	{
	    itmp = Wi [Wi [i]] ;
	    Wi [Wi [i]] = Wi [i] ;
	    Wi [i] = itmp ;
	    npiv++ ;
	}
    }

    
    d_sign = (npiv % 2) ? -1. : 1. ;

    
    
    

    (void) LU_free ((void *) Wi) ;

    
    
    

    if (Ex == (double *) NULL)
    {
	
	SCALE (d_mantissa, pow (10.0, d_exponent)) ;
    }
    else
    {
	Ex [0] = d_exponent ;
    }

    Mx [0] = d_sign * REAL_COMPONENT (d_mantissa) ;

    
    if (d_exponent + 1.0 > log10 (DBL_MAX))
    {
	Info [SparseLU_STATUS] = SparseLU_WARNING_determinant_overflow ;
    }
    else if (d_exponent - 1.0 < log10 (DBL_MIN))
    {
	Info [SparseLU_STATUS] = SparseLU_WARNING_determinant_underflow ;
    }

    return (SparseLU_OK) ;
}



GLOBAL Int SparseLU_get_lunz
(
    Int *lnz,
    Int *unz,
    Int *n_row,
    Int *n_col,
    Int *nz_udiag,
    void *NumericHandle
)
{
    NumericType *Numeric ;

    Numeric = (NumericType *) NumericHandle ;

    if (!LU_valid_numeric (Numeric))
    {
	return (SparseLU_ERROR_invalid_Numeric_object) ;
    }
    if (!lnz || !unz || !n_row || !n_col || !nz_udiag)
    {
	return (SparseLU_ERROR_argument_missing) ;
    }

    *n_row = Numeric->n_row ;
    *n_col = Numeric->n_col ;

    
    *lnz = Numeric->lnz + MIN (Numeric->n_row, Numeric->n_col) ;

    
    *unz = Numeric->unz + Numeric->nnzpiv ;

    
    *nz_udiag = Numeric->nnzpiv ;

    return (SparseLU_OK) ;
}


GLOBAL Int SparseLU_get_symbolic
(
    Int *p_n_row,
    Int *p_n_col,
    Int *p_n1,			
    Int *p_nz,
    Int *p_nfr,
    Int *p_nchains,
    Int P [ ],
    Int Q [ ],
    Int Front_npivcol [ ],
    Int Front_parent [ ],
    Int Front_1strow [ ],
    Int Front_leftmostdesc [ ],
    Int Chain_start [ ],
    Int Chain_maxrows [ ],
    Int Chain_maxcols [ ],
    void *SymbolicHandle
)
{
    SymbolicType *Symbolic ;
    Int k, n_row, n_col, n1, nfr, nchains, *p ;

    
    
    

    Symbolic = (SymbolicType *) SymbolicHandle ;
    if (!LU_valid_symbolic (Symbolic))
    {
	return (SparseLU_ERROR_invalid_Symbolic_object) ;
    }

    
    
    

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n1 = Symbolic->n1 ;
    nfr = Symbolic->nfr ;
    nchains = Symbolic->nchains ;

    if (p_n_row)
    {
	*p_n_row = n_row ;
    }

    if (p_n_col)
    {
	*p_n_col = n_col ;
    }

    if (p_n1)
    {
	*p_n1 = n1 ;
    }

    if (p_nz)
    {
	*p_nz = Symbolic->nz ;
    }

    if (p_nfr)
    {
	*p_nfr = nfr ;
    }

    if (p_nchains)
    {
	*p_nchains = nchains ;
    }

    if (P != (Int *) NULL)
    {
	Int *Rperm_init, *Diagonal_map ;
	Rperm_init = Symbolic->Rperm_init ;
	Diagonal_map = Symbolic->Diagonal_map ;
	if (Diagonal_map != (Int *) NULL)
	{
	    ASSERT (n_row == n_col) ;
	    
	    for (k = 0 ; k < n_row ; k++)
	    {
		P [k] = Rperm_init [Diagonal_map [k]] ;
	    }
	}
	else
	{
	    
	    for (k = 0 ; k < n_row ; k++)
	    {
		P [k] = Rperm_init [k] ;
	    }
	}
    }

    if (Q != (Int *) NULL)
    {
	p = Symbolic->Cperm_init ;
	for (k = 0 ; k < n_col ; k++)
	{
	    Q [k] = p [k] ;
	}
    }

    if (Front_npivcol != (Int *) NULL)
    {
	p = Symbolic->Front_npivcol ;
	for (k = 0 ; k <= nfr ; k++)
	{
	    Front_npivcol [k] = p [k] ;
	}
    }

    if (Front_parent != (Int *) NULL)
    {
	p = Symbolic->Front_parent ;
	for (k = 0 ; k <= nfr ; k++)
	{
	    Front_parent [k] = p [k] ;
	}
    }

    if (Front_1strow != (Int *) NULL)
    {
	p = Symbolic->Front_1strow ;
	for (k = 0 ; k <= nfr ; k++)
	{
	    Front_1strow [k] = p [k] ;
	}
    }

    if (Front_leftmostdesc != (Int *) NULL)
    {
	p = Symbolic->Front_leftmostdesc ;
	for (k = 0 ; k <= nfr ; k++)
	{
	    Front_leftmostdesc [k] = p [k] ;
	}
    }

    if (Chain_start != (Int *) NULL)
    {
	p = Symbolic->Chain_start ;
	for (k = 0 ; k <= nchains ; k++)
	{
	    Chain_start [k] = p [k] ;
	}
    }

    if (Chain_maxrows != (Int *) NULL)
    {
	p = Symbolic->Chain_maxrows ;
	for (k = 0 ; k < nchains ; k++)
	{
	    Chain_maxrows [k] = p [k] ;
	}
	Chain_maxrows [nchains] = 0 ;
    }

    if (Chain_maxcols != (Int *) NULL)
    {
	p = Symbolic->Chain_maxcols ;
	for (k = 0 ; k < nchains ; k++)
	{
	    Chain_maxcols [k] = p [k] ;
	}
	Chain_maxcols [nchains] = 0 ;
    }

    return (SparseLU_OK) ;
}



PRIVATE void get_L
(
    Int Lp [ ],
    Int Lj [ ],
    double Lx [ ],
    NumericType *Numeric,
    Int Pattern [ ],
    Int Wi [ ]
) ;

PRIVATE void get_U
(
    Int Up [ ],
    Int Ui [ ],
    double Ux [ ],
    NumericType *Numeric,
    Int Pattern [ ],
    Int Wi [ ]
) ;





GLOBAL Int SparseLU_get_numeric
(
    Int Lp [ ],
    Int Lj [ ],
    double Lx [ ],
    Int Up [ ],
    Int Ui [ ],
    double Ux [ ],
    Int P [ ],
    Int Q [ ],
    double Dx [ ],
    Int *p_do_recip,
    double Rs [ ],
    void *NumericHandle
)
{

    
    
    

    NumericType *Numeric ;
    Int getL, getU, *Rperm, *Cperm, k, nn, n_row, n_col, *Wi, *Pattern,
	n_inner ;
    double *Rs1 ;
    Entry *D ;

    Wi = (Int *) NULL ;
    Pattern = (Int *) NULL ;

    
    
    

    Numeric = (NumericType *) NumericHandle ;
    if (!LU_valid_numeric (Numeric))
    {
	return (SparseLU_ERROR_invalid_Numeric_object) ;
    }

    n_row = Numeric->n_row ;
    n_col = Numeric->n_col ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;

    
    
    

    getL = Lp && Lj && Lx ;
    getU = Up && Ui && Ux ;

    if (getL || getU)
    {
	Wi = (Int *) LU_malloc (nn, sizeof (Int)) ;
	Pattern = (Int *) LU_malloc (nn, sizeof (Int)) ;
	if (!Wi || !Pattern)
	{
	    (void) LU_free ((void *) Wi) ;
	    (void) LU_free ((void *) Pattern) ;
	    return (SparseLU_ERROR_out_of_memory) ;
	}
    }

    
    
    

    if (P != (Int *) NULL)
    {
	Rperm = Numeric->Rperm ;
	for (k = 0 ; k < n_row ; k++)
	{
	    P [k] = Rperm [k] ;
	}
    }

    if (Q != (Int *) NULL)
    {
	Cperm = Numeric->Cperm ;
	for (k = 0 ; k < n_col ; k++)
	{
	    Q [k] = Cperm [k] ;
	}
    }

    if (getL)
    {
	get_L (Lp, Lj, Lx,
	    Numeric, Pattern, Wi) ;
    }

    if (getU)
    {
	get_U (Up, Ui, Ux,
	    Numeric, Pattern, Wi) ;
    }

    if (Dx != (double *) NULL)
    {
	D = Numeric->D ;
	{
	    D = Numeric->D ;
	    for (k = 0 ; k < n_inner ; k++)
	    {
		Dx [k] = D [k] ;
	    }
	}
    }

    
    if (p_do_recip != (Int *) NULL)
    {
#ifndef NRECIPROCAL
	*p_do_recip = Numeric->do_recip ;
#else
	*p_do_recip = FALSE ;
#endif
    }

    if (Rs != (double *) NULL)
    {
	Rs1 = Numeric->Rs ;
	if (Rs1 == (double *) NULL)
	{
	    
	    for (k = 0 ; k < n_row ; k++)
	    {
		Rs [k] = 1.0 ;
	    }
	}
	else
	{
	    for (k = 0 ; k < n_row ; k++)
	    {
		Rs [k] = Rs1 [k] ;
	    }
	}
    }

    
    
    

    (void) LU_free ((void *) Wi) ;
    (void) LU_free ((void *) Pattern) ;
    ASSERT (LU_malloc_count == init_count) ;

    return (SparseLU_OK) ;
}









PRIVATE void get_L
(
    Int Lp [ ],		
    Int Lj [ ],		
    double Lx [ ],	
    NumericType *Numeric,
    Int Pattern [ ],	
    Int Wi [ ]		
)
{
    
    
    

    Entry value ;
    Entry *xp, *Lval ;
    Int deg, *ip, j, row, n_row, n_col, n_inner, *Lpos, *Lilen, *Lip, p, llen,
        lnz2, lp, newLchain, k, pos, npiv, *Li, n1 ;

    
    
    

    n_row = Numeric->n_row ;
    n_col = Numeric->n_col ;
    n_inner = MIN (n_row, n_col) ;
    npiv = Numeric->npiv ;
    n1 = Numeric->n1 ;
    Lpos = Numeric->Lpos ;
    Lilen = Numeric->Lilen ;
    Lip = Numeric->Lip ;
    deg = 0 ;

    
    
    

#pragma ivdep
    for (row = 0 ; row < n_inner ; row++)
    {
	
	Wi [row] = 1 ;
    }
#pragma ivdep
    for (row = n_inner ; row < n_row ; row++)
    {
	Wi [row] = 0 ;
    }

    
    for (k = 0 ; k < n1 ; k++)
    {
	deg = Lilen [k] ;
	if (deg > 0)
	{
	    lp = Lip [k] ;
	    Li = (Int *) (Numeric->Memory + lp) ;
	    lp += UNITS (Int, deg) ;
	    Lval = (Entry *) (Numeric->Memory + lp) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		row = Li [j] ;
		value = Lval [j] ;
		if (IS_NONZERO (value))
		{
		    Wi [row]++ ;
		}
	    }
	}
    }

    
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

	xp = (Entry *) (Numeric->Memory + lp + UNITS (Int, llen)) ;

	for (j = 0 ; j < deg ; j++)
	{
	    row = Pattern [j] ;
	    value = *xp++ ;
	    if (IS_NONZERO (value))
	    {
		Wi [row]++ ;
	    }
	}
    }

    
    
    

    
    lnz2 = 0 ;
    for (row = 0 ; row < n_row ; row++)
    {
	Lp [row] = lnz2 ;
	lnz2 += Wi [row] ;
	Wi [row] = Lp [row] ;
    }
    Lp [n_row] = lnz2 ;
    ASSERT (Numeric->lnz + n_inner == lnz2) ;

    
    for (k = 0 ; k < n1 ; k++)
    {
	deg = Lilen [k] ;
	if (deg > 0)
	{
	    lp = Lip [k] ;
	    Li = (Int *) (Numeric->Memory + lp) ;
	    lp += UNITS (Int, deg) ;
	    Lval = (Entry *) (Numeric->Memory + lp) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		row = Li [j] ;
		value = Lval [j] ;
		if (IS_NONZERO (value))
		{
		    p = Wi [row]++ ;
		    Lj [p] = k ;
		    Lx [p] = value ;
		}
	    }
	}
    }

    
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

	xp = (Entry *) (Numeric->Memory + lp + UNITS (Int, llen)) ;

	for (j = 0 ; j < deg ; j++)
	{
	    row = Pattern [j] ;
	    value = *xp++ ;
	    if (IS_NONZERO (value))
	    {
		p = Wi [row]++ ;
		Lj [p] = k ;
		Lx [p] = value ;
	    }
	}
    }

    
    for (row = 0 ; row < n_inner ; row++)
    {
	p = Wi [row]++ ;
	Lj [p] = row ;
	Lx [p] = 1. ;

    }

}









PRIVATE void get_U
(
    Int Up [ ],		
    Int Ui [ ],		
    double Ux [ ],	
    NumericType *Numeric,
    Int Pattern [ ],	
    Int Wi [ ]		
)
{
    
    
    

    Entry value ;
    Entry *xp, *D, *Uval ;
    Int deg, j, *ip, col, *Upos, *Uilen, *Uip, n_col, ulen, *Usi,
        unz2, p, k, up, newUchain, pos, npiv, n1 ;

    
    
    

    n_col = Numeric->n_col ;
    n1 = Numeric->n1 ;
    npiv = Numeric->npiv ;
    Upos = Numeric->Upos ;
    Uilen = Numeric->Uilen ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;

    
    
    

    for (col = 0 ; col < npiv ; col++)
    {
	
	Wi [col] = IS_NONZERO (D [col]) ;
    }
    for (col = npiv ; col < n_col ; col++)
    {
	
	Wi [col] = 0 ;
    }

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

	for (j = 0 ; j < deg ; j++)
	{
	    col = Pattern [j] ;
	    value = *xp++ ;
	    if (IS_NONZERO (value))
	    {
		Wi [col]++ ;
	    }
	}

	
	
	

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
	if (deg > 0)
	{
	    up = Uip [k] ;
	    Usi = (Int *) (Numeric->Memory + up) ;
	    up += UNITS (Int, deg) ;
	    Uval = (Entry *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		col = Usi [j] ;
		value = Uval [j] ;
		if (IS_NONZERO (value))
		{
		    Wi [col]++ ;
		}
	    }
	}
    }

    
    
    

    
    unz2 = 0 ;
    for (col = 0 ; col < n_col ; col++)
    {
	Up [col] = unz2 ;
	unz2 += Wi [col] ;
    }
    Up [n_col] = unz2 ;

    for (col = 0 ; col < n_col ; col++)
    {
	Wi [col] = Up [col+1] ;
    }

    
    for (col = 0 ; col < npiv ; col++)
    {
	if (IS_NONZERO (D [col]))
	{
	    p = --(Wi [col]) ;
	    Ui [p] = col ;
	    Ux [p] = D [col] ;
	}
    }

    

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

	xp += deg ;
	for (j = deg-1 ; j >= 0 ; j--)
	{
	    col = Pattern [j] ;
	    value = *(--xp) ;
	    if (IS_NONZERO (value))
	    {
		p = --(Wi [col]) ;
		Ui [p] = k ;
		Ux [p] = value ;

	    }
	}

	
	
	

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
	if (deg > 0)
	{
	    up = Uip [k] ;
	    Usi = (Int *) (Numeric->Memory + up) ;
	    up += UNITS (Int, deg) ;
	    Uval = (Entry *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		col = Usi [j] ;
		value = Uval [j] ;
		if (IS_NONZERO (value))
		{
		    p = --(Wi [col]) ;
		    Ui [p] = k ;
		    Ux [p] = value ;
		}
	    }
	}
    }

}



GLOBAL Int SparseLU_scale
(
    double Xx [ ],
    const double Bx [ ],
    void *NumericHandle
)
{
    
    
    

    NumericType *Numeric ;
    Int n, i ;
    double *Rs ;

    Numeric = (NumericType *) NumericHandle ;
    if (!LU_valid_numeric (Numeric))
    {
	return (SparseLU_ERROR_invalid_Numeric_object) ;
    }

    n = Numeric->n_row ;
    Rs = Numeric->Rs ;

    if (!Xx || !Bx)
    {
	return (SparseLU_ERROR_argument_missing) ;
    }

    
    
    

    if (Rs != (double *) NULL)
    {
#ifndef NRECIPROCAL
	if (Numeric->do_recip)
	{
	    
	    for (i = 0 ; i < n ; i++)
	    {
		Xx [i] = Bx [i] * Rs [i] ;
	    }
	}
	else
#endif
	{
	    
	    for (i = 0 ; i < n ; i++)
	    {
		Xx [i] = Bx [i] / Rs [i] ;
	    }
	}
    }
    else
    {
	
	for (i = 0 ; i < n ; i++)
	{
	    Xx [i] = Bx [i] ;
	}
    }

    return (SparseLU_OK) ;
}





GLOBAL Int SparseLU_transpose
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],	
    const Int Ai [ ],	
    const double Ax [ ], 

    const Int P [ ],	
			
			

    const Int Q [ ],	
			
			

    Int Rp [ ],		
    Int Ri [ ],		
    double Rx [ ]	
)
{

    
    
    

    Int status, *W, nn ;

    
    
    

    nn = MAX (n_row, n_col) ;
    nn = MAX (nn, 1) ;
    W = (Int *) LU_malloc (nn, sizeof (Int)) ;
    if (!W)
    {
	return (SparseLU_ERROR_out_of_memory) ;
    }

    
    
    

    status = LU_transpose (n_row, n_col, Ap, Ai, Ax, P, Q, n_col, Rp, Ri, Rx,
	W, TRUE
	) ;

    
    
    

    (void) LU_free ((void *) W) ;

    return (status) ;
}





GLOBAL Int SparseLU_triplet_to_col
(
    Int n_row,
    Int n_col,
    Int nz,
    const Int Ti [ ],		
    const Int Tj [ ],		
    const double Tx [ ],	
    Int Ap [ ],			
    Int Ai [ ],			
    double Ax [ ]		
    , Int Map [ ]		
)
{

    
    
    

    Int *RowCount, *Rp, *Rj, *W, nn, do_values, do_map, *Map2, status ;
    double *Rx ;

    
    
    

    if (!Ai || !Ap || !Ti || !Tj)
    {
	return (SparseLU_ERROR_argument_missing) ;
    }

    if (n_row <= 0 || n_col <= 0)		
    {
	return (SparseLU_ERROR_n_nonpositive) ;
    }

    if (nz < 0)		
    {
	return (SparseLU_ERROR_invalid_matrix) ;
    }

    nn = MAX (n_row, n_col) ;

    
    
    

    Rx = (double *) NULL ;

    do_values = Ax && Tx ;

    if (do_values)
    {
	Rx = (double *) LU_malloc (nz+1, sizeof (double)) ;
	if (!Rx)
	{
	    return (SparseLU_ERROR_out_of_memory) ;
	}
    }

    do_map = (Map != (Int *) NULL) ;
    Map2 = (Int *) NULL ;
    if (do_map)
    {
	Map2 = (Int *) LU_malloc (nz+1, sizeof (Int)) ;
	if (!Map2)
	{
	    (void) LU_free ((void *) Rx) ;
	    return (SparseLU_ERROR_out_of_memory) ;
	}
    }

    Rj = (Int *) LU_malloc (nz+1, sizeof (Int)) ;
    Rp = (Int *) LU_malloc (n_row+1, sizeof (Int)) ;
    RowCount = (Int *) LU_malloc (n_row, sizeof (Int)) ;
    W = (Int *) LU_malloc (nn, sizeof (Int)) ;
    if (!Rj || !Rp || !RowCount || !W)
    {
	(void) LU_free ((void *) Rx) ;
	(void) LU_free ((void *) Map2) ;
	(void) LU_free ((void *) Rp) ;
	(void) LU_free ((void *) Rj) ;
	(void) LU_free ((void *) RowCount) ;
	(void) LU_free ((void *) W) ;
	return (SparseLU_ERROR_out_of_memory) ;
    }

    
    
    

    if (do_map)
    {
	if (do_values)
	{
	    status = LU_triplet_map_x (n_row, n_col, nz, Ti, Tj, Ap, Ai, Rp,
		Rj, W, RowCount, Tx, Ax, Rx
		, Map, Map2) ;
	}
	else
	{
	    status = LU_triplet_map_nox (n_row, n_col, nz, Ti, Tj, Ap, Ai, Rp,
		Rj, W, RowCount, Map, Map2) ;
	}
    }
    else
    {
	if (do_values)
	{
	    status = LU_triplet_nomap_x (n_row, n_col, nz, Ti, Tj, Ap, Ai, Rp,
		Rj, W, RowCount , Tx, Ax, Rx
		) ;
	}
	else
	{
	    status = LU_triplet_nomap_nox (n_row, n_col, nz, Ti, Tj, Ap, Ai,
		Rp, Rj, W, RowCount) ;
	}
    }

    
    
    

    (void) LU_free ((void *) Rx) ;
    (void) LU_free ((void *) Map2) ;
    (void) LU_free ((void *) Rp) ;
    (void) LU_free ((void *) Rj) ;
    (void) LU_free ((void *) RowCount) ;
    (void) LU_free ((void *) W) ;
    ASSERT (LU_malloc_count == init_count) ;

    return (status) ;
}



#define WRITE(object,type,n) \
{ \
    ASSERT (object != (type *) NULL) ; \
    if (fwrite (object, sizeof (type), n, f) != (size_t) n) \
    { \
	fclose (f) ; \
	return (SparseLU_ERROR_file_IO) ; \
    } \
    fflush (f) ; \
}





GLOBAL Int SparseLU_save_symbolic
(
    void *SymbolicHandle,
    char *user_filename
)
{
    SymbolicType *Symbolic ;
    char *filename ;
    FILE *f ;

    
    Symbolic = (SymbolicType *) SymbolicHandle ;

    
    if (!LU_valid_symbolic (Symbolic))
    {
	return (SparseLU_ERROR_invalid_Symbolic_object) ;
    }

    
    if (user_filename == (char *) NULL)
    {
	filename = "symbolic.lu" ;
    }
    else
    {
	filename = user_filename ;
    }
    f = fopen (filename, "wb") ;
    if (!f)
    {
	return (SparseLU_ERROR_file_IO) ;
    }

    
    WRITE (Symbolic,                     SymbolicType, 1) ;
    WRITE (Symbolic->Cperm_init,         Int, Symbolic->n_col+1) ;
    WRITE (Symbolic->Rperm_init,         Int, Symbolic->n_row+1) ;
    WRITE (Symbolic->Front_npivcol,      Int, Symbolic->nfr+1) ;
    WRITE (Symbolic->Front_parent,       Int, Symbolic->nfr+1) ;
    WRITE (Symbolic->Front_1strow,       Int, Symbolic->nfr+1) ;
    WRITE (Symbolic->Front_leftmostdesc, Int, Symbolic->nfr+1) ;
    WRITE (Symbolic->Chain_start,        Int, Symbolic->nchains+1) ;
    WRITE (Symbolic->Chain_maxrows,      Int, Symbolic->nchains+1) ;
    WRITE (Symbolic->Chain_maxcols,      Int, Symbolic->nchains+1) ;
    WRITE (Symbolic->Cdeg,               Int, Symbolic->n_col+1) ;
    WRITE (Symbolic->Rdeg,               Int, Symbolic->n_row+1) ;
    if (Symbolic->esize > 0)
    {
	
	WRITE (Symbolic->Esize, Int, Symbolic->esize) ;
    }
    if (Symbolic->prefer_diagonal)
    {
	
	WRITE (Symbolic->Diagonal_map, Int, Symbolic->n_col+1) ;
    }

    
    fclose (f) ;

    return (SparseLU_OK) ;
}






GLOBAL Int SparseLU_save_numeric
(
    void *NumericHandle,
    char *user_filename
)
{
    NumericType *Numeric ;
    char *filename ;
    FILE *f ;

    
    Numeric = (NumericType *) NumericHandle ;

    
    if (!LU_valid_numeric (Numeric))
    {
	return (SparseLU_ERROR_invalid_Numeric_object) ;
    }

    
    if (user_filename == (char *) NULL)
    {
	filename = "numeric.lu" ;
    }
    else
    {
	filename = user_filename ;
    }
    f = fopen (filename, "wb") ;
    if (!f)
    {
	return (SparseLU_ERROR_file_IO) ;
    }

    
    WRITE (Numeric,        NumericType, 1) ;
    WRITE (Numeric->D,     Entry, MIN (Numeric->n_row, Numeric->n_col)+1) ;
    WRITE (Numeric->Rperm, Int,   Numeric->n_row+1) ;
    WRITE (Numeric->Cperm, Int,   Numeric->n_col+1) ;
    WRITE (Numeric->Lpos,  Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Lilen, Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Lip,   Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Upos,  Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Uilen, Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Uip,   Int,   Numeric->npiv+1) ;
    if (Numeric->scale != SparseLU_SCALE_NONE)
    {
	WRITE (Numeric->Rs, double, Numeric->n_row) ;
    }
    if (Numeric->ulen > 0)
    {
	WRITE (Numeric->Upattern, Int, Numeric->ulen+1) ;
    }
    
    WRITE (Numeric->Memory, Unit, Numeric->size) ;

    
    fclose (f) ;

    return (SparseLU_OK) ;
}



#define READ_SYMBOLIC(object,type,n) \
{ \
    object = (type *) LU_malloc (n, sizeof (type)) ; \
    if (object == (type *) NULL) \
    { \
	SparseLU_free_symbolic ((void **) &Symbolic) ; \
	fclose (f) ; \
	return (SparseLU_ERROR_out_of_memory) ; \
    } \
    if (fread (object, sizeof (type), n, f) != (size_t) n) \
    { \
	SparseLU_free_symbolic ((void **) &Symbolic) ; \
	fclose (f) ; \
	return (SparseLU_ERROR_file_IO) ; \
    } \
    if (ferror (f)) \
    { \
	SparseLU_free_symbolic ((void **) &Symbolic) ; \
	fclose (f) ; \
	return (SparseLU_ERROR_file_IO) ; \
    } \
}





GLOBAL Int SparseLU_load_symbolic
(
    void **SymbolicHandle,
    char *user_filename
)
{
    SymbolicType *Symbolic ;
    char *filename ;
    FILE *f ;

    *SymbolicHandle = (void *) NULL ;

    
    
    

    if (user_filename == (char *) NULL)
    {
	filename = "symbolic.lu" ;
    }
    else
    {
	filename = user_filename ;
    }
    f = fopen (filename, "rb") ;
    if (!f)
    {
	return (SparseLU_ERROR_file_IO) ;
    }

    
    
    

    Symbolic = (SymbolicType *) LU_malloc (1, sizeof (SymbolicType)) ;
    if (Symbolic == (SymbolicType *) NULL)
    {
	fclose (f) ;
	return (SparseLU_ERROR_out_of_memory) ;
    }
    if (fread (Symbolic, sizeof (SymbolicType), 1, f) != 1)
    {
	(void) LU_free ((void *) Symbolic) ;
	fclose (f) ;
	return (SparseLU_ERROR_file_IO) ;
    }
    if (ferror (f))
    {
	(void) LU_free ((void *) Symbolic) ;
	fclose (f) ;
	return (SparseLU_ERROR_file_IO) ;
    }

    if (Symbolic->valid != SYMBOLIC_VALID || Symbolic->n_row <= 0 ||
	Symbolic->n_col <= 0 || Symbolic->nfr < 0 || Symbolic->nchains < 0 ||
	Symbolic->esize < 0)
    {
	
	(void) LU_free ((void *) Symbolic) ;
	fclose (f) ;
	return (SparseLU_ERROR_invalid_Symbolic_object) ;
    }

    Symbolic->Cperm_init         = (Int *) NULL ;
    Symbolic->Rperm_init         = (Int *) NULL ;
    Symbolic->Front_npivcol      = (Int *) NULL ;
    Symbolic->Front_parent       = (Int *) NULL ;
    Symbolic->Front_1strow       = (Int *) NULL ;
    Symbolic->Front_leftmostdesc = (Int *) NULL ;
    Symbolic->Chain_start        = (Int *) NULL ;
    Symbolic->Chain_maxrows      = (Int *) NULL ;
    Symbolic->Chain_maxcols      = (Int *) NULL ;
    Symbolic->Cdeg               = (Int *) NULL ;
    Symbolic->Rdeg               = (Int *) NULL ;
    Symbolic->Esize              = (Int *) NULL ;
    Symbolic->Diagonal_map       = (Int *) NULL ;

    

    
    
    

    READ_SYMBOLIC (Symbolic->Cperm_init,         Int, Symbolic->n_col+1) ;
    READ_SYMBOLIC (Symbolic->Rperm_init,         Int, Symbolic->n_row+1) ;
    READ_SYMBOLIC (Symbolic->Front_npivcol,      Int, Symbolic->nfr+1) ;
    READ_SYMBOLIC (Symbolic->Front_parent,       Int, Symbolic->nfr+1) ;
    READ_SYMBOLIC (Symbolic->Front_1strow,       Int, Symbolic->nfr+1) ;
    READ_SYMBOLIC (Symbolic->Front_leftmostdesc, Int, Symbolic->nfr+1) ;
    READ_SYMBOLIC (Symbolic->Chain_start,        Int, Symbolic->nchains+1) ;
    READ_SYMBOLIC (Symbolic->Chain_maxrows,      Int, Symbolic->nchains+1) ;
    READ_SYMBOLIC (Symbolic->Chain_maxcols,      Int, Symbolic->nchains+1) ;
    READ_SYMBOLIC (Symbolic->Cdeg,               Int, Symbolic->n_col+1) ;
    READ_SYMBOLIC (Symbolic->Rdeg,               Int, Symbolic->n_row+1) ;
    if (Symbolic->esize > 0)
    {
	
	READ_SYMBOLIC (Symbolic->Esize, Int, Symbolic->esize) ;
    }
    if (Symbolic->prefer_diagonal)
    {
	
	READ_SYMBOLIC (Symbolic->Diagonal_map, Int, Symbolic->n_col+1) ;
    }

    
    fclose (f) ;

    
    if (!LU_valid_symbolic (Symbolic))
    {
	SparseLU_free_symbolic ((void **) &Symbolic) ;
	return (SparseLU_ERROR_invalid_Symbolic_object) ;
    }

    *SymbolicHandle = (void *) Symbolic ;
    return (SparseLU_OK) ;
}



#define READ_NUMERIC(object,type,n) \
{ \
    object = (type *) LU_malloc (n, sizeof (type)) ; \
    if (object == (type *) NULL) \
    { \
	SparseLU_free_numeric ((void **) &Numeric) ; \
	fclose (f) ; \
	return (SparseLU_ERROR_out_of_memory) ; \
    } \
    if (fread (object, sizeof (type), n, f) != (size_t) n) \
    { \
	SparseLU_free_numeric ((void **) &Numeric) ; \
	fclose (f) ; \
	return (SparseLU_ERROR_file_IO) ; \
    } \
    if (ferror (f)) \
    { \
	SparseLU_free_numeric ((void **) &Numeric) ; \
	fclose (f) ; \
	return (SparseLU_ERROR_file_IO) ; \
    } \
}






GLOBAL Int SparseLU_load_numeric
(
    void **NumericHandle,
    char *user_filename
)
{
    NumericType *Numeric ;
    char *filename ;
    FILE *f ;

    *NumericHandle = (void *) NULL ;

    
    
    

    if (user_filename == (char *) NULL)
    {
	filename = "numeric.lu" ;
    }
    else
    {
	filename = user_filename ;
    }
    f = fopen (filename, "rb") ;
    if (!f)
    {
	return (SparseLU_ERROR_file_IO) ;
    }

    
    
    

    Numeric = (NumericType *) LU_malloc (1, sizeof (NumericType)) ;
    if (Numeric == (NumericType *) NULL)
    {
	fclose (f) ;
	return (SparseLU_ERROR_out_of_memory) ;
    }
    if (fread (Numeric, sizeof (NumericType), 1, f) != 1)
    {
	(void) LU_free ((void *) Numeric) ;
	fclose (f) ;
	return (SparseLU_ERROR_file_IO) ;
    }
    if (ferror (f))
    {
	(void) LU_free ((void *) Numeric) ;
	fclose (f) ;
	return (SparseLU_ERROR_file_IO) ;
    }

    if (Numeric->valid != NUMERIC_VALID || Numeric->n_row <= 0 ||
	Numeric->n_col <= 0 || Numeric->npiv < 0 || Numeric->ulen < 0 ||
	Numeric->size <= 0)
    {
	
	(void) LU_free ((void *) Numeric) ;
	fclose (f) ;
	return (SparseLU_ERROR_invalid_Numeric_object) ;
    }

    Numeric->D        = (Entry *) NULL ;
    Numeric->Rperm    = (Int *) NULL ;
    Numeric->Cperm    = (Int *) NULL ;
    Numeric->Lpos     = (Int *) NULL ;
    Numeric->Lilen    = (Int *) NULL ;
    Numeric->Lip      = (Int *) NULL ;
    Numeric->Upos     = (Int *) NULL ;
    Numeric->Uilen    = (Int *) NULL ;
    Numeric->Uip      = (Int *) NULL ;
    Numeric->Rs       = (double *) NULL ;
    Numeric->Memory   = (Unit *) NULL ;
    Numeric->Upattern = (Int *) NULL ;

    

    
    
    

    READ_NUMERIC (Numeric->D,     Entry, MIN (Numeric->n_row, Numeric->n_col)+1) ;
    READ_NUMERIC (Numeric->Rperm, Int,   Numeric->n_row+1) ;
    READ_NUMERIC (Numeric->Cperm, Int,   Numeric->n_col+1) ;
    READ_NUMERIC (Numeric->Lpos,  Int,   Numeric->npiv+1) ;
    READ_NUMERIC (Numeric->Lilen, Int,   Numeric->npiv+1) ;
    READ_NUMERIC (Numeric->Lip,   Int,   Numeric->npiv+1) ;
    READ_NUMERIC (Numeric->Upos,  Int,   Numeric->npiv+1) ;
    READ_NUMERIC (Numeric->Uilen, Int,   Numeric->npiv+1) ;
    READ_NUMERIC (Numeric->Uip,   Int,   Numeric->npiv+1) ;
    if (Numeric->scale != SparseLU_SCALE_NONE)
    {
	READ_NUMERIC (Numeric->Rs, double, Numeric->n_row) ;
    }
    if (Numeric->ulen > 0)
    {
	READ_NUMERIC (Numeric->Upattern, Int, Numeric->ulen+1) ;
    }
    READ_NUMERIC (Numeric->Memory, Unit, Numeric->size) ;

    
    fclose (f) ;

    
    if (!LU_valid_numeric (Numeric))
    {
	SparseLU_free_numeric ((void **) &Numeric) ;
	return (SparseLU_ERROR_invalid_Numeric_object) ;
    }

    *NumericHandle = (void *) Numeric ;
    return (SparseLU_OK) ;
}

