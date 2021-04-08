
/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    SparseLU_internal.h
 * BRIEF:   稀疏LU内部函数头文件
 *****************************************************************************/

#ifndef _LU_INTERNAL
#define _LU_INTERNAL

#include <float.h>
#include <string.h>

#if defined (__sun) || defined (MSOL2) || defined (ARCH_SOL2)
#define LU_SOL2
#define SparseLU_ARCHITECTURE "Sun Solaris"

#elif defined (__sgi) || defined (MSGI) || defined (ARCH_SGI)
#define LU_SGI
#define SparseLU_ARCHITECTURE "SGI Irix"

#elif defined (__linux) || defined (MGLNX86) || defined (ARCH_GLNX86)
#define LU_LINUX
#define SparseLU_ARCHITECTURE "Linux"

#elif defined (__APPLE__)
#define LU_MAC
#define SparseLU_ARCHITECTURE "Mac"

#elif defined (_AIX) || defined (MIBM_RS) || defined (ARCH_IBM_RS)
#define LU_AIX
#define SparseLU_ARCHITECTURE "IBM AIX"

#elif defined (__alpha) || defined (MALPHA) || defined (ARCH_ALPHA)
#define LU_ALPHA
#define SparseLU_ARCHITECTURE "Compaq Alpha"

#elif defined (_WIN32) || defined (WIN32)
#if defined (__MINGW32__)
#define LU_MINGW
#elif defined (__CYGWIN32__)
#define LU_CYGWIN
#else
#define LU_WINDOWS
#endif
#define SparseLU_ARCHITECTURE "Microsoft Windows"

#elif defined (__hppa) || defined (__hpux) || defined (MHPUX) || defined (ARCH_HPUX)
#define LU_HP
#define SparseLU_ARCHITECTURE "HP Unix"

#elif defined (__hp700) || defined (MHP700) || defined (ARCH_HP700)
#define LU_HP
#define SparseLU_ARCHITECTURE "HP 700 Unix"

#else


#define SparseLU_ARCHITECTURE "unknown"
#endif



#if defined (DLONG)
#define LONG_INTEGER
#endif

#if defined (LU_WINDOWS) && !defined (MATHWORKS)

#define SCALAR_IS_NAN(x)	(((x) != (x)) || (((x) < (x))))
#define SCALAR_IS_ZERO(x)	(((x) == 0.) && !SCALAR_IS_NAN(x))
#define SCALAR_IS_NONZERO(x)	(((x) != 0.) || SCALAR_IS_NAN(x))
#define SCALAR_IS_LTZERO(x)	(((x) < 0.) && !SCALAR_IS_NAN(x))

#else



#define SCALAR_IS_NAN(x)	((x) != (x))
#define SCALAR_IS_ZERO(x)	((x) == 0.)
#define SCALAR_IS_NONZERO(x)	((x) != 0.)
#define SCALAR_IS_LTZERO(x)	((x) < 0.)

#endif


#define SCALAR_ABS(x) ((SCALAR_IS_LTZERO (x)) ? -(x) : (x))


#define INT_OVERFLOW(x) ((!((x) * (1.0+1e-8) <= (double) Int_MAX)) \
			|| SCALAR_IS_NAN (x))


#define PRINT_SCALAR(a) \
{ \
    if (SCALAR_IS_NONZERO (a)) \
    { \
	PRINTF ((" (%g)", (a))) ; \
    } \
    else \
    { \
	PRINTF ((" (0)")) ; \
    } \
}

#define Entry double

#define SPLIT(s)    		    (1)
#define REAL_COMPONENT(c)	    (c)
#define IMAG_COMPONENT(c)	    (0.)
#define ASSIGN(c,s1,s2,p,split)	    { (c) = (s1)[p] ; }
#define CLEAR(c)		    { (c) = 0. ; }
#define CLEAR_AND_INCREMENT(p)	    { *p++ = 0. ; }
#define IS_NAN(a)		    SCALAR_IS_NAN (a)
#define IS_ZERO(a)		    SCALAR_IS_ZERO (a)
#define IS_NONZERO(a)		    SCALAR_IS_NONZERO (a)
#define SCALE_DIV(c,s)		    { (c) /= (s) ; }
#define SCALE(c,s)		    { (c) *= (s) ; }
#define ASSEMBLE(c,a)		    { (c) += (a) ; }
#define ASSEMBLE_AND_INCREMENT(c,p) { (c) += *p++ ; }
#define DECREMENT(c,a)		    { (c) -= (a) ; }
#define MULT(c,a,b)		    { (c) = (a) * (b) ; }
#define MULT_CONJ(c,a,b)	    { (c) = (a) * (b) ; }
#define MULT_SUB(c,a,b)		    { (c) -= (a) * (b) ; }
#define MULT_SUB_CONJ(c,a,b)	    { (c) -= (a) * (b) ; }
#define DIV(c,a,b)		    { (c) = (a) / (b) ; }
#define DIV_CONJ(c,a,b)		    { (c) = (a) / (b) ; }
#define APPROX_ABS(s,a)		    { (s) = SCALAR_ABS (a) ; }
#define ABS(s,a)		    { (s) = SCALAR_ABS (a) ; }
#define PRINT_ENTRY(a)		    PRINT_SCALAR (a)


#define MULTSUB_FLOPS	2.	
#define DIV_FLOPS	1.	
#define ABS_FLOPS	0.	
#define ASSEMBLE_FLOPS	1.	
#define DECREMENT_FLOPS	1.	
#define MULT_FLOPS	1.	
#define SCALE_FLOPS	1.	





#ifdef DLONG

#define LU_analyze		 lu_analyze
#define LU_apply_order		 lu_apply_order
#define LU_assemble		 lu_assemble
#define LU_assemble_fixq	 lu_assemble_fixq
#define LU_blas3_update	 lu_blas3_update
#define LU_build_tuples	 lu_build_tuples
#define LU_build_tuples_usage	 lu_build_tuples_usage
#define LU_colamd		 lu_colamd
#define LU_colamd_set_defaults	 lu_colamd_set_defaults
#define LU_create_element	 lu_create_element
#define LU_extend_front	 lu_extend_front
#define LU_free		 lu_free
#define LU_fsize		 lu_fsize
#define LU_garbage_collection	 lu_garbage_collection
#define LU_get_memory		 lu_get_memory
#define LU_grow_front		 lu_grow_front
#define LU_init_front		 lu_init_front
#define LU_is_permutation	 lu_is_permutation
#define LU_kernel		 lu_kernel
#define LU_kernel_init		 lu_kernel_init
#define LU_kernel_init_usage	 lu_kernel_init_usage
#define LU_kernel_wrapup	 lu_kernel_wrapup
#define LU_local_search	 lu_local_search
#define LU_lsolve		 lu_lsolve
#define LU_ltsolve		 lu_ltsolve
#define LU_lhsolve		 lu_lhsolve
#define LU_malloc		 lu_malloc
#define LU_mem_alloc_element	 lu_mem_alloc_element
#define LU_mem_alloc_head_block lu_mem_alloc_head_block
#define LU_mem_alloc_tail_block lu_mem_alloc_tail_block
#define LU_mem_free_tail_block	 lu_mem_free_tail_block
#define LU_mem_init_memoryspace lu_mem_init_memoryspace
#define LU_realloc		 lu_realloc
#define LU_row_search		 lu_row_search
#define LU_scale		 lu_scale
#define LU_scale_column	 lu_scale_column
#define LU_set_stats		 lu_set_stats
#define LU_singletons		 lu_singletons
#define LU_solve		 lu_solve
#define LU_start_front		 lu_start_front
#define LU_store_lu		 lu_store_lu
#define LU_store_lu_drop	 lu_store_lu_drop
#define LU_symbolic_usage	 lu_symbolic_usage
#define LU_transpose		 lu_transpose
#define LU_tuple_lengths	 lu_tuple_lengths
#define LU_usolve		 lu_usolve
#define LU_utsolve		 lu_utsolve
#define LU_uhsolve		 lu_uhsolve
#define LU_valid_numeric	 lu_valid_numeric
#define LU_valid_symbolic	 lu_valid_symbolic
#define LU_triplet_map_x	 lu_triplet_map_x
#define LU_triplet_map_nox	 lu_triplet_map_nox
#define LU_triplet_nomap_x	 lu_triplet_nomap_x
#define LU_triplet_nomap_nox	 lu_triplet_nomap_nox
#define LU_chol		    lu_chol

#define SparseLU_col_to_triplet	 sparselu_col_to_triplet
#define SparseLU_defaults	 sparselu_defaults
#define SparseLU_free_numeric	 sparselu_free_numeric
#define SparseLU_free_symbolic	 sparselu_free_symbolic
#define SparseLU_get_lunz	 sparselu_get_lunz
#define SparseLU_get_numeric	 sparselu_get_numeric
#define SparseLU_get_symbolic	 sparselu_get_symbolic
#define SparseLU_get_determinant	 sparselu_get_determinant
#define SparseLU_numeric		 sparselu_numeric
#define SparseLU_qsymbolic	 sparselu_qsymbolic
#define SparseLU_fsymbolic	 sparselu_fsymbolic
#define SparseLU_save_numeric	 sparselu_save_numeric
#define SparseLU_save_symbolic	 sparselu_save_symbolic
#define SparseLU_load_numeric	 sparselu_load_numeric
#define SparseLU_load_symbolic	 sparselu_load_symbolic
#define SparseLU_scale		 sparselu_scale
#define SparseLU_solve		 sparselu_solve
#define SparseLU_symbolic	 sparselu_symbolic
#define SparseLU_transpose	 sparselu_transpose
#define SparseLU_triplet_to_col	 sparselu_triplet_to_col
#define SparseLU_wsolve		 sparselu_wsolve

#endif





#define ONES_COMPLEMENT(r) (-(r)-1)






#include "amd_internal.h"





#include "SparseLU_config.h"





#include "SparseLU.h"





#define ESTIMATE (SparseLU_NUMERIC_SIZE_ESTIMATE - SparseLU_NUMERIC_SIZE)
#define ACTUAL 0





#define GET_CONTROL(i,default) \
    ((Control != (double *) NULL) ? \
	(SCALAR_IS_NAN (Control [i]) ? default : Control [i]) \
	: default)





#define MAX_MARK(n) Int_MAX - (2*(n)+1)





#define MBYTES(units) (((units) * sizeof (Unit)) / 1048576.0)










#define SparseLU_DENSE_DEGREE_THRESHOLD(alpha,n) \
    ((Int) MAX (16.0, (alpha) * 16.0 * sqrt ((double) (n))))





#define PRINTFk(k,params) { if (prl >= (k)) { PRINTF (params) ; } }
#define PRINTF1(params) PRINTFk (1, params)
#define PRINTF2(params) PRINTFk (2, params)
#define PRINTF3(params) PRINTFk (3, params)
#define PRINTF4(params) PRINTFk (4, params)
#define PRINTF5(params) PRINTFk (5, params)
#define PRINTF6(params) PRINTFk (6, params)






#define MAX_CANDIDATES 128


#define LU_REALLOC_REDUCTION (0.95)


#define LU_REALLOC_INCREASE (1.2)


#define LU_FRONTAL_GROWTH (1.2)


#define MAXNB 64


#define RECIPROCAL_TOLERANCE 1e-12






typedef double Align ;




#define BYTES(type,n) (sizeof (type) * (n))


#define CEILING(b,u) (((b) + (u) - 1) / (u))


#define UNITS(type,n) (CEILING (BYTES (type, n), sizeof (Unit)))


#define DUNITS(type,n) (ceil (BYTES (type, (double) n) / sizeof (Unit)))

union Unit_union
{	
    struct
    {
	Int
	    size,	
			
			
	    prevsize ;	
			
			
			
    } header ;		
    Align  xxxxxx ;	
} ;

typedef union Unit_union Unit ;


#define GET_BLOCK_SIZE(p) (((p)-1)->header.size)







#ifdef DLONG
#define NUMERIC_VALID  399789720
#define SYMBOLIC_VALID 399192713
#endif

typedef struct	
{
    double
	flops,		
	relpt,		
	relpt2,		
	droptol,
	alloc_init,	
	front_alloc_init, 
	rsmin,		
	rsmax,		
	min_udiag,	
	max_udiag,	
	rcond ;		

    Int
	scale ;

    Int valid ;		

    
    Unit
	*Memory ;	
    Int
	ihead,		
	itail,		
			
	ibig,		
	size ;		

    Int
	*Rperm,		
			
			
	*Cperm,		
			
			

	*Upos,		
	*Lpos,
	*Lip,
	*Lilen,
	*Uip,
	*Uilen,
	*Upattern ;	

    Int
	ulen,		
	npiv,		
	nnzpiv ;	

    Entry
	*D ;		

    Int do_recip ;
    double *Rs ;	
			
			

    Int
	n_row, n_col,	
	n1 ;		

    
    Int
	tail_usage,	
			
	init_usage,	
	max_usage,	
			
	ngarbage,	
	nrealloc,	
	ncostly,	
	isize,		
	nLentries,	
	nUentries,	
			
	lnz,		
	all_lnz,	
	unz,		
	all_unz,	
	maxfrsize ;	

    Int maxnrows, maxncols ;	

} NumericType ;







typedef struct	
{
    
    Int
	e,		
	f ;		

} Tuple ;

#define TUPLES(t) MAX (4, (t) + 1)


#define NON_PIVOTAL_COL(col) (Col_degree [col] >= 0)
#define NON_PIVOTAL_ROW(row) (Row_degree [row] >= 0)





typedef struct	
{
    Int

	cdeg,		
	rdeg,		
	nrowsleft,	
	ncolsleft,	
	nrows,		
	ncols,		
	next ;		

    

} Element ;



#define GET_ELEMENT_SIZE(nr,nc) \
(UNITS (Element, 1) + UNITS (Int, (nc) + (nr)) + UNITS (Entry, (nc) * (nr)))

#define DGET_ELEMENT_SIZE(nr,nc) \
(DUNITS (Element, 1) + DUNITS (Int, (nc) + (nr)) + DUNITS (Entry, (nc) * (nr)))

#define GET_ELEMENT_COLS(ep,p,Cols) { \
    ASSERT (p != (Unit *) NULL) ; \
    ASSERT (p >= Numeric->Memory + Numeric->itail) ; \
    ASSERT (p <= Numeric->Memory + Numeric->size) ; \
    ep = (Element *) p ; \
    p += UNITS (Element, 1) ; \
    Cols = (Int *) p ; \
}

#define GET_ELEMENT_PATTERN(ep,p,Cols,Rows,ncm) { \
    GET_ELEMENT_COLS (ep, p, Cols) ; \
    ncm = ep->ncols ; \
    Rows = Cols + ncm ; \
}

#define GET_ELEMENT(ep,p,Cols,Rows,ncm,nrm,C) { \
    GET_ELEMENT_PATTERN (ep, p, Cols, Rows, ncm) ; \
    nrm = ep->nrows ; \
    p += UNITS (Int, ncm + nrm) ; \
    C = (Entry *) p ; \
}







typedef struct	
{
    Int	*E ;		
    Entry *Wx, *Wy ;	
    Int			
	*Wp,		
	*Wrp,		
	*Wm,		
	*Wio,		
	*Woi,		
	*Woo,		
	*Wrow,		
	*NewRows,	
	*NewCols ;	
    Int
	*Lpattern,	
	*Upattern,	
	ulen, llen ;	
    Int
	*Diagonal_map,	
	*Diagonal_imap ;
    Int
	n_row, n_col,	
	nz,		
	n1,		
	elen,		
	npiv,		
	ndiscard,	
	Wrpflag,
	nel,		
	noff_diagonal,
	prior_element,
	rdeg0, cdeg0,
	rrdeg, ccdeg,
	Candidates [MAX_CANDIDATES],	 
	nCandidates,	
	ksuper,
	firstsuper,
	jsuper,
	ncand,		
	nextcand,	
	lo,
	hi,
	pivrow,		
	pivcol,		
	do_extend,	
	do_update,	
	nforced,	
	any_skip,
	do_scan2row,
	do_scan2col,
	do_grow,
	pivot_case,
	frontid,	
	nfr ;		
    Int
	*Front_new1strow ;
    Int Pivrow [MAXNB],
	Pivcol [MAXNB] ;

    Entry
	*Flublock,	
	*Flblock,	
	*Fublock,	
	*Fcblock ;	
    Int
	*Frows,		
	*Fcols,		
	*Frpos,		
	*Fcpos,		
	fnrows,		
	fncols,		
	fnr_curr,	
	fnc_curr,	
	fcurr_size,	
	fnrows_max,	
	fncols_max,	
	nb,
	fnpiv,		
	fnzeros,	
	fscan_row,	
	fscan_col,	
	fnrows_new,	
	fncols_new,	
	pivrow_in_front,	
	pivcol_in_front ;	
} WorkType ;


typedef struct	
{

    double
	num_mem_usage_est,	
	num_mem_size_est,	
	peak_sym_usage,		
	sym,			
	dnum_mem_init_usage,	
	amd_lunz,	
	lunz_bound ;	

    Int valid,		
	max_nchains,
	nchains,
	*Chain_start,
	*Chain_maxrows,
	*Chain_maxcols,
	maxnrows,		
	maxncols,		
	*Front_npivcol,		
	*Front_1strow,		
	*Front_leftmostdesc,	
	*Front_parent,		
	*Cperm_init,		
	*Rperm_init,		
	*Cdeg, *Rdeg,
	*Esize,
	dense_row_threshold,
	n1,			
	nempty,			
	*Diagonal_map,		
	esize,			
	nfr,
	n_row, n_col,		
	nz,			
	nb,			
	num_mem_init_usage,	
	nempty_row, nempty_col,

	strategy,
	ordering,
	fixQ,
	prefer_diagonal,
	nzaat,
	nzdiag,
	amd_dmax ;

} SymbolicType ;



#endif
