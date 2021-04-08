/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    SparseLU_symbolic.c
 * BRIEF:   稀疏LU符号分解
 *****************************************************************************/



#include "SparseLU_internal.h"
#include "SparseLU_function.h"


#define SYM_WORK_USAGE(n_col,n_row,Clen) \
    (DUNITS (Int, Clen) + \
     DUNITS (Int, nz) + \
     4 * DUNITS (Int, n_row) + \
     4 * DUNITS (Int, n_col) + \
     2 * DUNITS (Int, n_col + 1) + \
     DUNITS (double, n_row))


#define LU_ANALYZE_CLEN(nz,n_row,n_col,nn) \
    ((n_col) + MAX ((nz),(n_col)) + 3*(nn)+1 + (n_col))


#define ELEMENT_SIZE(r,c) \
    (DGET_ELEMENT_SIZE (r, c) + 1 + (r + c) * UNITS (Tuple, 1))

typedef struct	
{
    Int *Front_npivcol ;    
    Int *Front_nrows ;	    
    Int *Front_ncols ;	    
    Int *Front_parent ;	    
    Int *Front_cols ;	    
    Int *InFront ;	    
    Int *Ci ;		    
    Int *Cperm1 ;	    
    Int *Rperm1 ;	    
    Int *InvRperm1 ;	    
    Int *Si ;		    
    Int *Sp ;		    
    double *Rs ;	    

} SWType ;

PRIVATE Int symbolic_analysis
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    const Int Quser [ ],
    int (*user_ordering) ( Int, Int, Int, Int *, Int *, Int *, void *, double * ),
    void *user_params, 
    void **SymbolicHandle,
    const double Control [SparseLU_CONTROL],
    double User_Info [SparseLU_INFO]
);

PRIVATE int inverse_permutation
(
    Int *P,  
    Int *Pinv,  
    Int n      
);

PRIVATE int do_amd_1
(
    Int n,	
    Int Ap [ ],	     
    Int Ai [ ],	      
    Int P [ ],		
    Int Pinv [ ],	
    Int Len [ ],	
    Int slen,	
    Int S [ ],		
    Int ordering_option,
    Int print_level,
    int (*user_ordering) ( Int, Int, Int, Int *, Int *, Int *, void *, double * ),
    void *user_params, 
    Int *ordering_used,
    double amd_Control [ ],	
    double amd_Info [ ] 
);

PRIVATE int do_amd
(
    Int n,
    Int Ap [ ],		  
    Int Ai [ ],		     
    Int Q [ ],		
    Int Qinv [ ],		
    Int Sdeg [ ],	
    Int Clen,			
    Int Ci [ ],			
    double amd_Control [ ],	
    double amd_Info [ ],
    SymbolicType *Symbolic,	
    double Info [ ],		
    Int ordering_option,
    Int print_level,
    int (*user_ordering) ( Int, Int, Int, Int *, Int *, Int *, void *, double * ),
    void *user_params, 
    Int *ordering_used
);

PRIVATE Int prune_singletons
(
    Int n1,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    Int Cperm1 [ ],
    Int InvRperm1 [ ],
    Int Si [ ],
    Int Sp [ ]
);

PRIVATE void combine_ordering
(
    Int n1,
    Int nempty_col,
    Int n_col,
    Int Cperm_init [ ],	
    Int Cperm1 [ ],	   
    Int Qinv [ ]	
);

PRIVATE void free_work
(
    SWType *SW
) ;

PRIVATE void error
(
    SymbolicType **Symbolic,
    SWType *SW
) ;






GLOBAL Int SparseLU_symbolic
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    void **SymbolicHandle,
    const double Control [SparseLU_CONTROL],
    double Info [SparseLU_INFO]
)
{
    return (symbolic_analysis (n_row, n_col, Ap, Ai, Ax, (Int *) NULL, (void *) NULL,         
        (void *) NULL, SymbolicHandle, Control, Info)) ;
}







PRIVATE Int symbolic_analysis
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],

    
    const Int Quser [ ],

    
    int (*user_ordering)    
    (
        
        Int,            
        Int,            
        Int,            
        Int *,          
        Int *,          
        
        Int *,          
        
        void *,         
        double *        
    ),
    void *user_params,  

    void **SymbolicHandle,
    const double Control [SparseLU_CONTROL],
    double User_Info [SparseLU_INFO]
)
{

    
    
    

    double knobs [COLAMD_KNOBS], flops, f, r, c, force_fixQ,
	Info2 [SparseLU_INFO], drow, dcol, dtail_usage, dlf, duf, dmax_usage,
	dhead_usage, dlnz, dunz, dmaxfrsize, dClen, dClen_analyze, sym,
	amd_Info [AMD_INFO], dClen_amd, dr, dc, cr, cc, cp,
	amd_Control [AMD_CONTROL], stats [2] ;
    double *Info ;
    Int i, nz, j, newj, status, f1, f2, maxnrows, maxncols, nfr, col,
	nchains, maxrows, maxcols, p, nb, nn, *Chain_start, *Chain_maxrows,
	*Chain_maxcols, *Front_npivcol, *Ci, Clen, colamd_stats [COLAMD_STATS],
	fpiv, n_inner, child, parent, *Link, row, *Front_parent,
	analyze_compactions, k, chain, is_sym, *Si, *Sp, n2, do_LU_analyze,
	fpivcol, fallrows, fallcols, *InFront, *F1, snz, *Front_1strow, f1rows,
	kk, *Cperm_init, *Rperm_init, newrow, *InvRperm1, *Front_leftmostdesc,
	Clen_analyze, strategy, Clen_amd, fixQ, prefer_diagonal, nzdiag, nzaat,
	*Wq, *Sdeg, *Fr_npivcol, nempty, *Fr_nrows, *Fr_ncols, *Fr_parent,
	*Fr_cols, nempty_row, nempty_col, user_auto_strategy, fail, max_rdeg,
	head_usage, tail_usage, lnz, unz, esize, *Esize, rdeg, *Cdeg, *Rdeg,
	*Cperm1, *Rperm1, n1, oldcol, newcol, n1c, n1r, oldrow,
	dense_row_threshold, tlen, aggressive, *Rp, *Ri ;
    Int do_singletons, ordering_option, print_level ;
    int ok ;

    SymbolicType *Symbolic ;
    SWType SWspace, *SW ;

    
    
    

    drow = GET_CONTROL (SparseLU_DENSE_ROW, SparseLU_DEFAULT_DENSE_ROW) ;
    dcol = GET_CONTROL (SparseLU_DENSE_COL, SparseLU_DEFAULT_DENSE_COL) ;
    nb = GET_CONTROL (SparseLU_BLOCK_SIZE, SparseLU_DEFAULT_BLOCK_SIZE) ;
    strategy = GET_CONTROL (SparseLU_STRATEGY, SparseLU_DEFAULT_STRATEGY) ;
    force_fixQ = GET_CONTROL (SparseLU_FIXQ, SparseLU_DEFAULT_FIXQ) ;
    do_singletons = GET_CONTROL (SparseLU_SINGLETONS,SparseLU_DEFAULT_SINGLETONS);
    AMD_defaults (amd_Control) ;
    amd_Control [AMD_DENSE] = GET_CONTROL (SparseLU_AMD_DENSE, SparseLU_DEFAULT_AMD_DENSE) ;
    aggressive = (GET_CONTROL (SparseLU_AGGRESSIVE, SparseLU_DEFAULT_AGGRESSIVE) != 0) ;
    amd_Control [AMD_AGGRESSIVE] = aggressive ;
    print_level = GET_CONTROL (SparseLU_PRL, SparseLU_DEFAULT_PRL) ;

    
    ordering_option = GET_CONTROL (SparseLU_ORDERING, SparseLU_DEFAULT_ORDERING) ;
    if (ordering_option < 0 || ordering_option > SparseLU_ORDERING_USER)
    {
        ordering_option = SparseLU_DEFAULT_ORDERING ;
    }
    if (Quser == (Int *) NULL)
    {
        
        
        if (ordering_option == SparseLU_ORDERING_GIVEN ||
           (ordering_option == SparseLU_ORDERING_USER && !user_ordering))
        {
            ordering_option = SparseLU_ORDERING_NONE ;
        }
    }
    else
    {
        
        ordering_option = SparseLU_ORDERING_GIVEN ;
    }

    nb = MAX (2, nb) ;
    nb = MIN (nb, MAXNB) ;
    if (nb % 2 == 1) nb++ ;	

    if (User_Info != (double *) NULL)
    {
        
        Info = User_Info ;
    }
    else
    {
        
        Info = Info2 ;
    }
    
    for (i = 0 ; i < SparseLU_INFO ; i++)
    {
	    Info [i] = EMPTY ;
    }

    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;

    Info [SparseLU_STATUS] = SparseLU_OK ;
    Info [SparseLU_NROW] = n_row ;
    Info [SparseLU_NCOL] = n_col ;
    Info [SparseLU_SIZE_OF_UNIT] = (double) (sizeof (Unit)) ;
    Info [SparseLU_SIZE_OF_INT] = (double) (sizeof (int)) ;
    Info [SparseLU_SIZE_OF_LONG] = (double) (sizeof (Sparse_long)) ;
    Info [SparseLU_SIZE_OF_POINTER] = (double) (sizeof (void *)) ;
    Info [SparseLU_SIZE_OF_ENTRY] = (double) (sizeof (Entry)) ;
    Info [SparseLU_SYMBOLIC_DEFRAG] = 0 ;
    Info [SparseLU_ORDERING_USED] = EMPTY ;

    if (SymbolicHandle != NULL)
    {
        *SymbolicHandle = (void *) NULL ;
    }

    if (!Ai || !Ap || !SymbolicHandle)
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_argument_missing ;
        return (SparseLU_ERROR_argument_missing) ;
    }

    if (n_row <= 0 || n_col <= 0)	
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_n_nonpositive ;
        return (SparseLU_ERROR_n_nonpositive) ;
    }

    nz = Ap [n_col] ;
    Info [SparseLU_NZ] = nz ;
    if (nz < 0)
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_invalid_matrix ;
        return (SparseLU_ERROR_invalid_matrix) ;
    }

    
    
    

    if (n_row != n_col)
    {
        
        strategy = SparseLU_STRATEGY_UNSYMMETRIC ;
    }

    if (strategy < SparseLU_STRATEGY_AUTO
     || strategy > SparseLU_STRATEGY_SYMMETRIC)
    {
        
        strategy = SparseLU_STRATEGY_AUTO ;
    }

    if (Quser != (Int *) NULL)
    {
        
        if (strategy != SparseLU_STRATEGY_SYMMETRIC)
        {
            strategy = SparseLU_STRATEGY_UNSYMMETRIC ;
        }
    }

    user_auto_strategy = (strategy == SparseLU_STRATEGY_AUTO) ;

    
    
    

    
    
    

    
    dClen = LU_COLAMD_RECOMMENDED ((double) nz, (double) n_row,
	(double) n_col) ;

    
    dClen_analyze = LU_ANALYZE_CLEN ((double) nz, (double) n_row,
	(double) n_col, (double) nn) ;
    dClen = MAX (dClen, dClen_analyze) ;

    
    dClen_amd = 2.4 * (double) nz + 8 * (double) n_inner + 1 ;

    dClen = MAX (dClen, dClen_amd) ;

    
    Info [SparseLU_SYMBOLIC_PEAK_MEMORY] =
	SYM_WORK_USAGE (n_col, n_row, dClen) +
	LU_symbolic_usage (n_row, n_col, n_col, n_col, n_col, TRUE) ;

    if (INT_OVERFLOW (dClen * sizeof (Int)))
    {
	
	
	
	
        
	Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
	return (SparseLU_ERROR_out_of_memory) ;
    }

    
    Clen = LU_COLAMD_RECOMMENDED (nz, n_row, n_col) ;
    Clen_analyze = LU_ANALYZE_CLEN (nz, n_row, n_col, nn) ;
    Clen = MAX (Clen, Clen_analyze) ;
    Clen_amd = 2.4 * nz + 8 * n_inner + 1 ;
    Clen = MAX (Clen, Clen_amd) ;

    
    
    

    

    Symbolic = (SymbolicType *) LU_malloc (1, sizeof (SymbolicType)) ;

    if (!Symbolic)
    {
	
	
	Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
	error (&Symbolic, (SWType *) NULL) ;
	return (SparseLU_ERROR_out_of_memory) ;
    }

    
    Symbolic->valid = 0 ;
    Symbolic->Chain_start = (Int *) NULL ;
    Symbolic->Chain_maxrows = (Int *) NULL ;
    Symbolic->Chain_maxcols = (Int *) NULL ;
    Symbolic->Front_npivcol = (Int *) NULL ;
    Symbolic->Front_parent = (Int *) NULL ;
    Symbolic->Front_1strow = (Int *) NULL ;
    Symbolic->Front_leftmostdesc = (Int *) NULL ;
    Symbolic->Esize = (Int *) NULL ;
    Symbolic->esize = 0 ;
    Symbolic->ordering = EMPTY ;    
    Symbolic->amd_lunz = EMPTY ;
    Symbolic->max_nchains = EMPTY ;

    Symbolic->Cperm_init   = (Int *) LU_malloc (n_col+1, sizeof (Int)) ;
    Symbolic->Rperm_init   = (Int *) LU_malloc (n_row+1, sizeof (Int)) ;
    Symbolic->Cdeg	   = (Int *) LU_malloc (n_col+1, sizeof (Int)) ;
    Symbolic->Rdeg	   = (Int *) LU_malloc (n_row+1, sizeof (Int)) ;
    Symbolic->Diagonal_map = (Int *) NULL ;

    Cperm_init = Symbolic->Cperm_init ;
    Rperm_init = Symbolic->Rperm_init ;
    Cdeg = Symbolic->Cdeg ;
    Rdeg = Symbolic->Rdeg ;

    if (!Cperm_init || !Rperm_init || !Cdeg || !Rdeg)
    {
	Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
	error (&Symbolic, (SWType *) NULL) ;
	return (SparseLU_ERROR_out_of_memory) ;
    }

    Symbolic->n_row = n_row ;
    Symbolic->n_col = n_col ;
    Symbolic->nz = nz ;
    Symbolic->nb = nb ;
    Cdeg [n_col] = EMPTY ;      
    Rdeg [n_row] = EMPTY ;

    
    
    

    if (Quser != (Int *) NULL)
    {
        
        if (!LU_is_permutation (Quser, Cperm_init, n_col, n_col))
        {
            Info [SparseLU_STATUS] = SparseLU_ERROR_invalid_permutation ;
            error (&Symbolic, (SWType *) NULL) ;
            return (SparseLU_ERROR_invalid_permutation) ;
        }
    }

    
    
    

    

    SW = &SWspace ;	

    
    

    
    SW->Si	      = (Int *) LU_malloc (nz, sizeof (Int)) ;
    SW->Sp	      = (Int *) LU_malloc (n_col + 1, sizeof (Int)) ;
    SW->InvRperm1     = (Int *) LU_malloc (n_row, sizeof (Int)) ;
    SW->Cperm1	      = (Int *) LU_malloc (n_col, sizeof (Int)) ;

    
    SW->Ci	      = (Int *) LU_malloc (Clen, sizeof (Int)) ;
    SW->Front_npivcol = (Int *) LU_malloc (n_col + 1, sizeof (Int)) ;
    SW->Front_nrows   = (Int *) LU_malloc (n_col, sizeof (Int)) ;
    SW->Front_ncols   = (Int *) LU_malloc (n_col, sizeof (Int)) ;
    SW->Front_parent  = (Int *) LU_malloc (n_col, sizeof (Int)) ;
    SW->Front_cols    = (Int *) LU_malloc (n_col, sizeof (Int)) ;
    SW->Rperm1	      = (Int *) LU_malloc (n_row, sizeof (Int)) ;
    SW->InFront	      = (Int *) LU_malloc (n_row, sizeof (Int)) ;

    
    SW->Rs	      = (double *) NULL ;	
    Ci	       = SW->Ci ;
    Fr_npivcol = SW->Front_npivcol ;
    Fr_nrows   = SW->Front_nrows ;
    Fr_ncols   = SW->Front_ncols ;
    Fr_parent  = SW->Front_parent ;
    Fr_cols    = SW->Front_cols ;
    Cperm1     = SW->Cperm1 ;
    Rperm1     = SW->Rperm1 ;
    Si	       = SW->Si ;
    Sp	       = SW->Sp ;
    InvRperm1  = SW->InvRperm1 ;
    InFront    = SW->InFront ;

    if (!Ci || !Fr_npivcol || !Fr_nrows || !Fr_ncols || !Fr_parent || !Fr_cols
	|| !Cperm1 || !Rperm1 || !Si || !Sp || !InvRperm1 || !InFront)
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
        error (&Symbolic, SW) ;
        return (SparseLU_ERROR_out_of_memory) ;
    }

    
    
    

    

    status = LU_singletons (n_row, n_col, Ap, Ai, Quser, strategy,
        do_singletons, 
	Cdeg, Cperm1, Rdeg,
	Rperm1, InvRperm1, &n1, &n1c, &n1r, &nempty_col, &nempty_row, &is_sym,
	&max_rdeg,  Rperm_init, Ci, Ci + nz, Ci + nz + n_row) ;

    

    

    if (status != SparseLU_OK)
    {
        Info [SparseLU_STATUS] = status ;
        error (&Symbolic, SW) ;
        return (status) ;
    }
    Info [SparseLU_NEMPTY_COL] = nempty_col ;
    Info [SparseLU_NEMPTY_ROW] = nempty_row ;
    Info [SparseLU_NDENSE_COL] = 0 ;	
    Info [SparseLU_NDENSE_ROW] = 0 ;
    Info [SparseLU_COL_SINGLETONS] = n1c ;
    Info [SparseLU_ROW_SINGLETONS] = n1r ;
    Info [SparseLU_S_SYMMETRIC] = is_sym ;

    nempty = MIN (nempty_col, nempty_row) ;
    Symbolic->nempty_row = nempty_row ;
    Symbolic->nempty_col = nempty_col ;

    

    Symbolic->n1 = n1 ;
    Symbolic->nempty = nempty ;
    n2 = nn - n1 - nempty ;

    dense_row_threshold =
	SparseLU_DENSE_DEGREE_THRESHOLD (drow, n_col - n1 - nempty_col) ;
    Symbolic->dense_row_threshold = dense_row_threshold ;

    if (!is_sym)
    {
        
        strategy = SparseLU_STRATEGY_UNSYMMETRIC ;
    }

    
    
    

    

    Wq = Rperm_init ;	    
    Sdeg = Cperm_init ;	    
    sym = EMPTY ;
    nzaat = EMPTY ;
    nzdiag = EMPTY ;
    for (i = 0 ; i < AMD_INFO ; i++)
    {
	    amd_Info [i] = EMPTY ;
    }

    if (strategy != SparseLU_STRATEGY_UNSYMMETRIC)
    {
	
	ASSERT (n_row == n_col) ;
	ASSERT (nempty_row == nempty_col) ;

	

	nzdiag = prune_singletons (n1, nn, Ap, Ai, Ax,
	    Cperm1, InvRperm1, Si, Sp
	    ) ;

	
	if (Quser != (Int *) NULL)
	{
	    
	    Rp = Ci ;
	    Ri = Ci + (n_row) + 1 ;
	    (void) LU_transpose (n2, n2, Sp, Si, (double *) NULL,
		(Int *) NULL, (Int *) NULL, 0,
		Rp, Ri, (double *) NULL, Wq, FALSE
		) ;
	}
	else
	{
	    
	    Rp = Sp ;
	    Ri = Si ;
	}
	ASSERT (AMD_valid (n2, n2, Rp, Ri) == AMD_OK) ;

	nzaat = AMD_aat (n2, Rp, Ri, Sdeg, Wq, amd_Info) ;
	sym = amd_Info [AMD_SYMMETRY] ;
	Info [SparseLU_N2] = n2 ;
	

	

    }

    
    Symbolic->sym = sym ;
    Symbolic->nzaat = nzaat ;
    Symbolic->nzdiag = nzdiag ;
    Symbolic->amd_dmax = EMPTY ;

    Info [SparseLU_PATTERN_SYMMETRY] = sym ;
    Info [SparseLU_NZ_A_PLUS_AT] = nzaat ;
    Info [SparseLU_NZDIAG] = nzdiag ;

    
    
    

    if (strategy == SparseLU_STRATEGY_AUTO)
    {
        if (sym >= 0.5 && nzdiag >= 0.9 * n2)
        {
            
	        strategy = SparseLU_STRATEGY_SYMMETRIC ;
        }
        else
        {
            
	        strategy = SparseLU_STRATEGY_UNSYMMETRIC ;
        }
    }

    
    
    

    if (strategy == SparseLU_STRATEGY_SYMMETRIC)
    {
        
        fixQ = TRUE ;
        prefer_diagonal = TRUE ;
    }
    else
    {
        
        fixQ = FALSE ;
        prefer_diagonal = FALSE ;
    }

    if (force_fixQ > 0)
    {
	    fixQ = TRUE ;
    }
    else if (force_fixQ < 0)
    {
	    fixQ = FALSE ;
    }

    
    Symbolic->strategy = strategy ;
    Symbolic->fixQ = fixQ ;
    Symbolic->prefer_diagonal = prefer_diagonal ;

    Info [SparseLU_STRATEGY_USED] = strategy ;
    Info [SparseLU_QFIXED] = fixQ ;
    Info [SparseLU_DIAG_PREFERRED] = prefer_diagonal ;

    
    
    

    if (strategy == SparseLU_STRATEGY_SYMMETRIC && Quser == (Int *) NULL)
    {
        
        Int ordering_used ;
        Int *Qinv = Fr_npivcol ;
        ok = do_amd (n2, Sp, Si, Wq, Qinv, Sdeg, Clen, Ci,
                amd_Control, amd_Info, Symbolic, Info,
                ordering_option, print_level, user_ordering, user_params,
                &ordering_used) ;
        if (!ok)
        {
            status = SparseLU_ERROR_ordering_failed ;
            Info [SparseLU_STATUS] = status ;
            error (&Symbolic, SW) ;
            return (status) ;
        }
        
        Symbolic->ordering = ordering_used ;
        combine_ordering (n1, nempty, nn, Cperm_init, Cperm1, Qinv) ;
    }
    
    

    

    
    
    

    if (Quser != (Int *) NULL)
    {
        for (k = 0 ; k < n_col ; k++)
        {
            Cperm_init [k] = Cperm1 [k] ;
        }
        Symbolic->ordering = SparseLU_ORDERING_GIVEN ;
    }

    
    
    

    if (strategy == SparseLU_STRATEGY_UNSYMMETRIC && Quser == (Int *) NULL)
    {
        Int nrow2, ncol2 ;

	
	
	

	

	(void) prune_singletons (n1, n_col, Ap, Ai,
	    (double *) NULL,
	    Cperm1, InvRperm1, Ci, Cperm_init
	    ) ;

        
        nrow2 = n_row - n1 - nempty_row ;
        ncol2 = n_col - n1 - nempty_col ;

        if ((ordering_option == SparseLU_ORDERING_USER
            || ordering_option == SparseLU_ORDERING_NONE
            || ordering_option == SparseLU_ORDERING_BEST)
            && nrow2 > 0 && ncol2 > 0)
        {

            
            
            

            double user_info [3] ;    
            Int *Qinv = Fr_npivcol ;  
            Int *QQ = Fr_nrows ;      

            
            do_LU_analyze = TRUE ;

            if (ordering_option == SparseLU_ORDERING_USER)
            {
                ok = (*user_ordering) (
                    
                    nrow2,
                    ncol2,
                    FALSE,
                    Cperm_init, 
                    Ci, 
                    
                    QQ, 
                    
                    user_params,
                    user_info) ;
                Symbolic->ordering = SparseLU_ORDERING_USER ;
            }
            else
            {
                Int params [3] ;
                params [0] = ordering_option ;
                params [1] = print_level ;
                ok = LU_chol (
                    
                    nrow2,
                    ncol2,
                    FALSE,
                    Cperm_init, 
                    Ci, 
                    
                    QQ, 
                    
                    &params,
                    user_info) ;
                Symbolic->ordering = params [2] ;
            }

            
            if (!ok || !inverse_permutation (QQ, Qinv, ncol2))
            {
                
                status = SparseLU_ERROR_ordering_failed ;
                Info [SparseLU_STATUS] = status ;
                error (&Symbolic, SW) ;
                return (status) ;
            }

            
            
            combine_ordering (n1, nempty_col, n_col, Cperm_init, Cperm1, Qinv) ;

        }
        else
        {

            
            
            

            LU_colamd_set_defaults (knobs) ;
            knobs [COLAMD_DENSE_ROW] = drow ;
            knobs [COLAMD_DENSE_COL] = dcol ;
            knobs [COLAMD_AGGRESSIVE] = aggressive ;

            
            
            

            

            (void) LU_colamd (
                    n_row - n1 - nempty_row,
                    n_col - n1 - nempty_col,
                    Clen, Ci, Cperm_init, knobs, colamd_stats,
                    Fr_npivcol, Fr_nrows, Fr_ncols, Fr_parent, Fr_cols, &nfr,
                    InFront) ;
            Symbolic->ordering = SparseLU_ORDERING_AMD ;

            
            Info [SparseLU_NDENSE_ROW]  = colamd_stats [COLAMD_DENSE_ROW] ;
            Info [SparseLU_NDENSE_COL]  = colamd_stats [COLAMD_DENSE_COL] ;
            Info [SparseLU_SYMBOLIC_DEFRAG] = colamd_stats [COLAMD_DEFRAG_COUNT];

            
            do_LU_analyze =
                colamd_stats [COLAMD_DENSE_ROW] > 0 ||
                colamd_stats [COLAMD_DENSE_COL] > 0 ;

            
            
            combine_ordering (n1, nempty_col, n_col, Cperm_init, Cperm1, Ci) ;
        }

	

    }
    else
    {

	
	
	

	
	do_LU_analyze = TRUE ;

    }

    
    Info [SparseLU_ORDERING_USED] = Symbolic->ordering ;

    Cperm_init [n_col] = EMPTY ;	

    
    
    

    
    
    

    if (do_LU_analyze)
    {

	Int *W, *Bp, *Bi, *Cperm2, *P, Clen2, bsize, Clen0 ;

	
	
	

	
	(void) prune_singletons (n1, n_col, Ap, Ai,
	    (double *) NULL,
	    Cperm_init, InvRperm1, Si, Sp
	    ) ;

	
	

	Clen0 = Clen - (nn+1 + 2*nn + n_col) ;
	Bp = Ci + Clen0 ;
	Link = Bp + (nn+1) ;
	W = Link + nn ;
	Cperm2 = W + nn ;

	
	
	

	
	for (row = 0 ; row < n_row - n1 ; row++)
	{
	    W [row] = FALSE ;
	}
	P = Link ;

	k = 0 ;

	for (col = 0 ; col < n_col - n1 ; col++)
	{
	    
	    for (p = Sp [col] ; p < Sp [col+1] ; p++)
	    {
		row = Si [p] ;
		if (!W [row])
		{
		    
		    W [row] = TRUE ;
		    P [k++] = row ;
		}
	    }
	}

	
	
	nempty_row = n_row - n1 - k ;
	if (k < n_row - n1)
	{
	    
	    
	    for (row = 0 ; row < n_row - n1 ; row++)
	    {
		if (!W [row])
		{
		    
		    P [k++] = row ;
		}
	    }
	}

	

	
	
	

	

	Clen2 = Clen0 ;
	snz = Sp [n_col - n1] ;
	bsize = MAX (snz, 1) ;
	Clen2 -= bsize ;
	Bi = Ci + Clen2 ;

	(void) LU_transpose (n_row - n1, n_col - n1 - nempty_col,
	    Sp, Si, (double *) NULL,
	    P, (Int *) NULL, 0, Bp, Bi, (double *) NULL, W, FALSE
	    ) ;

	

	
	

	
	for (i = 0 ; i <= n_row - n1 ; i++)
	{
	    Bp [i] += Clen2 ;
	}

	

	
	
	

	
	ok = LU_analyze (
		n_row - n1 - nempty_row,
		n_col - n1 - nempty_col,
		Ci, Bp, Cperm2, fixQ, W, Link,
		Fr_ncols, Fr_nrows, Fr_npivcol,
		Fr_parent, &nfr, &analyze_compactions) ;
	if (!ok)
	{
	    
	    Info [SparseLU_STATUS] = SparseLU_ERROR_internal_error ;
	    error (&Symbolic, SW) ;
	    return (SparseLU_ERROR_internal_error) ;
	}
	Info [SparseLU_SYMBOLIC_DEFRAG] += analyze_compactions ;

	
	
	

	if (!fixQ)
	{
	    

	    
	    for (k = 0 ; k < n_col - n1 - nempty_col ; k++)
	    {
		W [k] = Cperm_init [n1 + Cperm2 [k]] ;
	    }

	    for (k = 0 ; k < n_col - n1 - nempty_col ; k++)
	    {
		Cperm_init [n1 + k] = W [k] ;
	    }
	}
    }

    
    
    

    

    SW->Si = (Int *) LU_free ((void *) SW->Si) ;
    SW->Sp = (Int *) LU_free ((void *) SW->Sp) ;
    SW->Cperm1 = (Int *) LU_free ((void *) SW->Cperm1) ;
    ASSERT (SW->Rs == (double *) NULL) ;

    
    
    

    nchains = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	if (Fr_parent [i] != i+1)
	{
	    nchains++ ;
	}
    }

    Symbolic->nchains = nchains ;
    Symbolic->nfr = nfr ;
    Symbolic->esize = (max_rdeg > dense_row_threshold) ? (n_col - n1 - nempty_col) : 0 ;

    
    Info [SparseLU_SYMBOLIC_SIZE] = LU_symbolic_usage (n_row, n_col, nchains,
	    nfr, Symbolic->esize, prefer_diagonal) ;

    
    Info [SparseLU_SYMBOLIC_PEAK_MEMORY] =
	SYM_WORK_USAGE (n_col, n_row, Clen) + Info [SparseLU_SYMBOLIC_SIZE] ;
    Symbolic->peak_sym_usage = Info [SparseLU_SYMBOLIC_PEAK_MEMORY] ;

    
    
    

    

    
    Symbolic->Front_npivcol = (Int *) LU_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_parent = (Int *) LU_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_1strow = (Int *) LU_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_leftmostdesc = (Int *) LU_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Chain_start = (Int *) LU_malloc (nchains+1, sizeof (Int)) ;
    Symbolic->Chain_maxrows = (Int *) LU_malloc (nchains+1, sizeof (Int)) ;
    Symbolic->Chain_maxcols = (Int *) LU_malloc (nchains+1, sizeof (Int)) ;

    fail = (!Symbolic->Front_npivcol || !Symbolic->Front_parent ||
	!Symbolic->Front_1strow || !Symbolic->Front_leftmostdesc ||
	!Symbolic->Chain_start || !Symbolic->Chain_maxrows ||
	!Symbolic->Chain_maxcols) ;

    if (Symbolic->esize > 0)
    {
        Symbolic->Esize = (Int *) LU_malloc (Symbolic->esize, sizeof (Int)) ;
        fail = fail || !Symbolic->Esize ;
    }

    if (fail)
    {
        Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
        error (&Symbolic, SW) ;
        return (SparseLU_ERROR_out_of_memory) ;
    }

    Front_npivcol = Symbolic->Front_npivcol ;
    Front_parent = Symbolic->Front_parent ;
    Front_1strow = Symbolic->Front_1strow ;
    Front_leftmostdesc = Symbolic->Front_leftmostdesc ;

    Chain_start = Symbolic->Chain_start ;
    Chain_maxrows = Symbolic->Chain_maxrows ;
    Chain_maxcols = Symbolic->Chain_maxcols ;

    Esize = Symbolic->Esize ;

    
    
    

    
    if (do_LU_analyze)
    {
        
        for (row = 0 ; row < n_row ; row++)
        {
            InFront [row] = nfr ;
        }
        
        for (k = 0 ; k < n1 ; k++)
        {
            row = Rperm1 [k] ;
            InFront [row] = EMPTY ;
        }
        newj = n1 ;
        for (i = 0 ; i < nfr ; i++)
        {
            fpivcol = Fr_npivcol [i] ;
            f1rows = 0 ;
            
            for (kk = 0 ; kk < fpivcol ; kk++, newj++)
            {
            j = Cperm_init [newj] ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                row = Ai [p] ;
                if (InFront [row] == nfr)
                {
                
                InFront [row] = i ;
                f1rows++ ;
                }
            }
            }
            Front_1strow [i] = f1rows ;
        }

    }
    else
    {

	
	for (i = 0 ; i <= nfr ; i++)
	{
	    Front_1strow [i] = 0 ;
	}
	
	for (k = 0 ; k < n1 ; k++)
	{
	    row = Rperm1 [k] ;
	    Ci [row] = EMPTY ;
	}
	
	for ( ; k < n_row - nempty_row ; k++)
	{
	    row = Rperm1 [k] ;
	    i = InFront [k - n1] ;
	    ASSERT (i >= EMPTY && i < nfr) ;
	    if (i != EMPTY)
	    {
		Front_1strow [i]++ ;
	    }
	    
	    Ci [row] = i ;
	}
	
	for ( ; k < n_row ; k++)
	{
	    row = Rperm1 [k] ;
	    Ci [row] = nfr ;
	}
	
	for (row = 0 ; row < n_row ; row++)
	{
	    InFront [row] = Ci [row] ;
	}
	
    }

    
    
    

    k = n1 ;
    for (i = 0 ; i < nfr ; i++)
    {
	fpivcol = Fr_npivcol [i] ;
	k += fpivcol ;
	
	Front_npivcol [i] = fpivcol ;
	Front_parent [i] = Fr_parent [i] ;
    }

    
    Front_npivcol [nfr] = n_col - k ;
    Front_parent [nfr] = EMPTY ;

    
    
    

    
    for (k = 0 ; k < n1 ; k++)
    {
	Rperm_init [k] = Rperm1 [k] ;
    }

    
    for (i = 0 ; i < nfr ; i++)
    {
	f1rows = Front_1strow [i] ;
	Front_1strow [i] = k ;
	k += f1rows ;
    }

    
    Front_1strow [nfr] = k ;

    
    F1 = Ci ;				

    for (i = 0 ; i <= nfr ; i++)
    {
	F1 [i] = Front_1strow [i] ;
    }

    for (row = 0 ; row < n_row ; row++)
    {
	i = InFront [row] ;
	if (i != EMPTY)
	{
	    newrow = F1 [i]++ ;
	    Rperm_init [newrow] = row ;
	}
    }
    Rperm_init [n_row] = EMPTY ;	

    

    
    
    

    

    if (prefer_diagonal)
    {
	Int *Diagonal_map ;

	
	Symbolic->Diagonal_map = (Int *) LU_malloc (n_col+1, sizeof (Int)) ;
	Diagonal_map = Symbolic->Diagonal_map ;
	if (Diagonal_map == (Int *) NULL)
	{
	    
	    Info [SparseLU_STATUS] = SparseLU_ERROR_out_of_memory ;
	    error (&Symbolic, SW) ;
	    return (SparseLU_ERROR_out_of_memory) ;
	}

	
	for (newrow = 0 ; newrow < nn ; newrow++)
	{
	    oldrow = Rperm_init [newrow] ;
	    Ci [oldrow] = newrow ;
	}

        for (newcol = 0 ; newcol < nn ; newcol++)
        {
            oldcol = Cperm_init [newcol] ;
            oldrow = oldcol ;
            newrow = Ci [oldrow] ;
            ASSERT (newrow >= 0 && newrow < nn) ;
            Diagonal_map [newcol] = newrow ;
        }
	

    }

    
    
    

    for (i = 0 ; i <= nfr ; i++)
    {
	    Front_leftmostdesc [i] = EMPTY ;
    }

    for (i = 0 ; i < nfr ; i++)
    {
        
        j = i ;
        while (j != EMPTY && Front_leftmostdesc [j] == EMPTY)
        {
            Front_leftmostdesc [j] = i ;
            j = Front_parent [j] ;
        }
    }

    
    
    

    maxnrows = 1 ;		
    maxncols = 1 ;		
    dmaxfrsize = 1 ;		

    
    nchains = 0 ;		
    Chain_start [0] = 0 ;	
    maxrows = 1 ;		
    maxcols = 1 ;		

    for (i = 0 ; i < nfr ; i++)
    {
	
	fpivcol  = Front_npivcol [i] ;	    
	fallrows = Fr_nrows [i] ;	    
	fallcols = Fr_ncols [i] ;	    
	parent = Front_parent [i] ;	    
	fpiv = MIN (fpivcol, fallrows) ;    
	maxrows = MAX (maxrows, fallrows) ;
	maxcols = MAX (maxcols, fallcols) ;

	if (parent != i+1)
	{
	    
	    double s ;

	    
	    if (maxrows % 2 == 0) maxrows++ ;

	    Chain_maxrows [nchains] = maxrows ;
	    Chain_maxcols [nchains] = maxcols ;

	    

	    
	    s = (double) maxrows * (double) maxcols ;
	    dmaxfrsize = MAX (dmaxfrsize, s) ;

	    
	    maxnrows = MAX (maxnrows, maxrows) ;
	    maxncols = MAX (maxncols, maxcols) ;

	    
	    nchains++ ;
	    Chain_start [nchains] = i+1 ;
	    maxrows = 1 ;
	    maxcols = 1 ;
	}
    }

    Chain_maxrows [nchains] = 0 ;
    Chain_maxcols [nchains] = 0 ;

    
    dmaxfrsize = ceil (dmaxfrsize) ;

    
    Symbolic->maxnrows = maxnrows ;
    Symbolic->maxncols = maxncols ;

    
    
    

    if (max_rdeg > dense_row_threshold)
    {
	
	
	
	for (newrow = 0 ; newrow < n_row ; newrow++)
	{
	    oldrow = Rperm_init [newrow] ;
	    Ci [oldrow] = newrow ;
	}
	for (col = n1 ; col < n_col - nempty_col ; col++)
	{
	    oldcol = Cperm_init [col] ;
	    esize = Cdeg [oldcol] ;
	    for (p = Ap [oldcol] ; p < Ap [oldcol+1] ; p++)
	    {
		oldrow = Ai [p] ;
		newrow = Ci [oldrow] ;
		if (newrow >= n1 && Rdeg [oldrow] > dense_row_threshold)
		{
		    esize-- ;
		}
	    }
	    Esize [col - n1] = esize ;
	}
	
    }

    

    
    
    

    
    for (k = 0 ; k < n_col ; k++)
    {
	Ci [k] = Cdeg [Cperm_init [k]] ;
    }
    for (k = 0 ; k < n_col ; k++)
    {
	Cdeg [k] = Ci [k] ;
    }
    for (k = 0 ; k < n_row ; k++)
    {
	Ci [k] = Rdeg [Rperm_init [k]] ;
    }
    for (k = 0 ; k < n_row ; k++)
    {
	Rdeg [k] = Ci [k] ;
    }
    

    
    
    

    

    dlnz = n_inner ;	
    dunz = dlnz ;	

    
    head_usage  = 1 ;
    dhead_usage = 1 ;

    
    tail_usage  = 2 ;
    dtail_usage = 2 ;

    
    tail_usage  +=  UNITS (Int *, n_row+1) +  UNITS (Entry *, n_row+1) + 2 ;
    dtail_usage += DUNITS (Int *, n_row+1) + DUNITS (Entry *, n_row+1) + 2 ;

    
    for (k = 0 ; k < n1 ; k++)
    {
        lnz = Cdeg [k] - 1 ;
        unz = Rdeg [k] - 1 ;
        dlnz += lnz ;
        dunz += unz ;
        head_usage  +=  UNITS (Int, lnz) +  UNITS (Entry, lnz)
                +   UNITS (Int, unz) +  UNITS (Entry, unz) ;
        dhead_usage += DUNITS (Int, lnz) + DUNITS (Entry, lnz)
                +  DUNITS (Int, unz) + DUNITS (Entry, unz) ;
    }

    
    for (k = n1 ; k < n_col - nempty_col; k++)
    {
        esize = Esize ? Esize [k-n1] : Cdeg [k] ;
        if (esize > 0)
        {
            tail_usage  +=  GET_ELEMENT_SIZE (esize, 1) + 1 ;
            dtail_usage += DGET_ELEMENT_SIZE (esize, 1) + 1 ;
        }
    }

    
    if (Esize)
    {
        Int nrow_elements = 0 ;
        for (k = n1 ; k < n_row - nempty_row ; k++)
        {
            rdeg = Rdeg [k] ;
            if (rdeg > dense_row_threshold)
            {
            tail_usage  += GET_ELEMENT_SIZE (1, rdeg) + 1 ;
            dtail_usage += GET_ELEMENT_SIZE (1, rdeg) + 1 ;
            nrow_elements++ ;
            }
        }
        Info [SparseLU_NDENSE_ROW] = nrow_elements ;
    }

    
    if (Esize)
    {
        
        for (row = n1 ; row < n_row ; row++)
        {
            rdeg = Rdeg [row] ;
            tlen = (rdeg > dense_row_threshold) ? 1 : rdeg ;
            tail_usage  += 1 +  UNITS (Tuple, TUPLES (tlen)) ;
            dtail_usage += 1 + DUNITS (Tuple, TUPLES (tlen)) ;
        }
        
        for (col = n1 ; col < n_col - nempty_col ; col++)
        {
            
            esize = Esize [col - n1] ;
            tlen = (esize > 0) + (Cdeg [col] - esize) ;
            tail_usage  += 1 +  UNITS (Tuple, TUPLES (tlen)) ;
            dtail_usage += 1 + DUNITS (Tuple, TUPLES (tlen)) ;
        }
        for ( ; col < n_col ; col++)
        {
            tail_usage  += 1 +  UNITS (Tuple, TUPLES (0)) ;
            dtail_usage += 1 + DUNITS (Tuple, TUPLES (0)) ;
        }
    }
    else
    {
        
        for (row = n1 ; row < n_row ; row++)
        {
            tlen = Rdeg [row] ;
            tail_usage  += 1 +  UNITS (Tuple, TUPLES (tlen)) ;
            dtail_usage += 1 + DUNITS (Tuple, TUPLES (tlen)) ;
        }
        
        for (col = n1 ; col < n_col ; col++)
        {
            tail_usage  += 1 +  UNITS (Tuple, TUPLES (1)) ;
            dtail_usage += 1 + DUNITS (Tuple, TUPLES (1)) ;
        }
    }

    Symbolic->num_mem_init_usage = head_usage + tail_usage ;

    
    dmax_usage = dhead_usage + dtail_usage ;
    dmax_usage = MAX (Symbolic->num_mem_init_usage, ceil (dmax_usage)) ;
    Info [SparseLU_VARIABLE_INIT_ESTIMATE] = dmax_usage ;

    
    Symbolic->dnum_mem_init_usage = dmax_usage ;

    
    tail_usage  -=  UNITS (Int *, n_row+1) +  UNITS (Entry *, n_row+1) ;
    dtail_usage -= DUNITS (Int *, n_row+1) + DUNITS (Entry *, n_row+1) ;

    
    
    

    
    Link = Ci ;
    for (i = 0 ; i < nfr ; i++)
    {
	Link [i] = EMPTY ;
    }

    flops = 0 ;			

    for (chain = 0 ; chain < nchains ; chain++)
    {
        double fsize ;
        f1 = Chain_start [chain] ;
        f2 = Chain_start [chain+1] - 1 ;

        
        dr = Chain_maxrows [chain] ;
        dc = Chain_maxcols [chain] ;
        fsize =
            nb*nb	    
            + dr*nb	    
            + nb*dc	    
            + dr*dc ;	    
        dtail_usage += DUNITS (Entry, fsize) ;
        dmax_usage = MAX (dmax_usage, dhead_usage + dtail_usage) ;

        for (i = f1 ; i <= f2 ; i++)
        {
            
            fpivcol  = Front_npivcol [i] ; 
            fallrows = Fr_nrows [i] ;   
            fallcols = Fr_ncols [i] ;   
            parent = Front_parent [i] ; 
            fpiv = MIN (fpivcol, fallrows) ;	
            f = (double) fpiv ;
            r = fallrows - fpiv ;		
            c = fallcols - fpiv ;		

            
            for (child = Link [i] ; child != EMPTY ; child = Link [child])
            {
            
            cp = MIN (Front_npivcol [child], Fr_nrows [child]) ;
            cr = Fr_nrows [child] - cp ;
            cc = Fr_ncols [child] - cp ;
            dtail_usage -= ELEMENT_SIZE (cr, cc) ;

            }

            

            
            flops += DIV_FLOPS * (f*r + (f-1)*f/2)  
            
            + MULTSUB_FLOPS * (f*r*c + (r+c)*(f-1)*f/2 + (f-1)*f*(2*f-1)/6);

            
            dlf = (f*f-f)/2 + f*r ;		
            duf = (f*f-f)/2 + f*c ;		
            dlnz += dlf ;
            dunz += duf ;

            
            dhead_usage +=
            DUNITS (Entry, dlf + duf)   
            + DUNITS (Int, r + c + f) ; 

            if (parent != EMPTY)
            {
            
            dtail_usage += ELEMENT_SIZE (r, c) ;

            
            Link [i] = Link [parent] ;
            Link [parent] = i ;
            }

            
            dmax_usage = MAX (dmax_usage, dhead_usage + dtail_usage) ;
        }
        
        dtail_usage -= DUNITS (Entry, fsize) ;
    }

    dhead_usage = ceil (dhead_usage) ;
    dmax_usage = ceil (dmax_usage) ;
    Symbolic->num_mem_size_est = dhead_usage ;
    Symbolic->num_mem_usage_est = dmax_usage ;
    Symbolic->lunz_bound = dlnz + dunz - n_inner ;

    

    
    
    

    LU_set_stats (
	Info,
	Symbolic,
	dmax_usage,		
	dhead_usage,		
	flops,			
	dlnz,			
	dunz,			
	dmaxfrsize,		
	(double) n_col,		
	(double) n_inner,	
	(double) maxnrows,	
	(double) maxncols,	
	TRUE,			
	prefer_diagonal,
	ESTIMATE) ;

    

    
    
    

    Symbolic->valid = SYMBOLIC_VALID ;
    *SymbolicHandle = (void *) Symbolic ;

    
    
    

    

    free_work (SW) ;

    return (SparseLU_OK) ;
}








PRIVATE int inverse_permutation
(
    Int *P,     
    Int *Pinv,  
    Int n       
)
{
    Int i, k ;
    for (i = 0 ; i < n ; i++)
    {
        Pinv [i] = EMPTY ;
    }
    for (k = 0 ; k < n ; k++)
    {
        i = P [k] ;
        if (i < 0 || i >= n || Pinv [i] != EMPTY)
        {
            
            return (FALSE) ;
        }
        Pinv [i] = k ;
    }
    return (TRUE) ;
}








PRIVATE int do_amd_1
(
    Int n,		
    Int Ap [ ],	        
    Int Ai [ ],	        
    Int P [ ],		
    Int Pinv [ ],	
    Int Len [ ],	
    Int slen,		
    Int S [ ],		
    Int ordering_option,
    Int print_level,

    
    int (*user_ordering)    
    (
        
        Int,            
        Int,            
        Int,            
        Int *,          
        Int *,          
        
        Int *,          
        
        void *,         
        double *        
    ),
    void *user_params,  

    Int *ordering_used,

    double amd_Control [ ],	
    double amd_Info [ ] 	
)
{
    Int i, j, k, p, pfree, iwlen, pj, p1, p2, pj2, anz, *Iw, *Pe, *Nv, *Head,
	*Elen, *Degree, *s, *W, *Sp, *Tp ;

    
    
    

    ASSERT (n > 0) ;

    s = S ;
    Pe = s ;	    s += (n+1) ;    slen -= (n+1) ;
    Nv = s ;	    s += n ;        slen -= n ;

    if (user_ordering == NULL)
    {
        
        Head = s ;      s += n ;    slen -= n ;
        Elen = s ;      s += n ;    slen -= n ;
        Degree = s ;    s += n ;    slen -= n ;
    }
    else
    {
        
        Head = NULL ;
        Elen = NULL ;
        Degree = NULL ;
    }

    W = s ;	    
    s += n ;        
    slen -= n ;

    iwlen = slen ;
    Iw = s ;	    
    s += iwlen ;
    anz = Ap [n] ;

    
    Sp = Nv ;			
    Tp = W ;
    pfree = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	Pe [j] = pfree ;
	Sp [j] = pfree ;
	pfree += Len [j] ;
    }
    Pe [n] = pfree ;

    for (k = 0 ; k < n ; k++)
    {
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;

	
	for (p = p1 ; p < p2 ; )
	{
	    
	    j = Ai [p] ;
	    if (j < k)
	    {
		
		Iw [Sp [j]++] = k ;
		Iw [Sp [k]++] = j ;
		p++ ;
	    }
	    else if (j == k)
	    {
		
		p++ ;
		break ;
	    }
	    else 
	    {
		
		break ;
	    }
	    
	    ASSERT (Ap [j] <= Tp [j] && Tp [j] <= Ap [j+1]) ;
	    pj2 = Ap [j+1] ;
	    for (pj = Tp [j] ; pj < pj2 ; )
	    {
		i = Ai [pj] ;
		ASSERT (i >= 0 && i < n) ;
		if (i < k)
		{
		    
		    ASSERT (Sp [i] < (i == n-1 ? pfree : Pe [i+1])) ;
		    ASSERT (Sp [j] < (j == n-1 ? pfree : Pe [j+1])) ;
		    Iw [Sp [i]++] = j ;
		    Iw [Sp [j]++] = i ;
		    pj++ ;
		}
		else if (i == k)
		{
		    
		    pj++ ;
		    break ;
		}
		else 
		{
		    
		    break ;
		}
	    }
	    Tp [j] = pj ;
	}
	Tp [k] = p ;
    }

    
    for (j = 0 ; j < n ; j++)
    {
	for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
	{
	    i = Ai [pj] ;
	    ASSERT (i >= 0 && i < n) ;
	    
	    ASSERT (Sp [i] < (i == n-1 ? pfree : Pe [i+1])) ;
	    ASSERT (Sp [j] < (j == n-1 ? pfree : Pe [j+1])) ;
	    Iw [Sp [i]++] = j ;
	    Iw [Sp [j]++] = i ;
	}
    }

    

    
    
    

    if (ordering_option == SparseLU_ORDERING_AMD)
    {

        
        AMD_2 (n, Pe, Iw, Len, iwlen, pfree,
            Nv, Pinv, P, Head, Elen, Degree, W, amd_Control, amd_Info) ;
        *ordering_used = SparseLU_ORDERING_AMD ;
        return (TRUE) ;

    }
    else
    {

        
        double user_info [3], dmax, lnz, flops ;
        int ok ;
        user_info [0] = EMPTY ;
        user_info [1] = EMPTY ;
        user_info [2] = EMPTY ;

        if (ordering_option == SparseLU_ORDERING_USER)
        {
            ok = (*user_ordering) (n, n, TRUE, Pe, Iw, P,
                user_params, user_info) ;
            *ordering_used = SparseLU_ORDERING_USER ;
        }
        else
        {
            Int params [3] ;
            params [0] = ordering_option ;
            params [1] = print_level ;
            ok = LU_chol (n, n, TRUE, Pe, Iw, P, &params, user_info) ;
            *ordering_used = params [2] ;
        }

        if (!ok)
        {
            
            amd_Info [AMD_STATUS] = AMD_INVALID ;
            return (FALSE) ;
        }

        
        dmax  = user_info [0] ;
        lnz   = user_info [1] ;
        flops = user_info [2] ;

        
        amd_Info [AMD_STATUS] = AMD_OK ;
        amd_Info [AMD_N] = n ;
        amd_Info [AMD_NZ] = anz ;
        
        
        amd_Info [AMD_NZ_A_PLUS_AT] = pfree ;
        amd_Info [AMD_NDENSE] = 0 ;
        
        amd_Info [AMD_NCMPA] = 0 ;
        amd_Info [AMD_LNZ] = lnz ;
        amd_Info [AMD_NDIV] = lnz ;
        if (flops >= 0)
        {
            amd_Info [AMD_NMULTSUBS_LDL] = (flops - n) / 2 ;
            amd_Info [AMD_NMULTSUBS_LU]  = (flops - n) ;
        }
        else
        {
            amd_Info [AMD_NMULTSUBS_LDL] = EMPTY ;
            amd_Info [AMD_NMULTSUBS_LU]  = EMPTY ;
        }
        amd_Info [AMD_DMAX] = dmax ;

        
        return (inverse_permutation (P, Pinv, n)) ;
    }
}






PRIVATE int do_amd
(
    Int n,
    Int Ap [ ],		        
    Int Ai [ ],		        
    Int Q [ ],			
    Int Qinv [ ],		
    Int Sdeg [ ],		
    Int Clen,			
    Int Ci [ ],			
    double amd_Control [ ],	
    double amd_Info [ ],	
    SymbolicType *Symbolic,	
    double Info [ ],		
    Int ordering_option,
    Int print_level,

    
    int (*user_ordering)    
    (
        
        Int,            
        Int,            
        Int,            
        Int *,          
        Int *,          
        
        Int *,          
        
        void *,         
        double *        
    ),
    void *user_params,  
    Int *ordering_used
)
{
    int ok = TRUE ;
    *ordering_used = SparseLU_ORDERING_NONE ;

    if (n == 0)
    {
	Symbolic->amd_dmax = 0 ;
	Symbolic->amd_lunz = 0 ;
	Info [SparseLU_SYMMETRIC_LUNZ] = 0 ;
	Info [SparseLU_SYMMETRIC_FLOPS] = 0 ;
	Info [SparseLU_SYMMETRIC_DMAX] = 0 ;
	Info [SparseLU_SYMMETRIC_NDENSE] = 0 ;
    }
    else
    {
	ok = do_amd_1 (n, Ap, Ai, Q, Qinv, Sdeg, Clen,
            Ci, ordering_option, print_level, user_ordering, user_params,
            ordering_used, amd_Control, amd_Info) ;

        
        if (ok)
        {
            Symbolic->amd_dmax = amd_Info [AMD_DMAX] ;
            Symbolic->amd_lunz = 2 * amd_Info [AMD_LNZ] + n ;
            Info [SparseLU_SYMMETRIC_LUNZ] = Symbolic->amd_lunz ;
            Info [SparseLU_SYMMETRIC_FLOPS] = DIV_FLOPS * amd_Info [AMD_NDIV] +
                MULTSUB_FLOPS * amd_Info [AMD_NMULTSUBS_LU] ;
            Info [SparseLU_SYMMETRIC_DMAX] = Symbolic->amd_dmax ;
            Info [SparseLU_SYMMETRIC_NDENSE] = amd_Info [AMD_NDENSE] ;
            Info [SparseLU_SYMBOLIC_DEFRAG] += amd_Info [AMD_NCMPA] ;
        }
    }
    return (ok) ;
}







PRIVATE Int prune_singletons
(
    Int n1,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    Int Cperm1 [ ],
    Int InvRperm1 [ ],
    Int Si [ ],
    Int Sp [ ]
)
{
    Int row, k, pp, p, oldcol, newcol, newrow, nzdiag, do_nzdiag ;

    nzdiag = 0 ;
    do_nzdiag = (Ax != (double *) NULL) ;

    

    pp = 0 ;
    for (k = n1 ; k < n_col ; k++)
    {
	oldcol = Cperm1 [k] ;
	newcol = k - n1 ;
	Sp [newcol] = pp ;  
	for (p = Ap [oldcol] ; p < Ap [oldcol+1] ; p++)
	{
	    row = Ai [p] ;
	    newrow = InvRperm1 [row] - n1 ;
	    if (newrow >= 0)
	    {
		Si [pp++] = newrow ;
		if (do_nzdiag)
		{
		    
		    if (newrow == newcol)
		    {
			
			if (SCALAR_IS_NONZERO (Ax [p]))
			{
			    nzdiag++ ;
			}
		    }
		}
	    }
	}
    }
    Sp [n_col - n1] = pp ;

    return (nzdiag) ;
}





PRIVATE void combine_ordering
(
    Int n1,
    Int nempty_col,
    Int n_col,
    Int Cperm_init [ ],	    
    Int Cperm1 [ ],	    
    Int Qinv [ ]	    
)
{
    Int k, oldcol, newcol, knew ;

    
    for (k = 0 ; k < n1 ; k++)
    {
	Cperm_init [k] = Cperm1 [k] ;
    }
    for (k = n1 ; k < n_col - nempty_col ; k++)
    {
	
	oldcol = Cperm1 [k] ;	
	newcol = k - n1 ;	
	knew = Qinv [newcol] ;	
	knew += n1 ;		
	Cperm_init [knew] = oldcol ;
    }
    for (k = n_col - nempty_col ; k < n_col ; k++)
    {
	Cperm_init [k] = Cperm1 [k] ;
    }

}





PRIVATE void free_work
(
    SWType *SW
)
{
    if (SW)
    {
        SW->InvRperm1 = (Int *) LU_free ((void *) SW->InvRperm1) ;
        SW->Rs = (double *) LU_free ((void *) SW->Rs) ;
        SW->Si = (Int *) LU_free ((void *) SW->Si) ;
        SW->Sp = (Int *) LU_free ((void *) SW->Sp) ;
        SW->Ci = (Int *) LU_free ((void *) SW->Ci) ;
        SW->Front_npivcol = (Int *) LU_free ((void *) SW->Front_npivcol);
        SW->Front_nrows = (Int *) LU_free ((void *) SW->Front_nrows) ;
        SW->Front_ncols = (Int *) LU_free ((void *) SW->Front_ncols) ;
        SW->Front_parent = (Int *) LU_free ((void *) SW->Front_parent) ;
        SW->Front_cols = (Int *) LU_free ((void *) SW->Front_cols) ;
        SW->Cperm1 = (Int *) LU_free ((void *) SW->Cperm1) ;
        SW->Rperm1 = (Int *) LU_free ((void *) SW->Rperm1) ;
        SW->InFront = (Int *) LU_free ((void *) SW->InFront) ;
    }
}








PRIVATE void error
(
    SymbolicType **Symbolic,
    SWType *SW
)
{
    free_work (SW) ;
    SparseLU_free_symbolic ((void **) Symbolic) ;
    ASSERT (LU_malloc_count == init_count) ;
}



