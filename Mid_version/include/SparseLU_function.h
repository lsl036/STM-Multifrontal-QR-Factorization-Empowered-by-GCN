/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    SparseLU_function.h
 * BRIEF:   稀疏LU函数头文件
 *****************************************************************************/


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
) ;

GLOBAL void LU_apply_order
(
    Int Front [ ],
    const Int Order [ ],
    Int Temp [ ],
    Int n_col,
    Int nfr
) ;

GLOBAL void LU_assemble
(
    NumericType *Numeric,
    WorkType *Work
) ;

GLOBAL void LU_assemble_fixq
(
    NumericType *Numeric,
    WorkType *Work
) ;

GLOBAL void LU_blas3_update
(
    WorkType *Work
) ;

GLOBAL Int LU_build_tuples
(
    NumericType *Numeric,
    WorkType *Work
) ;

GLOBAL Int LU_create_element
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
) ;

GLOBAL Int LU_extend_front
(
    NumericType *Numeric,
    WorkType *Work
) ;

GLOBAL void *LU_free
(
    void *p
) ;

GLOBAL void LU_fsize
(
    Int nn,
    Int MaxFsize [ ],
    Int Fnrows [ ],
    Int Fncols [ ],
    Int Parent [ ],
    Int Npiv [ ]
) ;

GLOBAL void LU_garbage_collection
(
    NumericType *Numeric,
    WorkType *Work,
    Int drnew,
    Int dcnew,
    Int do_Fcpos
) ;

GLOBAL Int LU_get_memory
(
    NumericType *Numeric,
    WorkType *Work,
    Int needunits,
    Int r2,
    Int c2,
    Int do_Fcpos
) ;

GLOBAL Int LU_grow_front
(
    NumericType *Numeric,
    Int fnr2,
    Int fnc2,
    WorkType *Work,
    Int do_what
) ;

GLOBAL Int LU_init_front
(
    NumericType *Numeric,
    WorkType *Work
) ;

GLOBAL Int LU_is_permutation
(
    const Int P [ ],
    Int W [ ],
    Int n,
    Int r
) ;

GLOBAL Int LU_kernel
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
) ;

GLOBAL Int LU_kernel_init
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
) ;

GLOBAL void LU_kernel_wrapup
(
    NumericType *Numeric,
    SymbolicType *Symbolic,
    WorkType *Work
) ;

GLOBAL Int LU_local_search
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
) ;

#ifndef _LU_MALLOC
#define _LU_MALLOC

#if defined (LU_MALLOC_COUNT) || !defined (NDEBUG)

#ifndef EXTERN
#define EXTERN extern
#endif

GLOBAL EXTERN Int LU_malloc_count ;
#endif

GLOBAL void *LU_malloc
(
    Int n_objects,
    size_t size_of_object
) ;

#endif

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
) ;

GLOBAL Int LU_mem_alloc_head_block
(
    NumericType *Numeric,
    Int nunits
) ;

GLOBAL Int LU_mem_alloc_tail_block
(
    NumericType *Numeric,
    Int nunits
) ;

GLOBAL void LU_mem_free_tail_block
(
    NumericType *Numeric,
    Int i
) ;

GLOBAL void LU_mem_init_memoryspace
(
    NumericType *Numeric
) ;

GLOBAL void *LU_realloc
(
    void *p,
    Int n_objects,
    size_t size_of_object
) ;

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
    Int freebie [2]
) ;

#define IN 0
#define OUT 1

#define IN_IN 0
#define IN_OUT 1
#define OUT_IN 2
#define OUT_OUT 3

GLOBAL void LU_scale
(
    Int n,
    Entry alpha,
    Entry X [ ]
) ;

GLOBAL void LU_scale_column
(
    NumericType *Numeric,
    WorkType *Work
) ;

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
) ;

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
    Int *n1,
    Int *n1c,
    Int *n1r,
    Int *nempty_col,
    Int *nempty_row,
    Int *is_sym,
    Int *max_rdeg,
    Int Rp [ ],
    Int Ri [ ],
    Int W [ ],
    Int Next [ ]
) ;

GLOBAL Int LU_start_front
(
    Int chain,
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
) ;

GLOBAL Int LU_store_lu
(
    NumericType *Numeric,
    WorkType *Work
) ;

GLOBAL Int LU_store_lu_drop
(
    NumericType *Numeric,
    WorkType *Work
) ;

GLOBAL double LU_symbolic_usage
(
    Int n_row,
    Int n_col,
    Int nchains,
    Int nfr,
    Int esize,
    Int prefer_diagonal
) ;

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
) ;

GLOBAL Int LU_triplet_map_x
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
    , const double Tx [ ]
    , double Ax [ ]
    , double Rx [ ]
    , Int Map [ ]
    , Int Map2 [ ]
) ;

GLOBAL Int LU_triplet_map_nox
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
    , Int Map [ ]
    , Int Map2 [ ]
) ;

GLOBAL Int LU_triplet_nomap_x
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
    , const double Tx [ ]
    , double Ax [ ]
    , double Rx [ ]
) ;

GLOBAL Int LU_triplet_nomap_nox
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
) ;

GLOBAL Int LU_tuple_lengths
(
    NumericType *Numeric,
    WorkType *Work,
    double *dusage
) ;

GLOBAL Int LU_valid_numeric
(
    NumericType *Numeric
) ;

GLOBAL Int LU_valid_symbolic
(
    SymbolicType *Symbolic
) ;


GLOBAL Int LU_solve
(
    Int sys,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    double Xx [ ],
    const double Bx [ ],
    NumericType *Numeric,
    Int irstep,
    double Info [SparseLU_INFO],
    Int Pattern [ ],
    double SolveWork [ ]
) ; // lu_solve.c

GLOBAL double LU_lsolve
(
    NumericType *Numeric,
    Entry X [ ],
    Int Pattern [ ]
) ;   // core.c

GLOBAL double LU_usolve
(
    NumericType *Numeric,
    Entry X [ ],
    Int Pattern [ ]
) ;  //core.c

// solvebase.c
GLOBAL double LU_ltsolve
(
    NumericType *Numeric,
    Entry X [ ],
    Int Pattern [ ]
) ;

GLOBAL double LU_lhsolve
(
    NumericType *Numeric,
    Entry X [ ],
    Int Pattern [ ]
) ;

GLOBAL double LU_utsolve
(
    NumericType *Numeric,
    Entry X [ ],
    Int Pattern [ ]
) ;

GLOBAL double LU_uhsolve
(
    NumericType *Numeric,
    Entry X [ ],
    Int Pattern [ ]
) ;

int LU_chol
(
    
    Sparse_long nrow,      
    Sparse_long ncol,      
    Sparse_long symmetric, 
    Sparse_long Ap [ ],    
    Sparse_long Ai [ ],    
    
    Sparse_long Perm [ ],  
    
    void *ignore,           
    double user_info [3]    
) ;





#ifndef COLAMD_H
#define COLAMD_H





#include <stdlib.h>






#define COLAMD_KNOBS 20


#define COLAMD_STATS 20


#define COLAMD_DENSE_ROW 0


#define COLAMD_DENSE_COL 1


#define COLAMD_AGGRESSIVE 2


#define COLAMD_DEFRAG_COUNT 2


#define COLAMD_STATUS 3


#define COLAMD_INFO1 4
#define COLAMD_INFO2 5
#define COLAMD_INFO3 6




#define COLAMD_EMPTY_ROW 7

#define COLAMD_EMPTY_COL 8

#define COLAMD_NEWLY_EMPTY_ROW 9

#define COLAMD_NEWLY_EMPTY_COL 10



#define COLAMD_OK				(0)
#define COLAMD_ERROR_jumbled_matrix		(-11)
#define COLAMD_ERROR_A_not_present		(-1)
#define COLAMD_ERROR_p_not_present		(-2)
#define COLAMD_ERROR_nrow_negative		(-3)
#define COLAMD_ERROR_ncol_negative		(-4)
#define COLAMD_ERROR_nnz_negative		(-5)
#define COLAMD_ERROR_p0_nonzero			(-6)
#define COLAMD_ERROR_A_too_small		(-7)
#define COLAMD_ERROR_col_length_negative	(-8)
#define COLAMD_ERROR_row_index_out_of_bounds	(-9)
#define COLAMD_ERROR_out_of_memory		(-10)
#define COLAMD_ERROR_internal_error		(-999)









typedef struct Colamd_Col_struct
{
    Int start ;		
			
    Int length ;	
    union
    {
	Int thickness ;	
			
	Int parent ;	
			
    } shared1 ;
    union
    {
	Int score ;	
	Int order ;	
    } shared2 ;
    union
    {
	Int headhash ;	
			
	Int hash ;	
	Int prev ;	
			
    } shared3 ;
    union
    {
	Int degree_next ;	
	Int hash_next ;		
    } shared4 ;

    
    
    Int nextcol ;	
    Int lastcol ;	
    

} Colamd_Col ;

typedef struct Colamd_Row_struct
{
    Int start ;		
    Int length ;	
    union
    {
	Int degree ;	
	Int p ;		
    } shared1 ;
    union
    {
	Int mark ;	
	Int first_column ;
    } shared2 ;

    
    
    Int thickness ;	
			
    Int front ;		
			
			
    

} Colamd_Row ;










#define LU_COLAMD_C(n_col) ((n_col + 1) * sizeof (Colamd_Col) / sizeof (Int))


#define LU_COLAMD_R(n_row) ((n_row + 1) * sizeof (Colamd_Row) / sizeof (Int))


#define LU_COLAMD_RECOMMENDED(nnz, n_row, n_col)	\
(							\
((nnz) < 0 || (n_row) < 0 || (n_col) < 0)		\
?							\
    (-1)						\
:							\
    (MAX (2 * (nnz), 4 * (n_col)) +			\
    (Int) LU_COLAMD_C (n_col) +			\
    (Int) LU_COLAMD_R (n_row) + (n_col) + ((nnz) / 5))	\
)







void LU_colamd_set_defaults	
(				
    double knobs [COLAMD_KNOBS]	
) ;

Int LU_colamd			
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
    
) ;

#endif 

