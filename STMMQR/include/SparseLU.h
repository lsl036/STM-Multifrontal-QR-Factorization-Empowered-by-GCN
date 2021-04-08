/******************************************************************************
 * VERSION: 1.0
 * DATE:    2020年9月24日
 * FILE:    SparseLU.h
 * BRIEF:   稀疏LU头文件
 *****************************************************************************/

#ifndef SparseLU_H
#define SparseLU_H

#ifdef __cplusplus
extern "C" {
#endif

#include "SparseBase_config.h"

#define SparseLU_INFO 90
#define SparseLU_CONTROL 20


#include "amd.h"

#define SparseLU_STATUS 0	
#define SparseLU_NROW 1		
#define SparseLU_NCOL 16		
#define SparseLU_NZ 2		


#define SparseLU_SIZE_OF_UNIT 3		


#define SparseLU_SIZE_OF_INT 4		
#define SparseLU_SIZE_OF_LONG 5		
#define SparseLU_SIZE_OF_POINTER 6	
#define SparseLU_SIZE_OF_ENTRY 7		
#define SparseLU_NDENSE_ROW 8		
#define SparseLU_NEMPTY_ROW 9		
#define SparseLU_NDENSE_COL 10		
#define SparseLU_NEMPTY_COL 11		
#define SparseLU_SYMBOLIC_DEFRAG 12	
#define SparseLU_SYMBOLIC_PEAK_MEMORY 13	
#define SparseLU_SYMBOLIC_SIZE 14	
#define SparseLU_SYMBOLIC_TIME 15	
#define SparseLU_SYMBOLIC_WALLTIME 17	
#define SparseLU_STRATEGY_USED 18	
#define SparseLU_ORDERING_USED 19	
#define SparseLU_QFIXED 31		
#define SparseLU_DIAG_PREFERRED 32	
#define SparseLU_PATTERN_SYMMETRY 33	
#define SparseLU_NZ_A_PLUS_AT 34		
#define SparseLU_NZDIAG 35		


#define SparseLU_SYMMETRIC_LUNZ 36	
#define SparseLU_SYMMETRIC_FLOPS 37	
#define SparseLU_SYMMETRIC_NDENSE 38	
#define SparseLU_SYMMETRIC_DMAX 39	




#define SparseLU_COL_SINGLETONS 56	
#define SparseLU_ROW_SINGLETONS 57	
#define SparseLU_N2 58			
#define SparseLU_S_SYMMETRIC 59		


#define SparseLU_NUMERIC_SIZE_ESTIMATE 20    
#define SparseLU_PEAK_MEMORY_ESTIMATE 21	    
#define SparseLU_FLOPS_ESTIMATE 22	    
#define SparseLU_LNZ_ESTIMATE 23		    
#define SparseLU_UNZ_ESTIMATE 24		    
#define SparseLU_VARIABLE_INIT_ESTIMATE 25   
#define SparseLU_VARIABLE_PEAK_ESTIMATE 26   
#define SparseLU_VARIABLE_FINAL_ESTIMATE 27  
#define SparseLU_MAX_FRONT_SIZE_ESTIMATE 28  
#define SparseLU_MAX_FRONT_NROWS_ESTIMATE 29 
#define SparseLU_MAX_FRONT_NCOLS_ESTIMATE 30 


#define SparseLU_NUMERIC_SIZE 40		    
#define SparseLU_PEAK_MEMORY 41		    
#define SparseLU_FLOPS 42		    
#define SparseLU_LNZ 43			    
#define SparseLU_UNZ 44			    
#define SparseLU_VARIABLE_INIT 45	    
#define SparseLU_VARIABLE_PEAK 46	    
#define SparseLU_VARIABLE_FINAL 47	    
#define SparseLU_MAX_FRONT_SIZE 48	    
#define SparseLU_MAX_FRONT_NROWS 49	    
#define SparseLU_MAX_FRONT_NCOLS 50	    


#define SparseLU_NUMERIC_DEFRAG 60	    
#define SparseLU_NUMERIC_REALLOC 61	    
#define SparseLU_NUMERIC_COSTLY_REALLOC 62   
#define SparseLU_COMPRESSED_PATTERN 63	    
#define SparseLU_LU_ENTRIES 64		    
#define SparseLU_NUMERIC_TIME 65		    
#define SparseLU_UDIAG_NZ 66		    
#define SparseLU_RCOND 67		    
#define SparseLU_WAS_SCALED 68		    
#define SparseLU_RSMIN 69		    
#define SparseLU_RSMAX 70		    
#define SparseLU_UMIN 71			    
#define SparseLU_UMAX 72			    
#define SparseLU_ALLOC_INIT_USED 73	    
#define SparseLU_FORCED_UPDATES 74	    
#define SparseLU_NUMERIC_WALLTIME 75	    
#define SparseLU_NOFF_DIAG 76		    

#define SparseLU_ALL_LNZ 77		    
#define SparseLU_ALL_UNZ 78		    
#define SparseLU_NZDROPPED 79		    


#define SparseLU_IR_TAKEN 80	    
#define SparseLU_IR_ATTEMPTED 81	    
#define SparseLU_OMEGA1 82	    
#define SparseLU_OMEGA2 83	    
#define SparseLU_SOLVE_FLOPS 84	    
#define SparseLU_SOLVE_TIME 85	    
#define SparseLU_SOLVE_WALLTIME 86   

#define SparseLU_PRL 0			


#define SparseLU_DENSE_ROW 1		
#define SparseLU_DENSE_COL 2		
#define SparseLU_BLOCK_SIZE 4		
#define SparseLU_STRATEGY 5		
#define SparseLU_ORDERING 10             
#define SparseLU_FIXQ 13			
#define SparseLU_AMD_DENSE 14		
#define SparseLU_AGGRESSIVE 19		
#define SparseLU_SINGLETONS 11           


#define SparseLU_PIVOT_TOLERANCE 3	
#define SparseLU_ALLOC_INIT 6		
#define SparseLU_SYM_PIVOT_TOLERANCE 15	
#define SparseLU_SCALE 16		
#define SparseLU_FRONT_ALLOC_INIT 17	
#define SparseLU_DROPTOL 18		


#define SparseLU_IRSTEP 7		


#define SparseLU_COMPILED_WITH_BLAS 8	    

#define SparseLU_STRATEGY_AUTO 0		
#define SparseLU_STRATEGY_UNSYMMETRIC 1	
#define SparseLU_STRATEGY_SYMMETRIC 2	


#define SparseLU_SCALE_NONE 0	
#define SparseLU_SCALE_SUM 1	
#define SparseLU_SCALE_MAX 2	


#define SparseLU_ORDERING_AMD 1          
#define SparseLU_ORDERING_GIVEN 2        
#define SparseLU_ORDERING_BEST 4         
#define SparseLU_ORDERING_NONE 5         
#define SparseLU_ORDERING_USER 6         


#define SparseLU_DEFAULT_PRL 1
#define SparseLU_DEFAULT_DENSE_ROW 0.2
#define SparseLU_DEFAULT_DENSE_COL 0.2
#define SparseLU_DEFAULT_PIVOT_TOLERANCE 0.1
#define SparseLU_DEFAULT_SYM_PIVOT_TOLERANCE 0.001
#define SparseLU_DEFAULT_BLOCK_SIZE 32
#define SparseLU_DEFAULT_ALLOC_INIT 0.7
#define SparseLU_DEFAULT_FRONT_ALLOC_INIT 0.5
#define SparseLU_DEFAULT_IRSTEP 2
#define SparseLU_DEFAULT_SCALE SparseLU_SCALE_SUM
#define SparseLU_DEFAULT_STRATEGY SparseLU_STRATEGY_AUTO
#define SparseLU_DEFAULT_AMD_DENSE AMD_DEFAULT_DENSE
#define SparseLU_DEFAULT_FIXQ 0
#define SparseLU_DEFAULT_AGGRESSIVE 1
#define SparseLU_DEFAULT_DROPTOL 0
#define SparseLU_DEFAULT_ORDERING SparseLU_ORDERING_AMD
#define SparseLU_DEFAULT_SINGLETONS TRUE

#define SparseLU_OK (0)

#define SparseLU_WARNING_singular_matrix (1)

#define SparseLU_WARNING_determinant_underflow (2)
#define SparseLU_WARNING_determinant_overflow (3)

#define SparseLU_ERROR_out_of_memory (-1)
#define SparseLU_ERROR_invalid_Numeric_object (-3)
#define SparseLU_ERROR_invalid_Symbolic_object (-4)
#define SparseLU_ERROR_argument_missing (-5)
#define SparseLU_ERROR_n_nonpositive (-6)
#define SparseLU_ERROR_invalid_matrix (-8)
#define SparseLU_ERROR_different_pattern (-11)
#define SparseLU_ERROR_invalid_system (-13)
#define SparseLU_ERROR_invalid_permutation (-15)
#define SparseLU_ERROR_internal_error (-911) 
#define SparseLU_ERROR_file_IO (-17)

#define SparseLU_ERROR_ordering_failed (-18)

#define SparseLU_A	(0)	
#define SparseLU_At	(1)	
#define SparseLU_Aat	(2)	

#define SparseLU_Pt_L	(3)	
#define SparseLU_L	(4)	
#define SparseLU_Lt_P	(5)	
#define SparseLU_Lat_P	(6)	
#define SparseLU_Lt	(7)	
#define SparseLU_Lat	(8)	

#define SparseLU_U_Qt	(9)	
#define SparseLU_U	(10)	
#define SparseLU_Q_Ut	(11)	
#define SparseLU_Q_Uat	(12)	
#define SparseLU_Ut	(13)	
#define SparseLU_Uat	(14)	

Sparse_long sparselu_symbolic
(
    Sparse_long n_row,
    Sparse_long n_col,
    const Sparse_long Ap [ ],
    const Sparse_long Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [SparseLU_CONTROL],
    double Info [SparseLU_INFO]
) ;

Sparse_long sparselu_numeric
(
    const Sparse_long Ap [ ],
    const Sparse_long Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [SparseLU_CONTROL],
    double Info [SparseLU_INFO]
) ;

void sparselu_free_symbolic
(
    void **Symbolic
) ;

void sparselu_free_numeric
(
    void **Numeric
) ;

void sparselu_defaults
(
    double Control [SparseLU_CONTROL]
) ;

Sparse_long sparselu_solve
(
    Sparse_long sys,
    const Sparse_long Ap [ ],
    const Sparse_long Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [SparseLU_CONTROL],
    double Info [SparseLU_INFO]
) ;

Sparse_long sparselu_scale
(
    double X [ ],
    const double B [ ],
    void *Numeric
) ;

Sparse_long sparselu_transpose
(
    Sparse_long n_row,
    Sparse_long n_col,
    const Sparse_long Ap [ ],
    const Sparse_long Ai [ ],
    const double Ax [ ],
    const Sparse_long P [ ],
    const Sparse_long Q [ ],
    Sparse_long Rp [ ],
    Sparse_long Ri [ ],
    double Rx [ ]
) ;

Sparse_long sparselu_triplet_to_col
(
    Sparse_long n_row,
    Sparse_long n_col,
    Sparse_long nz,
    const Sparse_long Ti [ ],
    const Sparse_long Tj [ ],
    const double Tx [ ],
    Sparse_long Ap [ ],
    Sparse_long Ai [ ],
    double Ax [ ],
    Sparse_long Map [ ]
) ;

Sparse_long sparselu_col_to_triplet
(
    Sparse_long n_col,
    const Sparse_long Ap [ ],
    Sparse_long Tj [ ]
) ;

Sparse_long sparselu_get_determinant
(
    double *Mx,
    double *Ex,
    void *NumericHandle,
    double User_Info [SparseLU_INFO]
) ;

Sparse_long sparselu_get_lunz
(
    Sparse_long *lnz,
    Sparse_long *unz,
    Sparse_long *n_row,
    Sparse_long *n_col,
    Sparse_long *nz_udiag,
    void *Numeric
) ;

Sparse_long sparselu_get_symbolic
(
    Sparse_long *n_row,
    Sparse_long *n_col,
    Sparse_long *n1,
    Sparse_long *nz,
    Sparse_long *nfr,
    Sparse_long *nchains,
    Sparse_long P [ ],
    Sparse_long Q [ ],
    Sparse_long Front_npivcol [ ],
    Sparse_long Front_parent [ ],
    Sparse_long Front_1strow [ ],
    Sparse_long Front_leftmostdesc [ ],
    Sparse_long Chain_start [ ],
    Sparse_long Chain_maxrows [ ],
    Sparse_long Chain_maxcols [ ],
    void *Symbolic
) ;

Sparse_long sparselu_get_numeric
(
    Sparse_long Lp [ ],
    Sparse_long Lj [ ],
    double Lx [ ],
    Sparse_long Up [ ],
    Sparse_long Ui [ ],
    double Ux [ ],
    Sparse_long P [ ],
    Sparse_long Q [ ],
    double Dx [ ],
    Sparse_long *do_recip,
    double Rs [ ],
    void *Numeric
) ;

Sparse_long sparselu_save_symbolic
(
    void *Symbolic,
    char *filename
) ;

Sparse_long sparselu_save_numeric
(
    void *Numeric,
    char *filename
) ;

Sparse_long sparselu_load_symbolic
(
    void **Symbolic,
    char *filename
) ;

Sparse_long sparselu_load_numeric
(
    void **Numeric,
    char *filename
) ;

#ifdef __cplusplus
}
#endif

#endif 
