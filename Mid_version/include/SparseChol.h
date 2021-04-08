
#ifndef SPARSE_CHOL_H
#define SPARSE_CHOL_H

#include "SparseCore.h"
#include "tpsm.h"

#define SPARSE_A    0		/* 求解 Ax=b */
#define SPARSE_LDLt 1		/* 求解 LDL'x=b */
#define SPARSE_LD   2		/* 求解 LDx=b */
#define SPARSE_DLt  3		/* 求解 DL'x=b */
#define SPARSE_L    4		/* 求解 Lx=b */
#define SPARSE_Lt   5		/* 求解 L'x=b */
#define SPARSE_D    6		/* 求解 Dx=b */
#define SPARSE_P    7		/* 置换 x=Px */
#define SPARSE_Pt   8		/* 置换 x=P'x */

sparse_factor *SparseChol_analyze (sparse_csc *, sparse_common *) ;

sparse_factor *SparseChol_analyze_p (sparse_csc *, Sparse_long *,
    Sparse_long *, size_t, sparse_common *) ;

sparse_factor *SparseChol_analyze_p2 (int, sparse_csc *, Sparse_long *,
    Sparse_long *, size_t, sparse_common *) ;

int SparseChol_factorize (TPSM_t *, sparse_csc *, sparse_factor *, sparse_common *) ;

int SparseChol_factorize_p (TPSM_t *, sparse_csc *, double *, Sparse_long *,
    size_t, sparse_factor *, sparse_common *) ;

int SparseChol_etree (sparse_csc *, Sparse_long *, sparse_common *) ;

int SparseChol_rowcolcounts (sparse_csc *, Sparse_long *, size_t,
    Sparse_long *, Sparse_long *, Sparse_long *,
    Sparse_long *, Sparse_long *, Sparse_long *,
    sparse_common *) ;

int SparseChol_analyze_ordering (sparse_csc *, int, Sparse_long *,
    Sparse_long *, size_t, Sparse_long *, Sparse_long *,
    Sparse_long *, Sparse_long *, Sparse_long *,
    sparse_common *) ;

int SparseChol_amd (sparse_csc *, Sparse_long *, size_t,
    Sparse_long *, sparse_common *) ;

int SparseChol_colamd (sparse_csc *, Sparse_long *, size_t, int,
    Sparse_long *, sparse_common *) ;

int SparseChol_rowfac (sparse_csc *, sparse_csc *, double *, size_t,
    size_t, sparse_factor *, sparse_common *) ;

int SparseChol_rowfac_mask (sparse_csc *, sparse_csc *, double *, size_t,
    size_t, Sparse_long *, Sparse_long *, sparse_factor *,
    sparse_common *) ;

int SparseChol_rowfac_mask2 (sparse_csc *, sparse_csc *, double *,
    size_t, size_t, Sparse_long *, Sparse_long, Sparse_long *,
    sparse_factor *, sparse_common *) ;

int SparseChol_row_subtree (sparse_csc *, sparse_csc *, size_t,
    Sparse_long *, sparse_csc *, sparse_common *) ;

int SparseChol_lsolve_pattern (sparse_csc *B, sparse_factor *L, 
    sparse_csc *Yset, sparse_common *Common );

int SparseChol_row_lsubtree (sparse_csc *, Sparse_long *, size_t,
    size_t, sparse_factor *, sparse_csc *, sparse_common *) ;

int SparseChol_resymbol (sparse_csc *, Sparse_long *, size_t, int,
    sparse_factor *, sparse_common *) ;

int SparseChol_resymbol_noperm (sparse_csc *, Sparse_long *, size_t, int,
    sparse_factor *, sparse_common *) ;

double SparseChol_rcond (sparse_factor *, sparse_common *) ;

Sparse_long SparseChol_postorder (Sparse_long *, size_t,
    Sparse_long *, Sparse_long *, sparse_common *) ;

dense_array *SparseChol_solve (int sys, sparse_factor *L, 
        dense_array *B, sparse_common *Common);

int SparseChol_solve2( int sys, sparse_factor *L, dense_array *B,
    sparse_csc *Bset, dense_array **X_Handle, sparse_csc **Xset_Handle, 
    dense_array **Y_Handle, dense_array **E_Handle, sparse_common *Common);

int SparseChol_super_symbolic (sparse_csc *, sparse_csc *,
    Sparse_long *, sparse_factor *, sparse_common *) ;

int SparseChol_super_symbolic2 (int, sparse_csc *, sparse_csc *,
    Sparse_long *, sparse_factor *, sparse_common *) ;

int SparseChol_super_numeric (TPSM_t *, sparse_csc *, sparse_csc *, double *,
    sparse_factor *, sparse_common *) ;

int SparseChol_super_lsolve (sparse_factor *L, dense_array *X, 
        dense_array *E, sparse_common *Common );

int SparseChol_super_ltsolve (sparse_factor *L, dense_array *X,
        dense_array *E, sparse_common *Common );

#endif
