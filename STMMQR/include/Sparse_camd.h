#ifndef SPARSE_CAMD_H
#define SPARSE_CAMD_H

#include "SparseCore.h"
#define Int Sparse_long
/* -------------------------------------------------------------------------- */
/* SparseCore_ccolamd */
/* -------------------------------------------------------------------------- */

/* Order AA' or A(:,f)*A(:,f)' using CCOLAMD. */

int SparseCore_ccolamd
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    Int *Cmember,	/* size A->nrow.  Cmember [i] = c if row i is in the
			 * constraint set c.  c must be >= 0.  The # of
			 * constraint sets is max (Cmember) + 1.  If Cmember is
			 * NULL, then it is interpretted as Cmember [i] = 0 for
			 * all i */
    /* ---- output --- */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    sparse_common *Common
) ;

// int SparseCore_l_ccolamd (sparse_csc *, SuiteSparse_long *, size_t,
//     SuiteSparse_long *, SuiteSparse_long *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_csymamd */
/* -------------------------------------------------------------------------- */

/* Order A using CSYMAMD. */
// 在 SparseCore_csymamd.c 中
int SparseCore_csymamd
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to order */
    /* ---- output --- */
    Int *Cmember,	/* size nrow.  see SparseCore_ccolamd above */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    sparse_common *Common
) ;

// int SparseCore_l_csymamd (sparse_csc *, SuiteSparse_long *,
//     SuiteSparse_long *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_camd */
/* -------------------------------------------------------------------------- */

/* Order A using CAMD. */
// 在 SparseCore_camd.c 中
int SparseCore_camd
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    Int *Cmember,	/* size nrow.  see SparseCore_ccolamd above */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    sparse_common *Common
) ;

// int SparseCore_l_camd (sparse_csc *, SuiteSparse_long *, size_t,
//     SuiteSparse_long *, SuiteSparse_long *, sparse_common *) ;

#endif