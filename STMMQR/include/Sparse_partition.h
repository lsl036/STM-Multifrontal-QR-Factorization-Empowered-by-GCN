
#ifndef CHOLMOD_PARTITION_H
#define CHOLMOD_PARTITION_H

#include "SparseCore.h"
// #include "Sparse_internal.h"
#include "Sparse_camd.h"
#define Int Sparse_long
/* -------------------------------------------------------------------------- */
/* cholmod_nested_dissection */
/* -------------------------------------------------------------------------- */

/* Order A, AA', or A(:,f)*A(:,f)' using CHOLMOD's nested dissection method
 * (METIS's node bisector applied recursively to compute the separator tree
 * and constraint sets, followed by CCOLAMD using the constraints).  Usually
 * finds better orderings than METIS_NodeND, but takes longer.
 */

Sparse_long SparseCore_nested_dissection	/* returns # of components */
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    Int *Perm,		/* size A->nrow, output permutation */
    Int *CParent,	/* size A->nrow.  On output, CParent [c] is the parent
			 * of component c, or EMPTY if c is a root, and where
			 * c is in the range 0 to # of components minus 1 */
    Int *Cmember,	/* size A->nrow.  Cmember [j] = c if node j of A is
			 * in component c */
    /* --------------- */
    sparse_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_metis */
/* -------------------------------------------------------------------------- */

/* Order A, AA', or A(:,f)*A(:,f)' using METIS_NodeND. */

int SparseCore_metis
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int postorder,	/* if TRUE, follow with etree or coletree postorder */
    /* ---- output --- */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    sparse_common *Common
) ;


/* -------------------------------------------------------------------------- */
/* cholmod_bisect */
/* -------------------------------------------------------------------------- */

/* Finds a node bisector of A, A*A', A(:,f)*A(:,f)'. */

Sparse_long SparseCore_bisect	/* returns # of nodes in separator */
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to bisect */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int compress,	/* if TRUE, compress the graph first */
    /* ---- output --- */
    Int *Partition,	/* size A->nrow.  Node i is in the left graph if
			 * Partition [i] = 0, the right graph if 1, and in the
			 * separator if 2. */
    /* --------------- */
    sparse_common *Common
) ;


/* -------------------------------------------------------------------------- */
/* cholmod_metis_bisector */
/* -------------------------------------------------------------------------- */

/* Find a set of nodes that bisects the graph of A or AA' (direct interface
 * to METIS_ComputeVertexSeperator). */

Sparse_long SparseCore_metis_bisector	/* returns separator size */
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to bisect */
    Int *Anw,		/* size A->nrow, node weights, can be NULL, */
                        /* which means the graph is unweighted. */ 
    Int *Aew,		/* size nz, edge weights (silently ignored). */
                        /* This option was available with METIS 4, but not */
                        /* in METIS 5.  This argument is now unused, but */
                        /* it remains for backward compatibilty, so as not */
                        /* to change the API for cholmod_metis_bisector. */
    /* ---- output --- */
    Int *Partition,	/* size A->nrow */
    /* --------------- */
    sparse_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_collapse_septree */
/* -------------------------------------------------------------------------- */

/* Collapse nodes in a separator tree. */

Sparse_long SparseCore_collapse_septree
(
    /* ---- input ---- */
    size_t n,		/* # of nodes in the graph */
    size_t ncomponents,	/* # of nodes in the separator tree (must be <= n) */
    double nd_oksep,    /* collapse if #sep >= nd_oksep * #nodes in subtree */
    size_t nd_small,    /* collapse if #nodes in subtree < nd_small */
    /* ---- in/out --- */
    Int *CParent,	/* size ncomponents; from cholmod_nested_dissection */
    Int *Cmember,	/* size n; from cholmod_nested_dissection */
    /* --------------- */
    sparse_common *Common
) ;

// SuiteSparse_long cholmod_l_collapse_septree (size_t, size_t, double, size_t,
//     SuiteSparse_long *, SuiteSparse_long *, cholmod_common *) ;

#endif
