#ifndef NCAMD

#include "Sparse_internal.h"
#include "ccolamd.h"
#include "Sparse_camd.h"

#if (CCOLAMD_VERSION < CCOLAMD_VERSION_CODE (2,5))
#error "CCOLAMD v2.0 or later is required"
#endif

int SparseCore_csymamd
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to order */
    /* ---- output --- */
    Int *Cmember,	/* size nrow.  see cholmod_ccolamd.c for description */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    sparse_common *Common
)
{
    double knobs [CCOLAMD_KNOBS] ;
    Int *perm, *Head ;
    Int ok, i, nrow, stats [CCOLAMD_STATS] ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Perm, FALSE) ;
    //RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = SPARSE_OK ;

    if (A->nrow != A->ncol || !(A->packed))
    {
	ERROR (SPARSE_INVALID, "matrix must be square and packed") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    SparseCore_allocate_work (nrow, 0, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* order the matrix (does not affect A->p or A->i) */
    /* ---------------------------------------------------------------------- */

    perm = Common->Head ;	/* size nrow+1 (i/l/l) */

    /* get parameters */
#ifdef LONG
    ccolamd_l_set_defaults (knobs) ;
#else
    ccolamd_set_defaults (knobs) ;
#endif
    if (Common->current >= 0 && Common->current < SPARSE_MAXMETHODS)
    {
	/* get the knobs from the Common parameters */
	knobs [CCOLAMD_DENSE_ROW] =Common->method[Common->current].prune_dense ;
	knobs [CCOLAMD_AGGRESSIVE]=Common->method[Common->current].aggressive ;
    }
    {
#ifdef LONG
	csymamd_l (nrow, A->i, A->p, perm, knobs, stats,
                SuiteSparse_config.calloc_func,
                SuiteSparse_config.free_func,
                Cmember, A->stype) ;
#else
	csymamd (nrow, A->i, A->p, perm, knobs, stats,
                SuiteSparse_config.calloc_func,
                SuiteSparse_config.free_func,
                Cmember, A->stype) ;
#endif
	ok = stats [CCOLAMD_STATUS] ;
    }

    if (ok == CCOLAMD_ERROR_out_of_memory)
    {
	ERROR (SPARSE_OUT_OF_MEMORY, "out of memory") ; 
    }
    ok = (ok == CCOLAMD_OK || ok == CCOLAMD_OK_BUT_JUMBLED) ;

    /* ---------------------------------------------------------------------- */
    /* free the workspace and return result */
    /* ---------------------------------------------------------------------- */

    /* permutation returned in perm [0..n-1] */
    for (i = 0 ; i < nrow ; i++)
    {
	Perm [i] = perm [i] ;
    }

    /* clear Head workspace (used for perm, in csymamd): */
    Head = Common->Head ;
    for (i = 0 ; i <= nrow ; i++)
    {
	Head [i] = EMPTY ;
    }

    return (ok) ;
}
#endif