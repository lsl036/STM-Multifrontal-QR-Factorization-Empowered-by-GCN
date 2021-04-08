#ifndef NCAMD

#include "Sparse_internal.h"
#include "camd.h"
#include "Sparse_camd.h"

#if (CAMD_VERSION < CAMD_VERSION_CODE (2,0))
#error "CAMD v2.0 or later is required"
#endif

/* ========================================================================== */
/* === SparseCore_camd ========================================================= */
/* ========================================================================== */

int SparseCore_camd
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    Int *Cmember,	/* size nrow.  see SparseCore_ccolamd.c for description.*/
    /* ---- output ---- */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    sparse_common *Common
)
{
    double Info [CAMD_INFO], Control2 [CAMD_CONTROL], *Control ;
    Int *Cp, *Len, *Nv, *Head, *Elen, *Degree, *Wi, *Next, *BucketSet,
	*Work3n, *p ;
    sparse_csc *C ;
    Int j, n, cnz ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    n = A->nrow ;

    /* s = 4*n */
    s = SparseCore_mult_size_t (n, 4, &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    RETURN_IF_NULL (Perm, FALSE) ;
    //RETURN_IF_XTYPE_INVALID (A, SparseCore_PATTERN, SparseCore_ZOMPLEX, FALSE) ;
    Common->status = SPARSE_OK ;
    if (n == 0)
    {
	/* nothing to do */
	Common->fl = 0 ;
	Common->lnz = 0 ;
	Common->anz = 0 ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    /* SparseCore_analyze has allocated Cmember at Common->Iwork + 5*n+uncol, and
     * CParent at Common->Iwork + 4*n+uncol, where uncol is 0 if A is symmetric
     * or A->ncol otherwise.  Thus, only the first 4n integers in Common->Iwork
     * can be used here. */

    SparseCore_allocate_work (n, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    p = Common->Iwork ;
    Degree = p ; p += n ;	/* size n */
    Elen   = p ; p += n ;	/* size n */
    Len    = p ; p += n ;	/* size n */
    Nv     = p ; p += n ;	/* size n */

    Work3n = SparseCore_malloc (n+1, 3*sizeof (Int), Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }
    p = Work3n ;
    Next = p ; p += n ;		/* size n */
    Wi   = p ; p += (n+1) ;	/* size n+1 */
    BucketSet = p ;		/* size n */

    Head = Common->Head ;	/* size n+1 */

    /* ---------------------------------------------------------------------- */
    /* construct the input matrix for CAMD */
    /* ---------------------------------------------------------------------- */

    if (A->stype == 0)
    {
	/* C = A*A' or A(:,f)*A(:,f)', add extra space of nnz(C)/2+n to C */
	C = SparseCore_aat (A, fset, fsize, -2, Common) ;
    }
    else
    {
	/* C = A+A', but use only the upper triangular part of A if A->stype = 1
	 * and only the lower part of A if A->stype = -1.  Add extra space of
	 * nnz(C)/2+n to C. */
	C = SparseCore_copy (A, 0, -2, Common) ;
    }

    if (Common->status < SPARSE_OK)
    {
	/* out of memory, fset invalid, or other error */
	SparseCore_free (n+1, 3*sizeof (Int), Work3n, Common) ;
	return (FALSE) ;
    }

    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
	Len [j] = Cp [j+1] - Cp [j] ;
    }

    /* C does not include the diagonal, and both upper and lower parts.
     * Common->anz includes the diagonal, and just the lower part of C */
    cnz = Cp [n] ;
    Common->anz = cnz / 2 + n ;

    /* ---------------------------------------------------------------------- */
    /* order C using CAMD */
    /* ---------------------------------------------------------------------- */

    /* get parameters */
    if (Common->current < 0 || Common->current >= SPARSE_MAXMETHODS)
    {
	/* use CAMD defaults */
	Control = NULL ;
    }
    else
    {
	Control = Control2 ;
	Control [CAMD_DENSE] = Common->method [Common->current].prune_dense ;
	Control [CAMD_AGGRESSIVE] = Common->method [Common->current].aggressive;
    }

#ifdef LONG
    /* DEBUG (camd_l_debug_init ("SparseCore_l_camd")) ; */
    camd_l2 (n, C->p,  C->i, Len, C->nzmax, cnz, Nv, Next, Perm, Head, Elen,
	    Degree, Wi, Control, Info, Cmember, BucketSet) ;
#else
    /* DEBUG (camd_debug_init ("SparseCore_camd")) ; */
    camd_2 (n, C->p,  C->i, Len, C->nzmax, cnz, Nv, Next, Perm, Head, Elen,
	    Degree, Wi, Control, Info, Cmember, BucketSet) ;
#endif

    /* LL' flop count.  Need to subtract n for LL' flop count.  Note that this
     * is a slight upper bound which is often exact (see CAMD/Source/camd_2.c
     * for details).  SparseCore_analyze computes an exact flop count and
     * fill-in. */
    Common->fl = Info [CAMD_NDIV] + 2 * Info [CAMD_NMULTSUBS_LDL] + n ;

    /* Info [CAMD_LNZ] excludes the diagonal */
    Common->lnz = n + Info [CAMD_LNZ] ;

    /* ---------------------------------------------------------------------- */
    /* free the CAMD workspace and clear the persistent workspace in Common */
    /* ---------------------------------------------------------------------- */

    // ASSERT (IMPLIES (Common->status == SPARSE_OK,
	// 	CHOLMOD(dump_perm) (Perm, n, n, "CAMD2 perm", Common))) ;
    SparseCore_free_sparse (&C, Common) ;
    for (j = 0 ; j <= n ; j++)
    {
	Head [j] = EMPTY ;
    }
    SparseCore_free (n+1, 3*sizeof (Int), Work3n, Common) ;
    return (TRUE) ;
}
#endif
