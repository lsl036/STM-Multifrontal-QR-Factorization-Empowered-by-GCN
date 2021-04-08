/* ========================================================================== */
/* === colamd/symamd prototypes and definitions ============================= */
/* ========================================================================== */

#ifndef COLAMD_H
#define COLAMD_H

/* make it easy for C++ programs to include COLAMD */
#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include <stdlib.h>

/* ========================================================================== */
/* === COLAMD version ======================================================= */
/* ========================================================================== */

/* COLAMD Version 2.4 and later will include the following definitions.
 * As an example, to test if the version you are using is 2.4 or later:
 *
 * #ifdef COLAMD_VERSION
 *	if (COLAMD_VERSION >= COLAMD_VERSION_CODE (2,4)) ...
 * #endif
 *
 * This also works during compile-time:
 *
 *  #if defined(COLAMD_VERSION) && (COLAMD_VERSION >= COLAMD_VERSION_CODE (2,4))
 *    printf ("This is version 2.4 or later\n") ;
 *  #else
 *    printf ("This is an early version\n") ;
 *  #endif
 *
 * Versions 2.3 and earlier of COLAMD do not include a #define'd version number.
 */

#define COLAMD_DATE "May 4, 2016"
#define COLAMD_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define COLAMD_MAIN_VERSION 2
#define COLAMD_SUB_VERSION 9
#define COLAMD_SUBSUB_VERSION 6
#define COLAMD_VERSION \
	COLAMD_VERSION_CODE(COLAMD_MAIN_VERSION,COLAMD_SUB_VERSION)

/* ========================================================================== */
/* === Knob and statistics definitions ====================================== */
/* ========================================================================== */

/* size of the knobs [ ] array.  Only knobs [0..1] are currently used. */
#define COLAMD_KNOBS 20

/* number of output statistics.  Only stats [0..6] are currently used. */
#define COLAMD_STATS 20

/* knobs [0] and stats [0]: dense row knob and output statistic. */
#define COLAMD_DENSE_ROW 0

/* knobs [1] and stats [1]: dense column knob and output statistic. */
#define COLAMD_DENSE_COL 1

/* knobs [2]: aggressive absorption */
#define COLAMD_AGGRESSIVE 2

/* stats [2]: memory defragmentation count output statistic */
#define COLAMD_DEFRAG_COUNT 2

/* stats [3]: colamd status:  zero OK, > 0 warning or notice, < 0 error */
#define COLAMD_STATUS 3

/* stats [4..6]: error info, or info on jumbled columns */ 
#define COLAMD_INFO1 4
#define COLAMD_INFO2 5
#define COLAMD_INFO3 6

/* error codes returned in stats [3]: */
#define COLAMD_OK				(0)
#define COLAMD_OK_BUT_JUMBLED			(1)
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


/* ========================================================================== */
/* === Prototypes of user-callable routines ================================= */
/* ========================================================================== */

#include "SparseBase_config.h"

size_t colamd_recommended	/* returns recommended value of Alen, */
				/* or 0 if input arguments are erroneous */
(
    Sparse_long nnz,       /* nonzeros in A */
    Sparse_long n_row,     /* number of rows in A */
    Sparse_long n_col      /* number of columns in A */
) ;

void colamd_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [COLAMD_KNOBS]	/* parameter settings for colamd */
) ;

Sparse_long colamd       /* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    Sparse_long n_row,     /* number of rows in A */
    Sparse_long n_col,     /* number of columns in A */
    Sparse_long Alen,      /* size of the array A */
    Sparse_long A [],      /* row indices of A, of size Alen */
    Sparse_long p [],      /* column pointers of A, of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameter settings for colamd */
    Sparse_long stats [COLAMD_STATS]   /* colamd output statistics
                                             * and error codes */
) ;

Sparse_long symamd               /* return (1) if OK, (0) otherwise */
(
    Sparse_long n,                 /* number of rows and columns of A */
    Sparse_long A [],              /* row indices of A */
    Sparse_long p [],              /* column pointers of A */
    Sparse_long perm [],           /* output permutation, size n_col+1 */
    double knobs [COLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    Sparse_long stats [COLAMD_STATS],  /* output stats and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
) ;

void colamd_report
(
    Sparse_long stats [COLAMD_STATS]
) ;

void symamd_report
(
    Sparse_long stats [COLAMD_STATS]
) ;

#ifdef __cplusplus
}
#endif

#endif /* COLAMD_H */
