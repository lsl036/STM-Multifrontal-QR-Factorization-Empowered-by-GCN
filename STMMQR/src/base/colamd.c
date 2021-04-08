/* ========================================================================== */
/* === colamd/symamd - 一种稀疏矩阵列排序算法 ================================= */
/* ========================================================================== */

/* COLAMD / SYMAMD

    colamd:  一种近似最小度列排序算法，用于对称或非对称矩阵的LU分解，QR分解，
	最小二乘，线性规划问题的内点方法，以及其他相关问题。

    symamd:  对称矩阵Cholesky分解的一种近似最小度排序算法。

    用途:

	Colamd计算一个排列Q，这样cholesky分解(AQ)'(AQ)与A'A相比具有更少的填充并且只
	需要更少的浮点运算次数。这也为稀疏局部旋转方法提供了一个良好的排序，P(AQ) = LU，
	其中Q是在数值分解之前计算的，而P是在数值分解过程通过常规行交换的部分旋转计算的。
	Colamd是SuperLU中使用的列排序方法，SuperLU是ScaLAPACK库的一部分。它也可以在
	MATLAB版本6中作为内置函数使用，可以从MathWorks(http://www.mathworks.com)获得。
	这个例程可以在MATLAB中代替colmmd。

	Symamd计算对称矩阵A的排列P使得cholesky分解PAP'比A具有更少的填充并且只需要更少的
	浮点运算次数。 Symamd构造一个矩阵M, 使得M'M与A有相同的非零模式,然后使用colmmd对
	M的列进行排序。然后将M的列顺序返回为A的行和列顺序P。

*/

/* ========================================================================== */
/* === 代码调试 ============================================================== */
/* ========================================================================== */

/*
	当启用调试时，代码会变得异常的慢。要控制调试输出的级别，可以将环境变量D设置为0 (little)、
	1 (some)、2、3或4 (lots)。调试时，您应该会在标准输出中看到以下消息:

    	colamd: debug version, D = 1 (THIS WILL BE SLOW!)

    symamd也提供了类似的消息。如果没有，那么调试也没有被启用。
*/

/* ========================================================================== */
/* === 头文件 =============================================================== */
/* ========================================================================== */

#include "colamd.h"
#include <limits.h>
#include <math.h>

#if !defined (NPRINT)
#include <stdio.h>
#endif

#ifndef NULL
#define NULL ((void *) 0)
#endif

/* ========================================================================== */
/* === int or Sparse_long ============================================== */
/* ========================================================================== */

#ifdef DLONG

#define Int Sparse_long
#define ID  Sparse_long_id
#define Int_MAX Sparse_long_max

#define COLAMD_recommended colamd_recommended
#define COLAMD_set_defaults colamd_set_defaults
#define COLAMD_MAIN colamd
#define SYMAMD_MAIN symamd
#define COLAMD_report colamd_report
#define SYMAMD_report symamd_report

#endif

/* ========================================================================== */
/* === 行列结构 ============================================================== */
/* ========================================================================== */

/* 不被直接使用，仅用于 colamd_recommended. */

typedef struct Colamd_Col_struct
{
    Int start ;		/* index for A of first row in this column, or DEAD */
			/* if column is dead */
    Int length ;	/* number of rows in this column */
    union
    {
	Int thickness ;	/* number of original columns represented by this */
			/* col, if the column is alive */
	Int parent ;	/* parent in parent tree super-column structure, if */
			/* the column is dead */
    } shared1 ;
    union
    {
	Int score ;	/* the score used to maintain heap, if col is alive */
	Int order ;	/* pivot ordering of this column, if col is dead */
    } shared2 ;
    union
    {
	Int headhash ;	/* head of a hash bucket, if col is at the head of */
			/* a degree list */
	Int hash ;	/* hash value, if col is not in a degree list */
	Int prev ;	/* previous column in degree list, if col is in a */
			/* degree list (but not at the head of a degree list) */
    } shared3 ;
    union
    {
	Int degree_next ;	/* next column, if col is in a degree list */
	Int hash_next ;		/* next column, if col is in a hash list */
    } shared4 ;

} Colamd_Col ;

typedef struct Colamd_Row_struct
{
    Int start ;		/* index for A of first col in this row */
    Int length ;	/* number of principal columns in this row */
    union
    {
	Int degree ;	/* number of principal & non-principal columns in row */
	Int p ;		/* used as a row pointer in init_rows_cols () */
    } shared1 ;
    union
    {
	Int mark ;	/* for computing set differences and marking dead rows*/
	Int first_column ;/* first column in row (used in garbage collection) */
    } shared2 ;

} Colamd_Row ;

/* ========================================================================== */
/* === Definitions ========================================================== */
/* ========================================================================== */

/* Routines are either PUBLIC (user-callable) or PRIVATE (not user-callable) */
#define PUBLIC
#define PRIVATE static

#define DENSE_DEGREE(alpha,n) \
    ((Int) MAX (16.0, (alpha) * sqrt ((double) (n))))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define ONES_COMPLEMENT(r) (-(r)-1)

/* -------------------------------------------------------------------------- */
/* Change for version 2.1:  define TRUE and FALSE only if not yet defined */  
/* -------------------------------------------------------------------------- */

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

/* -------------------------------------------------------------------------- */

#define EMPTY	(-1)

/* Row and column status */
#define ALIVE	(0)
#define DEAD	(-1)

/* Column status */
#define DEAD_PRINCIPAL		(-1)
#define DEAD_NON_PRINCIPAL	(-2)

/* Macros for row and column status update and checking. */
#define ROW_IS_DEAD(r)			ROW_IS_MARKED_DEAD (Row[r].shared2.mark)
#define ROW_IS_MARKED_DEAD(row_mark)	(row_mark < ALIVE)
#define ROW_IS_ALIVE(r)			(Row [r].shared2.mark >= ALIVE)
#define COL_IS_DEAD(c)			(Col [c].start < ALIVE)
#define COL_IS_ALIVE(c)			(Col [c].start >= ALIVE)
#define COL_IS_DEAD_PRINCIPAL(c)	(Col [c].start == DEAD_PRINCIPAL)
#define KILL_ROW(r)			{ Row [r].shared2.mark = DEAD ; }
#define KILL_PRINCIPAL_COL(c)		{ Col [c].start = DEAD_PRINCIPAL ; }
#define KILL_NON_PRINCIPAL_COL(c)	{ Col [c].start = DEAD_NON_PRINCIPAL ; }

/* ========================================================================== */
/* === Colamd reporting mechanism =========================================== */
/* ========================================================================== */

#if defined (MATLAB_MEX_FILE) || defined (MATHWORKS)
/* In MATLAB, matrices are 1-based to the user, but 0-based internally */
#define INDEX(i) ((i)+1)
#else
/* In C, matrices are 0-based and indices are reported as such in *_report */
#define INDEX(i) (i)
#endif

/* ========================================================================== */
/* === Prototypes of PRIVATE routines ======================================= */
/* ========================================================================== */

PRIVATE Int init_rows_cols
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int p [],
    Int stats [COLAMD_STATS]
) ;

PRIVATE void init_scoring
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    double knobs [COLAMD_KNOBS],
    Int *p_n_row2,
    Int *p_n_col2,
    Int *p_max_deg
) ;

PRIVATE Int find_ordering
(
    Int n_row,
    Int n_col,
    Int Alen,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int n_col2,
    Int max_deg,
    Int pfree,
    Int aggressive
) ;

PRIVATE void order_children
(
    Int n_col,
    Colamd_Col Col [],
    Int p []
) ;

PRIVATE void detect_super_cols
(
    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int row_start,
    Int row_length
) ;

PRIVATE Int garbage_collection
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int *pfree
) ;

PRIVATE Int clear_mark
(
    Int tag_mark,
    Int max_mark,
    Int n_row,
    Colamd_Row Row []
) ;

PRIVATE void print_report
(
    char *method,
    Int stats [COLAMD_STATS]
) ;

/* === No debugging ========================================================= */

#define DEBUG0(params) ;
#define DEBUG1(params) ;
#define DEBUG2(params) ;
#define DEBUG3(params) ;
#define DEBUG4(params) ;

#define ASSERT(expression)

/* ========================================================================== */
/* === USER-CALLABLE ROUTINES: ============================================== */
/* ========================================================================== */

/* ========================================================================== */
/* === colamd_recommended =================================================== */
/* ========================================================================== */

/*
    The colamd_recommended routine returns the suggested size for Alen.  This
    value has been determined to provide good balance between the number of
    garbage collections and the memory requirements for colamd.  If any
    argument is negative, or if integer overflow occurs, a 0 is returned as an
    error condition.  2*nnz space is required for the row and column
    indices of the matrix. COLAMD_C (n_col) + COLAMD_R (n_row) space is
    required for the Col and Row arrays, respectively, which are internal to
    colamd (roughly 6*n_col + 4*n_row).  An additional n_col space is the
    minimal amount of "elbow room", and nnz/5 more space is recommended for
    run time efficiency.

    Alen is approximately 2.2*nnz + 7*n_col + 4*n_row + 10.

    This function is not needed when using symamd.
*/

/* add two values of type size_t, and check for integer overflow */
static size_t t_add (size_t a, size_t b, int *ok)
{
    (*ok) = (*ok) && ((a + b) >= MAX (a,b)) ;
    return ((*ok) ? (a + b) : 0) ;
}

/* compute a*k where k is a small integer, and check for integer overflow */
static size_t t_mult (size_t a, size_t k, int *ok)
{
    size_t i, s = 0 ;
    for (i = 0 ; i < k ; i++)
    {
	s = t_add (s, a, ok) ;
    }
    return (s) ;
}

/* size of the Col and Row structures */
#define COLAMD_C(n_col,ok) \
    ((t_mult (t_add (n_col, 1, ok), sizeof (Colamd_Col), ok) / sizeof (Int)))

#define COLAMD_R(n_row,ok) \
    ((t_mult (t_add (n_row, 1, ok), sizeof (Colamd_Row), ok) / sizeof (Int)))


PUBLIC size_t COLAMD_recommended	/* returns recommended value of Alen. */
(
    /* === Parameters ======================================================= */

    Int nnz,			/* number of nonzeros in A */
    Int n_row,			/* number of rows in A */
    Int n_col			/* number of columns in A */
)
{
    size_t s, c, r ;
    int ok = TRUE ;
    if (nnz < 0 || n_row < 0 || n_col < 0)
    {
	return (0) ;
    }
    s = t_mult (nnz, 2, &ok) ;	    /* 2*nnz */
    c = COLAMD_C (n_col, &ok) ;	    /* size of column structures */
    r = COLAMD_R (n_row, &ok) ;	    /* size of row structures */
    s = t_add (s, c, &ok) ;
    s = t_add (s, r, &ok) ;
    s = t_add (s, n_col, &ok) ;	    /* elbow room */
    s = t_add (s, nnz/5, &ok) ;	    /* elbow room */
    ok = ok && (s < Int_MAX) ;
    return (ok ? s : 0) ;
}


/* ========================================================================== */
/* === colamd_set_defaults ================================================== */
/* ========================================================================== */

/*
    The colamd_set_defaults routine sets the default values of the user-
    controllable parameters for colamd and symamd:

	Colamd: rows with more than max (16, knobs [0] * sqrt (n_col))
	entries are removed prior to ordering.  Columns with more than
	max (16, knobs [1] * sqrt (MIN (n_row,n_col))) entries are removed
	prior to ordering, and placed last in the output column ordering. 

	Symamd: Rows and columns with more than max (16, knobs [0] * sqrt (n))
	entries are removed prior to ordering, and placed last in the
	output ordering.

	knobs [0]	dense row control

	knobs [1]	dense column control

	knobs [2]	if nonzero, do aggresive absorption

	knobs [3..19]	unused, but future versions might use this

*/

PUBLIC void COLAMD_set_defaults
(
    /* === Parameters ======================================================= */

    double knobs [COLAMD_KNOBS]		/* knob array */
)
{
    /* === Local variables ================================================== */

    Int i ;

    if (!knobs)
    {
	return ;			/* no knobs to initialize */
    }
    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
	knobs [i] = 0 ;
    }
    knobs [COLAMD_DENSE_ROW] = 10 ;
    knobs [COLAMD_DENSE_COL] = 10 ;
    knobs [COLAMD_AGGRESSIVE] = TRUE ;	/* default: do aggressive absorption*/
}


/* ========================================================================== */
/* === symamd =============================================================== */
/* ========================================================================== */

PUBLIC Int SYMAMD_MAIN			/* return TRUE if OK, FALSE otherwise */
(
    /* === Parameters ======================================================= */

    Int n,				/* number of rows and columns of A */
    Int A [],				/* row indices of A */
    Int p [],				/* column pointers of A */
    Int perm [],			/* output permutation, size n+1 */
    double knobs [COLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    Int stats [COLAMD_STATS],		/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
)
{
    /* === Local variables ================================================== */

    Int *count ;		/* length of each column of M, and col pointer*/
    Int *mark ;			/* mark array for finding duplicate entries */
    Int *M ;			/* row indices of matrix M */
    size_t Mlen ;		/* length of M */
    Int n_row ;			/* number of rows in M */
    Int nnz ;			/* number of entries in A */
    Int i ;			/* row index of A */
    Int j ;			/* column index of A */
    Int k ;			/* row index of M */ 
    Int mnz ;			/* number of nonzeros in M */
    Int pp ;			/* index into a column of A */
    Int last_row ;		/* last row seen in the current column */
    Int length ;		/* number of nonzeros in a column */

    double cknobs [COLAMD_KNOBS] ;		/* knobs for colamd */
    double default_knobs [COLAMD_KNOBS] ;	/* default knobs for colamd */

    /* === Check the input arguments ======================================== */

    if (!stats)
    {
	DEBUG0 (("symamd: stats not present\n")) ;
	return (FALSE) ;
    }
    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

    if (!A)
    {
    	stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
	DEBUG0 (("symamd: A not present\n")) ;
	return (FALSE) ;
    }

    if (!p)		/* p is not present */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
	DEBUG0 (("symamd: p not present\n")) ;
    	return (FALSE) ;
    }

    if (n < 0)		/* n must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
	stats [COLAMD_INFO1] = n ;
	DEBUG0 (("symamd: n negative %d\n", n)) ;
    	return (FALSE) ;
    }

    nnz = p [n] ;
    if (nnz < 0)	/* nnz must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
	stats [COLAMD_INFO1] = nnz ;
	DEBUG0 (("symamd: number of entries negative %d\n", nnz)) ;
	return (FALSE) ;
    }

    if (p [0] != 0)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero ;
	stats [COLAMD_INFO1] = p [0] ;
	DEBUG0 (("symamd: p[0] not zero %d\n", p [0])) ;
	return (FALSE) ;
    }

    /* === If no knobs, set default knobs =================================== */

    if (!knobs)
    {
	COLAMD_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }

    /* === Allocate count and mark ========================================== */

    count = (Int *) ((*allocate) (n+1, sizeof (Int))) ;
    if (!count)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	DEBUG0 (("symamd: allocate count (size %d) failed\n", n+1)) ;
	return (FALSE) ;
    }

    mark = (Int *) ((*allocate) (n+1, sizeof (Int))) ;
    if (!mark)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	(*release) ((void *) count) ;
	DEBUG0 (("symamd: allocate mark (size %d) failed\n", n+1)) ;
	return (FALSE) ;
    }

    /* === Compute column counts of M, check if A is valid ================== */

    stats [COLAMD_INFO3] = 0 ;  /* number of duplicate or unsorted row indices*/

    for (i = 0 ; i < n ; i++)
    {
    	mark [i] = -1 ;
    }

    for (j = 0 ; j < n ; j++)
    {
	last_row = -1 ;

	length = p [j+1] - p [j] ;
	if (length < 0)
	{
	    /* column pointers must be non-decreasing */
	    stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
	    stats [COLAMD_INFO1] = j ;
	    stats [COLAMD_INFO2] = length ;
	    (*release) ((void *) count) ;
	    (*release) ((void *) mark) ;
	    DEBUG0 (("symamd: col %d negative length %d\n", j, length)) ;
	    return (FALSE) ;
	}

	for (pp = p [j] ; pp < p [j+1] ; pp++)
	{
	    i = A [pp] ;
	    if (i < 0 || i >= n)
	    {
		/* row index i, in column j, is out of bounds */
		stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
		stats [COLAMD_INFO1] = j ;
		stats [COLAMD_INFO2] = i ;
		stats [COLAMD_INFO3] = n ;
		(*release) ((void *) count) ;
		(*release) ((void *) mark) ;
		DEBUG0 (("symamd: row %d col %d out of bounds\n", i, j)) ;
		return (FALSE) ;
	    }

	    if (i <= last_row || mark [i] == j)
	    {
		/* row index is unsorted or repeated (or both), thus col */
		/* is jumbled.  This is a notice, not an error condition. */
		stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
		stats [COLAMD_INFO1] = j ;
		stats [COLAMD_INFO2] = i ;
		(stats [COLAMD_INFO3]) ++ ;
		DEBUG1 (("symamd: row %d col %d unsorted/duplicate\n", i, j)) ;
	    }

	    if (i > j && mark [i] != j)
	    {
		/* row k of M will contain column indices i and j */
		count [i]++ ;
		count [j]++ ;
	    }

	    /* mark the row as having been seen in this column */
	    mark [i] = j ;

	    last_row = i ;
	}
    }

    /* v2.4: removed free(mark) */

    /* === Compute column pointers of M ===================================== */

    /* use output permutation, perm, for column pointers of M */
    perm [0] = 0 ;
    for (j = 1 ; j <= n ; j++)
    {
	perm [j] = perm [j-1] + count [j-1] ;
    }
    for (j = 0 ; j < n ; j++)
    {
	count [j] = perm [j] ;
    }

    /* === Construct M ====================================================== */

    mnz = perm [n] ;
    n_row = mnz / 2 ;
    Mlen = COLAMD_recommended (mnz, n_row, n) ;
    M = (Int *) ((*allocate) (Mlen, sizeof (Int))) ;
    DEBUG0 (("symamd: M is %d-by-%d with %d entries, Mlen = %g\n",
    	n_row, n, mnz, (double) Mlen)) ;

    if (!M)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	(*release) ((void *) count) ;
	(*release) ((void *) mark) ;
	DEBUG0 (("symamd: allocate M (size %g) failed\n", (double) Mlen)) ;
	return (FALSE) ;
    }

    k = 0 ;

    if (stats [COLAMD_STATUS] == COLAMD_OK)
    {
	/* Matrix is OK */
	for (j = 0 ; j < n ; j++)
	{
	    ASSERT (p [j+1] - p [j] >= 0) ;
	    for (pp = p [j] ; pp < p [j+1] ; pp++)
	    {
		i = A [pp] ;
		ASSERT (i >= 0 && i < n) ;
		if (i > j)
		{
		    /* row k of M contains column indices i and j */
		    M [count [i]++] = k ;
		    M [count [j]++] = k ;
		    k++ ;
		}
	    }
	}
    }
    else
    {
	/* Matrix is jumbled.  Do not add duplicates to M.  Unsorted cols OK. */
	DEBUG0 (("symamd: Duplicates in A.\n")) ;
	for (i = 0 ; i < n ; i++)
	{
	    mark [i] = -1 ;
	}
	for (j = 0 ; j < n ; j++)
	{
	    ASSERT (p [j+1] - p [j] >= 0) ;
	    for (pp = p [j] ; pp < p [j+1] ; pp++)
	    {
		i = A [pp] ;
		ASSERT (i >= 0 && i < n) ;
		if (i > j && mark [i] != j)
		{
		    /* row k of M contains column indices i and j */
		    M [count [i]++] = k ;
		    M [count [j]++] = k ;
		    k++ ;
		    mark [i] = j ;
		}
	    }
	}
	/* v2.4: free(mark) moved below */
    }

    /* count and mark no longer needed */
    (*release) ((void *) count) ;
    (*release) ((void *) mark) ;	/* v2.4: free (mark) moved here */
    ASSERT (k == n_row) ;

    /* === Adjust the knobs for M =========================================== */

    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
	cknobs [i] = knobs [i] ;
    }

    /* there are no dense rows in M */
    cknobs [COLAMD_DENSE_ROW] = -1 ;
    cknobs [COLAMD_DENSE_COL] = knobs [COLAMD_DENSE_ROW] ;

    /* === Order the columns of M =========================================== */

    /* v2.4: colamd cannot fail here, so the error check is removed */
    (void) COLAMD_MAIN (n_row, n, (Int) Mlen, M, perm, cknobs, stats) ;

    /* Note that the output permutation is now in perm */

    /* === get the statistics for symamd from colamd ======================== */

    /* a dense column in colamd means a dense row and col in symamd */
    stats [COLAMD_DENSE_ROW] = stats [COLAMD_DENSE_COL] ;

    /* === Free M =========================================================== */

    (*release) ((void *) M) ;
    DEBUG0 (("symamd: done.\n")) ;
    return (TRUE) ;

}

/* ========================================================================== */
/* === colamd =============================================================== */
/* ========================================================================== */

/*
    The colamd routine computes a column ordering Q of a sparse matrix
    A such that the LU factorization P(AQ) = LU remains sparse, where P is
    selected via partial pivoting.   The routine can also be viewed as
    providing a permutation Q such that the Cholesky factorization
    (AQ)'(AQ) = LL' remains sparse.
*/

PUBLIC Int COLAMD_MAIN		/* returns TRUE if successful, FALSE otherwise*/
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows in A */
    Int n_col,			/* number of columns in A */
    Int Alen,			/* length of A */
    Int A [],			/* row indices of A */
    Int p [],			/* pointers to columns in A */
    double knobs [COLAMD_KNOBS],/* parameters (uses defaults if NULL) */
    Int stats [COLAMD_STATS]	/* output statistics and error codes */
)
{
    /* === Local variables ================================================== */

    Int i ;			/* loop index */
    Int nnz ;			/* nonzeros in A */
    size_t Row_size ;		/* size of Row [], in integers */
    size_t Col_size ;		/* size of Col [], in integers */
    size_t need ;		/* minimum required length of A */
    Colamd_Row *Row ;		/* pointer into A of Row [0..n_row] array */
    Colamd_Col *Col ;		/* pointer into A of Col [0..n_col] array */
    Int n_col2 ;		/* number of non-dense, non-empty columns */
    Int n_row2 ;		/* number of non-dense, non-empty rows */
    Int ngarbage ;		/* number of garbage collections performed */
    Int max_deg ;		/* maximum row degree */
    double default_knobs [COLAMD_KNOBS] ;	/* default knobs array */
    Int aggressive ;		/* do aggressive absorption */
    int ok ;

    /* === Check the input arguments ======================================== */

    if (!stats)
    {
	DEBUG0 (("colamd: stats not present\n")) ;
	return (FALSE) ;
    }
    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

    if (!A)		/* A is not present */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
	DEBUG0 (("colamd: A not present\n")) ;
	return (FALSE) ;
    }

    if (!p)		/* p is not present */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
	DEBUG0 (("colamd: p not present\n")) ;
    	return (FALSE) ;
    }

    if (n_row < 0)	/* n_row must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_nrow_negative ;
	stats [COLAMD_INFO1] = n_row ;
	DEBUG0 (("colamd: nrow negative %d\n", n_row)) ;
    	return (FALSE) ;
    }

    if (n_col < 0)	/* n_col must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
	stats [COLAMD_INFO1] = n_col ;
	DEBUG0 (("colamd: ncol negative %d\n", n_col)) ;
    	return (FALSE) ;
    }

    nnz = p [n_col] ;
    if (nnz < 0)	/* nnz must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
	stats [COLAMD_INFO1] = nnz ;
	DEBUG0 (("colamd: number of entries negative %d\n", nnz)) ;
	return (FALSE) ;
    }

    if (p [0] != 0)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero	;
	stats [COLAMD_INFO1] = p [0] ;
	DEBUG0 (("colamd: p[0] not zero %d\n", p [0])) ;
	return (FALSE) ;
    }

    /* === If no knobs, set default knobs =================================== */

    if (!knobs)
    {
	COLAMD_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }

    aggressive = (knobs [COLAMD_AGGRESSIVE] != FALSE) ;

    /* === Allocate the Row and Col arrays from array A ===================== */

    ok = TRUE ;
    Col_size = COLAMD_C (n_col, &ok) ;	    /* size of Col array of structs */
    Row_size = COLAMD_R (n_row, &ok) ;	    /* size of Row array of structs */

    /* need = 2*nnz + n_col + Col_size + Row_size ; */
    need = t_mult (nnz, 2, &ok) ;
    need = t_add (need, n_col, &ok) ;
    need = t_add (need, Col_size, &ok) ;
    need = t_add (need, Row_size, &ok) ;

    if (!ok || need > (size_t) Alen || need > Int_MAX)
    {
	/* not enough space in array A to perform the ordering */
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_too_small ;
	stats [COLAMD_INFO1] = need ;
	stats [COLAMD_INFO2] = Alen ;
	DEBUG0 (("colamd: Need Alen >= %d, given only Alen = %d\n", need,Alen));
	return (FALSE) ;
    }

    Alen -= Col_size + Row_size ;
    Col = (Colamd_Col *) &A [Alen] ;
    Row = (Colamd_Row *) &A [Alen + Col_size] ;

    /* === Construct the row and column data structures ===================== */

    if (!init_rows_cols (n_row, n_col, Row, Col, A, p, stats))
    {
	/* input matrix is invalid */
	DEBUG0 (("colamd: Matrix invalid\n")) ;
	return (FALSE) ;
    }

    /* === Initialize scores, kill dense rows/columns ======================= */

    init_scoring (n_row, n_col, Row, Col, A, p, knobs,
	&n_row2, &n_col2, &max_deg) ;

    /* === Order the supercolumns =========================================== */

    ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
	n_col2, max_deg, 2*nnz, aggressive) ;

    /* === Order the non-principal columns ================================== */

    order_children (n_col, Col, p) ;

    /* === Return statistics in stats ======================================= */

    stats [COLAMD_DENSE_ROW] = n_row - n_row2 ;
    stats [COLAMD_DENSE_COL] = n_col - n_col2 ;
    stats [COLAMD_DEFRAG_COUNT] = ngarbage ;
    DEBUG0 (("colamd: done.\n")) ; 
    return (TRUE) ;
}


/* ========================================================================== */
/* === colamd_report ======================================================== */
/* ========================================================================== */

PUBLIC void COLAMD_report
(
    Int stats [COLAMD_STATS]
)
{
    print_report ("colamd", stats) ;
}


/* ========================================================================== */
/* === symamd_report ======================================================== */
/* ========================================================================== */

PUBLIC void SYMAMD_report
(
    Int stats [COLAMD_STATS]
)
{
    print_report ("symamd", stats) ;
}



/* ========================================================================== */
/* === NON-USER-CALLABLE ROUTINES: ========================================== */
/* ========================================================================== */

/* There are no user-callable routines beyond this point in the file */


/* ========================================================================== */
/* === init_rows_cols ======================================================= */
/* ========================================================================== */

/*
    Takes the column form of the matrix in A and creates the row form of the
    matrix.  Also, row and column attributes are stored in the Col and Row
    structs.  If the columns are un-sorted or contain duplicate row indices,
    this routine will also sort and remove duplicate row indices from the
    column form of the matrix.  Returns FALSE if the matrix is invalid,
    TRUE otherwise.  Not user-callable.
*/

PRIVATE Int init_rows_cols	/* returns TRUE if OK, or FALSE otherwise */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* row indices of A, of size Alen */
    Int p [],			/* pointers to columns in A, of size n_col+1 */
    Int stats [COLAMD_STATS]	/* colamd statistics */ 
)
{
    /* === Local variables ================================================== */

    Int col ;			/* a column index */
    Int row ;			/* a row index */
    Int *cp ;			/* a column pointer */
    Int *cp_end ;		/* a pointer to the end of a column */
    Int *rp ;			/* a row pointer */
    Int *rp_end ;		/* a pointer to the end of a row */
    Int last_row ;		/* previous row */

    /* === Initialize columns, and check column pointers ==================== */

    for (col = 0 ; col < n_col ; col++)
    {
	Col [col].start = p [col] ;
	Col [col].length = p [col+1] - p [col] ;

	if (Col [col].length < 0)
	{
	    /* column pointers must be non-decreasing */
	    stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
	    stats [COLAMD_INFO1] = col ;
	    stats [COLAMD_INFO2] = Col [col].length ;
	    DEBUG0 (("colamd: col %d length %d < 0\n", col, Col [col].length)) ;
	    return (FALSE) ;
	}

	Col [col].shared1.thickness = 1 ;
	Col [col].shared2.score = 0 ;
	Col [col].shared3.prev = EMPTY ;
	Col [col].shared4.degree_next = EMPTY ;
    }

    /* p [0..n_col] no longer needed, used as "head" in subsequent routines */

    /* === Scan columns, compute row degrees, and check row indices ========= */

    stats [COLAMD_INFO3] = 0 ;	/* number of duplicate or unsorted row indices*/

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].length = 0 ;
	Row [row].shared2.mark = -1 ;
    }

    for (col = 0 ; col < n_col ; col++)
    {
	last_row = -1 ;

	cp = &A [p [col]] ;
	cp_end = &A [p [col+1]] ;

	while (cp < cp_end)
	{
	    row = *cp++ ;

	    /* make sure row indices within range */
	    if (row < 0 || row >= n_row)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		stats [COLAMD_INFO3] = n_row ;
		DEBUG0 (("colamd: row %d col %d out of bounds\n", row, col)) ;
		return (FALSE) ;
	    }

	    if (row <= last_row || Row [row].shared2.mark == col)
	    {
		/* row index are unsorted or repeated (or both), thus col */
		/* is jumbled.  This is a notice, not an error condition. */
		stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		(stats [COLAMD_INFO3]) ++ ;
		DEBUG1 (("colamd: row %d col %d unsorted/duplicate\n",row,col));
	    }

	    if (Row [row].shared2.mark != col)
	    {
		Row [row].length++ ;
	    }
	    else
	    {
		/* this is a repeated entry in the column, */
		/* it will be removed */
		Col [col].length-- ;
	    }

	    /* mark the row as having been seen in this column */
	    Row [row].shared2.mark = col ;

	    last_row = row ;
	}
    }

    /* === Compute row pointers ============================================= */

    /* row form of the matrix starts directly after the column */
    /* form of matrix in A */
    Row [0].start = p [n_col] ;
    Row [0].shared1.p = Row [0].start ;
    Row [0].shared2.mark = -1 ;
    for (row = 1 ; row < n_row ; row++)
    {
	Row [row].start = Row [row-1].start + Row [row-1].length ;
	Row [row].shared1.p = Row [row].start ;
	Row [row].shared2.mark = -1 ;
    }

    /* === Create row form ================================================== */

    if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
    {
	/* if cols jumbled, watch for repeated row indices */
	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		row = *cp++ ;
		if (Row [row].shared2.mark != col)
		{
		    A [(Row [row].shared1.p)++] = col ;
		    Row [row].shared2.mark = col ;
		}
	    }
	}
    }
    else
    {
	/* if cols not jumbled, we don't need the mark (this is faster) */
	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		A [(Row [*cp++].shared1.p)++] = col ;
	    }
	}
    }

    /* === Clear the row marks and set row degrees ========================== */

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].shared2.mark = 0 ;
	Row [row].shared1.degree = Row [row].length ;
    }

    /* === See if we need to re-create columns ============================== */

    if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
    {
    	DEBUG0 (("colamd: reconstructing column form, matrix jumbled\n")) ;

	/* === Compute col pointers ========================================= */

	/* col form of the matrix starts at A [0]. */
	/* Note, we may have a gap between the col form and the row */
	/* form if there were duplicate entries, if so, it will be */
	/* removed upon the first garbage collection */
	Col [0].start = 0 ;
	p [0] = Col [0].start ;
	for (col = 1 ; col < n_col ; col++)
	{
	    /* note that the lengths here are for pruned columns, i.e. */
	    /* no duplicate row indices will exist for these columns */
	    Col [col].start = Col [col-1].start + Col [col-1].length ;
	    p [col] = Col [col].start ;
	}

	/* === Re-create col form =========================================== */

	for (row = 0 ; row < n_row ; row++)
	{
	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		A [(p [*rp++])++] = row ;
	    }
	}
    }

    /* === Done.  Matrix is not (or no longer) jumbled ====================== */

    return (TRUE) ;
}


/* ========================================================================== */
/* === init_scoring ========================================================= */
/* ========================================================================== */

/*
    Kills dense or empty columns and rows, calculates an initial score for
    each column, and places all columns in the degree lists.  Not user-callable.
*/

PRIVATE void init_scoring
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* column form and row form of A */
    Int head [],		/* of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameters */
    Int *p_n_row2,		/* number of non-dense, non-empty rows */
    Int *p_n_col2,		/* number of non-dense, non-empty columns */
    Int *p_max_deg		/* maximum row degree */
)
{
    /* === Local variables ================================================== */

    Int c ;			/* a column index */
    Int r, row ;		/* a row index */
    Int *cp ;			/* a column pointer */
    Int deg ;			/* degree of a row or column */
    Int *cp_end ;		/* a pointer to the end of a column */
    Int *new_cp ;		/* new column pointer */
    Int col_length ;		/* length of pruned column */
    Int score ;			/* current column score */
    Int n_col2 ;		/* number of non-dense, non-empty columns */
    Int n_row2 ;		/* number of non-dense, non-empty rows */
    Int dense_row_count ;	/* remove rows with more entries than this */
    Int dense_col_count ;	/* remove cols with more entries than this */
    Int min_score ;		/* smallest column score */
    Int max_deg ;		/* maximum row degree */
    Int next_col ;		/* Used to add to degree list.*/

    /* === Extract knobs ==================================================== */

    /* Note: if knobs contains a NaN, this is undefined: */
    if (knobs [COLAMD_DENSE_ROW] < 0)
    {
	/* only remove completely dense rows */
	dense_row_count = n_col-1 ;
    }
    else
    {
	dense_row_count = DENSE_DEGREE (knobs [COLAMD_DENSE_ROW], n_col) ;
    }
    if (knobs [COLAMD_DENSE_COL] < 0)
    {
	/* only remove completely dense columns */
	dense_col_count = n_row-1 ;
    }
    else
    {
	dense_col_count =
	    DENSE_DEGREE (knobs [COLAMD_DENSE_COL], MIN (n_row, n_col)) ;
    }

    DEBUG1 (("colamd: densecount: %d %d\n", dense_row_count, dense_col_count)) ;
    max_deg = 0 ;
    n_col2 = n_col ;
    n_row2 = n_row ;

    /* === Kill empty columns =============================================== */

    /* Put the empty columns at the end in their natural order, so that LU */
    /* factorization can proceed as far as possible. */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	deg = Col [c].length ;
	if (deg == 0)
	{
	    /* this is a empty column, kill and order it last */
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("colamd: null columns killed: %d\n", n_col - n_col2)) ;

    /* === Kill dense columns =============================================== */

    /* Put the dense columns at the end, in their natural order */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* skip any dead columns */
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	deg = Col [c].length ;
	if (deg > dense_col_count)
	{
	    /* this is a dense column, kill and order it last */
	    Col [c].shared2.order = --n_col2 ;
	    /* decrement the row degrees */
	    cp = &A [Col [c].start] ;
	    cp_end = cp + Col [c].length ;
	    while (cp < cp_end)
	    {
		Row [*cp++].shared1.degree-- ;
	    }
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("colamd: Dense and null columns killed: %d\n", n_col - n_col2)) ;

    /* === Kill dense and empty rows ======================================== */

    for (r = 0 ; r < n_row ; r++)
    {
	deg = Row [r].shared1.degree ;
	ASSERT (deg >= 0 && deg <= n_col) ;
	if (deg > dense_row_count || deg == 0)
	{
	    /* kill a dense or empty row */
	    KILL_ROW (r) ;
	    --n_row2 ;
	}
	else
	{
	    /* keep track of max degree of remaining rows */
	    max_deg = MAX (max_deg, deg) ;
	}
    }
    DEBUG1 (("colamd: Dense and null rows killed: %d\n", n_row - n_row2)) ;

    /* === Compute initial column scores ==================================== */

    /* At this point the row degrees are accurate.  They reflect the number */
    /* of "live" (non-dense) columns in each row.  No empty rows exist. */
    /* Some "live" columns may contain only dead rows, however.  These are */
    /* pruned in the code below. */

    /* now find the initial matlab score for each column */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* skip dead column */
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	score = 0 ;
	cp = &A [Col [c].start] ;
	new_cp = cp ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    /* get a row */
	    row = *cp++ ;
	    /* skip if dead */
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }
	    /* compact the column */
	    *new_cp++ = row ;
	    /* add row's external degree */
	    score += Row [row].shared1.degree - 1 ;
	    /* guard against integer overflow */
	    score = MIN (score, n_col) ;
	}
	/* determine pruned column length */
	col_length = (Int) (new_cp - &A [Col [c].start]) ;
	if (col_length == 0)
	{
	    /* a newly-made null column (all rows in this col are "dense" */
	    /* and have already been killed) */
	    DEBUG2 (("Newly null killed: %d\n", c)) ;
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	}
	else
	{
	    /* set column length and set score */
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    Col [c].length = col_length ;
	    Col [c].shared2.score = score ;
	}
    }
    DEBUG1 (("colamd: Dense, null, and newly-null columns killed: %d\n",
    	n_col-n_col2)) ;

    /* At this point, all empty rows and columns are dead.  All live columns */
    /* are "clean" (containing no dead rows) and simplicial (no supercolumns */
    /* yet).  Rows may contain dead columns, but all live rows contain at */
    /* least one live column. */

    /* === Initialize degree lists ========================================== */

    /* clear the hash buckets */
    for (c = 0 ; c <= n_col ; c++)
    {
	head [c] = EMPTY ;
    }
    min_score = n_col ;
    /* place in reverse order, so low column indices are at the front */
    /* of the lists.  This is to encourage natural tie-breaking */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* only add principal columns to degree lists */
	if (COL_IS_ALIVE (c))
	{
	    DEBUG4 (("place %d score %d minscore %d ncol %d\n",
		c, Col [c].shared2.score, min_score, n_col)) ;

	    /* === Add columns score to DList =============================== */

	    score = Col [c].shared2.score ;

	    ASSERT (min_score >= 0) ;
	    ASSERT (min_score <= n_col) ;
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    ASSERT (head [score] >= EMPTY) ;

	    /* now add this column to dList at proper score location */
	    next_col = head [score] ;
	    Col [c].shared3.prev = EMPTY ;
	    Col [c].shared4.degree_next = next_col ;

	    /* if there already was a column with the same score, set its */
	    /* previous pointer to this new column */
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = c ;
	    }
	    head [score] = c ;

	    /* see if this score is less than current min */
	    min_score = MIN (min_score, score) ;

	}
    }

    /* === Return number of remaining columns, and max row degree =========== */

    *p_n_col2 = n_col2 ;
    *p_n_row2 = n_row2 ;
    *p_max_deg = max_deg ;
}


/* ========================================================================== */
/* === find_ordering ======================================================== */
/* ========================================================================== */

/*
    Order the principal columns of the supercolumn form of the matrix
    (no supercolumns on input).  Uses a minimum approximate column minimum
    degree ordering method.  Not user-callable.
*/

PRIVATE Int find_ordering	/* return the number of garbage collections */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Int Alen,			/* size of A, 2*nnz + n_col or larger */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* column form and row form of A */
    Int head [],		/* of size n_col+1 */
    Int n_col2,			/* Remaining columns to order */
    Int max_deg,		/* Maximum row degree */
    Int pfree,			/* index of first free slot (2*nnz on entry) */
    Int aggressive
)
{
    /* === Local variables ================================================== */

    Int k ;			/* current pivot ordering step */
    Int pivot_col ;		/* current pivot column */
    Int *cp ;			/* a column pointer */
    Int *rp ;			/* a row pointer */
    Int pivot_row ;		/* current pivot row */
    Int *new_cp ;		/* modified column pointer */
    Int *new_rp ;		/* modified row pointer */
    Int pivot_row_start ;	/* pointer to start of pivot row */
    Int pivot_row_degree ;	/* number of columns in pivot row */
    Int pivot_row_length ;	/* number of supercolumns in pivot row */
    Int pivot_col_score ;	/* score of pivot column */
    Int needed_memory ;		/* free space needed for pivot row */
    Int *cp_end ;		/* pointer to the end of a column */
    Int *rp_end ;		/* pointer to the end of a row */
    Int row ;			/* a row index */
    Int col ;			/* a column index */
    Int max_score ;		/* maximum possible score */
    Int cur_score ;		/* score of current column */
    unsigned Int hash ;		/* hash value for supernode detection */
    Int head_column ;		/* head of hash bucket */
    Int first_col ;		/* first column in hash bucket */
    Int tag_mark ;		/* marker value for mark array */
    Int row_mark ;		/* Row [row].shared2.mark */
    Int set_difference ;	/* set difference size of row with pivot row */
    Int min_score ;		/* smallest column score */
    Int col_thickness ;		/* "thickness" (no. of columns in a supercol) */
    Int max_mark ;		/* maximum value of tag_mark */
    Int pivot_col_thickness ;	/* number of columns represented by pivot col */
    Int prev_col ;		/* Used by Dlist operations. */
    Int next_col ;		/* Used by Dlist operations. */
    Int ngarbage ;		/* number of garbage collections performed */

    /* === Initialization and clear mark ==================================== */

    max_mark = INT_MAX - n_col ;	/* INT_MAX defined in <limits.h> */
    tag_mark = clear_mark (0, max_mark, n_row, Row) ;
    min_score = 0 ;
    ngarbage = 0 ;
    DEBUG1 (("colamd: Ordering, n_col2=%d\n", n_col2)) ;

    /* === Order the columns ================================================ */

    for (k = 0 ; k < n_col2 ; /* 'k' is incremented below */)
    {

	/* === Select pivot column, and order it ============================ */

	/* make sure degree list isn't empty */
	ASSERT (min_score >= 0) ;
	ASSERT (min_score <= n_col) ;
	ASSERT (head [min_score] >= EMPTY) ;

	/* get pivot column from head of minimum degree list */
	while (head [min_score] == EMPTY && min_score < n_col)
	{
	    min_score++ ;
	}
	pivot_col = head [min_score] ;
	ASSERT (pivot_col >= 0 && pivot_col <= n_col) ;
	next_col = Col [pivot_col].shared4.degree_next ;
	head [min_score] = next_col ;
	if (next_col != EMPTY)
	{
	    Col [next_col].shared3.prev = EMPTY ;
	}

	ASSERT (COL_IS_ALIVE (pivot_col)) ;

	/* remember score for defrag check */
	pivot_col_score = Col [pivot_col].shared2.score ;

	/* the pivot column is the kth column in the pivot order */
	Col [pivot_col].shared2.order = k ;

	/* increment order count by column thickness */
	pivot_col_thickness = Col [pivot_col].shared1.thickness ;
	k += pivot_col_thickness ;
	ASSERT (pivot_col_thickness > 0) ;
	DEBUG3 (("Pivot col: %d thick %d\n", pivot_col, pivot_col_thickness)) ;

	/* === Garbage_collection, if necessary ============================= */

	needed_memory = MIN (pivot_col_score, n_col - k) ;
	if (pfree + needed_memory >= Alen)
	{
	    pfree = garbage_collection (n_row, n_col, Row, Col, A, &A [pfree]) ;
	    ngarbage++ ;
	    /* after garbage collection we will have enough */
	    ASSERT (pfree + needed_memory < Alen) ;
	    /* garbage collection has wiped out the Row[].shared2.mark array */
	    tag_mark = clear_mark (0, max_mark, n_row, Row) ;
	}

	/* === Compute pivot row pattern ==================================== */

	/* get starting location for this new merged row */
	pivot_row_start = pfree ;

	/* initialize new row counts to zero */
	pivot_row_degree = 0 ;

	/* tag pivot column as having been visited so it isn't included */
	/* in merged pivot row */
	Col [pivot_col].shared1.thickness = -pivot_col_thickness ;

	/* pivot row is the union of all rows in the pivot column pattern */
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    /* get a row */
	    row = *cp++ ;
	    DEBUG4 (("Pivot col pattern %d %d\n", ROW_IS_ALIVE (row), row)) ;
	    /* skip if row is dead */
	    if (ROW_IS_ALIVE (row))
	    {
		rp = &A [Row [row].start] ;
		rp_end = rp + Row [row].length ;
		while (rp < rp_end)
		{
		    /* get a column */
		    col = *rp++ ;
		    /* add the column, if alive and untagged */
		    col_thickness = Col [col].shared1.thickness ;
		    if (col_thickness > 0 && COL_IS_ALIVE (col))
		    {
			/* tag column in pivot row */
			Col [col].shared1.thickness = -col_thickness ;
			ASSERT (pfree < Alen) ;
			/* place column in pivot row */
			A [pfree++] = col ;
			pivot_row_degree += col_thickness ;
		    }
		}
	    }
	}

	/* clear tag on pivot column */
	Col [pivot_col].shared1.thickness = pivot_col_thickness ;
	max_deg = MAX (max_deg, pivot_row_degree) ;

	/* === Kill all rows used to construct pivot row ==================== */

	/* also kill pivot row, temporarily */
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    /* may be killing an already dead row */
	    row = *cp++ ;
	    DEBUG3 (("Kill row in pivot col: %d\n", row)) ;
	    KILL_ROW (row) ;
	}

	/* === Select a row index to use as the new pivot row =============== */

	pivot_row_length = pfree - pivot_row_start ;
	if (pivot_row_length > 0)
	{
	    /* pick the "pivot" row arbitrarily (first row in col) */
	    pivot_row = A [Col [pivot_col].start] ;
	    DEBUG3 (("Pivotal row is %d\n", pivot_row)) ;
	}
	else
	{
	    /* there is no pivot row, since it is of zero length */
	    pivot_row = EMPTY ;
	    ASSERT (pivot_row_length == 0) ;
	}
	ASSERT (Col [pivot_col].length > 0 || pivot_row_length == 0) ;

	/* === Approximate degree computation =============================== */

	/* Here begins the computation of the approximate degree.  The column */
	/* score is the sum of the pivot row "length", plus the size of the */
	/* set differences of each row in the column minus the pattern of the */
	/* pivot row itself.  The column ("thickness") itself is also */
	/* excluded from the column score (we thus use an approximate */
	/* external degree). */

	/* The time taken by the following code (compute set differences, and */
	/* add them up) is proportional to the size of the data structure */
	/* being scanned - that is, the sum of the sizes of each column in */
	/* the pivot row.  Thus, the amortized time to compute a column score */
	/* is proportional to the size of that column (where size, in this */
	/* context, is the column "length", or the number of row indices */
	/* in that column).  The number of row indices in a column is */
	/* monotonically non-decreasing, from the length of the original */
	/* column on input to colamd. */

	/* === Compute set differences ====================================== */

	DEBUG3 (("** Computing set differences phase. **\n")) ;

	/* pivot row is currently dead - it will be revived later. */

	DEBUG3 (("Pivot row: ")) ;
	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    DEBUG3 (("Col: %d\n", col)) ;

	    /* clear tags used to construct pivot row pattern */
	    col_thickness = -Col [col].shared1.thickness ;
	    ASSERT (col_thickness > 0) ;
	    Col [col].shared1.thickness = col_thickness ;

	    /* === Remove column from degree list =========================== */

	    cur_score = Col [col].shared2.score ;
	    prev_col = Col [col].shared3.prev ;
	    next_col = Col [col].shared4.degree_next ;
	    ASSERT (cur_score >= 0) ;
	    ASSERT (cur_score <= n_col) ;
	    ASSERT (cur_score >= EMPTY) ;
	    if (prev_col == EMPTY)
	    {
		head [cur_score] = next_col ;
	    }
	    else
	    {
		Col [prev_col].shared4.degree_next = next_col ;
	    }
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = prev_col ;
	    }

	    /* === Scan the column ========================================== */

	    cp = &A [Col [col].start] ;
	    cp_end = cp + Col [col].length ;
	    while (cp < cp_end)
	    {
		/* get a row */
		row = *cp++ ;
		row_mark = Row [row].shared2.mark ;
		/* skip if dead */
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    continue ;
		}
		ASSERT (row != pivot_row) ;
		set_difference = row_mark - tag_mark ;
		/* check if the row has been seen yet */
		if (set_difference < 0)
		{
		    ASSERT (Row [row].shared1.degree <= max_deg) ;
		    set_difference = Row [row].shared1.degree ;
		}
		/* subtract column thickness from this row's set difference */
		set_difference -= col_thickness ;
		ASSERT (set_difference >= 0) ;
		/* absorb this row if the set difference becomes zero */
		if (set_difference == 0 && aggressive)
		{
		    DEBUG3 (("aggressive absorption. Row: %d\n", row)) ;
		    KILL_ROW (row) ;
		}
		else
		{
		    /* save the new mark */
		    Row [row].shared2.mark = set_difference + tag_mark ;
		}
	    }
	}

	/* === Add up set differences for each column ======================= */

	DEBUG3 (("** Adding set differences phase. **\n")) ;

	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    /* get a column */
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    hash = 0 ;
	    cur_score = 0 ;
	    cp = &A [Col [col].start] ;
	    /* compact the column */
	    new_cp = cp ;
	    cp_end = cp + Col [col].length ;

	    DEBUG4 (("Adding set diffs for Col: %d.\n", col)) ;

	    while (cp < cp_end)
	    {
		/* get a row */
		row = *cp++ ;
		ASSERT(row >= 0 && row < n_row) ;
		row_mark = Row [row].shared2.mark ;
		/* skip if dead */
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    DEBUG4 ((" Row %d, dead\n", row)) ;
		    continue ;
		}
		DEBUG4 ((" Row %d, set diff %d\n", row, row_mark-tag_mark));
		ASSERT (row_mark >= tag_mark) ;
		/* compact the column */
		*new_cp++ = row ;
		/* compute hash function */
		hash += row ;
		/* add set difference */
		cur_score += row_mark - tag_mark ;
		/* integer overflow... */
		cur_score = MIN (cur_score, n_col) ;
	    }

	    /* recompute the column's length */
	    Col [col].length = (Int) (new_cp - &A [Col [col].start]) ;

	    /* === Further mass elimination ================================= */

	    if (Col [col].length == 0)
	    {
		DEBUG4 (("further mass elimination. Col: %d\n", col)) ;
		/* nothing left but the pivot row in this column */
		KILL_PRINCIPAL_COL (col) ;
		pivot_row_degree -= Col [col].shared1.thickness ;
		ASSERT (pivot_row_degree >= 0) ;
		/* order it */
		Col [col].shared2.order = k ;
		/* increment order count by column thickness */
		k += Col [col].shared1.thickness ;
	    }
	    else
	    {
		/* === Prepare for supercolumn detection ==================== */

		DEBUG4 (("Preparing supercol detection for Col: %d.\n", col)) ;

		/* save score so far */
		Col [col].shared2.score = cur_score ;

		/* add column to hash table, for supercolumn detection */
		hash %= n_col + 1 ;

		DEBUG4 ((" Hash = %d, n_col = %d.\n", hash, n_col)) ;
		ASSERT (((Int) hash) <= n_col) ;

		head_column = head [hash] ;
		if (head_column > EMPTY)
		{
		    /* degree list "hash" is non-empty, use prev (shared3) of */
		    /* first column in degree list as head of hash bucket */
		    first_col = Col [head_column].shared3.headhash ;
		    Col [head_column].shared3.headhash = col ;
		}
		else
		{
		    /* degree list "hash" is empty, use head as hash bucket */
		    first_col = - (head_column + 2) ;
		    head [hash] = - (col + 2) ;
		}
		Col [col].shared4.hash_next = first_col ;

		/* save hash function in Col [col].shared3.hash */
		Col [col].shared3.hash = (Int) hash ;
		ASSERT (COL_IS_ALIVE (col)) ;
	    }
	}

	/* The approximate external column degree is now computed.  */

	/* === Supercolumn detection ======================================== */

	DEBUG3 (("** Supercolumn detection phase. **\n")) ;

	detect_super_cols (
		Col, A, head, pivot_row_start, pivot_row_length) ;

	/* === Kill the pivotal column ====================================== */

	KILL_PRINCIPAL_COL (pivot_col) ;

	/* === Clear mark =================================================== */

	tag_mark = clear_mark (tag_mark+max_deg+1, max_mark, n_row, Row) ;

	/* === Finalize the new pivot row, and column scores ================ */

	DEBUG3 (("** Finalize scores phase. **\n")) ;

	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	/* compact the pivot row */
	new_rp = rp ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    /* skip dead columns */
	    if (COL_IS_DEAD (col))
	    {
		continue ;
	    }
	    *new_rp++ = col ;
	    /* add new pivot row to column */
	    A [Col [col].start + (Col [col].length++)] = pivot_row ;

	    /* retrieve score so far and add on pivot row's degree. */
	    /* (we wait until here for this in case the pivot */
	    /* row's degree was reduced due to mass elimination). */
	    cur_score = Col [col].shared2.score + pivot_row_degree ;

	    /* calculate the max possible score as the number of */
	    /* external columns minus the 'k' value minus the */
	    /* columns thickness */
	    max_score = n_col - k - Col [col].shared1.thickness ;

	    /* make the score the external degree of the union-of-rows */
	    cur_score -= Col [col].shared1.thickness ;

	    /* make sure score is less or equal than the max score */
	    cur_score = MIN (cur_score, max_score) ;
	    ASSERT (cur_score >= 0) ;

	    /* store updated score */
	    Col [col].shared2.score = cur_score ;

	    /* === Place column back in degree list ========================= */

	    ASSERT (min_score >= 0) ;
	    ASSERT (min_score <= n_col) ;
	    ASSERT (cur_score >= 0) ;
	    ASSERT (cur_score <= n_col) ;
	    ASSERT (head [cur_score] >= EMPTY) ;
	    next_col = head [cur_score] ;
	    Col [col].shared4.degree_next = next_col ;
	    Col [col].shared3.prev = EMPTY ;
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = col ;
	    }
	    head [cur_score] = col ;

	    /* see if this score is less than current min */
	    min_score = MIN (min_score, cur_score) ;

	}

	/* === Resurrect the new pivot row ================================== */

	if (pivot_row_degree > 0)
	{
	    /* update pivot row length to reflect any cols that were killed */
	    /* during super-col detection and mass elimination */
	    Row [pivot_row].start  = pivot_row_start ;
	    Row [pivot_row].length = (Int) (new_rp - &A[pivot_row_start]) ;
	    ASSERT (Row [pivot_row].length > 0) ;
	    Row [pivot_row].shared1.degree = pivot_row_degree ;
	    Row [pivot_row].shared2.mark = 0 ;
	    /* pivot row is no longer dead */

	    DEBUG1 (("Resurrect Pivot_row %d deg: %d\n",
			pivot_row, pivot_row_degree)) ;
	}
    }

    /* === All principal columns have now been ordered ====================== */

    return (ngarbage) ;
}


/* ========================================================================== */
/* === order_children ======================================================= */
/* ========================================================================== */

/*
    The find_ordering routine has ordered all of the principal columns (the
    representatives of the supercolumns).  The non-principal columns have not
    yet been ordered.  This routine orders those columns by walking up the
    parent tree (a column is a child of the column which absorbed it).  The
    final permutation vector is then placed in p [0 ... n_col-1], with p [0]
    being the first column, and p [n_col-1] being the last.  It doesn't look
    like it at first glance, but be assured that this routine takes time linear
    in the number of columns.  Although not immediately obvious, the time
    taken by this routine is O (n_col), that is, linear in the number of
    columns.  Not user-callable.
*/

PRIVATE void order_children
(
    /* === Parameters ======================================================= */

    Int n_col,			/* number of columns of A */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int p []			/* p [0 ... n_col-1] is the column permutation*/
)
{
    /* === Local variables ================================================== */

    Int i ;			/* loop counter for all columns */
    Int c ;			/* column index */
    Int parent ;		/* index of column's parent */
    Int order ;			/* column's order */

    /* === Order each non-principal column ================================== */

    for (i = 0 ; i < n_col ; i++)
    {
	/* find an un-ordered non-principal column */
	ASSERT (COL_IS_DEAD (i)) ;
	if (!COL_IS_DEAD_PRINCIPAL (i) && Col [i].shared2.order == EMPTY)
	{
	    parent = i ;
	    /* once found, find its principal parent */
	    do
	    {
		parent = Col [parent].shared1.parent ;
	    } while (!COL_IS_DEAD_PRINCIPAL (parent)) ;

	    /* now, order all un-ordered non-principal columns along path */
	    /* to this parent.  collapse tree at the same time */
	    c = i ;
	    /* get order of parent */
	    order = Col [parent].shared2.order ;

	    do
	    {
		ASSERT (Col [c].shared2.order == EMPTY) ;

		/* order this column */
		Col [c].shared2.order = order++ ;
		/* collaps tree */
		Col [c].shared1.parent = parent ;

		/* get immediate parent of this column */
		c = Col [c].shared1.parent ;

		/* continue until we hit an ordered column.  There are */
		/* guarranteed not to be anymore unordered columns */
		/* above an ordered column */
	    } while (Col [c].shared2.order == EMPTY) ;

	    /* re-order the super_col parent to largest order for this group */
	    Col [parent].shared2.order = order ;
	}
    }

    /* === Generate the permutation ========================================= */

    for (c = 0 ; c < n_col ; c++)
    {
	p [Col [c].shared2.order] = c ;
    }
}


/* ========================================================================== */
/* === detect_super_cols ==================================================== */
/* ========================================================================== */

/*
    Detects supercolumns by finding matches between columns in the hash buckets.
    Check amongst columns in the set A [row_start ... row_start + row_length-1].
    The columns under consideration are currently *not* in the degree lists,
    and have already been placed in the hash buckets.

    The hash bucket for columns whose hash function is equal to h is stored
    as follows:

	if head [h] is >= 0, then head [h] contains a degree list, so:

		head [h] is the first column in degree bucket h.
		Col [head [h]].headhash gives the first column in hash bucket h.

	otherwise, the degree list is empty, and:

		-(head [h] + 2) is the first column in hash bucket h.

    For a column c in a hash bucket, Col [c].shared3.prev is NOT a "previous
    column" pointer.  Col [c].shared3.hash is used instead as the hash number
    for that column.  The value of Col [c].shared4.hash_next is the next column
    in the same hash bucket.

    Assuming no, or "few" hash collisions, the time taken by this routine is
    linear in the sum of the sizes (lengths) of each column whose score has
    just been computed in the approximate degree computation.
    Not user-callable.
*/

PRIVATE void detect_super_cols
(
    /* === Parameters ======================================================= */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* row indices of A */
    Int head [],		/* head of degree lists and hash buckets */
    Int row_start,		/* pointer to set of columns to check */
    Int row_length		/* number of columns to check */
)
{
    /* === Local variables ================================================== */

    Int hash ;			/* hash value for a column */
    Int *rp ;			/* pointer to a row */
    Int c ;			/* a column index */
    Int super_c ;		/* column index of the column to absorb into */
    Int *cp1 ;			/* column pointer for column super_c */
    Int *cp2 ;			/* column pointer for column c */
    Int length ;		/* length of column super_c */
    Int prev_c ;		/* column preceding c in hash bucket */
    Int i ;			/* loop counter */
    Int *rp_end ;		/* pointer to the end of the row */
    Int col ;			/* a column index in the row to check */
    Int head_column ;		/* first column in hash bucket or degree list */
    Int first_col ;		/* first column in hash bucket */

    /* === Consider each column in the row ================================== */

    rp = &A [row_start] ;
    rp_end = rp + row_length ;
    while (rp < rp_end)
    {
	col = *rp++ ;
	if (COL_IS_DEAD (col))
	{
	    continue ;
	}

	/* get hash number for this column */
	hash = Col [col].shared3.hash ;
	ASSERT (hash <= n_col) ;

	/* === Get the first column in this hash bucket ===================== */

	head_column = head [hash] ;
	if (head_column > EMPTY)
	{
	    first_col = Col [head_column].shared3.headhash ;
	}
	else
	{
	    first_col = - (head_column + 2) ;
	}

	/* === Consider each column in the hash bucket ====================== */

	for (super_c = first_col ; super_c != EMPTY ;
	    super_c = Col [super_c].shared4.hash_next)
	{
	    ASSERT (COL_IS_ALIVE (super_c)) ;
	    ASSERT (Col [super_c].shared3.hash == hash) ;
	    length = Col [super_c].length ;

	    /* prev_c is the column preceding column c in the hash bucket */
	    prev_c = super_c ;

	    /* === Compare super_c with all columns after it ================ */

	    for (c = Col [super_c].shared4.hash_next ;
		 c != EMPTY ; c = Col [c].shared4.hash_next)
	    {
		ASSERT (c != super_c) ;
		ASSERT (COL_IS_ALIVE (c)) ;
		ASSERT (Col [c].shared3.hash == hash) ;

		/* not identical if lengths or scores are different */
		if (Col [c].length != length ||
		    Col [c].shared2.score != Col [super_c].shared2.score)
		{
		    prev_c = c ;
		    continue ;
		}

		/* compare the two columns */
		cp1 = &A [Col [super_c].start] ;
		cp2 = &A [Col [c].start] ;

		for (i = 0 ; i < length ; i++)
		{
		    /* the columns are "clean" (no dead rows) */
		    ASSERT (ROW_IS_ALIVE (*cp1))  ;
		    ASSERT (ROW_IS_ALIVE (*cp2))  ;
		    /* row indices will same order for both supercols, */
		    /* no gather scatter nessasary */
		    if (*cp1++ != *cp2++)
		    {
			break ;
		    }
		}

		/* the two columns are different if the for-loop "broke" */
		if (i != length)
		{
		    prev_c = c ;
		    continue ;
		}

		/* === Got it!  two columns are identical =================== */

		ASSERT (Col [c].shared2.score == Col [super_c].shared2.score) ;

		Col [super_c].shared1.thickness += Col [c].shared1.thickness ;
		Col [c].shared1.parent = super_c ;
		KILL_NON_PRINCIPAL_COL (c) ;
		/* order c later, in order_children() */
		Col [c].shared2.order = EMPTY ;
		/* remove c from hash bucket */
		Col [prev_c].shared4.hash_next = Col [c].shared4.hash_next ;
	    }
	}

	/* === Empty this hash bucket ======================================= */

	if (head_column > EMPTY)
	{
	    /* corresponding degree list "hash" is not empty */
	    Col [head_column].shared3.headhash = EMPTY ;
	}
	else
	{
	    /* corresponding degree list "hash" is empty */
	    head [hash] = EMPTY ;
	}
    }
}


/* ========================================================================== */
/* === garbage_collection =================================================== */
/* ========================================================================== */

/*
    Defragments and compacts columns and rows in the workspace A.  Used when
    all avaliable memory has been used while performing row merging.  Returns
    the index of the first free position in A, after garbage collection.  The
    time taken by this routine is linear is the size of the array A, which is
    itself linear in the number of nonzeros in the input matrix.
    Not user-callable.
*/

PRIVATE Int garbage_collection  /* returns the new value of pfree */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows */
    Int n_col,			/* number of columns */
    Colamd_Row Row [],		/* row info */
    Colamd_Col Col [],		/* column info */
    Int A [],			/* A [0 ... Alen-1] holds the matrix */
    Int *pfree			/* &A [0] ... pfree is in use */
)
{
    /* === Local variables ================================================== */

    Int *psrc ;			/* source pointer */
    Int *pdest ;		/* destination pointer */
    Int j ;			/* counter */
    Int r ;			/* a row index */
    Int c ;			/* a column index */
    Int length ;		/* length of a row or column */

    /* === Defragment the columns =========================================== */

    pdest = &A[0] ;
    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    psrc = &A [Col [c].start] ;

	    /* move and compact the column */
	    ASSERT (pdest <= psrc) ;
	    Col [c].start = (Int) (pdest - &A [0]) ;
	    length = Col [c].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		r = *psrc++ ;
		if (ROW_IS_ALIVE (r))
		{
		    *pdest++ = r ;
		}
	    }
	    Col [c].length = (Int) (pdest - &A [Col [c].start]) ;
	}
    }

    /* === Prepare to defragment the rows =================================== */

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_DEAD (r) || (Row [r].length == 0))
	{
	    /* This row is already dead, or is of zero length.  Cannot compact
	     * a row of zero length, so kill it.  NOTE: in the current version,
	     * there are no zero-length live rows.  Kill the row (for the first
	     * time, or again) just to be safe. */
	    KILL_ROW (r) ;
	}
	else
	{
	    /* save first column index in Row [r].shared2.first_column */
	    psrc = &A [Row [r].start] ;
	    Row [r].shared2.first_column = *psrc ;
	    ASSERT (ROW_IS_ALIVE (r)) ;
	    /* flag the start of the row with the one's complement of row */
	    *psrc = ONES_COMPLEMENT (r) ;
	}
    }

    /* === Defragment the rows ============================================== */

    psrc = pdest ;
    while (psrc < pfree)
    {
	/* find a negative number ... the start of a row */
	if (*psrc++ < 0)
	{
	    psrc-- ;
	    /* get the row index */
	    r = ONES_COMPLEMENT (*psrc) ;
	    ASSERT (r >= 0 && r < n_row) ;
	    /* restore first column index */
	    *psrc = Row [r].shared2.first_column ;
	    ASSERT (ROW_IS_ALIVE (r)) ;
	    ASSERT (Row [r].length > 0) ;
	    /* move and compact the row */
	    ASSERT (pdest <= psrc) ;
	    Row [r].start = (Int) (pdest - &A [0]) ;
	    length = Row [r].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		c = *psrc++ ;
		if (COL_IS_ALIVE (c))
		{
		    *pdest++ = c ;
		}
	    }
	    Row [r].length = (Int) (pdest - &A [Row [r].start]) ;
	    ASSERT (Row [r].length > 0) ;
	}
    }
    /* ensure we found all the rows */
    ASSERT (debug_rows == 0) ;

    /* === Return the new value of pfree ==================================== */

    return ((Int) (pdest - &A [0])) ;
}


/* ========================================================================== */
/* === clear_mark =========================================================== */
/* ========================================================================== */

/*
    Clears the Row [].shared2.mark array, and returns the new tag_mark.
    Return value is the new tag_mark.  Not user-callable.
*/

PRIVATE Int clear_mark	/* return the new value for tag_mark */
(
    /* === Parameters ======================================================= */

    Int tag_mark,	/* new value of tag_mark */
    Int max_mark,	/* max allowed value of tag_mark */

    Int n_row,		/* number of rows in A */
    Colamd_Row Row []	/* Row [0 ... n_row-1].shared2.mark is set to zero */
)
{
    /* === Local variables ================================================== */

    Int r ;

    if (tag_mark <= 0 || tag_mark >= max_mark)
    {
	for (r = 0 ; r < n_row ; r++)
	{
	    if (ROW_IS_ALIVE (r))
	    {
		Row [r].shared2.mark = 0 ;
	    }
	}
	tag_mark = 1 ;
    }

    return (tag_mark) ;
}


/* ========================================================================== */
/* === print_report ========================================================= */
/* ========================================================================== */

PRIVATE void print_report
(
    char *method,
    Int stats [COLAMD_STATS]
)
{

    Int i1, i2, i3 ;

    SPARSE_PRINTF (("\n%s version %d.%d, %s: ", method,
            COLAMD_MAIN_VERSION, COLAMD_SUB_VERSION, COLAMD_DATE)) ;

    if (!stats)
    {
        SPARSE_PRINTF (("No statistics available.\n")) ;
	return ;
    }

    i1 = stats [COLAMD_INFO1] ;
    i2 = stats [COLAMD_INFO2] ;
    i3 = stats [COLAMD_INFO3] ;

    if (stats [COLAMD_STATUS] >= 0)
    {
        SPARSE_PRINTF (("OK.  ")) ;
    }
    else
    {
        SPARSE_PRINTF (("ERROR.  ")) ;
    }

    switch (stats [COLAMD_STATUS])
    {

	case COLAMD_OK_BUT_JUMBLED:

            SPARSE_PRINTF((
                    "Matrix has unsorted or duplicate row indices.\n")) ;

            SPARSE_PRINTF((
                    "%s: number of duplicate or out-of-order row indices: %d\n",
                    method, i3)) ;

            SPARSE_PRINTF((
                    "%s: last seen duplicate or out-of-order row index:   %d\n",
                    method, INDEX (i2))) ;

            SPARSE_PRINTF((
                    "%s: last seen in column:                             %d",
                    method, INDEX (i1))) ;

	    /* no break - fall through to next case instead */

	case COLAMD_OK:

            SPARSE_PRINTF(("\n")) ;

            SPARSE_PRINTF((
                    "%s: number of dense or empty rows ignored:           %d\n",
                    method, stats [COLAMD_DENSE_ROW])) ;

            SPARSE_PRINTF((
                    "%s: number of dense or empty columns ignored:        %d\n",
                    method, stats [COLAMD_DENSE_COL])) ;

            SPARSE_PRINTF((
                    "%s: number of garbage collections performed:         %d\n",
                    method, stats [COLAMD_DEFRAG_COUNT])) ;
	    break ;

	case COLAMD_ERROR_A_not_present:

	    SPARSE_PRINTF((
                    "Array A (row indices of matrix) not present.\n")) ;
	    break ;

	case COLAMD_ERROR_p_not_present:

            SPARSE_PRINTF((
                    "Array p (column pointers for matrix) not present.\n")) ;
	    break ;

	case COLAMD_ERROR_nrow_negative:

            SPARSE_PRINTF(("Invalid number of rows (%d).\n", i1)) ;
	    break ;

	case COLAMD_ERROR_ncol_negative:

            SPARSE_PRINTF(("Invalid number of columns (%d).\n", i1)) ;
	    break ;

	case COLAMD_ERROR_nnz_negative:

            SPARSE_PRINTF((
                    "Invalid number of nonzero entries (%d).\n", i1)) ;
	    break ;

	case COLAMD_ERROR_p0_nonzero:

            SPARSE_PRINTF((
                    "Invalid column pointer, p [0] = %d, must be zero.\n", i1));
	    break ;

	case COLAMD_ERROR_A_too_small:

            SPARSE_PRINTF(("Array A too small.\n")) ;
            SPARSE_PRINTF((
                    "        Need Alen >= %d, but given only Alen = %d.\n",
                    i1, i2)) ;
	    break ;

	case COLAMD_ERROR_col_length_negative:

            SPARSE_PRINTF
            (("Column %d has a negative number of nonzero entries (%d).\n",
            INDEX (i1), i2)) ;
	    break ;

	case COLAMD_ERROR_row_index_out_of_bounds:

            SPARSE_PRINTF
            (("Row index (row %d) out of bounds (%d to %d) in column %d.\n",
            INDEX (i2), INDEX (0), INDEX (i3-1), INDEX (i1))) ;
	    break ;

	case COLAMD_ERROR_out_of_memory:

            SPARSE_PRINTF(("Out of memory.\n")) ;
	    break ;

	/* v2.4: internal-error case deleted */
    }
}
