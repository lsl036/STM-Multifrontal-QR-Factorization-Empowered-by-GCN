/**
 * @file SparseCore_check.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef NCHECK

#include "Sparse_internal.h"
#include "SparseCore.h"

/* ========================================================================== */
/* === 输出定义 ================================================= */
/* ========================================================================== */

#ifdef LONG
#define I8 "%8ld"
#define I_8 "%-8ld"
#else
#define I8 "%8d"
#define I_8 "%-8d"
#endif

#define PR(i,format,arg) \
{ \
    if (print >= i && SparseBase_config.printf_func != NULL) \
    { \
	SparseBase_config.printf_func (format, arg) ; \
    } \
}

#define P1(format,arg) PR(1,format,arg)
#define P2(format,arg) PR(2,format,arg)
#define P3(format,arg) PR(3,format,arg)
#define P4(format,arg) PR(4,format,arg)

#define ERR(msg) \
{ \
    P1 ("\nHNUCHOL ERROR: %s: ", type) ; \
    if (name != NULL) \
    { \
	P1 ("%s", name) ; \
    } \
    P1 (": %s\n", msg) ; \
    ERROR (SPARSE_INVALID, "invalid") ; \
    return (FALSE) ; \
}

/* 输出数值 */
#define PRINTVALUE(value) \
{ \
    if (Common->precise) \
    { \
	P4 (" %23.15e", value) ; \
    } \
    else \
    { \
	P4 (" %.5g", value) ; \
    } \
}

/* 开始输出 */
#define ETC_START(count,limit) \
{ \
    count = (init_print == 4) ? (limit) : (-1) ; \
}

/* 如果满足条件，则重新启用输出 */
#define ETC_ENABLE(condition,count,limit) \
{ \
    if ((condition) && init_print == 4) \
    { \
	count = limit ; \
	print = 4 ; \
    } \
}

/* 如果达到限制，请关闭输出 */
#define ETC_DISABLE(count) \
{ \
    if ((count >= 0) && (count-- == 0) && print == 4) \
    { \
	P4 ("%s", "    ...\n")  ; \
	print = 3 ; \
    } \
}

/* 重新启用输出，或在达到限制后关闭 */
#define ETC(condition,count,limit) \
{ \
    ETC_ENABLE (condition, count, limit) ; \
    ETC_DISABLE (count) ; \
}

#define BOOLSTR(x) ((x) ? "true " : "false")

/* ========================================================================== */
/* === 输出值 ========================================================== */
/* ========================================================================== */

static void print_value
(
    Int print,
    Int xtype,
    double *Xx,
    double *Xz,
    Int p,
    sparse_common *Common)
{
    if (xtype == SPARSE_REAL)
    {
	PRINTVALUE (Xx [p]) ;
    }
}

/**
 * @brief 	输出并验证Common中的值
 * 
 * @param print 
 * @param name 
 * @param Common 
 * @return int 
 */
static int check_common
(
    Int print,
    const char *name,
    sparse_common *Common
)
{
    double fl, lnz ;
    double *Xwork ;
    Int *Flag, *Head ;
    Sparse_long mark ;
    Int i, nrow, nmethods, ordering, xworksize, amd_backup, init_print ;
    const char *type = "common" ;

    /* ---------------------------------------------------------------------- */
    /* 输出控制参数和统计数据 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    init_print = print ;

    if (name != NULL)
    {
	P1 ("%s: ", name) ;
    }
    switch (Common->status)
    {

	case SPARSE_OK:
	    P1 ("%s", "status: OK\n") ;
	    break ;

	case SPARSE_OUT_OF_MEMORY:
	    P1 ("%s", "status: ERROR, 内存溢出\n") ;
	    break ;

	case SPARSE_INVALID:
	    P1 ("%s", "status: ERROR, invalid parameter\n") ;
	    break ;

	case SPARSE_TOO_LARGE:
	    P1 ("%s", "status: ERROR, problem too large\n") ;
	    break ;

	case SPARSE_NOT_INSTALLED:
	    P1 ("%s", "status: ERROR, method not installed\n") ;
	    break ;

	case SPARSE_NOT_POSDEF:
	    P1 ("%s", "status: warning, matrix not positive definite\n") ;
	    break ;

	case SPARSE_DSMALL:
	    P1 ("%s", "status: warning, diagonal entry has tiny abs. value\n") ;
	    break ;

	default:
	    ERR ("unknown status") ;
    }

    P2 ("  Architecture: %s\n", SPARSE_ARCHITECTURE) ;
    P3 ("    sizeof(int):      %d\n", (int) sizeof (int)) ;
    P3 ("    sizeof(Sparse_long):  %d\n", (int) sizeof (Sparse_long));
    P3 ("    sizeof(void *):   %d\n", (int) sizeof (void *)) ;
    P3 ("    sizeof(double):   %d\n", (int) sizeof (double)) ;
    P3 ("    sizeof(Int):      %d (HNUCHOL's basic integer)\n", (int) sizeof (Int)) ;
    P3 ("    sizeof(BLAS_INT): %d (integer used in the BLAS)\n",
	    (int) sizeof (BLAS_INT)) ;

    if (Common->fl != EMPTY)
    {
	P2 ("%s", "  Results from most recent analysis:\n") ;
	P2 ("    Cholesky flop count: %.5g\n", Common->fl) ;
	P2 ("    Nonzeros in L:       %.5g\n", Common->lnz) ;
    }

    P2 ("  memory blocks in use:    %8.0f\n", (double) (Common->malloc_count)) ;
    P2 ("  memory in use (MB):      %8.1f\n", 
	(double) (Common->memory_inuse) / 1048576.) ;
    P2 ("  peak memory usage (MB):  %8.1f\n", 
	(double) (Common->memory_usage) / 1048576.) ;

    /* ---------------------------------------------------------------------- */
    /* 主要控制参数及相关排序统计 */
    /* ---------------------------------------------------------------------- */

    P3 ("  maxrank:    update/downdate rank:   "ID"\n",
	    (Int) CORE(maxrank) (0, Common)) ;
    P3 ("  supernodal control: %d", Common->supernodal) ;
    P3 (" %g ", Common->supernodal_switch) ;
    if (Common->supernodal <= SPARSE_SIMPLICIAL)
    {
	P3 ("%s", "(always do simplicial)\n") ;
    }
    else if (Common->supernodal == SPARSE_AUTO)
    {
	P3 ("(supernodal if flops/lnz >= %g)\n", Common->supernodal_switch) ;
    }
    else
    {
	P3 ("%s", "(always do supernodal)\n") ;
    }

    nmethods = MIN (Common->nmethods, SPARSE_MAXMETHODS) ;
    nmethods = MAX (0, nmethods) ;

    if (nmethods > 0)
    {
	P3 ("%s", "  nmethods:   number of ordering methods to try: ") ;
	P3 (""ID"\n", nmethods) ;
        amd_backup = (nmethods > 1) ;
    }
    else
    {
	P3 ("%s", "  nmethods=0: default strategy:  Try user permutation if "
		"given.  Try AMD.\n") ;
	P3 ("%s", "    Select best ordering tried.\n") ;
	Common->method [0].ordering = SPARSE_GIVEN ;
	Common->method [1].ordering = SPARSE_AMD ;
	amd_backup = FALSE ;
	nmethods = 2 ;
    }

    for (i = 0 ; i < nmethods ; i++)
    {
	P3 ("    method "ID": ", i) ;
	ordering = Common->method [i].ordering ;
	fl = Common->method [i].fl ;
	lnz = Common->method [i].lnz ;
	switch (ordering)
	{

	    case SPARSE_NATURAL:
		P3 ("%s", "natural\n") ;
		break ;

	    case SPARSE_GIVEN:
		P3 ("%s", "user permutation (if given)\n") ;
		break ;

	    case SPARSE_AMD:
		P3 ("%s", "AMD (or COLAMD if factorizing AA')\n") ;
		amd_backup = FALSE ;
		break ;

	    case SPARSE_COLAMD:
		P3 ("%s", "AMD if factorizing A, COLAMD if factorizing AA')\n");
		amd_backup = FALSE ;
		break ;

	    default:
		P3 (ID, ordering) ;
		ERR ("unknown ordering method") ;
		break ;

	}

	if (!(ordering == SPARSE_NATURAL || ordering == SPARSE_GIVEN))
	{
	    if (Common->method [i].prune_dense < 0)
	    {
		P3 ("        prune_dense: for pruning dense nodes:   %s\n",
			" none pruned") ;
	    }
	    else
	    {
		P3 ("        prune_dense: for pruning dense nodes:   "
		    "%.5g\n",
		    Common->method [i].prune_dense) ;
		P3 ("        a dense node has degree "
			">= max(16,(%.5g)*sqrt(n))\n",
		    Common->method [i].prune_dense) ;
	    }
	}

	if (ordering == SPARSE_COLAMD)
	{
	    if (Common->method [i].prune_dense2 < 0)
	    {
		P3 ("        prune_dense2: for pruning dense rows for AA':"
			"  %s\n", " none pruned") ;
	    }
	    else
	    {
		P3 ("        prune_dense2: for pruning dense rows for AA':"
		    " %.5g\n", Common->method [i].prune_dense2) ;
		P3 ("        a dense row has degree "
			">= max(16,(%.5g)*sqrt(ncol))\n",
		    Common->method [i].prune_dense2) ;
	    }
	}

	if (fl  != EMPTY) P3 ("        flop count: %.5g\n", fl) ;
	if (lnz != EMPTY) P3 ("        nnz(L):     %.5g\n", lnz) ;
    }

    /* 备份AMD结果，如果有 */
    if (amd_backup)
    {
	P3 ("%s", "    backup method: ") ;
	P3 ("%s", "AMD (or COLAMD if factorizing AA')\n") ;
	fl = Common->method [nmethods].fl ;
	lnz = Common->method [nmethods].lnz ;
	if (fl  != EMPTY) P3 ("        AMD flop count: %.5g\n", fl) ;
	if (lnz != EMPTY) P3 ("        AMD nnz(L):     %.5g\n", lnz) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 神秘的控制参数 */
    /* ---------------------------------------------------------------------- */

    if (Common->final_asis)
    {
	P4 ("%s", "  final_asis: TRUE, leave as is\n") ;
    }
    else
    {
	P4 ("%s", "  final_asis: FALSE, convert when done\n") ;
	if (Common->final_super)
	{
	    P4 ("%s", "  final_super: TRUE, leave in supernodal form\n") ;
	}
	else
	{
	    P4 ("%s", "  final_super: FALSE, convert to simplicial form\n") ;
	}
	if (Common->final_ll)
	{
	    P4 ("%s", "  final_ll: TRUE, convert to LL' form\n") ;
	}
	else
	{
	    P4 ("%s", "  final_ll: FALSE, convert to LDL' form\n") ;
	}
	if (Common->final_pack)
	{
	    P4 ("%s", "  final_pack: TRUE, pack when done\n") ;
	}
	else
	{
	    P4 ("%s", "  final_pack: FALSE, do not pack when done\n") ;
	}
	if (Common->final_monotonic)
	{
	    P4 ("%s", "  final_monotonic: TRUE, ensure L is monotonic\n") ;
	}
	else
	{
	    P4 ("%s",
		"  final_monotonic: FALSE, do not ensure L is monotonic\n") ;
	}
	P4 ("  final_resymbol: remove zeros from amalgamation: %s\n",
		BOOLSTR (Common->final_resymbol)) ;
    }

    P4 ("  dbound:  LDL' diagonal threshold: % .5g\n    Entries with abs. value"
	    " less than dbound are replaced with +/- dbound.\n",
	    Common->dbound) ;

    P4 ("  grow0: memory reallocation: % .5g\n", Common->grow0) ;
    P4 ("  grow1: memory reallocation: % .5g\n", Common->grow1) ;
    P4 ("  grow2: memory reallocation: %g\n", (double) (Common->grow2)) ;

    P4 ("%s", "  nrelax, zrelax:  supernodal amalgamation rule:\n") ;
    P4 ("%s", "    s = # columns in two adjacent supernodes\n") ;
    P4 ("%s", "    z = % of zeros in new supernode if they are merged.\n") ;
    P4 ("%s", "    Two supernodes are merged if") ;
    P4 (" (s <= %g) or (no new zero entries) or\n",
	    (double) (Common->nrelax [0])) ;
    P4 ("    (s <= %g and ",  (double) (Common->nrelax [1])) ;
    P4 ("z < %.5g%%) or",      Common->zrelax [0] * 100) ;
    P4 (" (s <= %g and ",     (double) (Common->nrelax [2])) ;
    P4 ("z < %.5g%%) or",      Common->zrelax [1] * 100) ;
    P4 (" (z < %.5g%%)\n",     Common->zrelax [2] * 100) ;

    /* ---------------------------------------------------------------------- */
    /* 检查工作区 */
    /* ---------------------------------------------------------------------- */

    mark = Common->mark ;
    nrow = Common->nrow ;
    Flag = Common->Flag ;
    Head = Common->Head ;
    if (nrow > 0)
    {
	if (mark < 0 || Flag == NULL || Head == NULL)
	{
	    ERR ("workspace corrupted (Flag and/or Head missing)") ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    if (Flag [i] >= mark)
	    {
		ERR ("workspace corrupted (Flag)") ;
	    }
	}
	for (i = 0 ; i <= nrow ; i++)
	{
	    if (Head [i] != EMPTY)
	    {
		ERR ("workspace corrupted (Head)") ;
	    }
	}
    }
    xworksize = Common->xworksize ;
    Xwork = Common->Xwork ;
    if (xworksize > 0)
    {
	if (Xwork == NULL)
	{
	    ERR ("workspace corrupted (Xwork missing)") ;
	}
	for (i = 0 ; i < xworksize ; i++)
	{
	    if (Xwork [i] != 0.)
	    {
		ERR ("workspace corrupted (Xwork)") ;
	    }
	}
    }

    /* 工作空间和参数是有效的 */
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}


int CORE(check_common)
(
    sparse_common *Common
)
{
    return (check_common (0, NULL, Common)) ;
}


int CORE(print_common)
(
    /* ---- input ---- */
    const char *name,		/* 输出的Common对象的昵称 */
    /* --------------- */
    sparse_common *Common
)
{
    Int print = (Common == NULL) ? 3 : (Common->print) ;
    return (check_common (print, name, Common)) ;
}

/**
 * @brief 	打印CPU统计数据。如果没有安装计时器，
 * 			时间将被报告为零，但是这个函数仍然可以工作。
 * 
 */
// int CORE(cpu_stats)
// (
//     sparse_common *Common      /* input */
// )
// {
//     double cpu_time ;
//     int print ;

//     RETURN_IF_NULL_COMMON (FALSE) ;
//     print = Common->print ;

//     P2 ("%s", "\nHNUCHOL CPU statistics:\n") ;
//     P2 ("SYRK  CPU calls %12.0f", (double) Common->SPARSE_CPU_SYRK_CALLS) ;
//     P2 (" time %12.4e\n", Common->SPARSE_CPU_SYRK_TIME) ;
//     P2 ("GEMM  CPU calls %12.0f", (double) Common->SPARSE_CPU_GEMM_CALLS) ;
//     P2 (" time %12.4e\n", Common->SPARSE_CPU_GEMM_TIME) ;
//     P2 ("POTRF CPU calls %12.0f", (double) Common->SPARSE_CPU_POTRF_CALLS) ;
//     P2 (" time %12.4e\n", Common->SPARSE_CPU_POTRF_TIME) ;
//     P2 ("TRSM  CPU calls %12.0f", (double) Common->SPARSE_CPU_TRSM_CALLS) ;
//     P2 (" time %12.4e\n", Common->SPARSE_CPU_TRSM_TIME) ;

//     cpu_time = Common->SPARSE_CPU_SYRK_TIME + Common->SPARSE_CPU_TRSM_TIME +
//                Common->SPARSE_CPU_GEMM_TIME + Common->SPARSE_CPU_POTRF_TIME ;

//     P2 ("time in the BLAS: CPU %12.4e", cpu_time) ;

//     P2 ("assembly time %12.4e", Common->SPARSE_ASSEMBLE_TIME) ;
//     P2 ("  %12.4e\n", Common->SPARSE_ASSEMBLE_TIME2) ;
//     return (TRUE) ;
// }

/**
 * @brief 	确保一个面向列的稀疏矩阵是有效的，并可选地打印它。
 * 			返回对角线上的条目数，如果错误则返回-1。
 * 			工作空间：Iwork (nrow)
 * 
 * @param Wi 
 * @param print 
 * @param name 
 * @param A 
 * @param nnzdiag 
 * @param Common 
 * @return Sparse_long 
 */
static Sparse_long check_sparse
(
    Int *Wi,
    Int print,
    const char *name,
    sparse_csc *A,
    Sparse_long *nnzdiag,
    sparse_common *Common
)
{
    double *Ax, *Az ;
    Int *Ap, *Ai, *Anz ;
    Int nrow, ncol, nzmax, sorted, packed, j, p, pend, i, nz, ilast,
	init_print, dnz, count, xtype ;
    const char *type = "sparse" ;

    /* ---------------------------------------------------------------------- */
    /* 打印头信息 */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "HNUCHOL sparse:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (A == NULL)
    {
	ERR ("null") ;
    }

    nrow = A->nrow ;
    ncol = A->ncol ;
    nzmax = A->nzmax ;
    sorted = A->sorted ;
    packed = A->packed ;
    xtype = A->xtype ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    nz = CORE(nnz) (A, Common) ;

    P3 (" "ID"", nrow) ;
    P3 ("-by-"ID", ", ncol) ;
    P3 ("nz "ID",", nz) ;
    if (A->stype > 0)
    {
	P3 ("%s", " upper.") ;
    }
    else if (A->stype < 0)
    {
	P3 ("%s", " lower.") ;
    }
    else
    {
	P3 ("%s", " up/lo.") ;
    }

    P4 ("\n  nzmax "ID", ", nzmax) ;
    if (nz > nzmax)
    {
	ERR ("nzmax too small") ;
    }
    if (!sorted)
    {
	P4 ("%s", "un") ;
    }
    P4 ("%s", "sorted, ") ;
    if (!packed)
    {
	P4 ("%s", "un") ;
    }
    P4 ("%s", "packed, ") ;

    switch (A->itype)
    {
	case SPARSE_INT:     P4 ("%s", "\n  scalar types: int, ") ; break ;
	
	case SPARSE_LONG:    P4 ("%s", "\n  scalar types: Sparse_long, ");
        break ;
	default:	      ERR ("unknown itype") ;
    }

    switch (A->xtype)
    {
	case SPARSE_PATTERN: P4 ("%s", "pattern") ;	break ;
	case SPARSE_REAL:    P4 ("%s", "real") ;	break ;
	default:	      ERR ("unknown xtype") ;
    }

    switch (A->dtype)
    {
	case SPARSE_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	default:	      ERR ("unknown dtype") ;
    }

    if (A->itype != ITYPE || A->dtype != DTYPE)
    {
	ERR ("integer and real type must match routine") ;
    }

    if (A->stype && nrow != ncol)
    {
	ERR ("symmetric but not square") ;
    }

    /* 检查是否存在Ap, Ai, Anz, Ax和Az数组 */
    if (Ap == NULL)
    {
	ERR ("p array not present") ;
    }
    if (Ai == NULL)
    {
	ERR ("i array not present") ;
    }
    if (!packed && Anz == NULL)
    {
	ERR ("nz array not present") ;
    }
    if (xtype != SPARSE_PATTERN && Ax == NULL)
    {
	ERR ("x array not present") ;
    }
    /* 填充矩阵必须从Ap[0] = 0开始 */
    if (packed && Ap [0] != 0)
    {
	ERR ("p [0] must be zero") ;
    }
    if (packed && (Ap [ncol] < Ap [0] || Ap [ncol] > nzmax))
    {
	ERR ("p [ncol] invalid") ;
    }

    /* ---------------------------------------------------------------------- */
    /* 如果需要，分配工作空间 */
    /* ---------------------------------------------------------------------- */

    if (!sorted)
    {
	if (Wi == NULL)
	{
	    CORE(allocate_work) (0, nrow, 0, Common) ;
	    Wi = Common->Iwork ;	/* 大小为nrow, (i/i/l) */
	}
	if (Common->status < SPARSE_OK)
	{
	    return (FALSE) ;	    /* 内存溢出 */
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = EMPTY ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 检查并打印每一列 */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    dnz = 0 ;
    ETC_START (count, 8) ;

    for (j = 0 ; j < ncol ; j++)
    {
	ETC (j == ncol-1, count, 4) ;
	p = Ap [j] ;
	if (packed)
	{
	    pend = Ap [j+1] ;
	    nz = pend - p ;
	}
	else
	{
	    /* 注意Anz [j] < 0 视为0 */
	    nz = MAX (0, Anz [j]) ;
	    pend = p + nz ;
	}
	P4 ("  col "ID":", j) ;
	P4 (" nz "ID"", nz) ;
	P4 (" start "ID"", p) ;
	P4 (" end "ID"", pend) ;
	P4 ("%s", ":\n") ;
	if (p < 0 || pend > nzmax)
	{
	    ERR ("pointer invalid") ;
	}
	if (nz < 0 || nz > nrow)
	{
	    ERR ("nz invalid") ;
	}
	ilast = EMPTY ;

	for ( ; p < pend ; p++)
	{
	    ETC (j == ncol-1 && p >= pend-4, count, -1) ;
	    i = Ai [p] ;
	    P4 ("  "I8":", i) ;

	    print_value (print, xtype, Ax, Az, p, Common) ;

	    if (i == j)
	    {
		dnz++ ;
	    }
	    if (i < 0 || i >= nrow)
	    {
		ERR ("row index out of range") ;
	    }
	    if (sorted && i <= ilast)
	    {
		ERR ("row indices out of order") ;
	    }
	    if (!sorted && Wi [i] == j)
	    {
		ERR ("duplicate row index") ;
	    }
	    P4 ("%s", "\n") ;
	    ilast = i ;
	    if (!sorted)
	    {
		Wi [i] = j ;
	    }
	}
    }

    /* 矩阵有效 */
    P4 ("  nnz on diagonal: "ID"\n", dnz) ;
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    *nnzdiag = dnz ;
    return (TRUE) ;
}

/**
 * @brief 
 * 
 */
int CORE(check_sparse)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要检查的稀疏矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    Sparse_long nnzdiag ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_sparse (NULL, 0, NULL, A, &nnzdiag, Common)) ;
}

/**
 * @brief 
 * 
 */
int CORE(print_sparse)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要打印的稀疏矩阵 */
    const char *name,	/* 打印的稀疏矩阵名称 */
    /* --------------- */
    sparse_common *Common
)
{
    Sparse_long nnzdiag ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_sparse (NULL, Common->print, name, A, &nnzdiag, Common)) ;
}

/**
 * @brief 	确保稠密矩阵是有效的，并可选地打印它。
 * 
 * @param print 
 * @param name 
 * @param X 
 * @param Common 
 * @return int 
 */
static int check_dense
(
    Int print,
    const char *name,
    dense_array *X,
    sparse_common *Common
)
{
    double *Xx, *Xz ;
    Int i, j, d, nrow, ncol, nzmax, nz, init_print, count, xtype ;
    const char *type = "dense" ;

    /* ---------------------------------------------------------------------- */
    /* 打印头信息 */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "HNUCHOL dense:   ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (X == NULL)
    {
	ERR ("null") ;
    }

    nrow = X->nrow ;
    ncol = X->ncol ;
    nzmax = X->nzmax ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;
    xtype = X->xtype ;

    P3 (" "ID"", nrow) ;
    P3 ("-by-"ID", ", ncol) ;
    P4 ("\n  leading dimension "ID", ", d) ;
    P4 ("nzmax "ID", ", nzmax) ;
    if (d * ncol > nzmax)
    {
	ERR ("nzmax too small") ;
    }
    if (d < nrow)
    {
	ERR ("leading dimension must be >= # of rows") ;
    }
    if (Xx == NULL)
    {
	ERR ("null") ;
    }

    switch (X->xtype)
    {
	case SPARSE_PATTERN: ERR ("pattern unsupported") ;  break ;
	case SPARSE_REAL:    P4 ("%s", "real") ;	break ;
	default:	      ERR ("unknown xtype") ;
    }

    switch (X->dtype)
    {
	case SPARSE_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	default:	      ERR ("unknown dtype") ;
    }

    /* ---------------------------------------------------------------------- */
    /* 检查并打印每一项 */
    /* ---------------------------------------------------------------------- */

    if (print >= 4)
    {
	init_print = print ;
	ETC_START (count, 9) ;
	nz = nrow * ncol ;
	for (j = 0 ; j < ncol ; j++)
	{
	    ETC (j == ncol-1, count, 5) ;
	    P4 ("  col "ID":\n", j) ;
	    for (i = 0 ; i < nrow ; i++)
	    {
		ETC (j == ncol-1 && i >= nrow-4, count, -1) ;
		P4 ("  "I8":", i) ;

		print_value (print, xtype, Xx, Xz, i+j*d, Common) ;

		P4 ("%s", "\n") ;
	    }
	}
    }

    /* 稠密有效 */
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}

/**
 * @brief 
 * 
 */
int CORE(check_dense)
(
    /* ---- input ---- */
    dense_array *X,	/* 需要检查的稠密矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_dense (0, NULL, X, Common)) ;
}

/**
 * @brief 
 * 
 */
int CORE(print_dense)
(
    /* ---- input ---- */
    dense_array *X,	/* 需要打印的稠密矩阵 */
    const char *name,	/* 需要打印的稠密矩阵的名称 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_dense (Common->print, name, X, Common)) ;
}

/**
 * @brief 	确保S (0:len-1)是0:n-1的子集。
 * 			允许重复。S可以为空。负的len表示集合0:n-1。
 * 
 * @param S 
 * @param len 
 * @param n 
 * @param print 
 * @param name 
 * @param Common 
 * @return int 
 */
static int check_subset
(
    Int *S,
    Sparse_long len,
    size_t n,
    Int print,
    const char *name,
    sparse_common *Common
)
{
    Int i, k, init_print, count ;
    const char *type = "subset" ;

    init_print = print ;

    if (S == NULL)
    {
	/* len=0表示S =[]，负len表示S = 0:n-1 */
	len = (len < 0) ? (-1) : 0 ;
    }

    P4 ("%s", "\n") ;
    P3 ("%s", "HNUCHOL subset:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    P3 (" len: %ld ", len) ;
    if (len < 0)
    {
	P3 ("%s", "(denotes 0:n-1) ") ;
    }
    P3 ("n: "ID"", (Int) n) ;
    P4 ("%s", "\n") ;

    if (len <= 0 || S == NULL)
    {
	P3 ("%s", "  OK\n") ;
	P4 ("%s", "\n") ;
	return (TRUE) ;
    }

    if (print >= 4)
    {
	ETC_START (count, 8) ;
	for (k = 0 ; k < ((Int) len) ; k++)
	{
	    ETC (k == ((Int) len) - 4, count, -1) ;
	    i = S [k] ;
	    P4 ("  "I8":", k) ;
	    P4 (" "ID"\n", i) ;
	    if (i < 0 || i >= ((Int) n))
	    {
		ERR ("entry out range") ;
	    }
	}
    }
    else
    {
	for (k = 0 ; k < ((Int) len) ; k++)
	{
	    i = S [k] ;
	    if (i < 0 || i >= ((Int) n))
	    {
		ERR ("entry out range") ;
	    }
	}
    }
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}


int CORE(check_subset)
(
    /* ---- input ---- */
    Int *Set,			/* 集合[0:len-1]是0:n-1的子集。可以重复 */
    Sparse_long len, /* Set的大小(一个整数数组)，0:n - 1或者<0 */
    size_t n,			/* 0:n-1是有效的范围 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_subset (Set, len, n, 0, NULL, Common)) ;
}

/**
 * @brief 
 * 
 */
int CORE(print_subset)
(
    /* ---- input ---- */
    Int *Set,			/* 集合[0:len-1]是0:n-1的子集。可以重复 */
    Sparse_long len, /* 集合的大小(一个整数数组)，0:n - 1或者<0 */
    size_t n,			/* 0:n-1是有效的范围 */
    const char *name,	/* 打印的集合的名称 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_subset (Set, len, n, Common->print, name, Common)) ;
}

/* Ensure that Perm [0..len-1] is a permutation of a subset of 0:n-1.  Perm
 * may be NULL, which is interpreted as the identity permutation.  There can
 * be no duplicate entries (len must be <= n).
 */

/**
 * @brief 	确保Perm[0..len-1]是0:n-1子集的排列。
 * 			Perm可以为NULL，这被解释为恒等置换。
 * 			不能有重复的条目(len必须是<= n)。
 * 
 * @param Wi 
 * @param print 
 * @param name 
 * @param Perm 
 * @param len 
 * @param n 
 * @param Common 
 * @return int 
 */
static int check_perm
(
    Int *Wi,
    Int print,
    const char *name,
    Int *Perm,
    size_t len,
    size_t n,
    sparse_common *Common
)
{
    Int *Flag ;
    Int i, k, mark, init_print, count ;
    const char *type = "perm" ;

    /* ---------------------------------------------------------------------- */
    /* 检查是否花费O(1)时间 */
    /* ---------------------------------------------------------------------- */

    if (Perm == NULL || n == 0)
    {
	/* Perm是有效的隐式标识，或空 */
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 检查是否花费O(n)时间或需要内存分配 */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    ETC_START (count, 8) ;

    if (Wi == NULL && n <= Common->nrow)
    {
	/* 如果数组足够大，则使用Common->Flag数组 */
	mark = CORE(clear_flag) (Common) ;
	Flag = Common->Flag ;
	if (print >= 4)
	{
	    for (k = 0 ; k < ((Int) len) ; k++)
	    {
		ETC (k >= ((Int) len) - 4, count, -1) ;
		i = Perm [k] ;
		P4 ("  "I8":", k) ;
		P4 (""ID"\n", i) ;
		if (i < 0 || i >= ((Int) n) || Flag [i] == mark)
		{
		    CORE(clear_flag) (Common) ;
		    ERR ("invalid permutation") ;
		}
		Flag [i] = mark ;
	    }
	}
	else
	{
	    for (k = 0 ; k < ((Int) len) ; k++)
	    {
		i = Perm [k] ;
		if (i < 0 || i >= ((Int) n) || Flag [i] == mark)
		{
		    CORE(clear_flag) (Common) ;
		    ERR ("invalid permutation") ;
		}
		Flag [i] = mark ;
	    }
	}
	CORE(clear_flag) (Common) ;
    }
    else
    {
	if (Wi == NULL)
	{
	    /* 使用Common->Iwork替代前先进行初始化 */
	    CORE(allocate_work) (0, n, 0, Common) ;
	    Wi = Common->Iwork ;		    /* size n, (i/i/i) is OK */
	}
	if (Common->status < SPARSE_OK)
	{
	    return (FALSE) ;	    /* 内存溢出 */
	}
	for (i = 0 ; i < ((Int) n) ; i++)
	{
	    Wi [i] = FALSE ;
	}
	if (print >= 4)
	{
	    for (k = 0 ; k < ((Int) len) ; k++)
	    {
		ETC (k >= ((Int) len) - 4, count, -1) ;
		i = Perm [k] ;
		P4 ("  "I8":", k) ;
		P4 (""ID"\n", i) ;
		if (i < 0 || i >= ((Int) n) || Wi [i])
		{
		    ERR ("invalid permutation") ;
		}
		Wi [i] = TRUE ;
	    }
	}
	else
	{
	    for (k = 0 ; k < ((Int) len) ; k++)
	    {
		i = Perm [k] ;
		if (i < 0 || i >= ((Int) n) || Wi [i])
		{
		    ERR ("invalid permutation") ;
		}
		Wi [i] = TRUE ;
	    }
	}
    }

    /* perm有效 */
    return (TRUE) ;
}

/**
 * @brief 
 * 
 */
int CORE(check_perm)
(
    /* ---- input ---- */
    Int *Perm,		/* Perm [0:len-1]是0:n-1的子集 */
    size_t len,		/* Perm的大小 */
    size_t n,		/* 0:n-1是有效范围 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_perm (NULL, 0, NULL, Perm, len, n, Common)) ;
}

/**
 * @brief 
 * 
 */
int CORE(print_perm)
(
	/* ---- input ---- */
    Int *Perm,			/* Perm [0:len-1]是0:n-1的子集 */
    size_t len,			/* Perm的大小 */
    size_t n,			/* 0:n-1是有效范围 */
	const char *name,	/* 要打印的Perm的名称 */
    /* --------------- */
    sparse_common *Common    
)
{
    Int ok, print ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    print = Common->print ;
    P4 ("%s", "\n") ;
    P3 ("%s", "HNUCHOL perm:    ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }
    P3 (" len: "ID"", (Int) len) ;
    P3 (" n: "ID"", (Int) n) ;
    P4 ("%s", "\n") ;
    ok = check_perm (NULL, print, name, Perm, len, n, Common) ;
    if (ok)
    {
	P3 ("%s", "  OK\n") ;
	P4 ("%s", "\n") ;
    }
    return (ok) ;
}

/**
 * @brief 	确保父节点是节点0到n-1的有效消除树。如果j是树的根结点，那么父结点[j]为空(-1)。
 * 			不需要额外工作区
 * 
 * @param Parent 
 * @param n 
 * @param print 
 * @param name 
 * @param Common 
 * @return int 
 */
static int check_parent
(
    Int *Parent,
    size_t n,
    Int print,
    const char *name,
    sparse_common *Common
)
{
    Int j, p, init_print, count ;
    const char *type = "parent" ;

    init_print = print ;

    P4 ("%s", "\n") ;
    P3 ("%s", "HNUCHOL parent:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    P3 (" n: "ID"", (Int) n) ;
    P4 ("%s", "\n") ;

    if (Parent == NULL)
    {
	ERR ("null") ;
    }

    /* ---------------------------------------------------------------------- */
    /* 检查是否华为O(n)时间 */
    /* ---------------------------------------------------------------------- */

    ETC_START (count, 8) ;
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	ETC (j == ((Int) n) - 4, count, -1) ;
	p = Parent [j] ;
	P4 ("  "I8":", j) ;
	P4 (" "ID"\n", p) ;
	if (!(p == EMPTY || p > j))
	{
	    ERR ("invalid") ;
	}
    }
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}

/**
 * @brief 
 * 
 */
int CORE(check_parent)
(
    /* ---- input ---- */
    Int *Parent,	/* Parent [0:n-1] 是一个消去树 */
    size_t n,		/* 双亲节点数量 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_parent (Parent, n, 0, NULL, Common)) ;
}

/**
 * @brief 
 * 
 */
int CORE(print_parent)
(
    /* ---- input ---- */
    Int *Parent,		/* Parent [0:n-1] 是一棵消去树 */
    size_t n,			/* 双亲节点的数目 */
    const char *name,	/* 打印的双亲节点名称 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_parent (Parent, n, Common->print, name, Common)) ;
}

/**
 * @brief 
 * 
 * @param Wi 
 * @param print 
 * @param name 
 * @param L 
 * @param Common 
 * @return int 
 */
static int check_factor
(
    Int *Wi,
    Int print,
    const char *name,
    sparse_factor *L,
    sparse_common *Common
)
{
    double *Lx, *Lz ;
    Int *Lp, *Li, *Lnz, *Lnext, *Lprev, *Perm, *ColCount, *Lpi, *Lpx, *Super,
	*Ls ;
    Int n, nzmax, j, p, pend, i, nz, ordering, space, is_monotonic, minor,
	count, precise, init_print, ilast, lnz, head, tail, jprev, plast,
	jnext, examine_super, nsuper, s, k1, k2, psi, psend, psx, nsrow, nscol,
	ps2, psxend, ssize, xsize, maxcsize, maxesize, nsrow2, jj, ii, xtype ;
    Int check_Lpx ;
    const char *type = "factor" ;

    /* ---------------------------------------------------------------------- */
    /* 打印头信息 */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "HNUCHOL factor:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (L == NULL)
    {
	ERR ("null") ;
    }

    n = L->n ;
    minor = L->minor ;
    ordering = L->ordering ;
    xtype = L->xtype ;

    Perm = L->Perm ;
    ColCount = L->ColCount ;
    lnz = 0 ;

    precise = Common->precise ;

    P3 (" "ID"", n) ;
    P3 ("-by-"ID"", n) ;

    if (minor < n)
    {
	P3 (" not positive definite (column "ID")", minor) ;
    }

    switch (L->itype)
    {
	case SPARSE_INT:     P4 ("%s", "\n  scalar types: int, ") ; break ;
	
	case SPARSE_LONG:    P4 ("%s", "\n  scalar types: Sparse_long, ");
        break ;
	default:	      ERR ("unknown itype") ;
    }

    switch (L->xtype)
    {
	case SPARSE_PATTERN: P4 ("%s", "pattern") ;	break ;
	case SPARSE_REAL:    P4 ("%s", "real") ;	break ;
	default:	      ERR ("unknown xtype") ;
    }

    switch (L->dtype)
    {
	case SPARSE_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	default:	      ERR ("unknown dtype") ;
    }

    if (L->itype != ITYPE || L->dtype != DTYPE)
    {
	ERR ("integer and real type must match routine") ;
    }

    if (L->is_super)
    {
	P3 ("%s", "  supernodal") ;
    }
    else
    {
	P3 ("%s", "  simplicial") ;
    }

    if (L->is_ll)
    {
	P3 ("%s", ", LL'.") ;
    }
    else
    {
	P3 ("%s", ", LDL'.") ;
    }

    P4 ("%s", "\n  ordering method used: ") ;
    switch (L->ordering)
    {
	case SPARSE_POSTORDERED:P4 ("%s", "natural (postordered)") ;	 break ;
	case SPARSE_NATURAL:	P4 ("%s", "natural") ;			 break ;
	case SPARSE_GIVEN:	P4 ("%s", "user-provided") ;		 break ;
	case SPARSE_AMD:	P4 ("%s", "AMD") ;			 break ;
	case SPARSE_COLAMD:	P4 ("%s", "AMD for A, COLAMD for A*A'") ;break ;
	default:		ERR ("unknown ordering") ;
    }

    P4 ("%s", "\n") ;

    init_print = print ;

    /* ---------------------------------------------------------------------- */
    /* 检查L->Perm */
    /* ---------------------------------------------------------------------- */

    if (!check_perm (Wi, print, name, Perm, n, n, Common))
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 检查L->ColCount */
    /* ---------------------------------------------------------------------- */

    if (ColCount == NULL)
    {
	ERR ("ColCount vector invalid") ;
    }

    ETC_START (count, 8) ;
    for (j = 0 ; j < n ; j++)
    {
	ETC (j >= n-4, count, -1) ;
	P4 ("  col: "ID" ", j) ;
	nz = ColCount [j] ;
	P4 ("colcount: "ID"\n", nz) ;
	if (nz < 0 || nz > n-j)
	{
	    ERR ("ColCount out of range") ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 检查factor */
    /* ---------------------------------------------------------------------- */

    if (L->xtype == SPARSE_PATTERN && !(L->is_super))
    {

	/* ------------------------------------------------------------------ */
	/* 检查单节点符号分解 */
	/* ------------------------------------------------------------------ */

	;

    }
    else if (L->xtype != SPARSE_PATTERN && !(L->is_super))
    {

	/* ------------------------------------------------------------------ */
	/* 检查单节点数值分解 */
	/* ------------------------------------------------------------------ */

	P4 ("monotonic: %d\n", L->is_monotonic) ;
	nzmax = L->nzmax ;
	P3 (" nzmax "ID".", nzmax) ;
	P4 ("%s", "\n") ;
	Lp = L->p ;
	Li = L->i ;
	Lx = L->x ;
	Lz = L->z ;
	Lnz = L->nz ;
	Lnext = L->next ;
	Lprev = L->prev ;

	/* 检查Lp, Li, Lnz, Lnext, Lprev, Lx数组 */
	if (Lp == NULL)
	{
	    ERR ("p array not present") ;
	}
	if (Li == NULL)
	{
	    ERR ("i array not present") ;
	}
	if (Lnz == NULL)
	{
	    ERR ("nz array not present") ;
	}
	if (Lx == NULL)
	{
	    ERR ("x array not present") ;
	}
	if (Lnext == NULL)
	{
	    ERR ("next array not present") ;
	}
	if (Lprev == NULL)
	{
	    ERR ("prev array not present") ;
	}

	ETC_START (count, 8) ;

	/* 检查L的每一列 */
	plast = 0 ;
	is_monotonic = TRUE ;
	for (j = 0 ; j < n ; j++)
	{
	    ETC (j >= n-3, count, -1) ;
	    p = Lp [j] ;
	    nz = Lnz [j] ;
	    pend = p + nz ;
	    lnz += nz ;

	    P4 ("  col "ID":", j) ;
	    P4 (" nz "ID"", nz) ;
	    P4 (" start "ID"", p) ;
	    P4 (" end "ID"", pend) ;

	    if (Lnext [j] < 0 || Lnext [j] > n)
	    {
		ERR ("invalid link list")  ;
	    }
	    space = Lp [Lnext [j]] - p ;

	    P4 (" space "ID"", space) ;
	    P4 (" free "ID":\n", space - nz) ;

	    if (p < 0 || pend > nzmax || space < 1)
	    {
		ERR ("pointer invalid") ;
	    }
	    if (nz < 1 || nz > (n-j) || nz > space)
	    {
		ERR ("nz invalid") ;
	    }
	    ilast = j-1 ;

	    if (p < plast)
	    {
		is_monotonic = FALSE ;
	    }
	    plast = p ;

	    i = Li [p] ;
	    P4 ("  "I8":", i) ;
	    if (i != j)
	    {
		ERR ("diagonal missing") ;
	    }

	    print_value (print, xtype, Lx, Lz, p, Common) ;

	    P4 ("%s", "\n") ;
	    ilast = j ;
	    for (p++ ; p < pend ; p++)
	    {
		ETC_DISABLE (count) ;
		i = Li [p] ;
		P4 ("  "I8":", i) ;
		if (i < j || i >= n)
		{
		    ERR ("row index out of range") ;
		}
		if (i <= ilast)
		{
		    ERR ("row indices out of order") ;
		}

		print_value (print, xtype, Lx, Lz, p, Common) ;

		P4 ("%s", "\n") ;
		ilast = i ;
	    }
	}

	if (L->is_monotonic && !is_monotonic)
	{
	    ERR ("columns not monotonic") ;
	}

	/* 检查链表 */
	head = n+1 ;
	tail = n ;
	j = head ;
	jprev = EMPTY ;
	count = 0 ;
	for ( ; ; )
	{
	    if (j < 0 || j > n+1 || count > n+2)
	    {
		ERR ("invalid link list") ;
	    }
	    jnext = Lnext [j] ;
	    if (j >= 0 && j < n)
	    {
		if (jprev != Lprev [j])
		{
		    ERR ("invalid link list") ;
		}
	    }
	    count++ ;
	    if (j == tail)
	    {
		break ;
	    }
	    jprev = j ;
	    j = jnext ;
	}
	if (Lnext [tail] != EMPTY || count != n+2)
	{
	    ERR ("invalid link list") ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* 检查超节点数值分解和符号分解 */
	/* ------------------------------------------------------------------ */

	nsuper = L->nsuper ;
	ssize = L->ssize ;
	xsize = L->xsize ;
	maxcsize = L->maxcsize ;
	maxesize = L->maxesize ;
	Ls = L->s ;
	Lpi = L->pi ;
	Lpx = L->px ;
	Super = L->super ;
	Lx = L->x ;
	ETC_START (count, 8) ;

	P4 ("  ssize "ID" ", ssize) ;
	P4 ("xsize "ID" ", xsize) ;
	P4 ("maxcsize "ID" ", maxcsize) ;
	P4 ("maxesize "ID"\n", maxesize) ;

	if (Ls == NULL)
	{
	    ERR ("invalid: L->s missing") ;
	}
	if (Lpi == NULL)
	{
	    ERR ("invalid: L->pi missing") ;
	}
	if (Lpx == NULL)
	{
	    ERR ("invalid: L->px missing") ;
	}
	if (Super == NULL)
	{
	    ERR ("invalid: L->super missing") ;
	}

	if (L->xtype != SPARSE_PATTERN)
	{
	    /* 数值超节点分解 */
	    if (Lx == NULL)
	    {
		ERR ("invalid: L->x missing") ;
	    }
	    if (Ls [0] == EMPTY)
	    {
		ERR ("invalid: L->s not defined") ;
	    }
	    examine_super = TRUE ;
	}
	else
	{
	    /* 超节点符号分解，但仅当它已经被计算 */
	    examine_super = (Ls [0] != EMPTY) ;
	}

	if (examine_super)
	{
	    if (Lpi [0] != 0 || MAX (1, Lpi [nsuper]) != ssize)
	    {
		ERR ("invalid: L->pi invalid") ;
	    }

            /* 如果Lpx [0] == 123456, 超节点就出现了，
             * 但是Lpx [0...nsuper]未定义, 所以不需要检查.  
			 * 这用于SPQR */
            check_Lpx = (Lpx [0] != 123456) ;
	    if (check_Lpx && (Lpx [0] != 0 || MAX (1, Lpx[nsuper]) != xsize))
	    {
		ERR ("invalid: L->px invalid") ;
	    }

	    /* 检查并打印每个超节点 */
	    for (s = 0 ; s < nsuper ; s++)
	    {
		k1 = Super [s] ;
		k2 = Super [s+1] ;
		psi = Lpi [s] ;
		psend = Lpi [s+1] ;
		nsrow = psend - psi ;
		nscol = k2 - k1 ;
		nsrow2 = nsrow - nscol ;
		ps2 = psi + nscol ;

                if (check_Lpx)
                {
                    psx = Lpx [s] ;
                    psxend = Lpx [s+1] ;
                }

		ETC (s == nsuper-1, count, 4) ;

		P4 ("  supernode "ID", ", s) ;
		P4 ("col "ID" ", k1) ;
		P4 ("to "ID". ", k2-1) ;
		P4 ("nz in first col: "ID".\n", nsrow) ;

                if (check_Lpx)
                {
                    P4 ("  values start "ID", ", psx) ;
                    P4 ("end "ID"\n", psxend) ;
                }

		if (k1 > k2 || k1 < 0 || k2 > n || nsrow < nscol || nsrow2 < 0
                    || (check_Lpx && psxend - psx != nsrow * nscol))
		{
		    ERR ("invalid supernode") ;
		}

		lnz += nscol * nsrow - (nscol*nscol - nscol)/2 ;

		if (L->xtype != SPARSE_PATTERN)
		{
		    /* 打印超节点的每一列 */
		    for (jj = 0 ; jj < nscol ; jj++)
		    {
			ETC_ENABLE (s == nsuper-1 && jj >= nscol-3, count, -1) ;
			j = k1 + jj ;
			P4 ("  col "ID"\n", j) ;
			ilast = j ;
			i = Ls [psi + jj] ;
			P4 ("  "I8":", i) ;
			if (i != j)
			{
			    ERR ("row index invalid") ;
			}

			/* PRINTVALUE (Lx [psx + jj + jj*nsrow]) ; */
			print_value (print, xtype, Lx, NULL,
				psx + jj + jj*nsrow, Common) ;

			P4 ("%s", "\n") ;
			for (ii = jj + 1 ; ii < nsrow ; ii++)
			{
			    ETC_DISABLE (count) ;
			    i = Ls [psi + ii] ;
			    P4 ("  "I8":", i) ;
			    if (i <= ilast || i > n)
			    {
				ERR ("row index out of range") ;
			    }

			    /* PRINTVALUE (Lx [psx + ii + jj*nsrow]) ; */
			    print_value (print, xtype, Lx, NULL,
				    psx + ii + jj*nsrow, Common) ;

			    P4 ("%s", "\n") ;
			    ilast = i ;
			}
		    }
		}
		else
		{
		    /* 只需打印超节点的前导列 */
		    P4 ("  col "ID"\n", k1) ;
		    for (jj = 0 ; jj < nscol ; jj++)
		    {
			ETC (s == nsuper-1 && jj >= nscol-3, count, -1) ;
			j = k1 + jj ;
			i = Ls [psi + jj] ;
			P4 ("  "I8"", i) ;
			if (i != j)
			{
			    ERR ("row index invalid") ;
			}
			P4 ("%s", "\n") ;
		    }
		    ilast = j ;
		    for (ii = nscol ; ii < nsrow ; ii++)
		    {
			ETC_DISABLE (count) ;
			i = Ls [psi + ii] ;
			P4 ("  "I8"", i) ;
			if (i <= ilast || i > n)
			{
			    ERR ("row index out of range") ;
			}
			P4 ("%s", "\n") ;
			ilast = i ;
		    }
		}
	    }
	}
    }

    /* 分解因子有效 */
    P3 ("  nz "ID"", lnz) ;
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}

/**
 * @brief 
 * 
 */
int CORE(check_factor)
(
    /* ---- input ---- */
    sparse_factor *L,	/* 待检查的分解因子 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_factor (NULL, 0, NULL, L, Common)) ;
}

/**
 * @brief 
 * 
 */
int CORE(print_factor)
(
    /* ---- input ---- */
    sparse_factor *L,	/* 待打印的分解因子 */
    const char *name,	/* 打印的分解因子的名称 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_factor (NULL, Common->print, name, L, Common)) ;
}

/**
 * @brief 	确保一个三元组矩阵是有效的，并可选地打印它。
 * 
 * @param print 
 * @param name 
 * @param T 
 * @param Common 
 * @return int 
 */
static int check_triplet
(
    Int print,
    const char *name,
    sparse_triplet *T,
    sparse_common *Common
)
{
    double *Tx, *Tz ;
    Int *Ti, *Tj ;
    Int i, j, p, nrow, ncol, nzmax, nz, xtype, init_print, count ;
    const char *type = "triplet" ;

    /* ---------------------------------------------------------------------- */
    /* 打印头信息 */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "HNUCHOL triplet: ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (T == NULL)
    {
	ERR ("null") ;
    }

    nrow = T->nrow ;
    ncol = T->ncol ;
    nzmax = T->nzmax ;
    nz = T->nnz ;
    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    Tz = T->z ;
    xtype = T->xtype ;


    P3 (" "ID"", nrow) ;
    P3 ("-by-"ID", ", ncol) ;
    P3 ("nz "ID",", nz) ;
    if (T->stype > 0)
    {
	P3 ("%s", " upper.") ;
    }
    else if (T->stype < 0)
    {
	P3 ("%s", " lower.") ;
    }
    else
    {
	P3 ("%s", " up/lo.") ;
    }

    P4 ("\n  nzmax "ID", ", nzmax) ;
    if (nz > nzmax)
    {
	ERR ("nzmax too small") ;
    }

    switch (T->itype)
    {
	case SPARSE_INT:     P4 ("%s", "\n  scalar types: int, ") ; break ;
	
	case SPARSE_LONG:    P4 ("%s", "\n  scalar types: Sparse_long, ");
        break ;
	default:	      ERR ("unknown itype") ;
    }

    switch (T->xtype)
    {
	case SPARSE_PATTERN: P4 ("%s", "pattern") ;	break ;
	case SPARSE_REAL:    P4 ("%s", "real") ;	break ;
	default:	      ERR ("unknown xtype") ;
    }

    switch (T->dtype)
    {
	case SPARSE_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	default:	      ERR ("unknown dtype") ;
    }

    if (T->itype != ITYPE || T->dtype != DTYPE)
    {
	ERR ("integer and real type must match routine") ;
    }

    if (T->stype && nrow != ncol)
    {
	ERR ("symmetric but not square") ;
    }

    /* 检查Ti, Tj, Tx数组 */
    if (Tj == NULL)
    {
	ERR ("j array not present") ;
    }
    if (Ti == NULL)
    {
	ERR ("i array not present") ;
    }

    if (xtype != SPARSE_PATTERN && Tx == NULL)
    {
	ERR ("x array not present") ;
    }

    /* ---------------------------------------------------------------------- */
    /* 检查并打印每一项 */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    ETC_START (count, 8) ;

    for (p = 0 ; p < nz ; p++)
    {
	ETC (p >= nz-4, count, -1) ;
	i = Ti [p] ;
	P4 ("  "I8":", p) ;
	P4 (" "I_8"", i) ;
	if (i < 0 || i >= nrow)
	{
	    ERR ("row index out of range") ;
	}
	j = Tj [p] ;
	P4 (" "I_8"", j) ;
	if (j < 0 || j >= ncol)
	{
	    ERR ("column index out of range") ;
	}

	print_value (print, xtype, Tx, Tz, p, Common) ;

	P4 ("%s", "\n") ;
    }

    /* 三元组矩阵有效 */
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}


/**
 * @brief 
 * 
 */
int CORE(check_triplet)
(
    /* ---- input ---- */
    sparse_triplet *T,	/* 待检查的三元组矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_triplet (0, NULL, T, Common)) ;
}

/**
 * @brief 
 * 
 */
int CORE(print_triplet)
(
    /* ---- input ---- */
    sparse_triplet *T,	/* 待打印的三元组矩阵 */
    const char *name,	/* 待打印的三元组矩阵的名称 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;
    return (check_triplet (Common->print, name, T, Common)) ;
}

#endif
