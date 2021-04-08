/**
 * @file SparseCore_common.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_internal.h"
#include "SparseCore.h"

/**
 * @brief   初始化常见的默认参数和统计信息。将工作区指针设置为空。
 */
int CORE(start)
(
    sparse_common *Common
)
{
    int k ;

    if (Common == NULL)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 用户错误处理程序 */
    /* ---------------------------------------------------------------------- */

    Common->error_handler = NULL ;

    /* ---------------------------------------------------------------------- */
    /* 整型和数值类型 */
    /* ---------------------------------------------------------------------- */

    Common->itype = ITYPE ;
    Common->dtype = DTYPE ;

    /* ---------------------------------------------------------------------- */
    /* 默认的控制参数 */
    /* ---------------------------------------------------------------------- */

    CORE(defaults) (Common) ;
    Common->try_catch = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* 工作空间 */
    /* ---------------------------------------------------------------------- */

    /* 此代码假设公共的工作区没有初始化。如果是，则会发生内存泄漏，因为指针被NULL覆盖了。 */

    Common->nrow = 0 ;
    Common->mark = EMPTY ;
    Common->xworksize = 0 ;
    Common->iworksize = 0 ;
    Common->Flag = NULL ;
    Common->Head = NULL ;
    Common->Iwork = NULL ;
    Common->Xwork = NULL ;
    Common->no_workspace_reallocate = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* 统计数据 */
    /* ---------------------------------------------------------------------- */

    /* fl和lnz在SparseCore_analyze和SparseCore_rowcolcounts中计算 */
    Common->fl = EMPTY ;
    Common->lnz = EMPTY ;


    /* 所有例程都使用status作为错误报告代码 */
    Common->status = SPARSE_OK ;

    Common->malloc_count = 0 ;	/* #调用malloc减去#调用free */
    Common->memory_usage = 0 ;	/* 内存峰值使用(以字节为单位) */
    Common->memory_inuse = 0 ;	/* 当前正在使用的内存(以字节为单位) */

    Common->nrealloc_col = 0 ;
    Common->nrealloc_factor = 0 ;
    Common->ndbounds_hit = 0 ;
    Common->rowfacfl = 0 ;
    Common->aatfl = EMPTY ;

    Common->blas_ok = TRUE ;    /* 如果发生BLAS int溢出，则为false */

    /* ---------------------------------------------------------------------- */
    /* 默认的HnuSparseQR旋钮和统计 */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < 10 ; k++) Common->SPQR_istat [k] = 0 ;

    Common->SPQR_flopcount_bound = 0 ;   /* 浮点运算次数的上限 */
    Common->SPQR_tol_used = 0 ;          /* 公差使用 */
    Common->SPQR_norm_E_fro = 0 ;        /* 删除条目的范数 */

    Common->SPQR_grain = 1 ;    /* 默认没有英特尔TBB多任务处理 */
    Common->SPQR_small = 1e6 ;  /* 为TBB设定最小任务大小 */
    Common->SPQR_shrink = 1 ;   /* 控制SPQR shrink回收 */

    Common->SPQR_flopcount = 0 ;         /* SPQR的浮点计算次数 */

    Common->thread_pool_time = 0;
    return (TRUE) ;
}


/**
 * @brief   设置常见的默认参数，除了函数指针。
 *          不需要额外工作空间
 */
int CORE(defaults)
(
    sparse_common *Common
)
{
    Int i ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    /* ---------------------------------------------------------------------- */
    /* 默认的控制参数 */
    /* ---------------------------------------------------------------------- */

    Common->dbound = 0.0 ;
    Common->grow0 = 1.2 ;
    Common->grow1 = 1.2 ;
    Common->grow2 = 5 ;
    Common->maxrank = 8 ;

    Common->final_asis = TRUE ;
    Common->final_super = TRUE ;
    Common->final_ll = FALSE ;
    Common->final_pack = TRUE ;
    Common->final_monotonic = TRUE ;
    Common->final_resymbol = FALSE ;

    /* 用简单因子分解处理flop/nnz(L)<40,否则为超节点 */
    Common->supernodal = SPARSE_AUTO ;
    Common->supernodal_switch = 40 ;

    Common->nrelax [0] = 4 ;
    Common->nrelax [1] = 16 ;
    Common->nrelax [2] = 48 ;
    Common->zrelax [0] = 0.8 ;
    Common->zrelax [1] = 0.1 ;
    Common->zrelax [2] = 0.05 ;

    Common->prefer_upper = TRUE ;
    Common->prefer_binary = FALSE ;
    Common->quick_return_if_not_posdef = FALSE ;

    Common->print = 3 ;
    Common->precise = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* 默认的排序方法 */
    /* ---------------------------------------------------------------------- */

    /* 默认策略.
     */
    Common->nmethods = 0 ;		/* 使用默认策略 */

    Common->current = 0 ;	/* 正在使用的方法 */
    Common->selected = 0 ;	/* 选择的最优方法 */

    /* 首先，用默认参数填充每个方法 */
    for (i = 0 ; i <= SPARSE_MAXMETHODS ; i++)
    {
	/* HNUCHOL的默认方法是AMD为A或AA' */
	Common->method [i].ordering = SPARSE_AMD ;

	/* HNUCHOL嵌套分离和敏度参数 */
	Common->method [i].prune_dense = 10.0 ;	/* 稠密row/col控制 */

	/* 敏度参数 (AMD, COLAMD, SYMAMD, CAMD, CCOLAMD, CSYMAMD)*/
	Common->method [i].prune_dense2 = -1 ;	/* COLAMD 稠密行控制 */
	Common->method [i].aggressive = TRUE ;	/* aggressive absorption */
	Common->method [i].order_for_lu = FALSE ;/* Cholesky的顺序，而不是LU */

	/* 还没有计算每种方法的统计信息 */
	Common->method [i].fl = EMPTY ;
	Common->method [i].lnz = EMPTY ;
    }

    Common->postorder = TRUE ;	/* 按照加权后序排序 */

    /* 接下来，定义一些方法。前五个使用默认参数。 */
    Common->method [0].ordering = SPARSE_GIVEN ;   /* 如果UserPerm为NULL，则跳过 */
    Common->method [1].ordering = SPARSE_AMD ;
    Common->method [4].ordering = SPARSE_NATURAL ;
    /* COLAMD for A*A', AMD for A */
    Common->method [8].ordering = SPARSE_COLAMD ;

    return (TRUE) ;
}


/**
 * @brief 对HNUCHOL的最后一次调用必须是SparseCore_finish。
 * 
 */
int CORE(finish)
(
    sparse_common *Common
)
{
    return (CORE(free_work) (Common)) ;
}

/**
 * @brief   为HNUCHOL例程分配和初始化工作区，或者增加已经分配的工作区的大小。
 *          如果已经分配了足够的工作区，那么什么也不会发生。
 * 
 *          工作空间：Flag (nrow), Head (nrow+1), Iwork (iworksize), Xwork (xworksize)
 */
int CORE(allocate_work)
(
    /* ---- input ---- */
    size_t nrow,	    /* 矩阵A的行数 */
    size_t iworksize,	/* Iwork的大小 */
    size_t xworksize,	/* Xwork的大小 */
    /* --------------- */
    sparse_common *Common
)
{
    double *W ;
    Int *Head ;
    Int i ;
    size_t nrow1 ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 指定Flag (nrow)和Head (nrow+1) */
    /* ---------------------------------------------------------------------- */

    nrow = MAX (1, nrow) ;

    /* nrow1 = nrow + 1 */
    nrow1 = CORE(add_size_t) (nrow, 1, &ok) ;
    if (!ok)
    {
	/* nrow+1导致size_t溢出 */
	Common->status = SPARSE_TOO_LARGE ;
	CORE(free_work) (Common) ;
	return (FALSE) ;
    }

    if (nrow > Common->nrow)
    {

	if (Common->no_workspace_reallocate)
	{
	    /* HNUCHOL不允许在这里更改工作区 */
	    Common->status = SPARSE_INVALID ;
	    return (FALSE) ;
	}

	/* 释放旧的工作空间(如果有的话)并分配新的空间 */
	Common->Flag = CORE(free) (Common->nrow,  sizeof (Int), Common->Flag,
		Common) ;
	Common->Head = CORE(free) (Common->nrow+1,sizeof (Int), Common->Head,
		Common) ;
	Common->Flag = CORE(malloc) (nrow,   sizeof (Int), Common) ;
	Common->Head = CORE(malloc) (nrow1, sizeof (Int), Common) ;

	/* 记录Flag和Head的新的大小 */
	Common->nrow = nrow ;

	if (Common->status < SPARSE_OK)
	{
	    CORE(free_work) (Common) ;
	    return (FALSE) ;
	}

	/* 初始化Flag和Head */
	Common->mark = EMPTY ;
	CORE(clear_flag) (Common) ;
	Head = Common->Head ;
	for (i = 0 ; i <= (Int) (nrow) ; i++)
	{
	    Head [i] = EMPTY ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 指定Iwork (iworksize) */
    /* ---------------------------------------------------------------------- */

    iworksize = MAX (1, iworksize) ;
    if (iworksize > Common->iworksize)
    {

	if (Common->no_workspace_reallocate)
	{
	    /* HNUCHOL不允许在这里更改工作区 */
	    Common->status = SPARSE_INVALID ;
	    return (FALSE) ;
	}

	/* 释放旧的工作空间(如果有的话)并分配新的空间。在SparseCore_malloc中安全检测到整数溢出 */
	CORE(free) (Common->iworksize, sizeof (Int), Common->Iwork, Common) ;
	Common->Iwork = CORE(malloc) (iworksize, sizeof (Int), Common) ;

	/* 记录Iwork新的大小 */
	Common->iworksize = iworksize ;

	if (Common->status < SPARSE_OK)
	{
	    CORE(free_work) (Common) ;
	    return (FALSE) ;
	}

	/* 注意，Iwork不需要初始化 */
    }

    /* ---------------------------------------------------------------------- */
    /* 指定Xwork (xworksize)并且设定它为((double) 0.) */
    /* ---------------------------------------------------------------------- */

    /* 确保xworksize>= 1 */
    xworksize = MAX (1, xworksize) ;
    if (xworksize > Common->xworksize)
    {

	if (Common->no_workspace_reallocate)
	{
	    /* HNUCHOL不允许在这里更改工作区 */
	    Common->status = SPARSE_INVALID ;
	    return (FALSE) ;
	}

	/* 释放旧的工作空间(如果有的话)并分配新的空间 */
	CORE(free) (Common->xworksize, sizeof (double), Common->Xwork,
		Common) ;
	Common->Xwork = CORE(malloc) (xworksize, sizeof (double), Common) ;

	/* 记录Xwork新的大小 */
	Common->xworksize = xworksize ;

	if (Common->status < SPARSE_OK)
	{
	    CORE(free_work) (Common) ;
	    return (FALSE) ;
	}

	/* 初始化 Xwork */
	W = Common->Xwork ;
	for (i = 0 ; i < (Int) xworksize ; i++)
	{
	    W [i] = 0. ;
	}
    }

    return (TRUE) ;
}


/**
 * @brief   释放HNUCHOL工作区。
 *          工作区：共同释放所有工作区
 */
int CORE(free_work)
(
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->Flag  = CORE(free) (Common->nrow, sizeof (Int),
	    Common->Flag, Common) ;
    Common->Head  = CORE(free) (Common->nrow+1, sizeof (Int),
	    Common->Head, Common) ;
    Common->Iwork = CORE(free) (Common->iworksize, sizeof (Int),
	    Common->Iwork, Common) ;
    Common->Xwork = CORE(free) (Common->xworksize, sizeof (double),
	    Common->Xwork, Common) ;
    Common->nrow = 0 ;
    Common->iworksize = 0 ;
    Common->xworksize = 0 ;

    return (TRUE) ;
}


/**
 * @brief   递增mark以确保Flag [0..nrow-1] < mark。
 *          如果发生整数溢出，或者mark最初为负，则重置整个数组。
 *          这不是一个错误条件，而是标志工作区的一个预期功能。
 *          工作区：Flag (nrow)。如果nrow为零，则不修改标记。
 */
Sparse_long CORE(clear_flag)
(
    sparse_common *Common
)
{
    Int i, nrow, *Flag ;

    RETURN_IF_NULL_COMMON (-1) ;

    Common->mark++ ;
    if (Common->mark <= 0)
    {
	nrow = Common->nrow ;
	Flag = Common->Flag ;
	for (i = 0 ; i < nrow ; i++)
	{
	    Flag [i] = EMPTY ;
	}
	Common->mark = 0 ;
    }
    return (Common->mark) ;
}


/**
 * @brief   求公共最大值的有效值。如果错误返回0，如果成功返回2、4或8。
 * 
 * @return size_t   返回Common->maxrank的验证值
 */
size_t CORE(maxrank)	
(
    /* ---- input ---- */
    size_t n,		/* A和L有n行 */
    /* --------------- */
    sparse_common *Common
)
{
    size_t maxrank ;
    RETURN_IF_NULL_COMMON (0) ;
    maxrank = Common->maxrank ;
    if (n > 0)
    {
	/* 确保maxrank*n*sizeof(double)不会导致整数溢出。如果n太大，
     * 以至于2*n*sizeof(double)导致整数溢出(如果Int是32位，则n = 268,435,455)，
     * 那么maxrank将为0或1，但maxrank将被设置为2。2*n不会导致整数溢出，
     * HNUCHOL将耗尽内存或在其他地方安全检测整数溢出。
	 */
	maxrank = MIN (maxrank, Size_max / (n * sizeof (double))) ;
    }
    if (maxrank <= 2)
    {
	maxrank = 2 ;
    }
    else if (maxrank <= 4)
    {
	maxrank = 4 ;
    }
    else
    {
	maxrank = 8 ;
    }
    return (maxrank) ;
}


/**
 * @brief   确保对角项D(j,j)的绝对值大于Common->dbound。这个例程不是用来让用户调用的。
 *          它被各种LDL的分解和 update/downdate 例程所使用。Common->dbound的默认值是零，
 *          在这种情况下，根本不调用这个例程。D(j,j)为NaN时不变。
 *          如果Common->dbound是NaN, HNUCHOL不调用这个例程。
 * @return double   返回修改后的对角项D
 */
double CORE(dbound)	
(
    /* ---- input ---- */
    double dj,		    /* D的对角项，用于LDL的分解 */
    /* --------------- */
    sparse_common *Common
)
{
    double dbound ;
    RETURN_IF_NULL_COMMON (0) ;
    if (!IS_NAN (dj))
    {
	dbound = Common->dbound ;
	if (dj < 0)
	{
	    if (dj > -dbound)
	    {
		dj = -dbound ;
		Common->ndbounds_hit++ ;
		if (Common->status == SPARSE_OK)
		{
		    ERROR (SPARSE_DSMALL, "diagonal below threshold") ;
		}
	    }
	}
	else
	{
	    if (dj < dbound)
	    {
		dj = dbound ;
		Common->ndbounds_hit++ ;
		if (Common->status == SPARSE_OK)
		{
		    ERROR (SPARSE_DSMALL, "diagonal below threshold") ;
		}
	    }
	}
    }
    return (dj) ;
}

/**
 * @brief 用于用qsort排序子代超节点
 * 
 */
int CORE(score_comp) (struct SparseCore_descendant_score_t *i, 
			       struct SparseCore_descendant_score_t *j)
{
  if ((*i).score < (*j).score)
    {
	return (1) ;
    }
    else
    {
	return (-1) ;
    }
}

/**
 * @brief   安全计算a+b，并检查整数溢出。如果发生溢出，则返回0并将OK设置为FALSE。
 *          如果OK在输入时为FALSE，也返回0。
 */
size_t CORE(add_size_t) (size_t a, size_t b, int *ok)
{
    size_t s = a + b ;
    (*ok) = (*ok) && (s >= a) ;
    return ((*ok) ? s : 0) ;
}

/**
 * @brief   安全地计算a*k，其中k应该很小，并检查整数溢出。如果发生溢出，
 *          则返回0并将OK设置为FALSE。如果OK在输入时为FALSE，也返回0。
 * 
 */
size_t CORE(mult_size_t) (size_t a, size_t k, int *ok)
{
    size_t p = 0, s ;
    while (*ok)
    {
	if (k % 2)
	{
	    p = p + a ;
	    (*ok) = (*ok) && (p >= a) ;
	}
	k = k / 2 ;
	if (!k) return (p) ;
	s = a + a ;
	(*ok) = (*ok) && (s >= a) ;
	a = s ;
    }
    return (0) ;
}


/**
 * @brief malloc例程的包装器。分配大小为MAX(1,n)*size的空间，其中size通常为sizeof(…)。
 * 
 * @return void*    返回指向新malloc的块的指针
 */
void *CORE(malloc)	
(
    /* ---- input ---- */
    size_t n,		/* 条目的数量 */
    size_t size,	/* 每个条目的大小 */
    /* --------------- */
    sparse_common *Common
)
{
    void *p ;
    size_t s ;
    /*
    int ok = TRUE ;
    */

    RETURN_IF_NULL_COMMON (NULL) ;
    if (size == 0)
    {
	ERROR (SPARSE_INVALID, "sizeof(item) must be > 0")  ;
	p = NULL ;
    }
    else if (n >= (Size_max / size) || n >= Int_max)
    {
	/* 对象太大，无法在不导致整型溢出的情况下分配 */
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	p = NULL ;
    }
    else
    {
	/* 调用malloc或它的等效程序 */
	p = SparseBase_malloc (n, size) ;

	if (p == NULL)
	{
	    /* 失败：内存溢出 */
	    ERROR (SPARSE_OUT_OF_MEMORY, "out of memory") ;
	}
	else
	{
	    /* 成功:增加已分配对象的计数 */
	    Common->malloc_count++ ;
	    Common->memory_inuse += (n * size) ;
	    Common->memory_usage =
		MAX (Common->memory_usage, Common->memory_inuse) ;
	}
    }
    return (p) ;
}

/**
 * @brief   包装自由例行公事。返回NULL，它可以被分配给被释放的指针，如:
 * 
 * @return void*    总是返回空
 */
void *CORE(free)	
(
    /* ---- input ---- */
    size_t n,		/* 条目数 */
    size_t size,	/* 每个条目大小 */
    /* ---- in/out --- */
    void *p,		/* 要释放的内存块 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (NULL) ;
    if (p != NULL)
    {
	/* 只有当指针不为空时才释放对象 */
	/* 调用free，或者它的等价 */
	SparseBase_free (p) ;

	Common->malloc_count-- ;
	Common->memory_inuse -= (n * size) ;
    }
    /* 返回NULL，调用者应该将这个赋值给p。这样就避免了两次释放相同的指针。 */
    return (NULL) ;
}

/* Wrapper around calloc routine.
 *
 * Uses a pointer to the calloc routine (or its equivalent) defined in Common.
 * This routine is identical to malloc, except that it zeros the newly allocated
 * block to zero.
 */

/**
 * @brief   封装calloc例程。
 * 
 *          使用一个指向通用定义的calloc例程(或等效的例程)的指针。
 *          这个例程与malloc相同，只是它将新分配的块归零。
 * 
 * @return void*    返回指向新调用的块的指针
 */
void *CORE(calloc)	
(
    /* ---- input ---- */
    size_t n,		/* 条目数 */
    size_t size,	/* 每个条目大小 */
    /* --------------- */
    sparse_common *Common
)
{
    void *p ;

    RETURN_IF_NULL_COMMON (NULL) ;
    if (size == 0)
    {
	ERROR (SPARSE_INVALID, "sizeof(item) must be > 0") ;
	p = NULL ;
    }
    else if (n >= (Size_max / size) || n >= Int_max)
    {
	/* 对象太大，无法在不导致整型溢出的情况下分配 */
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	p = NULL ;
    }
    else
    {
	/* 调用calloc或它的等效项 */
	p = SparseBase_calloc (n, size) ;

	if (p == NULL)
	{
	    /* 失败：内存溢出 */
	    ERROR (SPARSE_OUT_OF_MEMORY, "out of memory") ;
	}
	else
	{
	    /* 成功:增加已分配对象的计数 */
	    Common->malloc_count++ ;
	    Common->memory_inuse += (n * size) ;
	    Common->memory_usage =
		MAX (Common->memory_usage, Common->memory_inuse) ;
	}
    }
    return (p) ;
}

/* Wrapper around realloc routine.  Given a pointer p to a block of size
 * (*n)*size memory, it changes the size of the block pointed to by p to be
 * MAX(1,nnew)*size in size.  It may return a pointer different than p.  This
 * should be used as (for a pointer to int):
 */

/**
 * @brief   围绕realloc例程的包装。给定一个指针p指向一个大小为(*n)*size的内存块，
 *          它将p指向的内存块的大小更改为MAX(1,nnew)*size。它可能返回一个不同于p的指针。
 *          这应该被用作(一个指向int的指针):
 * 
 * @return void*    返回指向重新分配块的指针
 */
void *CORE(realloc)	
(
    /* ---- input ---- */
    size_t nnew,	/* 重新分配块中请求的项目# */
    size_t size,	/* 每个条目大小 */
    /* ---- in/out --- */
    void *p,		/* 块内存到realloc */
    size_t *n,		/* 输入当前的大小，如果成功，输出为nnew */
    /* --------------- */
    sparse_common *Common
)
{
    size_t nold = (*n) ;
    void *pnew ;
    size_t s ;
    int ok = TRUE ;

    RETURN_IF_NULL_COMMON (NULL) ;
    if (size == 0)
    {
	ERROR (SPARSE_INVALID, "sizeof(item) must be > 0") ;
	p = NULL ;
    }
    else if (p == NULL)
    {
	/* 正在分配一个新对象. */
	p = CORE(malloc) (nnew, size, Common) ;
	*n = (p == NULL) ? 0 : nnew ;
    }
    else if (nold == nnew)
    {
    }
    else if (nnew >= (Size_max / size) || nnew >= Int_max)
    {
	/* 失败：nnew太大 */
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
    }
    else
    {
	/* 对象存在，并且正在改变为其他的非零大小。 */
	/* 调用realloc或其等效函数 */
        pnew = SparseBase_realloc (nnew, nold, size, p, &ok) ;
        if (ok)
        {
	    /* 成功:返回修改后的p并改变块的大小 */
	    p = pnew ;
	    *n = nnew ;
	    Common->memory_inuse += ((nnew-nold) * size) ;
	}
        else
        {
            /* 增加块的大小失败了。 */
            ERROR (SPARSE_OUT_OF_MEMORY, "out of memory") ;
        }

	Common->memory_usage = MAX (Common->memory_usage, Common->memory_inuse);
    }

    return (p) ;
}

/**
 * @brief   重新分配多个内存块，所有块大小相同(最多两个整数和两个实际块)。
 *          要么重新分配成功，要么以原始大小返回(如果原始大小为零，则释放它们)。
 *          nnew块的大小为1或更多。
 */
int CORE(realloc_multiple)
(
    /* ---- input ---- */
    size_t nnew,	/* 重新分配块中请求的项目#s */
    int nint,		/* int/Sparse_long块的数目 */
    int xtype,		/* SPARSE_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
    /* ---- in/out --- */
    void **Iblock,	/* int/Sparse_long的块 */
    void **Jblock,	/* int/Sparse_long的块 */
    void **Xblock,	/* double块 */
    void **Zblock,	/* zomplex case only: double block */
    size_t *nold_p,	/* 输入时当前的块I,J,X,Z的大小，如果成功，输入为nnew */
    /* --------------- */
    sparse_common *Common
)
{
    double *xx, *zz ;
    size_t i, j, x, z, nold ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (xtype < SPARSE_PATTERN || xtype > SPARSE_REAL)
    {
	ERROR (SPARSE_INVALID, "invalid xtype") ;
	return (FALSE) ;
    }

    nold = *nold_p ;

    if (nint < 1 && xtype == SPARSE_PATTERN)
    {
	return (TRUE) ;
    }

    i = nold ;
    j = nold ;
    x = nold ;
    z = nold ;

    if (nint > 0)
    {
	*Iblock = CORE(realloc) (nnew, sizeof (Int), *Iblock, &i, Common) ;
    }
    if (nint > 1)
    {
	*Jblock = CORE(realloc) (nnew, sizeof (Int), *Jblock, &j, Common) ;
    }

    switch (xtype)
    {
	case SPARSE_REAL:
	    *Xblock = CORE(realloc) (nnew, sizeof (double), *Xblock, &x,
                    Common) ;
	    break ;
    }

    if (Common->status < SPARSE_OK)
    {
	/* 一个或多个realloc失败。将所有的大小调整为nold。 */

	if (nold == 0)
	{

	    if (nint > 0)
	    {
		*Iblock = CORE(free) (i, sizeof (Int), *Iblock, Common) ;
	    }
	    if (nint > 1)
	    {
		*Jblock = CORE(free) (j, sizeof (Int), *Jblock, Common) ;
	    }

	    switch (xtype)
	    {
		case SPARSE_REAL:
		    *Xblock = CORE(free) (x, sizeof (double), *Xblock,
                            Common) ;
		    break ;
	    }

	}
	else
	{
	    if (nint > 0)
	    {
		*Iblock = CORE(realloc) (nold, sizeof (Int), *Iblock, &i,
                            Common) ;
	    }
	    if (nint > 1)
	    {
		*Jblock = CORE(realloc) (nold, sizeof (Int), *Jblock, &j,
                            Common) ;
	    }

	    switch (xtype)
	    {
		case SPARSE_REAL:
		    *Xblock = CORE(realloc) (nold, sizeof (double),
                            *Xblock, &x, Common) ;
		    break ;
	    }

	}

	return (FALSE) ;
    }

    if (nold == 0)
    {
	/* 分配了新的空间。清除第一个条目，这样valgrind在访问change_complexity时有问题 */
	xx = *Xblock ;
	zz = *Zblock ;
	switch (xtype)
	{
	    case SPARSE_REAL:
		xx [0] = 0 ;
		break ;
	}
    }

    /* 所有的realloc成功，改变大小以反映realloc的大小。 */
    *nold_p = nnew ;
    return (TRUE) ;
}

/* An error has occurred.  Set the status, optionally print an error message,
 * and call the user error-handling routine (if it exists).  If
 * Common->try_catch is TRUE, then HNUCHOL is inside a try/catch block.
 * The status is set, but no message is printed and the user error handler
 * is not called.  This is not (yet) an error, since HNUCHOL may recover.
 *
 * In the current version, this try/catch mechanism is used internally only in
 * SparseCore_analyze, which tries multiple ordering methods and picks the best
 * one.  If one or more ordering method fails, it keeps going.  Only one
 * ordering needs to succeed for SparseCore_analyze to succeed.
 */

/**
 * @brief   发生错误。设置状态，可选地打印错误消息，并调用用户错误处理例程(如果存在)。
 *          如果Common->try_catch为真，则HNUCHOL在一个try/catch块中。
 *          设置了状态，但不打印任何消息，也不调用用户错误处理程序。
 *          这还不是错误，因为HNUCHOL可以恢复。
 */
int CORE(error)
(
    /* ---- input ---- */
    int status,		        /* 错误状态 */
    const char *file,	    /* 发生错误的源代码文件的名称 */ 
    int line,		        /* 源代码文件中发生错误的行号 */
    const char *message,    /* 错误信息 */
    /* --------------- */
    sparse_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;

    Common->status = status ;

    if (!(Common->try_catch))
    {

#ifndef NPRINT
	/* 打印错误信息 */
	if (SparseBase_config.printf_func != NULL)
	{
	    if (status > 0 && Common->print > 1)
	    {
                SparseBase_config.printf_func ("HNUCHOL warning:") ;
                if (message != NULL)
                {
                    SparseBase_config.printf_func (" %s.", message) ;
                }
                if (file != NULL)
                {
                    SparseBase_config.printf_func (" file: %s", file) ;
                    SparseBase_config.printf_func (" line: %d", line) ;
                }
                SparseBase_config.printf_func ("\n") ;
		fflush (stdout) ;
		fflush (stderr) ;
	    }
	    else if (Common->print > 0)
	    {
                SparseBase_config.printf_func ("HNUCHOL error:") ;
                if (message != NULL)
                {
                    SparseBase_config.printf_func (" %s.", message) ;
                }
                if (file != NULL)
                {
                    SparseBase_config.printf_func (" file: %s", file) ;
                    SparseBase_config.printf_func (" line: %d", line) ;
                }
                SparseBase_config.printf_func ("\n") ;
		fflush (stdout) ;
		fflush (stderr) ;
	    }
	}
#endif

	/* 调用用户错误处理程序(如果存在) */
	if (Common->error_handler != NULL)
	{
	    Common->error_handler (status, file, line, message) ;
	}
    }

    return (TRUE) ;
}
