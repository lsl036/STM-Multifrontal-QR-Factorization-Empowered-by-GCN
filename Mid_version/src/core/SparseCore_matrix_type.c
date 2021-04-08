/**
 * @file SparseCore_matrix_type.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-20
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_internal.h"
#include "SparseCore.h"


/**
 * @brief   为矩阵分配空间。A->i和A->x没有初始化。
 *          A-p(以及A-nz(如果A未被打包)被设为零，
 *          因此返回一个不包含任何项(全部为零)的矩阵。
 *          可见SparseCore_spzeros。
 * 
 *          工作空间：无
 * 
 */
sparse_csc *CORE(allocate_sparse)
(
    /* ---- input ---- */
    size_t nrow,	/* A的行数 */
    size_t ncol,	/* A的列数 */
    size_t nzmax,	/* A的非零元最大值 */
    int sorted,		/* 如果A的列已排序为真，否则为假 */
    int packed,		/* 如果A要填充，则为真，否则为假 */
    int stype,		/* A的数据类型 */
    int xtype,		/* SPARSE_PATTERN, _REAL */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *A ;
    Int *Ap, *Anz ;
    size_t nzmax0 ;
    Int j ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    if (stype != 0 && nrow != ncol)
    {
	ERROR (SPARSE_INVALID, "rectangular matrix with stype != 0 invalid") ;
	return (NULL) ;
    }
    if (xtype < SPARSE_PATTERN || xtype > SPARSE_REAL)
    {
	ERROR (SPARSE_INVALID, "xtype invalid") ;
	return (NULL) ;
    }
    /* 确保维度不会导致整型溢出 */
    (void) CORE(add_size_t) (ncol, 2, &ok) ;
    if (!ok || nrow > Int_max || ncol > Int_max || nzmax > Int_max)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 给抬头分配空间 */
    /* ---------------------------------------------------------------------- */

    A = CORE(malloc) (sizeof (sparse_csc), 1, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    nzmax = MAX (1, nzmax) ;

    A->nrow = nrow ;
    A->ncol = ncol ;
    A->nzmax = nzmax ;
    A->packed = packed ;    /* 默认是填充的(A->nz不存在) */
    A->stype = stype ;
    A->itype = ITYPE ;
    A->xtype = xtype ;
    A->dtype = DTYPE ;

    A->nz = NULL ;
    A->p = NULL ;
    A->i = NULL ;
    A->x = NULL ;
    A->z = NULL ;

    /* 1*m的矩阵A默认列有序 */
    A->sorted = (nrow <= 1) ? TRUE : sorted ;

    /* ---------------------------------------------------------------------- */
    /* 给矩阵本身分配内存空间 */
    /* ---------------------------------------------------------------------- */

    /* 分配大小为O(ncol)的空间 */
    A->p = CORE(malloc) (((size_t) ncol)+1, sizeof (Int), Common) ;
    if (!packed)
    {
	A->nz = CORE(malloc) (ncol, sizeof (Int), Common) ;
    }

    /* 分配大小为O(nz) 的空间*/
    nzmax0 = 0 ;
    CORE(realloc_multiple) (nzmax, 1, xtype, &(A->i), NULL, &(A->x), &(A->z),
	    &nzmax0, Common) ;

    if (Common->status < SPARSE_OK)
    {
	CORE(free_sparse) (&A, Common) ;
	return (NULL) ;	    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 初始化A->p和A->nz，使A为空矩阵 */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    for (j = 0 ; j <= (Int) ncol ; j++)
    {
	Ap [j] = 0 ;
    }
    if (!packed)
    {
	Anz = A->nz ;
	for (j = 0 ; j < (Int) ncol ; j++)
	{
	    Anz [j] = 0 ;
	}
    }
    return (A) ;
}


/**
 * @brief   释放稀疏矩阵存储空间
 *          不需要额外工作空间
 * 
 */
int CORE(free_sparse)
(
    /* ---- in/out --- */
    sparse_csc **AHandle,	/* 需要释放的矩阵，输入为空 */
    /* --------------- */
    sparse_common *Common
)
{
    Int n, nz ;
    sparse_csc *A ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (AHandle == NULL)
    {
	return (TRUE) ;
    }
    A = *AHandle ;
    if (A == NULL)
    {
	return (TRUE) ;
    }
    n = A->ncol ;
    nz = A->nzmax ;
    A->p  = CORE(free) (n+1, sizeof (Int), A->p,  Common) ;
    A->i  = CORE(free) (nz,  sizeof (Int), A->i,  Common) ;
    A->nz = CORE(free) (n,   sizeof (Int), A->nz, Common) ;

    switch (A->xtype)
    {
	case SPARSE_REAL:
	    A->x = CORE(free) (nz, sizeof (double), A->x,  Common) ;
	    break ;
    }

    *AHandle = CORE(free) (1, sizeof (sparse_csc), (*AHandle), Common) ;
    return (TRUE) ;
}


/**
 * @brief   改变A->i,A->x,A->z的大小，或者在它们当前大小为0时分配它们
 *          如果A->xtype==SPARSE_PATTERN，则A->x,A->z不会被修改。
 * 
 */
int CORE(reallocate_sparse)
(
    /* ---- input ---- */
    size_t nznew,	    /* A中条目的新# */
    /* ---- in/out --- */
    sparse_csc *A,	/* 需要重新分配的矩阵 */
    /* --------------- */
    sparse_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 调整矩阵 */
    /* ---------------------------------------------------------------------- */

    CORE(realloc_multiple) (MAX (1,nznew), 1, A->xtype, &(A->i), NULL,
	    &(A->x), &(A->z), &(A->nzmax), Common) ;

    return (Common->status == SPARSE_OK) ;
}


/**
 * @brief 返回稀疏单位矩阵
 * 
 */
sparse_csc *CORE(speye)
(
    /* ---- input ---- */
    size_t nrow,	/* A的行数 */
    size_t ncol,	/* A的列数 */
    int xtype,		/* SPARSE_PATTERN, _REAL */
    /* --------------- */
    sparse_common *Common
)
{
    double *Ax, *Az ;
    sparse_csc *A ;
    Int *Ap, *Ai ;
    Int j, n ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 给矩阵分配空间 */
    /* ---------------------------------------------------------------------- */

    n = MIN (nrow, ncol) ;
    A = CORE(allocate_sparse) (nrow, ncol, n, TRUE, TRUE, 0, xtype,
	    Common) ;

    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出或者非法输入 */
    }

    /* ---------------------------------------------------------------------- */
    /* 创建单位矩阵 */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;

    for (j = 0 ; j < n ; j++)
    {
	Ap [j] = j ;
    }
    for (j = n ; j <= ((Int) ncol) ; j++)
    {
	Ap [j] = n ;
    }
    for (j = 0 ; j < n ; j++)
    {
	Ai [j] = j ;
    }

    switch (xtype)
    {
	case SPARSE_REAL:
	    for (j = 0 ; j < n ; j++)
	    {
		Ax [j] = 1 ;
	    }
	    break ;
    }

    return (A) ;
}


/**
 * @brief 返回稀疏0矩阵
 * 
 */
sparse_csc *CORE(spzeros)
(
    /* ---- input ---- */
    size_t nrow,	/* A的行数 */
    size_t ncol,	/* A的列数 */
    size_t nzmax,	/* A的非零元最大值 */
    int xtype,		/* SPARSE_PATTERN, _REAL */
    /* --------------- */
    sparse_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 给矩阵分配内存空间 */
    /* ---------------------------------------------------------------------- */

    return (CORE(allocate_sparse) (nrow, ncol, nzmax, TRUE, TRUE, 0, xtype,
	    Common)) ;
}


/* Return the number of entries in a sparse matrix.
 *
 * workspace: none
 * integer overflow cannot occur, since the matrix is already allocated.
 */

/**
 * @brief   返回稀疏矩阵中的项数
 *          不需要额外的工作空间
 *          不能发生整型溢出，因为已经分配了矩阵。
 * 
 */
Sparse_long CORE(nnz)
(
    /* ---- input ---- */
    sparse_csc *A,
    /* --------------- */
    sparse_common *Common
)
{
    Int *Ap, *Anz ;
    size_t nz ;
    Int j, ncol ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, EMPTY) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 返回nnz (A) */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    if (A->packed)
    {
	Ap = A->p ;
	RETURN_IF_NULL (Ap, EMPTY) ;
	nz = Ap [ncol] ;
    }
    else
    {
	Anz = A->nz ;
	RETURN_IF_NULL (Anz, EMPTY) ;
	nz = 0 ;
	for (j = 0 ; j < ncol ; j++)
	{
	    nz += MAX (0, Anz [j]) ;
	}
    }
    return (nz) ;
}


/**
 * @brief   C=A.创建一个稀疏矩阵的精确副本，但有一个例外。
 *          未使用的空间中的条目没有被复制(它们可能没有被初始化，
 *          并且复制它们会导致程序检查程序，例如purify和valgrind抱怨)。
 *          得到的矩阵C的xtype与输入矩阵A的xtype相同。
 * 
 */
sparse_csc *CORE(copy_sparse)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要拷贝的矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Ax, *Cx, *Az, *Cz ;
    Int *Ap, *Ai, *Anz, *Cp, *Ci, *Cnz ;
    sparse_csc *C ;
    Int p, pend, j, ncol, packed, nzmax, nz, xtype ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    if (A->stype != 0 && A->nrow != A->ncol)
    {
	ERROR (SPARSE_INVALID, "rectangular matrix with stype != 0 invalid") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    nzmax = A->nzmax ;
    packed = A->packed ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    xtype = A->xtype ;

    /* ---------------------------------------------------------------------- */
    /* 为副本分配空间 */
    /* ---------------------------------------------------------------------- */

    C = CORE(allocate_sparse) (A->nrow, A->ncol, A->nzmax, A->sorted,
	    A->packed, A->stype, A->xtype, Common) ;

    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;
    Cz = C->z ;
    Cnz = C->nz ;

    /* ---------------------------------------------------------------------- */
    /* 拷贝矩阵 */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j <= ncol ; j++)
    {
	Cp [j] = Ap [j] ;
    }

    if (packed)
    {
	nz = Ap [ncol] ;
	for (p = 0 ; p < nz ; p++)
	{
	    Ci [p] = Ai [p] ;
	}

	switch (xtype)
	{
	    case SPARSE_REAL:
		for (p = 0 ; p < nz ; p++)
		{
		    Cx [p] = Ax [p] ;
		}
		break ;
	}

    }
    else
    {

	for (j = 0 ; j < ncol ; j++)
	{
	    Cnz [j] = Anz [j] ;
	}

	switch (xtype)
	{
	    case SPARSE_PATTERN:
		for (j = 0 ; j < ncol ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			Ci [p] = Ai [p] ;
		    }
		}
		break ;

	    case SPARSE_REAL:
		for (j = 0 ; j < ncol ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			Ci [p] = Ai [p] ;
			Cx [p] = Ax [p] ;
		    }
		}
		break ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */
    return (C) ;
}

/* ========================================================================== */
/* === 模板 ============================================================= */
/* ========================================================================== */

#define PATTERN
#include "SparseCore_t_dense.c"
#define REAL
#include "SparseCore_t_dense.c"


/**
 * @brief 分配一个前导维数为d的稠密矩阵。该空间未初始化。
 * 
 */
dense_array *CORE(allocate_dense)
(
    /* ---- input ---- */
    size_t nrow,	/* 矩阵的行数 */
    size_t ncol,	/* 矩阵的列数 */
    size_t d,		/* 前导维数 */
    int xtype,		/* SPARSE_REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *X ;
    size_t nzmax, nzmax0 ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    if (d < nrow)
    {
	ERROR (SPARSE_INVALID, "leading dimension invalid") ;
	return (NULL) ;
    }
    if (xtype != SPARSE_REAL)
    {
	ERROR (SPARSE_INVALID, "xtype invalid") ;
	return (NULL) ;
    }

    /* 确保维度不会导致整型溢出 */
    (void) CORE(add_size_t) (ncol, 2, &ok) ;

    /* nzmax = MAX (1, d*ncol) ; */
    nzmax = CORE(mult_size_t) (d, ncol, &ok) ;
    nzmax = MAX (1, nzmax) ;

    if (!ok || nrow > Int_max || ncol > Int_max || nzmax > Int_max)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 给抬头分配空间 */
    /* ---------------------------------------------------------------------- */

    X = CORE(malloc) (sizeof (dense_array), 1, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    X->nrow = nrow ;
    X->ncol = ncol ;
    X->nzmax = nzmax ;
    X->xtype = xtype ;
    X->dtype = DTYPE ;
    X->x = NULL ;
    X->z = NULL ;
    X->d = d ;

    /* ---------------------------------------------------------------------- */
    /* 给矩阵分配内存 */
    /* ---------------------------------------------------------------------- */

    nzmax0 = 0 ;
    CORE(realloc_multiple) (nzmax, 0, xtype, NULL, NULL, &(X->x), &(X->z),
	    &nzmax0, Common) ;

    if (Common->status < SPARSE_OK)
    {
	CORE(free_dense) (&X, Common) ;
	return (NULL) ;	    /* 内存溢出 */
    }

    return (X) ;
}


/**
 * @brief 给稠密矩阵分配空间并置为0
 * 
 */
dense_array *CORE(zeros)
(
    /* ---- input ---- */
    size_t nrow,	/* 矩阵的行数 */
    size_t ncol,	/* 矩阵的列数 */
    int xtype,		/* SPARSE_REAL */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *X ;
    double *Xx, *Xz ;
    Int i, nz ;

    /* ---------------------------------------------------------------------- */
    /* 给稠密矩阵分配空间并置为0 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    X = CORE(allocate_dense) (nrow, ncol, nrow, xtype, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 检查错误 */
    }

    Xx = X->x ;
    Xz = X->z ;
    nz = MAX (1, X->nzmax) ;

    switch (xtype)
    {
	case SPARSE_REAL:
	    for (i = 0 ; i < nz ; i++)
	    {
		Xx [i] = 0 ;
	    }
	    break ;
    }

    return (X) ;
}


/**
 * @brief 给稠密矩阵分配空间并置为0 
 * 
 */
dense_array *CORE(ones)
(
    /* ---- input ---- */
    size_t nrow,	/* 矩阵的行数 */
    size_t ncol,	/* 矩阵的列数 */
    int xtype,		/* SPARSE_REAL */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *X ;
    double *Xx, *Xz ;
    Int i, nz ;

    /* ---------------------------------------------------------------------- */
    /* 给一个稠密矩阵分配空间，并将其全部设置为1 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    X = CORE(allocate_dense) (nrow, ncol, nrow, xtype, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 检查错误 */
    }

    Xx = X->x ;
    Xz = X->z ;
    nz = MAX (1, X->nzmax) ;

    switch (xtype)
    {
	case SPARSE_REAL:
	    for (i = 0 ; i < nz ; i++)
	    {
		Xx [i] = 1 ;
	    }
	    break ;
    }

    return (X) ;
}


/**
 * @brief 给一个稠密矩阵分配空间，并将其设置为单位矩阵
 * 
 */
dense_array *CORE(eye)
(
    /* ---- input ---- */
    size_t nrow,	/* 矩阵的行数 */
    size_t ncol,	/* 矩阵的列数 */
    int xtype,		/* SPARSE_REAL */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *X ;
    double *Xx, *Xz ;
    Int i, n, nz ;

    /* ---------------------------------------------------------------------- */
    /* 给一个稠密矩阵分配空间，并将其设置为单位矩阵 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    X = CORE(zeros) (nrow, ncol, xtype, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 检查错误 */
    }

    nz = MAX (1, nrow*ncol) ;
    Xx = X->x ;
    Xz = X->z ;

    n = MIN (nrow, ncol) ;
    return (X) ;
}


/**
 * @brief   释放稠密矩阵存储空间
 *          不需要额外工作空间
 */
int CORE(free_dense)
(
    /* ---- in/out --- */
    dense_array **XHandle,	/* 需要释放的稠密矩阵，输出为空 */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *X ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (XHandle == NULL)
    {
	return (TRUE) ;
    }
    X = *XHandle ;
    if (X == NULL)
    {
	return (TRUE) ;
    }

    switch (X->xtype)
    {
	case SPARSE_REAL:
	    X->x = CORE(free) (X->nzmax, sizeof (double), X->x, Common) ;
	    break ;
    }

    *XHandle = CORE(free) (1, sizeof (dense_array), (*XHandle), Common) ;
    return (TRUE) ;
}


/* Ensure that the input matrix has a certain size and type.  If not, free
 * the existing matrix and reallocate one of the right size and type.
 * Returns a pointer to the dense_array matrix, possibly reallocated.
 * Also modifies the input matrix handle, XHandle, if necessary.
 */

/**
 * @brief   确保输入矩阵有一定的大小和类型。
 *          如果没有，释放现有矩阵并重新分配一个合适的大小和类型。
 *          返回一个指向SparseCore_dense矩阵的指针，可能是重新分配的。
 *          如果需要，还可以修改输入矩阵句柄XHandle。
 */
dense_array *CORE(ensure_dense)
(
    /* ---- input/output ---- */
    dense_array **XHandle,    /* 需要检查的矩阵句柄 */
    /* ---- input ---- */
    size_t nrow,	            /* 矩阵的行数 */
    size_t ncol,	            /* 矩阵的列数 */
    size_t d,		            /* 前导维度 */
    int xtype,		            /* SPARSE_REAL */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *X ;

    RETURN_IF_NULL_COMMON (NULL) ;
    if (XHandle == NULL)
    {
        ERROR (SPARSE_INVALID, "matrix invalid") ;
        return (NULL) ;
    }

    X = *XHandle ;
    if (X == NULL || X->nrow != nrow || X->ncol != ncol
        || X->d != d || X->xtype != xtype)
    {
        /* 矩阵X没有分配，或者大小不对。释放它并重新分配到正确的大小和形状。
         * 如果发生错误(内存不足或输入nrow无效等)，则在SparseCore_allocate_dense中
         * 设置错误，并将X返回为NULL。 */
        CORE(free_dense) (XHandle, Common) ;
        X = CORE(allocate_dense) (nrow, ncol, d, xtype, Common) ;
        *XHandle = X ;
    }
    return (X) ;
}


/**
 * @brief   将稀疏矩阵转换为稠密矩阵。输出的稠密矩阵的xtype与输入的稀疏矩阵相同，
 *          只是将一个模式稀疏矩阵A转换为一个实稠密矩阵X，其值分别为1和0。
 *          支持所有的xtype。
 */
dense_array *CORE(sparse_to_dense)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要拷贝的稀疏矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *X = NULL ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入*/
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    if (A->stype && A->nrow != A->ncol)
    {
	ERROR (SPARSE_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 使用模板历程转换矩阵 */
    /* ---------------------------------------------------------------------- */

    switch (A->xtype)
    {
	case SPARSE_PATTERN:
	    X = p_SparseCore_sparse_to_dense (A, Common) ;
	    break ;

	case SPARSE_REAL:
	    X = r_SparseCore_sparse_to_dense (A, Common) ;
	    break ;
    }
    return (X) ;
}

/**
 * @brief 将稠密矩阵转换为稀疏矩阵
 */
sparse_csc *CORE(dense_to_sparse)
(
    /* ---- input ---- */
    dense_array *X,	/* 需要拷贝的稠密矩阵 */
    int values,		    /* 如果要复制值为TRUE，否则为FALSE */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *C = NULL ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (X, NULL) ;
    RETURN_IF_XTYPE_INVALID (X, SPARSE_REAL, SPARSE_REAL, NULL) ;
    if (X->d < X->nrow)
    {
	ERROR (SPARSE_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 使用模板历程转换矩阵 */
    /* ---------------------------------------------------------------------- */

    switch (X->xtype)
    {
	case SPARSE_REAL:
	    C = r_SparseCore_dense_to_sparse (X, values, Common) ;
	    break ;
    }
    return (C) ;
}


/**
 * @brief   Y=X，其中X和Y都已经分配好了。X和Y的前导维可能不同，
 *          但都必须>=X和Y的行数。行nrow到d-1中的项不会从X复制，
 *          因为空间可能不会初始化。Y->nzmax是不变的。X->nzmax通常是
 *          (X->d)*(X->ncol)。
 *          但是用户可以在任何HNUCHOL例程之外修改该条件。
 * 
 *          两个密集矩阵X和Y必须具有相同的数值类型。
 */
int CORE(copy_dense2)
(
    /* ---- input ---- */
    dense_array *X,	/* 需要拷贝的矩阵 */
    /* ---- output --- */
    dense_array *Y,	/* 矩阵X的副本 */
    /* --------------- */
    sparse_common *Common
)
{
    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (X, FALSE) ;
    RETURN_IF_NULL (Y, FALSE) ;
    RETURN_IF_XTYPE_INVALID (X, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (Y, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    if (X->nrow != Y->nrow || X->ncol != Y->ncol || X->xtype != Y->xtype)
    {
	ERROR (SPARSE_INVALID, "X and Y must have same dimensions and xtype") ;
	return (FALSE) ;
    }
    if (X->d < X->nrow || Y->d < Y->nrow
	    || (X->d * X->ncol) > X->nzmax || (Y->d * Y->ncol) > Y->nzmax)
    {
	ERROR (SPARSE_INVALID, "X and/or Y invalid") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 使用模板例程复制矩阵 */
    /* ---------------------------------------------------------------------- */

    switch (X->xtype)
    {
	case SPARSE_REAL:
	    r_SparseCore_copy_dense2 (X, Y) ;
	    break ;
    }
    return (TRUE) ;
}


/**
 * @brief Y=X 复制稠密矩阵
 * 
 */
dense_array *CORE(copy_dense)
(
    /* ---- input ---- */
    dense_array *X,	/* 需要拷贝的矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *Y ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (X, NULL) ;
    RETURN_IF_XTYPE_INVALID (X, SPARSE_REAL, SPARSE_REAL, NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 给结果分配空间 */
    /* ---------------------------------------------------------------------- */

    Y = CORE(allocate_dense) (X->nrow, X->ncol, X->d, X->xtype, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出或者X非法 */
    }

    /* ---------------------------------------------------------------------- */
    /* Y = X */
    /* ---------------------------------------------------------------------- */

    /* 这不会失败(分配了X和Y，并且具有相同的nrow、ncol d和xtype) */
    CORE(copy_dense2) (X, Y, Common) ;

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */

    return (Y) ;
}


/* ========================================================================== */
/* === 模板 ============================================================= */
/* ========================================================================== */

#define PATTERN
#include "SparseCore_t_triplet.c"
#define REAL
#include "SparseCore_t_triplet.c"

/**
 * @brief   为一个三元组矩阵分配空间
 *          不需要额外工作空间
 */
sparse_triplet *CORE(allocate_triplet)
(
    /* ---- input ---- */
    size_t nrow,	/* T的行数 */
    size_t ncol,	/* T的列数 */
    size_t nzmax,	/* T的非零元最大值 */
    int stype,		/* T的数据类型 */
    int xtype,		/* SPARSE_PATTERN, _REAL */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_triplet *T ;
    size_t nzmax0 ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    if (xtype < SPARSE_PATTERN || xtype > SPARSE_REAL)
    {
	ERROR (SPARSE_INVALID, "xtype invalid") ;
	return (NULL) ;
    }
    /* 确保维度不会导致整数溢出 */
    (void) CORE(add_size_t) (ncol, 2, &ok) ;
    if (!ok || nrow > Int_max || ncol > Int_max || nzmax > Int_max)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }

    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配头 */
    /* ---------------------------------------------------------------------- */

    T = CORE(malloc) (sizeof (sparse_triplet), 1, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    nzmax = MAX (1, nzmax) ;

    T->nrow = nrow ;
    T->ncol = ncol ;
    T->nzmax = nzmax ;
    T->nnz = 0 ;
    T->stype = stype ;
    T->itype = ITYPE ;
    T->xtype = xtype ;
    T->dtype = DTYPE ;

    T->j = NULL ;
    T->i = NULL ;
    T->x = NULL ;
    T->z = NULL ;

    /* ---------------------------------------------------------------------- */
    /* 分配矩阵本身 */
    /* ---------------------------------------------------------------------- */

    nzmax0 = 0 ;
    CORE(realloc_multiple) (nzmax, 2, xtype, &(T->i), &(T->j),
		&(T->x), &(T->z), &nzmax0, Common) ;

    if (Common->status < SPARSE_OK)
    {
	CORE(free_triplet) (&T, Common) ;
	return (NULL) ;	    /* 内存溢出 */
    }

    return (T) ;
}

/**
 * @brief   释放三元组矩阵空间
 *          不需要额外存储空间
 */
int CORE(free_triplet)
(
    /* ---- in/out --- */
    sparse_triplet **THandle,    /* 需要释放的矩阵，输出为空 */
    /* --------------- */
    sparse_common *Common
)
{
    Int nz ;
    sparse_triplet *T ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (THandle == NULL)
    {
	return (TRUE) ;
    }
    T = *THandle ;
    if (T == NULL)
    {
	return (TRUE) ;
    }
    nz = T->nzmax ;
    T->j = CORE(free) (nz, sizeof (Int), T->j, Common) ;
    T->i = CORE(free) (nz, sizeof (Int), T->i, Common) ;
    if (T->xtype == SPARSE_REAL)
    {
	T->x = CORE(free) (nz, sizeof (double), T->x, Common) ;
    }
    *THandle = CORE(free) (1, sizeof (sparse_triplet), (*THandle), Common) ;
    return (TRUE) ;
}


/* Change the size of T->i, T->j, and T->x, or allocate them if their current
 * size is zero.  T->x is not modified if T->xtype is SPARSE_PATTERN.
 *
 * workspace: none
 */

/**
 * @brief   改变T->i T->j T->x的大小，如果它们的当前大小为零，就赋值给它们。
 *          如果T->xtype==SPARSE_PATTERN，则不修改T->x
 */
int CORE(reallocate_triplet)
(
    /* ---- input ---- */
    size_t nznew,	    /* T中的新条目# */
    /* ---- in/out --- */
    sparse_triplet *T,	/* 需要修改的三元组 */
    /* --------------- */
    sparse_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (T, FALSE) ;
    RETURN_IF_XTYPE_INVALID (T, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 调整矩阵 */
    /* ---------------------------------------------------------------------- */

    CORE(realloc_multiple) (MAX (1,nznew), 2, T->xtype, &(T->i), &(T->j),
	    &(T->x), &(T->z), &(T->nzmax), Common) ;

    return (Common->status == SPARSE_OK) ;
}


/**
 * @brief   将一系列三元组转换为SparseCore_sparse矩阵。
 * 
 */
sparse_csc *CORE(triplet_to_sparse)
(
    /* ---- input ---- */
    sparse_triplet *T,	/* 需要拷贝的矩阵 */
    size_t nzmax,	    /* 在输出矩阵中分配至少这么大的空间 */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *R, *A = NULL ;
    Int *Wj, *Rp, *Ri, *Rnz, *Ti, *Tj ;
    Int i, j, p, k, stype, nrow, ncol, nz, ok ;
    size_t anz = 0 ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (T, NULL) ;
    Ti = T->i ;
    Tj = T->j ;
    RETURN_IF_NULL (Ti, NULL) ;
    RETURN_IF_NULL (Tj, NULL) ;
    RETURN_IF_XTYPE_INVALID (T, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    stype = SIGN (T->stype) ;
    if (stype && T->nrow != T->ncol)
    {
	/* 非法输入 */
	ERROR (SPARSE_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    nrow = T->nrow ;
    ncol = T->ncol ;
    nz = T->nnz ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    CORE(allocate_work) (0, MAX (nrow, ncol), 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 指定临时矩阵R */
    /* ---------------------------------------------------------------------- */

    R = CORE(allocate_sparse) (ncol, nrow, nz, FALSE, FALSE, -stype,
	    T->xtype, Common) ;

    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    Rp = R->p ;
    Ri = R->i ;
    Rnz = R->nz ;

    /* ---------------------------------------------------------------------- */
    /* 计算A的每一行中的条目(也计算重复项) */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nrow ; i++)
    {
	Rnz [i] = 0 ;	
    }

    if (stype > 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i >= nrow || j < 0 || j >= ncol)
	    {
		ERROR (SPARSE_INVALID, "index out of range") ;
		break ;
	    }
	    /* A将是对称的，只存储上三角形部分。
         * 创建一个下三角形矩阵R。
         * R的上半部分的元素被转置到下半部分。 */
	    Rnz [MIN (i,j)]++ ;
	}
    }
    else if (stype < 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i >= nrow || j < 0 || j >= ncol)
	    {
		ERROR (SPARSE_INVALID, "index out of range") ;
		break ;
	    }
	    /* A将是对称的，只存储下三角形部分。
         * 创建一个上三角形矩阵R。
         * R的下三角部分的元素被转置到上三角部分。 */
	    Rnz [MAX (i,j)]++ ;
	}
    }
    else
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i >= nrow || j < 0 || j >= ncol)
	    {
		ERROR (SPARSE_INVALID, "index out of range") ;
		break ;
	    }
	    /* 构造一个不对称矩阵 */
	    Rnz [i]++ ;
	}
    }

    if (Common->status < SPARSE_OK)
    {
	/* 三元组矩阵非法 */
	CORE(free_sparse) (&R, Common) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 构造行指针 */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (i = 0 ; i < nrow ; i++)
    {
	Rp [i] = p ;
	p += Rnz [i] ;
    }
    Rp [nrow] = p ;

    /* 使用Wj (i/l/l)作为临时行指针 */
    Wj = Common->Iwork ;	/* 大小为MAX (nrow,ncol) FUTURE WORK: (i/l/l) */
    for (i = 0 ; i < nrow ; i++)
    {
	Wj [i] = Rp [i] ;
    }

    /* ---------------------------------------------------------------------- */
    /* 使用模板例程构造三重矩阵 */
    /* ---------------------------------------------------------------------- */

    switch (T->xtype)
    {
	case SPARSE_PATTERN:
	    anz = p_SparseCore_triplet_to_sparse (T, R, Common) ;
	    break ;

	case SPARSE_REAL:
	    anz = r_SparseCore_triplet_to_sparse (T, R, Common) ;
	    break ;
    }

    /* ---------------------------------------------------------------------- */
    /* A = R' (数组转置，共轭转置) */
    /* ---------------------------------------------------------------------- */

    /* 工作空间: Iwork (R->nrow), 大小为A->ncol */

    A = CORE(allocate_sparse) (nrow, ncol, MAX (anz, nzmax), TRUE, TRUE,
	stype, T->xtype, Common) ;

    if (stype)
    {
	ok = CORE(transpose_sym) (R, 1, NULL, A, Common) ;
    }
    else
    {
	ok = CORE(transpose_unsym) (R, 1, NULL, NULL, 0, A, Common) ; 
    }

    CORE(free_sparse) (&R, Common) ;
    if (Common->status < SPARSE_OK)
    {
	CORE(free_sparse) (&A, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */
    return (A) ;
}


/* Converts a sparse column-oriented matrix to triplet form.
 * The resulting triplet matrix has the same xtype as the sparse matrix.
 *
 * workspace: none
 */

/**
 * @brief   将一个稀疏的面向列的矩阵转换为三元组形式。
 *          得到的三元组矩阵具有与稀疏矩阵相同的xtype。
 *          不需要额外的工作空间
 */
sparse_triplet *CORE(sparse_to_triplet)
(
    /* ---- input ---- */
    sparse_csc *A,	/* 需要拷贝的矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Ax, *Az, *Tx, *Tz ;
    Int *Ap, *Ai, *Ti, *Tj, *Anz ;
    sparse_triplet *T ;
    Int i, xtype, p, pend, k, j, nrow, ncol, nz, stype, packed, up, lo,
	both ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    stype = SIGN (A->stype) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if (stype && nrow != ncol)
    {
	/* 非法输入 */
	ERROR (SPARSE_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Ax = A->x ;
    Az = A->z ;
    xtype = A->xtype ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 指定三元组矩阵 */
    /* ---------------------------------------------------------------------- */

    nz = CORE(nnz) (A, Common) ;
    T = CORE(allocate_triplet) (nrow, ncol, nz, A->stype, A->xtype, Common) ;

    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 转换为稀疏矩阵 */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;

    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    Tz = T->z ;
    T->stype = A->stype ;

    both = (A->stype == 0) ;
    up = (A->stype > 0) ;
    lo = (A->stype < 0) ;

    k = 0 ;

    for (j = 0 ; j < ncol ; j++)
    {
	p = Ap [j] ;
	pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	for ( ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if (both || (up && i <= j) || (lo && i >= j))
	    {
		Ti [k] = Ai [p] ;
		Tj [k] = j ;

		if (xtype == SPARSE_REAL)
		{
		    Tx [k] = Ax [p] ;
		}

		k++ ;
	    }
	}
    }

    T->nnz = k ;

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */
    return (T) ;
}


/* Create an exact copy of a triplet matrix, except that entries in unused
 * space are not copied (they might not be initialized, and copying them would
 * cause program checkers such as purify and valgrind to complain).
 * The output triplet matrix has the same xtype as the input triplet matrix.
 */

/**
 * @brief   创建一个精确的三元组矩阵的副本，除了未使用的空间中的条目
 *          没有被复制(它们可能没有被初始化，并且复制它们会开启检
 *          查程序，如purify和内存溢出)。
 *          输出的三元组矩阵具有与输入三元组矩阵相同的xtype。
 */
sparse_triplet *CORE(copy_triplet)
(
    /* ---- input ---- */
    sparse_triplet *T,	/* 需要拷贝的矩阵 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Tx, *Tz, *Cx, *Cz ;
    Int *Ci, *Cj, *Ti, *Tj ;
    sparse_triplet *C ;
    Int xtype, k, nz ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (T, NULL) ;
    RETURN_IF_XTYPE_INVALID (T, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    nz = T->nnz ;
    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    Tz = T->z ;
    xtype = T->xtype ;
    RETURN_IF_NULL (Ti, NULL) ;
    RETURN_IF_NULL (Tj, NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 指定副本 */
    /* ---------------------------------------------------------------------- */

    C = CORE(allocate_triplet) (T->nrow, T->ncol, T->nzmax, T->stype,
	    xtype, Common) ;

    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    /* ---------------------------------------------------------------------- */
    /* 拷贝三元组矩阵 */
    /* ---------------------------------------------------------------------- */

    Ci = C->i ;
    Cj = C->j ;
    Cx = C->x ;
    Cz = C->z ;
    C->nnz = nz ;

    for (k = 0 ; k < nz ; k++)
    {
	Ci [k] = Ti [k] ;
    }
    for (k = 0 ; k < nz ; k++)
    {
	Cj [k] = Tj [k] ;
    }

    if (xtype == SPARSE_REAL)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    Cx [k] = Tx [k] ;
	}
    }
    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */
    return (C) ;
}


/**
 * @brief   指定一个已经初始化L->Perm和L->ColCount
 *          为“空”(Perm [k] = k, and ColCount[k] = 1)
 *          的简单的符号因子。不分配L的整数部分和数值部分。
 *          L->xtype作为SPARSE_PATTERN返回，L->is_super作为FALSE返回。
 *          L->is_ll返回也为FALSE，但是当矩阵被分解时，这可能会被修改。
 */
sparse_factor *CORE(allocate_factor)
(
    /* ---- input ---- */
    size_t n,		/* L为n*n */
    /* --------------- */
    sparse_common *Common
)
{
    Int j ;
    Int *Perm, *ColCount ;
    sparse_factor *L ;
    int ok = TRUE ;

    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = SPARSE_OK ;

    /* 确保维度不会导致整型溢出 */
    (void) CORE(add_size_t) (n, 2, &ok) ;
    if (!ok || n > Int_max)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }

    L = CORE(malloc) (sizeof (sparse_factor), 1, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }
    L->n = n ;
    L->is_ll = FALSE ;
    L->is_super = FALSE ;
    L->is_monotonic = TRUE ;
    L->itype = ITYPE ;
    L->xtype = SPARSE_PATTERN ;
    L->dtype = DTYPE ;

    /* 指定L的纯符号部分 */
    L->ordering = SPARSE_NATURAL ;
    L->Perm = CORE(malloc) (n, sizeof (Int), Common) ;
    L->IPerm = NULL ;       /* 仅在需要时由cholmod_solve2创建 */
    L->ColCount = CORE(malloc) (n, sizeof (Int), Common) ;

    /* L的单纯部分是空的 */
    L->nzmax = 0 ;
    L->p = NULL ;
    L->i = NULL ;
    L->x = NULL ;
    L->z = NULL ;
    L->nz = NULL ;
    L->next = NULL ;
    L->prev = NULL ;

    /* L的超节点部分也是空的 */
    L->nsuper = 0 ;
    L->ssize = 0 ;
    L->xsize = 0 ;
    L->maxesize = 0 ;
    L->maxcsize = 0 ;
    L->super = NULL ;
    L->pi = NULL ;
    L->px = NULL ;
    L->s = NULL ;

    /* L还没有被分解 */
    L->minor = n ;

    if (Common->status < SPARSE_OK)
    {
	CORE(free_factor) (&L, Common) ;
	return (NULL) ;		/* 内存溢出 */
    }

    /* 初始化Perm和ColCount */
    Perm = L->Perm ;
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	Perm [j] = j ;
    }
    ColCount = L->ColCount ;
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	ColCount [j] = 1 ;
    }

    return (L) ;
}

/**
 * @brief   释放一个因子对象。
 *          不需要额外的工作空间
 */
int CORE(free_factor)
(
    /* ---- in/out --- */
    sparse_factor **LHandle,	/* 需要释放的因子，输出为空 */
    /* --------------- */
    sparse_common *Common
)
{
    Int n, lnz, xs, ss, s ;
    sparse_factor *L ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (LHandle == NULL)
    {
	return (TRUE) ;
    }
    L = *LHandle ;
    if (L == NULL)
    {
	return (TRUE) ;
    }

    n = L->n ;
    lnz = L->nzmax ;
    s = L->nsuper + 1 ;
    xs = (L->is_super) ? ((Int) (L->xsize)) : (lnz) ;
    ss = L->ssize ;

    /* L的符号部分 */
    CORE(free) (n,   sizeof (Int), L->Perm,     Common) ;
    CORE(free) (n,   sizeof (Int), L->IPerm,    Common) ;
    CORE(free) (n,   sizeof (Int), L->ColCount, Common) ;

    /* L的simplicial形式 */
    CORE(free) (n+1, sizeof (Int), L->p,        Common) ;
    CORE(free) (lnz, sizeof (Int), L->i,        Common) ;
    CORE(free) (n,   sizeof (Int), L->nz,       Common) ;
    CORE(free) (n+2, sizeof (Int), L->next,     Common) ;
    CORE(free) (n+2, sizeof (Int), L->prev,     Common) ;

    /* L的超节点形式 */
    CORE(free) (s,   sizeof (Int), L->pi,       Common) ;
    CORE(free) (s,   sizeof (Int), L->px,       Common) ;
    CORE(free) (s,   sizeof (Int), L->super,    Common) ;
    CORE(free) (ss,  sizeof (Int), L->s,        Common) ;

    /* simplicial和超节点L的数值 */
    if (L->xtype == SPARSE_REAL)
    {
	CORE(free) (xs, sizeof (double), L->x, Common) ;
    }
    *LHandle = CORE(free) (1, sizeof (sparse_factor), (*LHandle), Common) ;
    return (TRUE) ;
}


/* Change the size of L->i and L->x, or allocate them if their current size
 * is zero.  L must be simplicial.
 *
 * workspace: none
 */

/**
 * @brief   改变L->i和L->x的大小，或者它们的当前大小为零时指定它们。
 *          L必须为simplicial。
 *          不需要额外工作空间
 */
int CORE(reallocate_factor)
(
    /* ---- input ---- */
    size_t nznew,	    /* L中条目的新# */
    /* ---- in/out --- */
    sparse_factor *L,	/* 需要修改的因子 */
    /* --------------- */
    sparse_common *Common
)
{
    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    if (L->is_super)
    {
	/* L必须为simplicial而不是symbolic */
	ERROR (SPARSE_INVALID, "L invalid") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 调整（或者指定）该因子的L->i和L->x分量 */
    /* ---------------------------------------------------------------------- */

    CORE(realloc_multiple) (nznew, 1, L->xtype, &(L->i), NULL,
	    &(L->x), &(L->z), &(L->nzmax), Common) ;
    return (Common->status == SPARSE_OK) ;
}

/**
 * @brief   列j需要更多的空间，重新分配到L->i和L->x的末尾。
 *          如果重新分配失败，则将该因子转换为单纯符号因子(没
 *          有模式，只有L->Perm和L->ColCount)。
 *          不需要额外的工作空间
 */
int CORE(reallocate_column)
(
    /* ---- input ---- */
    size_t j,		    /* 要重新分配的列 */
    size_t need,	    /* 第j列所需的大小 */
    /* ---- in/out --- */
    sparse_factor *L,	/* 需要修改的因子 */
    /* --------------- */
    sparse_common *Common
)
{
    double xneed ;
    double *Lx, *Lz ;
    Int *Lp, *Lprev, *Lnext, *Li, *Lnz ;
    Int n, pold, pnew, len, k, tail ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    if (L->is_super)
    {
	ERROR (SPARSE_INVALID, "L must be simplicial") ;
	return (FALSE) ;
    }
    n = L->n ;
    if (j >= L->n || need == 0)
    {
	ERROR (SPARSE_INVALID, "j invalid") ;
	return (FALSE) ;	    /* j超过范围 */
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 如果需要，增加L的大小 */
    /* ---------------------------------------------------------------------- */

    /* head = n+1 ; */
    tail = n ;
    Lp = L->p ;
    Lnz = L->nz ;
    Lprev = L->prev ;
    Lnext = L->next ;

    /* 如果所有元素都存在，列j不能有超过n-j个元素 */
    need = MIN (need, n-j) ;

    /* 计算需要在双精度中避免整数溢出 */
    if (Common->grow1 >= 1.0)
    {
	xneed = (double) need ;
	xneed = Common->grow1 * xneed + Common->grow2 ;
	xneed = MIN (xneed, n-j) ;
	need = (Int) xneed ;
    }


    if (Lp [Lnext [j]] - Lp [j] >= (Int) need)
    {
	/* 不需要重新分配列，它已经足够大了 */
	return (TRUE) ;

    }

    if (Lp [tail] + need > L->nzmax)
    {
	/* 使用double来避免整数溢出 */
	xneed = (double) need ;
	if (Common->grow0 < 1.2)	    /* 比较，如果NaN为假 */
	{
	    /* 如果grow0小于1.2或NaN，不要使用它 */
	    xneed = 1.2 * (((double) L->nzmax) + xneed + 1) ;
	}
	else
	{
	    xneed = Common->grow0 * (((double) L->nzmax) + xneed + 1) ;
	}
	if (xneed > Size_max ||
		!CORE(reallocate_factor) ((Int) xneed, L, Common))
	{
	    /* 内存溢出, 转换为简单的符号 */
	    CORE(change_factor) (SPARSE_PATTERN, L->is_ll, FALSE, TRUE,
		    TRUE, L, Common) ;
	    ERROR (SPARSE_OUT_OF_MEMORY, "内存溢出; L now symbolic") ;
	    return (FALSE) ;	    /* 内存溢出 */
	}
	/* 填充所有列，以便每个列最多拥有grow2空闲空间 */
	CORE(pack_factor) (L, Common) ;
	Common->nrealloc_factor++ ;
    }

    /* ---------------------------------------------------------------------- */
    /* 重新分配列 */
    /* ---------------------------------------------------------------------- */

    Common->nrealloc_col++ ;

    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;

    /* 从列表中的当前位置删除j */
    Lnext [Lprev [j]] = Lnext [j] ;
    Lprev [Lnext [j]] = Lprev [j] ;

    /* 把j放在列表的最后 */
    Lnext [Lprev [tail]] = j ;
    Lprev [j] = Lprev [tail] ;
    Lnext [j] = n ;
    Lprev [tail] = j ;

    /* L不再是单调的;列是无序的 */
    L->is_monotonic = FALSE ;

    /* 为列j分配空间 */
    pold = Lp [j] ;
    pnew = Lp [tail] ;
    Lp [j] = pnew  ;
    Lp [tail] += need ;

    /* 将列j复制到新空间 */
    len = Lnz [j] ;
    for (k = 0 ; k < len ; k++)
    {
	Li [pnew + k] = Li [pold + k] ;
    }

    if (L->xtype == SPARSE_REAL)
    {
	for (k = 0 ; k < len ; k++)
	{
	    Lx [pnew + k] = Lx [pold + k] ;
	}
    }

    /* L的j列重新分配成功 */
    return (TRUE) ;
}


/**
 * @brief   填充简单LDL'或者LL'因子的列。
 *          接下来可以调用SparseCore_reallocate_factor来将L的大小减小
 *          到该因子所需的精确大小(如果需要的话)。
 *          或者，您可以保持L->i和L->x的大小不变，以便为将来的 更新/添加行 留出空间。
 * 
 *          每个列的大小都减小了，因此在列的末尾最多有Common->grow2的空闲空间。
 * 
 *          如果给定任何其他类型的因子，则不执行任何操作并静默返回。
 * 
 *          不会使L的列是单调的。因此，它不同于填充列并确保它们以单调的顺序出现的
 *          SparseCore_change_factor (xtype, -, FALSE, TRUE, TRUE, L, Common)。
 */
int CORE(pack_factor)
(
    /* ---- in/out --- */
    sparse_factor *L,	/* 需要修改的因子 */
    /* --------------- */
    sparse_common *Common
)
{
    double *Lx, *Lz ;
    Int *Lp, *Li, *Lnz, *Lnext ;
    Int pnew, j, k, pold, len, n, head, tail, grow2 ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    Common->status = SPARSE_OK ;

    if (L->xtype == SPARSE_PATTERN || L->is_super)
    {
	/* 什么也不做，除非L是单纯的数字 */
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 填充 */
    /* ---------------------------------------------------------------------- */

    grow2 = Common->grow2 ;

    pnew = 0 ;
    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;
    Lnz = L->nz ;
    Lnext = L->next ;

    head = n+1 ;
    tail = n ;

    for (j = Lnext [head] ; j != tail ; j = Lnext [j])
    {
	/* 填充第j列 */
	pold = Lp [j] ;
	len = Lnz [j] ;
	if (pnew < pold)
	{

	    for (k = 0 ; k < len ; k++)
	    {
		Li [pnew + k] = Li [pold + k] ;
	    }

	    if (L->xtype == SPARSE_REAL)
	    {
		for (k = 0 ; k < len ; k++)
		{
		    Lx [pnew + k] = Lx [pold + k] ;
		}
	    }
	    Lp [j] = pnew ;
	}
	len = MIN (len + grow2, n - j) ;
	pnew = MIN (Lp [j] + len, Lp [Lnext [j]]) ;
    }
    return (TRUE) ;
}


/**
 * @brief   构造含有simplicial或超节点数值因子的模式和值的列主稀疏矩阵，
 *          并将其转化为simplicial符号因子。
 *          如果L已经被填充、monotonic和simplicial(SparseCore_factorize使用
 *          simplicial Cholesky分解算法时就是这种情况)，
 *          那么这个例程只需要O(1)内存和O(1)时间。
 */
sparse_csc *CORE(factor_to_sparse)
(
    /* ---- in/out --- */
    sparse_factor *L,	/* 要复制的因子，在输出时转换为符号 */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *Lsparse ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (L, NULL) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_REAL, SPARSE_REAL, NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 转换为packed, monotonic, simplicial, numeric */
    /* ---------------------------------------------------------------------- */

    /* 作为LL或者LDL'留下 */
    if (!CORE(change_factor) (L->xtype, L->is_ll, FALSE, TRUE, TRUE, L,
		Common))
    {
	ERROR (SPARSE_INVALID, "cannot convert L") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 创建Lsparse */
    /* ---------------------------------------------------------------------- */

    /* 为Lsparse指定头，L的稀疏矩阵版本 */
    Lsparse = CORE(malloc) (sizeof (sparse_csc), 1, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;		/* 内存溢出 */
    }

    /* 将内容从L转移到Lsparse */
    Lsparse->nrow = L->n ;
    Lsparse->ncol = L->n ;
    Lsparse->p = L->p ;
    Lsparse->i = L->i ;
    Lsparse->x = L->x ;
    Lsparse->z = L->z ;
    Lsparse->nz = NULL ;
    Lsparse->stype = 0 ;
    Lsparse->itype = L->itype ;
    Lsparse->xtype = L->xtype ;
    Lsparse->dtype = L->dtype ;
    Lsparse->sorted = TRUE ;
    Lsparse->packed = TRUE ;
    Lsparse->nzmax = L->nzmax ;

    /* ---------------------------------------------------------------------- */
    /* 将L转换为符号，但不释放转移到Lsparse的内容 */
    /* ---------------------------------------------------------------------- */

    L->p = NULL ;
    L->i = NULL ;
    L->x = NULL ;
    L->z = NULL ;
    L->xtype = SPARSE_PATTERN ;
    CORE(change_factor) (SPARSE_PATTERN, FALSE, FALSE, TRUE, TRUE, L,
	    Common) ;

    return (Lsparse) ;
}


/**
 * @brief   创建一个因子的精确副本，但有一个例外:
 *          未使用空间中的条目没有被复制(它们可能没有被初始化，
 *          并且复制它们会开启检查程序，例如purify和valgrind)。
 */
sparse_factor *CORE(copy_factor)
(
    /* ---- input ---- */
    sparse_factor *L,	/* 要拷贝的因子 */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_factor *L2 ;
    double *Lx, *L2x, *Lz, *L2z ;
    Int *Perm, *ColCount, *Lp, *Li, *Lnz, *Lnext, *Lprev, *Lsuper, *Lpi, *Lpx,
	*Ls, *Perm2, *ColCount2, *L2p, *L2i, *L2nz, *L2next, *L2prev, *L2super,
	*L2pi, *L2px, *L2s ;
    Int n, j, p, pend, s, xsize, ssize, nsuper ;

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (L, NULL) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_PATTERN, SPARSE_REAL, NULL) ;
    Common->status = SPARSE_OK ;

    n = L->n ;

    /* ---------------------------------------------------------------------- */
    /* 指定一个简单的符号因子  */
    /* ---------------------------------------------------------------------- */

    /* 指定L2->Perm和L2->ColCount */
    L2 = CORE(allocate_factor) (n, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    Perm = L->Perm ;
    ColCount = L->ColCount ;
    Perm2 = L2->Perm ;
    ColCount2 = L2->ColCount ;
    L2->ordering = L->ordering ;

    for (j = 0 ; j < n ; j++)
    {
	Perm2 [j] = Perm [j] ;
    }
    for (j = 0 ; j < n ; j++)
    {
	ColCount2 [j] = ColCount [j] ;
    }
    L2->is_ll = L->is_ll ;

    /* ---------------------------------------------------------------------- */
    /* 复制因子的剩余部分 */
    /* ---------------------------------------------------------------------- */

    if (L->xtype != SPARSE_PATTERN && !(L->super))
    {

	/* ------------------------------------------------------------------ */
	/* 指定一个简单数值因子 */
	/* ------------------------------------------------------------------ */

	/* 指定 L2->p, L2->nz, L2->prev, L2->next, L2->i, and L2->x.
	 * packed = -1，因此SparseCore_change_factor分配的空间大小为L2->nzmax */
	L2->nzmax = L->nzmax ;
	if (!CORE(change_factor) (L->xtype, L->is_ll, FALSE, -1, TRUE,
		    L2, Common))
	{
	    CORE(free_factor) (&L2, Common) ;
	    return (NULL) ;	/* 内存溢出 */
	}

	/* ------------------------------------------------------------------ */
	/* 复制简单数值因子的内容 */
	/* ------------------------------------------------------------------ */

	Lp = L->p ;
	Li = L->i ;
	Lx = L->x ;
	Lz = L->z ;
	Lnz = L->nz ;
	Lnext = L->next ;
	Lprev = L->prev ;

	L2p = L2->p ;
	L2i = L2->i ;
	L2x = L2->x ;
	L2z = L2->z ;
	L2nz = L2->nz ;
	L2next = L2->next ;
	L2prev = L2->prev ;
	L2->xtype = L->xtype ;
	L2->dtype = L->dtype ;

	for (j = 0 ; j <= n ; j++)
	{
	    L2p [j] = Lp [j] ;
	}

	for (j = 0 ; j < n+2 ; j++)
	{
	    L2prev [j] = Lprev [j] ;
	}

	for (j = 0 ; j < n+2 ; j++)
	{
	    L2next [j] = Lnext [j] ;
	}

	for (j = 0 ; j < n ; j++)
	{
	    L2nz [j] = Lnz [j] ;
	}

	for (j = 0 ; j < n ; j++)
	{
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;
	    for ( ; p < pend ; p++)
	    {
		L2i [p] = Li [p] ;
	    }
	    p = Lp [j] ;

	    if (L->xtype == SPARSE_REAL)
	    {
		for ( ; p < pend ; p++)
		{
		    L2x [p] = Lx [p] ;
		}
	    }
	}

    }
    else if (L->is_super)
    {

	/* ------------------------------------------------------------------ */
	/* 复制一个超节点因子 */
	/* ------------------------------------------------------------------ */

	xsize = L->xsize ;
	ssize = L->ssize ;
	nsuper = L->nsuper ;

	L2->xsize = xsize ;
	L2->ssize = ssize ;
	L2->nsuper = nsuper ;

	/* 指定 L2->super, L2->pi, L2->px, and L2->s.  如果L是熟知的，指定 L2->x if */
	if (!CORE(change_factor) (L->xtype, TRUE, TRUE, TRUE, TRUE, L2,
		    Common))
	{
	    CORE(free_factor) (&L2, Common) ;
	    return (NULL) ;	/* 内存溢出 */
	}

	/* ------------------------------------------------------------------ */
	/* 复制超节点因子的内容 */
	/* ------------------------------------------------------------------ */

	Lsuper = L->super ;
	Lpi = L->pi ;
	Lpx = L->px ;
	Ls = L->s ;
	Lx = L->x ;

	L2super = L2->super ;
	L2pi = L2->pi ;
	L2px = L2->px ;
	L2s = L2->s ;
	L2x = L2->x ;

	L2->maxcsize = L->maxcsize ;
	L2->maxesize = L->maxesize ;

	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2super [s] = Lsuper [s] ;
	}
	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2pi [s] = Lpi [s] ;
	}
	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2px [s] = Lpx [s] ;
	}

	L2s [0] = 0 ;
	for (p = 0 ; p < ssize ; p++)
	{
	    L2s [p] = Ls [p] ;
	}

	if (L->xtype == SPARSE_REAL)
	{
	    for (p = 0 ; p < xsize ; p++)
	    {
		L2x [p] = Lx [p] ;
	    }
	}
    }

    L2->minor = L->minor ;
    L2->is_monotonic = L->is_monotonic ;
    return (L2) ;
}
