/**
 * @file SparseChol_solve_change.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-22
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef NCHOLESKY

#include "Sparse_internal.h"
#include "SparseChol.h"

/* ========================================================================== */
/* === 模板 ============================================================= */
/* ========================================================================== */
#define REAL
#include "SparseChol_t_solve.c"

#define P(k) ((Perm == NULL) ? (k) : Perm [k])

/**
 * @brief 
 * 
 */
static void perm
(
    /* ---- input ---- */
    dense_array *B,	/* 输入矩阵B */
    Int *Perm,		    /* 可选的输入排列(可以为空) */
    Int k1,		        /* 待复制B的第一列 */
    Int ncols,		    /* 待拷贝的最后一列为min(k1+ncols,B->ncol)-1 */
    /* ---- in/out --- */
    dense_array *Y	/* 输出matrix Y,已经分配好 */
)
{
    double *Yx, *Yz, *Bx, *Bz ;
    Int k2, nk, p, k, j, nrow, ncol, d, dual, dj, j2 ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    ncol = B->ncol ;
    nrow = B->nrow ;
    k2 = MIN (k1+ncols, ncol) ;
    nk = MAX (k2 - k1, 0) ;
    dual = (Y->xtype == SPARSE_REAL && B->xtype != SPARSE_REAL) ? 2 : 1 ;
    d = B->d ;
    Bx = B->x ;
    Bz = B->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    Y->nrow = nrow ;
    Y->ncol = dual*nk ;
    Y->d = nrow ;
    

    /* ---------------------------------------------------------------------- */
    /* Y = B (P (1:nrow), k1:k2-1) */
    /* ---------------------------------------------------------------------- */

    switch (Y->xtype)
    {

	case SPARSE_REAL:

	    switch (B->xtype)
	    {

		case SPARSE_REAL:
		    /* Y real, B real */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [k + j2] = Bx [p] ;		/* real */
			}
		    }
		    break ;

	    }
	    break ;

    }
}

/**
 * @brief   X (P (1:nrow), k1 : min (k1+ncols,ncol)-1) = Y 其中X大小为nrow*ncol。
 *          
 *          复制和永久保存Y到X的一组连续列中。X已经在输入上分配了。
 *          Y的大小必须足够大。设nk为在X中访问的列数。X->xtype决定结果的复杂性。
 *          
 *          如果X是实的，而Y是复数，则只将B的实部复制到X中，忽略Y的虚部。
 * 
 *          如果X是复数(或zomplex)，而Y是实数，那么Y的实、虚和部分都返回到X中。
 *          Y是nrow2 * 2*nk。Y的偶数列包含B的实部，奇列包含B的虚部。
 *          Y->nzmax必须>= 2*nrow*nk。另一方面，Y是n乘nk，前导维数为nrow。
 *          Y->nzmax必须是>= nrow*nk。
 * 
 *          不使用输入(Y)为复数而输出(X)为实数的情况
 *          以及输入(Y)为zomplex而输出(X)为实数的情况。
 */
static void iperm
(
    /* ---- input ---- */
    dense_array *Y,	/* 输入矩阵Y */
    Int *Perm,		    /* 可选的输入排列(可以为空) */
    Int k1,		        /* 待拷贝的B的第一列 */
    Int ncols,		    /* 待拷贝的最后一列为min(k1+ncols,B->ncol)-1 */
    /* ---- in/out --- */
    dense_array *X	/* 输出矩阵X, 已经分配好 */
)
{
    double *Yx, *Yz, *Xx, *Xz ;
    Int k2, nk, p, k, j, nrow, ncol, d, dj, j2 ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    ncol = X->ncol ;
    nrow = X->nrow ;
    k2 = MIN (k1+ncols, ncol) ;
    nk = MAX (k2 - k1, 0) ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    

    /* ---------------------------------------------------------------------- */
    /* X (P (1:nrow), k1:k2-1) = Y */
    /* ---------------------------------------------------------------------- */

    switch (Y->xtype)
    {

	case SPARSE_REAL:

	    switch (X->xtype)
	    {

		case SPARSE_REAL:
		    /* Y real, X real */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [k + j2] ;		/* real */
			}
		    }
		    break ;
	    }
	    break ;

    }
}

/**
 * @brief Y = B (P (1:nrow), k1 : min (k1+ncols,ncol)-1)' 其中B大小为nrow*ncol.
 * 
 */
static void ptrans
(
    /* ---- input ---- */
    dense_array *B,	/* 输入矩阵B */
    Int *Perm,		    /* 可选的输入排列(可以为空) */
    Int k1,		        /* 待拷贝的B的第一列 */
    Int ncols,		    /* 待拷贝的最后一列min(k1+ncols,B->ncol)-1 */
    /* ---- in/out --- */
    dense_array *Y	/* 输出矩阵Y,已经分配好 */
)
{
    double *Yx, *Yz, *Bx, *Bz ;
    Int k2, nk, p, k, j, nrow, ncol, d, dual, dj, j2 ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    ncol = B->ncol ;
    nrow = B->nrow ;
    k2 = MIN (k1+ncols, ncol) ;
    nk = MAX (k2 - k1, 0) ;
    dual = (Y->xtype == SPARSE_REAL && B->xtype != SPARSE_REAL) ? 2 : 1 ;
    d = B->d ;
    Bx = B->x ;
    Bz = B->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    Y->nrow = dual*nk ;
    Y->ncol = nrow ;
    Y->d = dual*nk ;

    /* ---------------------------------------------------------------------- */
    /* Y = B (P (1:nrow), k1:k2-1)' */
    /* ---------------------------------------------------------------------- */

    switch (Y->xtype)
    {

	case SPARSE_REAL:

	    switch (B->xtype)
	    {

		case SPARSE_REAL:
		    /* Y real, B real  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2 + k*nk] = Bx [p] ;		/* real */
			}
		    }
		    break ;
	    }
	    break ;
    }
}

/**
 * @brief X (P (1:nrow), k1 : min (k1+ncols,ncol)-1) = Y'其中X大小为nrow*ncol.
 * 
 */
static void iptrans
(
    /* ---- input ---- */
    dense_array *Y,	/* 输入矩阵Y */
    Int *Perm,		    /* 可选的输入排列(可以为空) */
    Int k1,		        /* 复制到X的第一列 */
    Int ncols,		    /* 待拷贝的最后一列min(k1+ncols,X->ncol)-1 */
    /* ---- in/out --- */
    dense_array *X	/* output matrix X, already allocated */
)
{
    double *Yx, *Yz, *Xx, *Xz ;
    Int k2, nk, p, k, j, nrow, ncol, d, dj, j2 ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    ncol = X->ncol ;
    nrow = X->nrow ;
    k2 = MIN (k1+ncols, ncol) ;
    nk = MAX (k2 - k1, 0) ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    

    /* ---------------------------------------------------------------------- */
    /* X (P (1:nrow), k1:k2-1) = Y' */
    /* ---------------------------------------------------------------------- */

    switch (Y->xtype)
    {

	case SPARSE_REAL:

	    switch (X->xtype)
	    {

		case SPARSE_REAL:
		    /* Y real, X real  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [j2 + k*nk] ;		/* real */
			}
		    }
		    break ;
	    }
	    break ;

    }
}

/**
 * @brief   求解一个线性系统。
 *          
 *          因子分解可以是单纯的LDL'、单纯的LL'或超节点LL'。
 *          对于LL'分解，Dx=b解直接返回(它是隐式恒等)。
 */
dense_array *SparseChol_solve
(
    /* ---- input ---- */
    int sys,		    /* 求解系统 */
    sparse_factor *L,	/* 使用的分解 */
    dense_array *B,	/* right-hand-side */
    /* --------------- */
    sparse_common *Common
)
{
    dense_array *Y = NULL, *X = NULL ;
    dense_array *E = NULL ;
    int ok ;

    /* 求解，根据需要分配工作空间  */
    ok = SparseChol_solve2 (sys, L, B, NULL, &X, NULL, &Y, &E, Common) ;

    /* 如果分配了工作区，则释放;如果发生错误，则释放结果 */
    SparseCore_free_dense (&Y, Common) ;
    SparseCore_free_dense (&E, Common) ;
    if (!ok)
    {
        SparseCore_free_dense (&X, Common) ;
    }
    return (X) ;
}

/**
 * @brief 
 * 
 * @return int  成功返回TRUE，失败返回FLASE
 */
int SparseChol_solve2   
(
    /* ---- input ---- */
    int sys,		                /* 求解系统 */
    sparse_factor *L,	            /* 使用的分解 */
    dense_array *B,               /* right-hand-side */
    sparse_csc *Bset,
    /* ---- output --- */
    dense_array **X_Handle,       /* 解决方案，如有需要分配 */
    sparse_csc **Xset_Handle,
    /* ---- workspace  */
    dense_array **Y_Handle,       /* workspace, or NULL */
    dense_array **E_Handle,       /* workspace, or NULL */
    /* --------------- */
    sparse_common *Common
)
{
    double *Yx, *Yz, *Bx, *Bz, *Xx, *Xz ;
    dense_array *Y = NULL, *X = NULL ;
    sparse_csc *C, *Yset, C_header, Yset_header, *Xset ;
    Int *Perm = NULL, *IPerm = NULL ;
    Int n, nrhs, ncols, xtype, k1, nr, ytype, k, blen, p, i, d, nrow ;
    Int Cp [2], Ysetp [2], *Ci, *Yseti, ysetlen ;
    Int *Bsetp, *Bseti, *Bsetnz, *Xseti, *Xsetp, *Iwork ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (sys < SPARSE_A || sys > SPARSE_Pt)
    {
	ERROR (SPARSE_INVALID, "invalid system") ;
	return (FALSE) ;
    }
    
    nrhs = B->ncol ;
    n = (Int) L->n ;
    d = (Int) B->d ;
    nrow = (Int) B->nrow ;
    if (d < n || nrow != n)
    {
	ERROR (SPARSE_INVALID, "dimensions of L and B do not match") ;
	return (FALSE) ;
    }
    if (Bset)
    {
        if (nrhs != 1)
        {
            ERROR (SPARSE_INVALID, "Bset requires a single right-hand side") ;
            return (FALSE) ;
        }
        if (L->xtype != B->xtype)
        {
            ERROR (SPARSE_INVALID, "Bset requires xtype of L and B to match") ;
            return (FALSE) ;
        }
        
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    if ((sys == SPARSE_P || sys == SPARSE_Pt || sys == SPARSE_A)
	    && L->ordering != SPARSE_NATURAL)
    {
        /* 否则，Perm为NULL，使用恒等置换 */
	Perm = L->Perm ;
    }

    /* ---------------------------------------------------------------------- */
    /* 分配结果X(或重新使用之前调用的空间) */
    /* ---------------------------------------------------------------------- */

    if (Bset)
    {
        xtype = L->xtype ;
    }
    else if (sys == SPARSE_P || sys == SPARSE_Pt)
    {
	/* 如果B是实数，则x=Pb和x=P'b返回x实数;
     * 如果B是complex或zomplex, X是首选的complex/zcomplex类型 */
	
    xtype = SPARSE_REAL ;
	}
    else if (L->xtype == SPARSE_REAL && B->xtype == SPARSE_REAL)
    {
	/* 若L和B为实数，则X为实数 */
	xtype = SPARSE_REAL ;
    }
    else
    {
	/* X是complex，使用首选的complex/zomplex类型 */
	
    xtype = SPARSE_REAL ;
	}

    /* 确保X有正确的大小和类型 */
    X = SparseCore_ensure_dense (X_Handle, n, nrhs, n, xtype, Common) ; //在hnucore_matrix_type.c 里面
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 使用L, D, L', P或者他们的组合求解 */
    /* ---------------------------------------------------------------------- */

    if (Bset)
    {

        /* ------------------------------------------------------------------ */
        /* 用一个稀疏的b，求解一个x的子集 */
        /* ------------------------------------------------------------------ */

        Int save_realloc_state ;

#ifndef NSUPERNODAL
        /* 使用Bset时，将超节点L转换为单纯形 */
        if (L->is_super)
        {
            /* 只能用Bset对一个单纯因子分解。超节点因子L被转换为单纯型，
             * 而x类型保持不变(实型、复型或复型)。因为超节点因数分解已经是LL'了，
             * 所以它就保留这种形式。该转换使用SparseCore_change_factor中
             * 的ll_super_to_simplicial_numeric函数。
             */
            SparseCore_change_factor (   // 在 hnucore_change_factor.c 里面
                SPARSE_REAL,   /* ignored, since L is already numeric */
                TRUE,           /* convert to LL' (no change to num. values) */
                FALSE,          /* convert to simplicial */
                FALSE,          /* do not pack the columns of L */
                FALSE,          /* (ignored) */
                L, Common) ;
            if (Common->status < SPARSE_OK)
            {
                /* 内存不足，L将原样返回 */
                return (FALSE) ;
            }
        }
#endif

        /* L, X, B是相同的数据类型*/
        /* 确保Y的大小正确 */
	Y = SparseCore_ensure_dense (Y_Handle, 1, n, 1, L->xtype, Common) ;
	if (Common->status < SPARSE_OK)
	{
	    /* 内存溢出 */
	    return (FALSE) ;
	}

        /* ------------------------------------------------------------------ */
        /* 求逆置换，如果需要的话，构造它 */
        /* ------------------------------------------------------------------ */


        if ((sys == SPARSE_A || sys == SPARSE_P) && Perm != NULL)
        {
            /* 在求解Ax=b或x=Pb时，用逆置换IPerm求解c=Pb步骤。
             * 没有其他步骤应该使用IPerm */
            if (L->IPerm == NULL)
            {
                /* 构造逆置换。这只执行一次，然后永久存储在L中。  */
                L->IPerm = SparseCore_malloc (n, sizeof (Int), Common) ;
                if (Common->status < SPARSE_OK)
                {
                    /* 内存溢出 */
                    return (FALSE) ;
                }
                IPerm = L->IPerm ;
                for (k = 0 ; k < n ; k++)
                {
                    IPerm [Perm [k]] = k ;
                }
            }
            /* x=A\b和x=Pb都需要IPerm */
            IPerm = L->IPerm ;
        }

        if (sys == SPARSE_P)
        {
            /* x=Pb需要关闭后续的x=P'b置换 */
            Perm = NULL ;
        }

        /* ------------------------------------------------------------------ */
        /* 确保Xset的类型和大小正确 */
        /* ------------------------------------------------------------------ */

        /* Xset 大小为 n*1, nzmax >= n, pattern-only, packed, unsorted */
        Xset = *Xset_Handle ;
        if (Xset == NULL || (Int) Xset->nrow != n || (Int) Xset->ncol != 1 ||
            (Int) Xset->nzmax < n || Xset->itype != SPARSE_PATTERN)
        {
            /* 这只在第一次调用SparseChol_solve时执行一次 */
            SparseCore_free_sparse (Xset_Handle, Common) ;
            Xset = SparseCore_allocate_sparse (n, 1, n, FALSE, TRUE, 0,
                SPARSE_PATTERN, Common) ;
            *Xset_Handle = Xset ;
        }
        Xset->sorted = FALSE ;
        Xset->stype = 0 ;
        if (Common->status < SPARSE_OK)
        {
            /* 内存溢出 */
            return (FALSE) ;
        }

        /* -------------------------------------------------------------- */
        /* 确保Flag大小为n,而且有3*n Int的可用工作空间 */
        /* -------------------------------------------------------------- */

        /* 如果之前的调用已经分配了足够的空间，是否不工作 */
        SparseCore_allocate_work (n, 3*n, 0, Common) ;
        if (Common->status < SPARSE_OK)
        {
            /* 内存溢出 */
            return (FALSE) ;
        }

        /* Ci和Yseti使用Iwork (n:3n-1) */
        Iwork = Common->Iwork ;
        /* Iwork (0:n-1)没有被使用，因为它被check_perm、print_perm、check_sparse和print_sparse使用 */
        Ci = Iwork + n ;
        Yseti = Ci + n ;

        /* 重新分配工作空间会破坏Ci和Yseti */
        save_realloc_state = Common->no_workspace_reallocate ;
        Common->no_workspace_reallocate = TRUE ;

        /* -------------------------------------------------------------- */
        /* C = permuted Bset,对应于L的排列 */
        /* -------------------------------------------------------------- */

        /* C = IPerm (Bset) */

        Bsetp = Bset->p ;
        Bseti = Bset->i ;
        Bsetnz = Bset->nz ;
        blen = (Bset->packed) ? Bsetp [1] : Bsetnz [0] ;

        /* C = spones (P*B) or C = spones (B) if IPerm is NULL */
        C = &C_header ;
        C->nrow = n ;
        C->ncol = 1 ;
        C->nzmax = n ;
        C->packed = TRUE ;
        C->stype = 0 ;
        C->itype = ITYPE ;
        C->xtype = SPARSE_PATTERN ;
        C->dtype = SPARSE_DOUBLE ;
        C->nz = NULL ;
        C->p = Cp ;
        C->i = Ci ;
        C->x = NULL ;
        C->z = NULL ;
        C->sorted = FALSE ;
        Cp [0] = 0 ;
        Cp [1] = blen ;
        for (p = 0 ; p < blen ; p++)
        {
            Int iold = Bseti [p] ;
            Ci [p] = IPerm ? IPerm [iold] : iold ;
        }

        /* 从Iwork(n:2n-1)创建一个稀疏列Yset  */
        Yset = &Yset_header ;
        Yset->nrow = n ;
        Yset->ncol = 1 ;
        Yset->nzmax = n ;
        Yset->packed = TRUE ;
        Yset->stype = 0 ;
        Yset->itype = ITYPE ;
        Yset->xtype = SPARSE_PATTERN ;
        Yset->dtype = SPARSE_DOUBLE ;
        Yset->nz = NULL ;
        Yset->p = Ysetp ;
        Yset->i = Yseti ;
        Yset->x = NULL ;
        Yset->z = NULL ;
        Yset->sorted = FALSE ;
        Ysetp [0] = 0 ;
        Ysetp [1] = 0 ;

        /* -------------------------------------------------------------- */
        /* Yset = nonzero pattern of L\C, or just C itself */
        /* -------------------------------------------------------------- */

        /* 这将花费O(ysetlen)的时间  */
        if (sys == SPARSE_P || sys == SPARSE_Pt || sys == SPARSE_D)
        {
            Ysetp [1] = blen ;
            for (p = 0 ; p < blen ; p++)
            {
                Yseti [p] = Ci [p] ;
            }
        }
        else
        {
            if (!SparseChol_lsolve_pattern (C, L, Yset, Common)) //SparseChol_analyze.c中
            {
                Common->no_workspace_reallocate = save_realloc_state ;
                return (FALSE) ;
            }
        }

        /* -------------------------------------------------------------- */
        /* 清除Y中用于求解的部分 */
        /* -------------------------------------------------------------- */

        Yx = Y->x ;
        Yz = Y->z ;
        ysetlen = Ysetp [1] ;

        switch (L->xtype)
        {

            case SPARSE_REAL:
                for (p = 0 ; p < ysetlen ; p++)
                {
                    i = Yseti [p] ;
                    Yx [i] = 0 ;
                }
                break ;
        }

        /* -------------------------------------------------------------- */
        /* 分散和排列B成Y */
        /* -------------------------------------------------------------- */

        /* Y (C) = B (Bset) */
        Bx = B->x ;
        Bz = B->z ;

        switch (L->xtype)
        {

            case SPARSE_REAL:
                for (p = 0 ; p < blen ; p++)
                {
                    Int iold = Bseti [p] ;
                    Int inew = Ci [p] ;
                    Yx [inew] = Bx [iold] ;
                }
                break ;
        }

        /* -------------------------------------------------------------- */
        /* 求解Y = (L' \ (L \ Y'))'或使用例程的其他分解 */
        /* -------------------------------------------------------------- */

        /* 这个解决方案在Yseti[0...ysetlen-1]中只迭代列 */

        if (! (sys == SPARSE_P || sys == SPARSE_Pt))
        {
            switch (L->xtype)
            {
                case SPARSE_REAL:
                    r_simplicial_solver (sys, L, Y, Yseti, ysetlen) ; //t_SparseChol_solve.c
                    break ;

            }
        }

        /* -------------------------------------------------------------- */
        /* X = P'*Y，但只适用于Yset中的行，并创建Xset */
        /* -------------------------------------------------------------- */

        /* X (Perm (Yset)) = Y (Yset) */
        Xx = X->x ;
        Xz = X->z ;
        Xseti = Xset->i ;
        Xsetp = Xset->p ;

        switch (L->xtype)
        {

            case SPARSE_REAL:
                for (p = 0 ; p < ysetlen ; p++)
                {
                    Int inew = Yseti [p] ;
                    Int iold = Perm ? Perm [inew] : inew ;
                    Xx [iold] = Yx [inew] ;
                    Xseti [p] = iold ;
                }
                break ;
        }

        Xsetp [0] = 0 ;
        Xsetp [1] = ysetlen ;

        
        Common->no_workspace_reallocate = save_realloc_state ;

    }
    else if (sys == SPARSE_P)
    {

	/* ------------------------------------------------------------------ */
	/* x = P*b */
	/* ------------------------------------------------------------------ */

	perm (B, Perm, 0, nrhs, X) ;

    }
    else if (sys == SPARSE_Pt)
    {

	/* ------------------------------------------------------------------ */
	/* x = P'*b */
	/* ------------------------------------------------------------------ */

	iperm (B, Perm, 0, nrhs, X) ;

    }
    else if (L->is_super)
    {

	/* ------------------------------------------------------------------ */
	/* 使用超节点LL'分解求解 */
	/* ------------------------------------------------------------------ */

#ifndef NSUPERNODAL
	/* 分配工作空间 */
	dense_array *E ;
	Int dual ;
        Common->blas_ok = TRUE ;
	dual = (L->xtype == SPARSE_REAL && B->xtype != SPARSE_REAL) ? 2 : 1 ;
	Y = SparseCore_ensure_dense (Y_Handle, n, dual*nrhs, n, L->xtype, Common);
	E = SparseCore_ensure_dense (E_Handle, dual*nrhs, L->maxesize, dual*nrhs,
		L->xtype, Common) ;

	if (Common->status < SPARSE_OK)
	{
	    /* 内存溢出 */
            return (FALSE) ;
	}

	perm (B, Perm, 0, nrhs, Y) ;			    /* Y = P*B */

	if (sys == SPARSE_A || sys == SPARSE_LDLt)
	{
	    SparseChol_super_lsolve (L, Y, E, Common) ;	    /* Y = L\Y */
            SparseChol_super_ltsolve (L, Y, E, Common) ;	    /* Y = L'\Y*/
	}
	else if (sys == SPARSE_L || sys == SPARSE_LD)
	{
	    SparseChol_super_lsolve (L, Y, E, Common) ;	    /* Y = L\Y */
	}
	else if (sys == SPARSE_Lt || sys == SPARSE_DLt)
	{
	    SparseChol_super_ltsolve (L, Y, E, Common) ;      /* Y = L'\Y*/
	}

	iperm (Y, Perm, 0, nrhs, X) ;			    /* X = P'*Y */

	if (CHECK_BLAS_INT && !Common->blas_ok)
	{
	    /* BLAS中的整型溢出。这可能是不可能的，因为BLAS被用来创建超节点因数分解。
         * 不过，对BLAS的调用可能会在分解和转发/回退之间有所不同。
         * 这种说法未经检验;如果CHECK_BLAS_INT为真(当在HNUCHOL和BLAS中
         * 使用相同的整数时)，则它不会出现在编译的代码中。 */
	    return (FALSE) ;
	}

#else
	/* HNUCHOL Supernodal模块没有安装 */
	ERROR (SPARSE_NOT_INSTALLED,"Supernodal module not installed") ;
#endif

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* 使用简单的LL'或LDL'分解来解决 */
	/* ------------------------------------------------------------------ */

        if (L->xtype == SPARSE_REAL && B->xtype == SPARSE_REAL)
	{
	    /* L, B, Y都是实数 */
	    /* 一次最多解决4列B */
            ncols = 4 ;
            nr = MAX (4, nrhs) ;
	    ytype = SPARSE_REAL ;
	}
	else if (L->xtype == SPARSE_REAL)
	{
            /* L是实数型，B是复数型 */
	    /* 用B的一列(real/imag)一次求解 */
	    ncols = 1 ;
	    nr = 2 ;
	    ytype = SPARSE_REAL ;
	}
	else
	{
	    /* L是复数型, B是real/complex/zomplex, 
         * Y和L一样.  一次解B的一列。 */
	    ncols = 1 ;
	    nr = 1 ;
	    ytype = L->xtype ;
	}

	Y = SparseCore_ensure_dense (Y_Handle, nr, n, nr, ytype, Common) ;
	if (Common->status < SPARSE_OK)
	{
	    /* 内存溢出 */
	    return (FALSE) ;
	}

        for (k1 = 0 ; k1 < nrhs ; k1 += ncols)
        {

            /* -------------------------------------------------------------- */
            /* Y = B (P, k1:k1+ncols-1)' = (P * B (:,...))' */
            /* -------------------------------------------------------------- */

            ptrans (B, Perm, k1, ncols, Y) ;

            /* -------------------------------------------------------------- */
            /* 求解Y = (L' \ (L \ Y'))',或者模板的其他分解 */
            /* -------------------------------------------------------------- */

            switch (L->xtype)
            {
                case SPARSE_REAL:
                    r_simplicial_solver (sys, L, Y, NULL, 0) ;
                    break ;
            }

            /* -------------------------------------------------------------- */
            /* X (P, k1:k2+ncols-1) = Y' */
            /* -------------------------------------------------------------- */

            iptrans (Y, Perm, k1, ncols, X) ;
        }
    }

    return (TRUE) ;
}

#endif

