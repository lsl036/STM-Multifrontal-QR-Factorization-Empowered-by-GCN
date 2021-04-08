/**
 * @file SparseChol_analyze.c
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
#include "amd.h"
#include "colamd.h"
#include "SparseChol.h"

/**
 * @brief   置换并转置矩阵。 如果需要，分配A1和A2矩阵，
 *          如果不需要，则将它们返回为NULL。
 */
static int permute_matrices
(
    /* ---- input ---- */
    sparse_csc *A,	/* 置换矩阵 */
    Int ordering,	    /* 使用的排序方法 */
    Int *Perm,		    /* 减少填充的排列 */
    Int *fset,		    /* 0:(A->ncol)-1的子集 */
    size_t fsize,	    /* fset的大小 */
    Int do_rowcolcounts,/* 如果为TRUE，则同时计算S和F。
                         * 如果为FALSE，则对称情况仅需要S，非对称情况仅需要F */
    /* ---- output --- */
    sparse_csc **A1_handle,	    /* 请参阅以下有关A1，A2，S，F的注释 */
    sparse_csc **A2_handle,
    sparse_csc **S_handle,
    sparse_csc **F_handle,
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *A1, *A2, *S, *F ;

    *A1_handle = NULL ;
    *A2_handle = NULL ;
    *S_handle = NULL ;
    *F_handle = NULL ;
    A1 = NULL ;
    A2 = NULL ;

    if (ordering == SPARSE_NATURAL)
    {

	/* ------------------------------------------------------------------ */
	/* A的自然排序 */
	/* ------------------------------------------------------------------ */

	if (A->stype < 0)
	{
	    /* 对称下三角：A已经是下三角形式，所以S = A' */
	    /* workspace: Iwork (nrow) */
	    A2 = SparseCore_ptranspose (A, 0, NULL, NULL, 0, Common) ;
	    F = A ;
	    S = A2 ;
	}
	else if (A->stype > 0)
	{
	    /* 对称上三角形式: F = pattern of triu (A)', S = A */
	    /* workspace: Iwork (nrow) */
	    if (do_rowcolcounts)
	    {
		/* 如果do_rowcolcounts为FALSE，则在对称情况下不需要F */
		A1 = SparseCore_ptranspose (A, 0, NULL, fset, fsize, Common) ;
	    }
	    F = A1 ;
	    S = A ;
	}
	else
	{
	    /* 非对称情况: F = pattern of A (:,f)',  S = A */
	    /* workspace: Iwork (nrow if no fset, MAX(nrow,ncol) if fset) */
	    A1 = SparseCore_ptranspose (A, 0, NULL, fset, fsize, Common) ;
	    F = A1 ;
	    S = A ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* A被置换 */
	/* ------------------------------------------------------------------ */

	if (A->stype < 0)
	{
	    /* 对称下三角形式: S = tril (A (p,p))' and F = S' */
	    /* workspace: Iwork (2*nrow) */
	    A2 = SparseCore_ptranspose (A, 0, Perm, NULL, 0, Common) ;
	    S = A2 ;
	    /* workspace: Iwork (nrow) */
	    if (do_rowcolcounts)
	    {
		/* 如果do_rowcolcounts为FALSE，则在对称情况下不需要F */
		A1 = SparseCore_ptranspose (A2, 0, NULL, NULL, 0, Common) ;
	    }
	    F = A1 ;
	}
	else if (A->stype > 0)
	{
	    /* 对称上三角情况: F = triu (A (p,p))' and S = F' */
	    /* workspace: Iwork (2*nrow) */
	    A1 = SparseCore_ptranspose (A, 0, Perm, NULL, 0, Common) ;
	    F = A1 ;
	    /* workspace: Iwork (nrow) */
	    A2 = SparseCore_ptranspose (A1, 0, NULL, NULL, 0, Common) ;
	    S = A2 ;
	}
	else
	{
	    /* 非对称情况:     F = A (p,f)'         and S = F' */
	    /* workspace: Iwork (nrow if no fset, MAX(nrow,ncol) if fset) */
	    A1 = SparseCore_ptranspose (A, 0, Perm, fset, fsize, Common) ;
	    F = A1 ;
	    if (do_rowcolcounts)
	    {
		/* 如果do_rowcolcounts为FALSE，则对于非对称情况不需要S */
		/* workspace: Iwork (nrow) */
		A2 = SparseCore_ptranspose (A1, 0, NULL, NULL, 0, Common) ;
	    }
	    S = A2 ;
	}
    }

    /* 如果任何SparseCore_*transpose失败，则一个或多个矩阵将为NULL */
    *A1_handle = A1 ;
    *A2_handle = A2 ;
    *S_handle = S ;
    *F_handle = F ;
    return (Common->status == SPARSE_OK) ;
}

/**
 * @brief   给定矩阵A及其填充减少的排列，计算消去树，（非加权）后序遍历以及
 *          L中每列的非零元素数目。还计算flop计数，L中的总非零元素和A中的非零元素
 *          (Common->fl, Common->lnz, and Common->anz).
 * 
 */
int SparseChol_analyze_ordering
(
    /* ---- input ---- */
    sparse_csc *A,	/* 分析的矩阵 */
    int ordering,	    /* 使用的排序方法 */
    Int *Perm,		    /* 大小为n, 减少填充的排列来分析 */
    Int *fset,		    /* 0:(A->ncol)-1的子集 */
    size_t fsize,	    /* fset的大小 */
    /* ---- output --- */
    Int *Parent,	    /* size n, 消去树 */
    Int *Post,		    /* size n, 消去树的后序遍历 */
    Int *ColCount,	    /* size n, L中每列的nnz */
    /* ---- workspace  */
    Int *First,		    /* SparseCore_postorder的大小为n的工作区 */
    Int *Level,		    /* SparseCore_postorder的大小为n的工作区 */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *A1, *A2, *S, *F ;
    Int n, ok, do_rowcolcounts ;

    /* 检查输入 */
    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;

    n = A->nrow ;

    do_rowcolcounts = (ColCount != NULL) ;

    /* 根据Perm和fset置换A */
    ok = permute_matrices (A, ordering, Perm, fset, fsize, do_rowcolcounts,
	    &A1, &A2, &S, &F, Common) ;

    /* 查找S（对称的上三角/下三角）或F（非对称）的消去树 */
    /* 工作空间：对称: Iwork (nrow), unsym: Iwork (nrow+ncol) 构造树*/
    ok = ok && SparseChol_etree (A->stype ? S:F, Parent, Common) ;

    /* 后序遍历消去树（由SparseCore_rowcolcounts要求） */
    /* workspace: Iwork (2*nrow) */
    ok = ok && (SparseChol_postorder (Parent, n, NULL, Post, Common) == n) ;

    /* 如果SparseCore_postorder的返回<n，则不会设置Common-> status */
    Common->status = (!ok && Common->status == SPARSE_OK) ?
	SPARSE_INVALID : Common->status ;

    /* 分析LL'=S or SS' or S(:,f)*S(:,f)' */
    /* workspace:
     *	对称:   Flag (nrow), Iwork (2*nrow)
     *	非对称: Flag (nrow), Iwork (2*nrow+ncol), Head (nrow+1)
     */
    if (do_rowcolcounts)
    {
	ok = ok && SparseChol_rowcolcounts (A->stype ? F:S, fset, fsize, Parent,
	    Post, NULL, ColCount, First, Level, Common) ;
    }

    /* 释放临时矩阵和返回结果 */
    SparseCore_free_sparse (&A1, Common) ;
    SparseCore_free_sparse (&A2, Common) ;
    return (ok) ;
}


/* ========================================================================== */
/* === 释放工作空间和返回结果L ========================================== */
/* ========================================================================== */

#define FREE_WORKSPACE_AND_RETURN \
{ \
    Common->no_workspace_reallocate = FALSE ; \
    SparseCore_free (n, sizeof (Int), Lparent,  Common) ; \
    SparseCore_free (n, sizeof (Int), Perm,     Common) ; \
    SparseCore_free (n, sizeof (Int), ColCount, Common) ; \
    if (Common->status < SPARSE_OK) \
    { \
	SparseCore_free_factor (&L, Common) ; \
    } \
    return (L) ; \
}

/**
 * @brief   默认使用AMD重排序，置换，生成消去树，符号分解，返回符号分解因子
 * 
 */
sparse_factor *SparseChol_analyze
(
    /* ---- input ---- */
    sparse_csc *A,	/* 进行排序和分析的矩阵 */
    /* --------------- */
    sparse_common *Common,
    char *result_name
)
{
    // Int * Perm;
    // Int n = A->nrow ;
    // Perm     = SparseCore_malloc (n, sizeof (Int), Common) ;
    // FILE *fp;
    // //fp = fopen("./Perm_array/s3rmt3m3.txt","r");
    // fp = fopen(result_name,"r");
    // if(fp == NULL) {
    //     printf("file is not exist!\n");
    //     exit(-1);
    // }
    // char buf[1024];
    // for (int k = 0 ; k < n ; k++)
    // {
    //     if ( fgets(buf, 1024, fp) == NULL )
    //     {
    //         break;
    //     }
    //     sscanf(buf, "%d", &(Perm[k]) );
    //     /* 在SparseCore_ptranspose中检查UserPerm */
    //     // Perm [k] = Useperm;
        
    // }
    // // for (k = 0 ; k < n ; k++)
    // // {
    // //     printf ("%d ", Perm[k]);
    // // }
    // printf("\n");
    // fclose(fp);
    // return (SparseChol_analyze_p2 (TRUE, A, Perm, NULL, 0, Common, result_name)) ;
    return (SparseChol_analyze_p2 (TRUE, A, NULL, NULL, 0, Common, result_name)) ;
}

/**
 * @brief   稀疏Cholesky或稀疏QR的排序和分析
 * 
 */
sparse_factor *SparseChol_analyze_p2
(
    /* ---- input ---- */
    int for_whom,       /* FOR_SPQR     (0): for SPQR
                           FOR_CHOLESKY (1): for Cholesky */
    sparse_csc *A,	/* 进行排序和分析的矩阵 */
    Int *UserPerm,	    /* 用户提供的排列size A->nrow */
    Int *fset,		    /* 0:(A->ncol)-1的子集 */
    size_t fsize,	    /* fset的大小 */
    /* --------------- */
    sparse_common *Common,
    char *result_name
)
{
    double lnz_best ;
    Int *First, *Level, *Work4n, *Cmember, *CParent, *ColCount, *Lperm, *Parent,
	*Post, *Perm, *Lparent, *Lcolcount ;
    sparse_factor *L ;
    Int k, n, ordering, method, nmethods, status, default_strategy, ncol, uncol,
	skip_analysis, skip_best ;
    Int amd_backup ;
    size_t s ;
    int ok = TRUE ;
    #ifdef all_methods_time
    FILE* fresult;
    fresult = fopen("./Results/Brute-force-fill.txt","a+");
    #endif
    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    Common->status = SPARSE_OK ;
    status = SPARSE_OK ;
    Common->selected = EMPTY ;
    Common->called_nd = FALSE ;
    /* ---------------------------------------------------------------------- */
    /* 获取输入参数*/
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    ncol = A->ncol ;
    uncol = (A->stype == 0) ? (A->ncol) : 0 ;

    /* ---------------------------------------------------------------------- */
    /* 设定默认策略 */
    /* ---------------------------------------------------------------------- */

    lnz_best = (double) EMPTY ;
    skip_best = FALSE ;
    nmethods = MIN (Common->nmethods, SPARSE_MAXMETHODS) ;  // Common->nmethods = 0
    nmethods = MAX (0, nmethods) ;

    default_strategy = (nmethods == 0) ;

    if (for_whom )  // Cholesky 直接选择这种方法
    {
        //Common->method [0].ordering = SPARSE_NATURAL;
        // Common->method [0].ordering = SPARSE_GIVEN;
        Common->method [0].ordering = SPARSE_AMD;
        //Common->method [0].ordering = SPARSE_NESDIS;
        // Common->method [0].ordering = SPARSE_COLAMD;
        //Common->method [0].ordering = SPARSE_METIS;
        #ifdef ONLY_METIS
        Common->method [0].ordering = SPARSE_METIS;
        #endif
        #ifdef ONLY_COLAMD
        Common->method [0].ordering = SPARSE_COLAMD;
        #endif
        #ifdef ONLY_AMD
        Common->method [0].ordering = SPARSE_AMD;
        #endif
        amd_backup = FALSE ;
        nmethods = 1 ;

    #ifdef all_methods_time
        //QR统计总时间
        nmethods = 4 ;
        Common->method [0].ordering = SPARSE_AMD ;
        Common->method [1].ordering = SPARSE_COLAMD;
        Common->method [2].ordering = SPARSE_METIS ;
        Common->method [3].ordering = SPARSE_NESDIS;
    #endif
    }
    else //QR default
    {
        /* 如果尝试2种或更多方法，请启用AMD备份 */
        amd_backup = (nmethods > 1) ;   // QR 为 0
        
    }

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间 */
    /* ---------------------------------------------------------------------- */

    /* 注意：此处需要分配足够的空间，以便SparseChol_analyze调用的例程不会重新分配空间 */

    /* s = 6*n + uncol */
    s = SparseCore_mult_size_t (n, 6, &ok) ;
    s = SparseCore_add_size_t (s, uncol, &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }

    SparseCore_allocate_work (n, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (NULL) ;	    /* 内存溢出 */
    }

    /* 确保由SparseChol_analyze调用的后续例程不
     * 重新分配任何工作区。 在FREE_WORKSPACE_AND_RETURN宏中将
     * 其设置回FALSE，这是此函数返回其调用方的唯一方法。 */
    Common->no_workspace_reallocate = TRUE ;

    /* 从Iwork中使用最后4 * n个Int作为Parent，First，Level和Post，因为
     * 其他例程将使用前2n + uncol空间。*/
    Work4n = Common->Iwork ;
    Work4n += 2*((size_t) n) + uncol ;
    Parent = Work4n ;
    First  = Work4n + n ;
    Level  = Work4n + 2*((size_t) n) ;
    Post   = Work4n + 3*((size_t) n) ;

    Cmember = Post ;
    CParent = Level ;

    /* ---------------------------------------------------------------------- */
    /* 分配更多的工作空间和一个空的简单符号因子 */
    /* ---------------------------------------------------------------------- */

    L = SparseCore_allocate_factor (n, Common) ;
    Lparent  = SparseCore_malloc (n, sizeof (Int), Common) ;
    Perm     = SparseCore_malloc (n, sizeof (Int), Common) ;
    ColCount = SparseCore_malloc (n, sizeof (Int), Common) ;
    if (Common->status < SPARSE_OK)
    {
	/* 内存溢出 */
	FREE_WORKSPACE_AND_RETURN ;
    }
    Lperm = L->Perm ;
    Lcolcount = L->ColCount ;
    Common->anz = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* 尝试所有要求的排序选项， 如果需要的话备份到AMD */
    /* ---------------------------------------------------------------------- */

    /* 关闭错误处理 */
    Common->try_catch = TRUE ;
    for (method = 0 ; method <= nmethods ; method++) // 开始尝试所有的方法
    {
        /* ------------------------------------------------------------------ */
        /* 确定尝试的方法 */
        /* ------------------------------------------------------------------ */
        //printf("nmethods = %d\n", nmethods);
        Common->fl = EMPTY ;
        Common->lnz = EMPTY ;
        skip_analysis = FALSE ;

        if (method == nmethods)
        {
            /* 所有方法均失败：备份到AMD */
            if (Common->selected == EMPTY && amd_backup)
            {
            ordering = SPARSE_AMD ;
            }
            else
            {
            break ;
            }
        }
        else
        {
            ordering = Common->method [method].ordering ;
        }
        Common->current = method ;

        /* ------------------------------------------------------------------ */
        /* 寻找填充简化排序的方法 */
        /* ------------------------------------------------------------------ */

        if (ordering == SPARSE_AMD)
        {
            printf("AMD USED\n");
            /* -------------------------------------------------------------- */
            /* A，A * A'或AN（：，f）* A（：，f）'的AMD顺序 */
            /* -------------------------------------------------------------- */
            // printf ("1: SparseChol_amd used \n");
            amd_backup = FALSE ;    /* 不需要尝试AMD两次 */
            SparseChol_amd (A, fset, fsize, Perm, Common) ;
            skip_analysis = TRUE ;
        }
        else if (ordering == SPARSE_NATURAL)
        {
            // printf("NATURAL USED\n");
            /* -------------------------------------------------------------- */
            /* 自然排序 */
            /* -------------------------------------------------------------- */

            for (k = 0 ; k < n ; k++)
            {
            Perm [k] = k ;
            }

        }
        else if (ordering == SPARSE_GIVEN)
        {
            // printf("GIVEN USED\n");
            /* -------------------------------------------------------------- */
            /* 如果提供的话，使用给定的A的排序 */
            /* -------------------------------------------------------------- */
            if (UserPerm == NULL)
            {
            /* 这不是错误情况 */
            continue ;
            }
            // printf("n = %ld \n", n);
            for (k = 0 ; k < n ; k++)
            {
            /* 在SparseCore_ptranspose中检查UserPerm */
            Perm [k] = UserPerm [k] ;
            }
        }
        else if (ordering == SPARSE_COLAMD)
        {
            printf("COLAMD USED\n");
            /* -------------------------------------------------------------- */
            /* AMD用于对称情况，COLAMD用于A * A'或A（：，f）* A（：，f）' */
            /* -------------------------------------------------------------- */
            if (A->stype)
            {
                // printf ("2: SparseChol_amd used \n");
                SparseChol_amd (A, fset, fsize, Perm, Common) ;
                skip_analysis = TRUE ;
            }
            else
            {
                /* 不要后序遍历，稍后再做 */
                /* workspace: Iwork (4*nrow+uncol), Flag (nrow), Head (nrow+1)*/
                // printf ("3: SparseChol_colamd used \n");
                SparseChol_colamd (A, fset, fsize, FALSE, Perm, Common) ;
            }
        }
        else if (ordering == SPARSE_METIS)
        {
            printf("METIS USED\n");
        #ifndef NMETIS
            Common->called_nd = TRUE;
            SparseCore_metis (A, fset, fsize, FALSE, Perm, Common);
        #else
            Common->status = SPARSE_NOT_INSTALLED ;
        #endif
        }
        else if (ordering == SPARSE_NESDIS)
        {
            printf("NESDIS USED\n");
        
            Common->called_nd = TRUE;
            SparseCore_nested_dissection (A, fset, fsize, Perm, CParent, Cmember,
		    Common) ;
        }
        else
        {
            /* -------------------------------------------------------------- */
            /* 无效的排序方法 */
            /* -------------------------------------------------------------- */
            Common->status = SPARSE_INVALID ;
        }

        if (Common->status < SPARSE_OK)
        {
            /* 内存溢出或方法失败 */
            status = MIN (status, Common->status) ;
            Common->status = SPARSE_OK ;
            continue ;
        }

        /* ------------------------------------------------------------------ */
        /* 分析排序 */
        /* ------------------------------------------------------------------ */

        if (!skip_analysis)
        {
            // printf ("   || analyze_ordering ... || \n");
            if (!SparseChol_analyze_ordering (A, ordering, Perm, fset, fsize,
                Parent, Post, ColCount, First, Level, Common))
            {
            /* 排序方法失败； 清除状态并尝试下一种方法 */
            status = MIN (status, Common->status) ;
            Common->status = SPARSE_OK ;
            continue ;
            }
        }

        // 保存当前重排序方法的 数据fl 和 lnz
        Common->method [method].fl  = Common->fl ;
        Common->method [method].lnz = Common->lnz ;
        printf ("lnz = %lf  \n", Common->lnz );
        #ifdef all_methods_time
        if (for_whom == 1)
            fprintf (fresult,"%lf ",Common->lnz);
        #endif
        /* ------------------------------------------------------------------ */
        /* 选择最好的方法 通过比较flop count ：fl */
        /* ------------------------------------------------------------------ */

        /* fl.pt. 比较，但lnz永远不能是NaN */
        if (Common->selected == EMPTY || Common->lnz < lnz_best)
        {
            //printf ("SAVE ORDERING %d ...\n", ordering);
            Common->selected = method ;
            L->ordering = ordering ;
            lnz_best = Common->lnz ;
            for (k = 0 ; k < n ; k++)
            {
                Lperm [k] = Perm [k] ;
            }
            /* 保存SparseCore_analyze_ordering的结果（如果已调用） */
            skip_best = skip_analysis ;
            if (!skip_analysis)
            {
                /* 保存列数； 成为L的永久部分 */
                for (k = 0 ; k < n ; k++)
                {
                    Lcolcount [k] = ColCount [k] ;
                }
                /* 对于加权后序遍历和超节点分析需要Parent。 不成为L的永久部分 */
                for (k = 0 ; k < n ; k++)
                {
                    Lparent [k] = Parent [k] ;
                }
            }
        }
    }
    #ifdef all_methods_time
    if (for_whom == 1)
            fprintf (fresult,"\n");
    fclose(fresult);
    #endif
    /* 重新打开错误打印 */
    Common->try_catch = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* 如果没有成功的排序方法，则返回 */
    /* ---------------------------------------------------------------------- */

    if (Common->selected == EMPTY)
    {
	/* 所有方法均失败。 
	 * 如果两个或更多个方法失败，则它们可能由于不同的原因而失败
	 * 两者都将清除Common-> status并跳至下一个方法。
	 *  此处需要将Common-> status恢复为使用任何一种方法所获得的最严重的错误。
	 * SPARSE_INVALID比SPARSE_OUT_OF_MEMORY更糟糕
	 * 因为前者暗示用户的输入可能有问题。
	 *  SPARSE_OUT_OF_MEMORY只是缺乏资源的指示。 */
        if (status >= SPARSE_OK)
        {
            /* 如果nmethods = 1，ordering= SPARSE_GIVEN，但UserPerm为NULL 则可能会发生这种情况 */
            status = SPARSE_INVALID ;
        }
        ERROR (status, "all methods failed") ;
        FREE_WORKSPACE_AND_RETURN ;
    }

    /* ---------------------------------------------------------------------- */
    /* 如果跳过，则对AND进行分析 */
    /* ---------------------------------------------------------------------- */


    Common->fl  = Common->method [Common->selected].fl  ;
    Common->lnz = Common->method [Common->selected].lnz ;

    if (skip_best)
    {
        if (!SparseChol_analyze_ordering (A, L->ordering, Lperm, fset, fsize,
            Lparent, Post, Lcolcount, First, Level, Common))  // QR,Chol都走了这里
        {
            /* 内存溢出或者方法失败 */
            FREE_WORKSPACE_AND_RETURN ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* 对消去树进行后序遍历，按列计数加权 */
    /* ---------------------------------------------------------------------- */

    if (Common->postorder)
    {
        /* 结合减少填充的排序和加权的后序遍历 */
        /* workspace: Iwork (2*nrow) */
        if (SparseChol_postorder (Lparent, n, Lcolcount, Post, Common) == n)
        {
            /* 使用First和Level作为工作空间 */
            Int *Wi = First, *InvPost = Level ;
            Int newchild, oldchild, newparent, oldparent ;

            for (k = 0 ; k < n ; k++)
            {
                Wi [k] = Lperm [Post [k]] ;
            }
            for (k = 0 ; k < n ; k++)
            {
                Lperm [k] = Wi [k] ;
            }

            for (k = 0 ; k < n ; k++)
            {
                Wi [k] = Lcolcount [Post [k]] ;
            }
            for (k = 0 ; k < n ; k++)
            {
                Lcolcount [k] = Wi [k] ;
            }
            for (k = 0 ; k < n ; k++)
            {
                InvPost [Post [k]] = k ;
            }

            /* 仅在超结点情况下需要更新的Lparent */
            for (newchild = 0 ; newchild < n ; newchild++)
            {
                oldchild = Post [newchild] ;
                oldparent = Lparent [oldchild] ;
                newparent = (oldparent == EMPTY) ? EMPTY : InvPost [oldparent] ;
                Wi [newchild] = newparent ;
            }
            for (k = 0 ; k < n ; k++)
            {
                Lparent [k] = Wi [k] ;
            }
            /* 使用Iwork作为工作空间完成 */

            /* L已被后序遍历，不再自然排序 */
            if (L->ordering == SPARSE_NATURAL)
            {
                L->ordering = SPARSE_POSTORDERED ;
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* 超结点分析，如果需要或自动选择 */
    /* ---------------------------------------------------------------------- */
    double switch_value = Common->fl / Common->lnz;
    // #ifdef PRINT_TIME
    // printf("||--- switch value = %lf ---|| \n", switch_value);
    // #endif
    if (!for_whom && switch_value >= 1000 ) // 单独对QR分解修改 CHUNK参数
    {
        Common->QR_CHUNK_FLAG = 1;
    }

    
#ifndef NSUPERNODAL
    if (Common->supernodal > SPARSE_AUTO
    || (Common->supernodal == SPARSE_AUTO && // 是AUTO
	Common->lnz > 0 &&  // 非零元数目 > 0
	(switch_value) >= Common->supernodal_switch)) //Common->supernodal_switch
    {
	sparse_csc *S, *F, *A2, *A1 ;
	permute_matrices (A, L->ordering, Lperm, fset, fsize, TRUE,
		&A1, &A2, &S, &F, Common) ;

	/* workspace: Flag (nrow), Head (nrow), Iwork (5*nrow) */
	SparseChol_super_symbolic2 (for_whom, S, F, Lparent, L, Common) ;
	SparseCore_free_sparse (&A1, Common) ;
	SparseCore_free_sparse (&A2, Common) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* 释放临时矩阵和工作空间，并返回结果L */
    /* ---------------------------------------------------------------------- */

    FREE_WORKSPACE_AND_RETURN ; // 这里return L
}

/* HNUCHOL对AMD排序例程的接口。 如果矩阵是对称的，则对A排序。
 *  在输出上，如果A的第i行/列是P * A * P'的第k行/列，则Perm [k] = i。
 */

#if (!defined (AMD_VERSION) || (AMD_VERSION < AMD_VERSION_CODE (2,0)))
#error "AMD v2.0 or later is required"
#endif

/**
 * @brief   
 * 
 */
int SparseChol_amd
(
    /* ---- input ---- */
    sparse_csc *A,	/* 要排序的矩阵 */
    Int *fset,		    /* 0:(A->ncol)-1的子集 */
    size_t fsize,	    /* fset的大小 */
    /* ---- output --- */
    Int *Perm,		    /* size A->nrow, output permutation */
    /* --------------- */
    sparse_common *Common
)
{
    double Info [AMD_INFO], Control2 [AMD_CONTROL], *Control ;
    Int *Cp, *Len, *Nv, *Head, *Elen, *Degree, *Wi, *Iwork, *Next ;
    sparse_csc *C ;
    Int j, n, cnz ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    // RETURN_IF_NULL_COMMON (FALSE) ;
    // RETURN_IF_NULL (A, FALSE) ;
    n = A->nrow ;

    // RETURN_IF_NULL (Perm, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    Common->status = SPARSE_OK ;
    if (n == 0)
    {
	Common->fl = 0 ;
	Common->lnz = 0 ;
	Common->anz = 0 ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 获取工作空间 */
    /* ---------------------------------------------------------------------- */

    /* 注意：此空间小于SparseChol_analyze中使用的空间，因此如果
     * 该例程正在调用SparseCore_amd，将不会分配空间
     */

    /* s = MAX (6*n, A->ncol) */
    s = SparseCore_mult_size_t (n, 6, &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }
    s = MAX (s, A->ncol) ;

    SparseCore_allocate_work (n, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    Iwork  = Common->Iwork ;
    Degree = Iwork ;			        /* size n */
    Wi     = Iwork + n ;		        /* size n */
    Len    = Iwork + 2*((size_t) n) ;	/* size n */
    Nv     = Iwork + 3*((size_t) n) ;   /* size n */
    Next   = Iwork + 4*((size_t) n) ;   /* size n */
    Elen   = Iwork + 5*((size_t) n) ;   /* size n */

    Head = Common->Head ;   /* size n+1, but only n is used */

    /* ---------------------------------------------------------------------- */
    /* 构造AMD的输入矩阵 */
    /* ---------------------------------------------------------------------- */

    if (A->stype == 0) // 非对称矩阵
    {
	/* C = A*A' or A(:,f)*A(:,f)',向C中添加nnz（C）/ 2 + n的额外空间 */
	C = SparseCore_aat (A, fset, fsize, -2, Common) ;
    }
    else // 对称情况
    {
	/* C = A + A'，但如果A-> stype = 1，则仅使用A的上三角部分
	 * 如果A-> stype = -1，则仅使用A的下三部分部。 
     * 向C添加nnz（C）/ 2 + n的额外空间 */
	C = SparseCore_copy (A, 0, -2, Common) ;
    }

    if (Common->status < SPARSE_OK)
    {
	/* 内存不足，fset无效或其他错误 */
	return (FALSE) ;
    }

    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
	Len [j] = Cp [j+1] - Cp [j] ;
    }

    /* C不包括对角线，并且不包括上部和下部。
     * Common-> anz包括对角线和C的下部 */
    cnz = Cp [n] ;
    Common->anz = cnz / 2 + n ;

    /* ---------------------------------------------------------------------- */
    /* order C using AMD */
    /* ---------------------------------------------------------------------- */

    /* 获取参数 */
    if (Common->current < 0 || Common->current >= SPARSE_MAXMETHODS)
    {
	/* 使用AMD默认值 */
	Control = NULL ;
    }
    else
    {
	Control = Control2 ;
	Control [AMD_DENSE] = Common->method [Common->current].prune_dense ;
	Control [AMD_AGGRESSIVE] = Common->method [Common->current].aggressive ;
    }

    amd_2 (n, C->p,  C->i, Len, C->nzmax, cnz, Nv, Next, Perm, Head, Elen,
	    Degree, Wi, Control, Info) ;

    /*  LL的flop计数。 需要为LL'flop计数减去n。 注意这个
     * 是一个很小的上限，通常很精确（有关详细信息，请参阅AMD / 
     * details). SparseChol_analyze计算准确的flop数和填入 */
    Common->fl = Info [AMD_NDIV] + 2 * Info [AMD_NMULTSUBS_LDL] + n ;

    /* Info[AMD_LNZ]不包括对角线 */
    Common->lnz = n + Info [AMD_LNZ] ;

    /* ---------------------------------------------------------------------- */
    /* 释放AMD工作区并清除“通用”中的持久性工作区 */
    /* ---------------------------------------------------------------------- */

    SparseCore_free_sparse (&C, Common) ;
    for (j = 0 ; j <= n ; j++)
    {
	Head [j] = EMPTY ;
    }
    return (TRUE) ;
}

/* HNUCHOL对于COLAMD排序例程（版本2.4或更高版本）。
 * 找到一个置换p，使得PAA'P'的Cholesky因式分解比使用colamd的AA'分解更稀疏。
 * 如果后序遍历的输入参数为TRUE，则找到etree列并对其进行后序遍历，然后将colamd排序与其后序遍历相结合。
 * A必须是不对称的。
 */

#if (!defined (COLAMD_VERSION) || (COLAMD_VERSION < COLAMD_VERSION_CODE (2,5)))
#error "COLAMD v2.5 or later is required"
#endif

int SparseChol_colamd
(
    /* ---- input ---- */
    sparse_csc *A,	/* 待排序矩阵 */               // QR: AT
    Int *fset,		    /* 0:(A->ncol)-1的子集 */        // QR:NULL
    size_t fsize,	    /* fset的大小 */                // QR:0
    int postorder,	    /* 如果为TRUE，请遵循coletree的后序遍历 */   // QR: TRUE
    /* ---- output --- */
    Int *Perm,		    /* size A->nrow, output permutation */    // QR:(Long *) (Q1fill + n1cols)
    /* --------------- */
    sparse_common *Common
)
{
    double knobs [COLAMD_KNOBS] ;
    sparse_csc *C ;
    Int *NewPerm, *Parent, *Post, *Work2n ;
    Int k, nrow, ncol ;
    size_t s, alen ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Perm, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    if (A->stype != 0)
    {
	ERROR (SPARSE_INVALID, "matrix must be unsymmetric") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作区 */
    /* ---------------------------------------------------------------------- */

    /* 注意：此空间小于SparseChol_analyze中使用的空间，因此如果
     * 该例程正在调用SparseCore_colamd，将不会有空间被分配。
     */

    /* s = 4*nrow + ncol */
    s = SparseCore_mult_size_t (nrow, 4, &ok) ;
    s = SparseCore_add_size_t (s, ncol, &ok) ;

    alen = colamd_recommended (A->nzmax, ncol, nrow) ;
    colamd_set_defaults (knobs) ;

    if (!ok || alen == 0)
    {
	ERROR (SPARSE_TOO_LARGE, "matrix invalid or too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (0, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 分配COLAMD工作区 */
    /* ---------------------------------------------------------------------- */

    C = SparseCore_allocate_sparse (ncol, nrow, alen, TRUE, TRUE, 0,
	    SPARSE_PATTERN, Common) ;

    /* ---------------------------------------------------------------------- */
    /* 将输入矩阵A复制（并转置）到Colamd工作区中 */
    /* ---------------------------------------------------------------------- */

    /* C = A (:,f)',如果需要可以打包A */
    /* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset) */
    ok = SparseCore_transpose_unsym (A, 0, NULL, fset, fsize, C, Common) ;

    /* ---------------------------------------------------------------------- */
    /* 对矩阵进行排序（破坏C-> i和C-> p的内容） */
    /* ---------------------------------------------------------------------- */

    /* 获取参数 */
    if (Common->current < 0 || Common->current >= SPARSE_MAXMETHODS)
    {
	/* 这是HNUCHOL默认值，而不是COLAMD默认值 */
	knobs [COLAMD_DENSE_ROW] = -1 ;
    }
    else
    {
	/* 从通用参数获取knobs */
	knobs [COLAMD_DENSE_COL] = Common->method[Common->current].prune_dense ;
	knobs [COLAMD_DENSE_ROW] = Common->method[Common->current].prune_dense2;
	knobs [COLAMD_AGGRESSIVE] = Common->method[Common->current].aggressive ;
    }

    if (ok)
    {
	Int *Cp ;
	Int stats [COLAMD_STATS] ;
	Cp = C->p ;

	colamd (ncol, nrow, alen, C->i, Cp, knobs, stats) ;

	ok = stats [COLAMD_STATUS] ;
	ok = (ok == COLAMD_OK || ok == COLAMD_OK_BUT_JUMBLED) ;
	/* 如果排序成功，则在C->p中返回置换 */
	for (k = 0 ; k < nrow ; k++)
	{
	    Perm [k] = Cp [k] ;
	}
    }

    SparseCore_free_sparse (&C, Common) ;

    /* ---------------------------------------------------------------------- */
    /* 列消去树后序遍历 */
    /* ---------------------------------------------------------------------- */

    if (postorder)
    {
	/* 将Iwork中的最后2 * n空间用于Parent 和 Post */
	Work2n = Common->Iwork ;
	Work2n += 2*((size_t) nrow) + ncol ;
	Parent = Work2n ;		/* size nrow (i/i/l) */
	Post   = Work2n + nrow ;	/* size nrow (i/i/l) */

	/* workspace: Iwork (2*nrow+ncol), Flag (nrow), Head (nrow+1) */
	ok = ok && SparseChol_analyze_ordering (A, SPARSE_COLAMD, Perm, fset,
		fsize, Parent, Post, NULL, NULL, NULL, Common) ;

	/* 将Colamd置换与其后序遍历结合 */
	if (ok)
	{
	    NewPerm = Common->Iwork ;		/* size nrow (i/i/l) */
	    for (k = 0 ; k < nrow ; k++)
	    {
		NewPerm [k] = Perm [Post [k]] ;
	    }
	    for (k = 0 ; k < nrow ; k++)
	    {
		Perm [k] = NewPerm [k] ;
	    }
	}
    }

    return (ok) ;
}


/* 计算A或A'* A的消除树 */

static void update_etree
(
    /* inputs, not modified */
    Int k,		        /* 处理输入图中的边（k，i） */
    Int i,
    /* inputs, modified on output */
    Int Parent [ ],	    /* 如果p是t的父节点Parent [t] = p */
    Int Ancestor [ ]	/* Ancestor[t]是节点t的祖先部分构造的消去树 */
)
{
    Int a ;
    for ( ; ; )		/* 遍历从k到树的根的路径 */
    {
        a = Ancestor [k] ;
        if (a == i)
        {
            /* 到达的最后一个祖先；tree没有改变 */
            return ;
        }
        /* 执行路径压缩 */
        Ancestor [k] = i ;
        if (a == EMPTY)
        {
            /* 未定义的最后一个祖先；这是树的新边 */
            Parent [k] = i ;
            return ;
        }
        /* 遍历到k的祖先 */
        k = a ;
    }
}

/**
 * @brief   找出A或A'* A的消除树
 * 
 */
int SparseChol_etree
(
    /* ---- input ---- */
    sparse_csc *A,
    /* ---- output --- */
    Int *Parent,	/* size ncol.  Parent [j] = p if p is the parent of j */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Ap, *Ai, *Anz, *Ancestor, *Prev, *Iwork ;
    Int i, j, jprev, p, pend, nrow, ncol, packed, stype ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作区 */
    /* ---------------------------------------------------------------------- */

    stype = A->stype ;

    /* s = A->nrow + (stype ? 0 : A->ncol) */
    s = SparseCore_add_size_t (A->nrow, (stype ? 0 : A->ncol), &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (0, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;	/* 内存溢出 */
    }

    Iwork = Common->Iwork ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;	/* A的列数 */
    nrow = A->nrow ;	/* A的行数 */
    Ap = A->p ;		    /* A的列指针 */
    Ai = A->i ;		    /* A的行下标 */
    Anz = A->nz ;	    /* A的每一列中的非零元素 */
    packed = A->packed ;
    Ancestor = Iwork ;	/* size ncol (i/i/l) */

    for (j = 0 ; j < ncol ; ++j)
    {
        Parent [j] = EMPTY ;
        Ancestor [j] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* 计算消去树 */
    /* ---------------------------------------------------------------------- */

    if (stype > 0)
    {
        /* ------------------------------------------------------------------ */
        /* 对称 (上半部分) 情况: 计算消去树 (A) */
        /* ------------------------------------------------------------------ */
        for (j = 0 ; j < ncol ; ++j)
        {
            /* 对于triu（A）的第j列中的每一行i，不包括对角线 */
            p = Ap [j] ;
            pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
            for ( ; p < pend ; p++)
            {
                i = Ai [p] ;
                if (i < j)
                {
                    update_etree (i, j, Parent, Ancestor) ;
                }
            }
        }

    }
    else if (stype == 0)
    {
        /* ------------------------------------------------------------------ */
        /* 非对称情况: 计算消去树 (A'*A) */
        /* ------------------------------------------------------------------ */
        Prev = Iwork + ncol ;	/* size nrow (i/i/l) */
        for (i = 0 ; i < nrow ; i++)
        {
            Prev [i] = EMPTY ;
        }
        for (j = 0 ; j < ncol ; j++)
        {
            /* 对于A列j中的每一行i */
            p = Ap [j] ;
            pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
            for ( ; p < pend ; p++)
            {
            /* 图是动态构造的，A的每行只有一条路径
            * A的第i行包含列下标（j1，j2，j3，j4），
            * 则新图具有边（j1，j2），（j2，j3），和（j3，j4）。 
            *  在此路径图的节点i处，所有边考虑（jprev，j），其中jprev <j */
            i = Ai [p] ;
            jprev = Prev [i] ;
            if (jprev != EMPTY)
            {
                update_etree (jprev, j, Parent, Ancestor) ;
            }
            Prev [i] = j ;
            }
        }
    }
    else
    {
        /* ------------------------------------------------------------------ */
        /* 不支持下三角部分的对称情况 */
        /* ------------------------------------------------------------------ */
        ERROR (SPARSE_INVALID, "symmetric lower not supported") ;
        return (FALSE) ;
    }

    return (TRUE) ;
}

/* 下面的代码包括对树的递归和非递归深度优先搜索。
 * 递归代码更简单，但是会导致栈溢出。
 * 此处留作参考，以了解非递归代码正在计算什么。
 * 要尝试递归版本，请取消注释下面的#define，或使用-DRECURSIVE编译代码。
 *请注意，可能会发生堆栈溢出。
#define RECURSIVE
 */

#ifdef RECURSIVE

/**
 * @brief   递归版本：工作代码仅供参考，不能实际使用
 * 
 */
static Int dfs		/* 返回新的k值 */
(
    Int p,		    /* 在节点p上启动DFS */
    Int k,		    /* 从k开始节点编号 */
    Int Post [ ],	/* 后序遍历，输出修改 */
    Int Head [ ],	/* Head [p] =p的最小的孩子；输出为空 */
    Int Next [ ],	/* Next[j] = j的兄弟；未修改 */
    Int Pstack [ ]	/* 未使用 */
)
{
    Int j ;
    /* 在节点p的每个子节点上启动DFS */
    for (j = Head [p] ; j != EMPTY ; j = Next [j])
    {
	/* 在子节点j上启动DFS */
	k = dfs (j, k, Post, Head, Next, Pstack) ;
    }
    Post [k++] = p ;	/* 将节点p作为第k个节点 */
    Head [p] = EMPTY ;	/* 不再需要链接列表p */
    return (k) ;	    /* 下一个节点将编号为k */
}

#else

/**
 * @brief   实际使用的非递归版本
 * 
 */
static Int dfs		/* 返回新的k值 */
(
    Int p,		    /* 在节点p上启动DFS */
    Int k,		    /* 从k开始节点编号 */
    Int Post [ ],	/* 后序遍历，输出修改 */
    Int Head [ ],	/* Head [p] =p的最小的孩子；输出为空 */
    Int Next [ ],	/* Next[j] = j的兄弟；未修改 */
    Int Pstack [ ]	/* 大小为n的工作空间，在输入或输出上未定义 */
)
{
    Int j, phead ;

    /* 将根节点放在栈上 */
    Pstack [0] = p ;
    phead = 0 ;

    /* 当堆栈不为空时，请执行以下操作 */
    while (phead >= 0)
    {
	/* 从栈顶获取节点p并获得其最小的孩子j */
	p = Pstack [phead] ;
	j = Head [p] ;
	if (j == EMPTY)
	{
	    /* p的所有子元素均已排序。 从栈中删除p并对其进行排序 */
	    phead-- ;
	    Post [k++] = p ;	/* 将节点p作为第k个节点 */
	}
	else
	{
	    /* 将p留在堆栈上。 通过将j放在堆栈上并从p的子节点列表中删除j，在子节点j上启动DFS。 */
	    Head [p] = Next [j] ;
	    Pstack [++phead] = j ;
	}
    }
    return (k) ;	/* 下一个节点将编号为k */
}

#endif

/**
 * @brief   后序遍历一棵树。 该树可以是消除树（来自SparseCore_etree的输出）或
 *          组件树（来自SparseCore_nested_dissection）
 * 
 * @return Sparse_long 
 */
Sparse_long SparseChol_postorder	/* 返回后序遍历节点的# */
(
    /* ---- input ---- */
    Int *Parent,	/* 如果p是j的父母，size n，Parent [j] = p */
    size_t n,
    Int *Weight,	/* size n,Weight[j]是节点j的权重 */
    /* ---- output --- */
    Int *Post,		/* size n.Post[k] = j是后序遍历树中的第k个 */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Head, *Next, *Pstack, *Iwork ;
    Int j, p, k, w, nextj ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (Parent, EMPTY) ;
    RETURN_IF_NULL (Post, EMPTY) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作区 */
    /* ---------------------------------------------------------------------- */

    /* s = 2*n */
    s = SparseCore_mult_size_t (n, 2, &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (EMPTY) ;
    }

    SparseCore_allocate_work (n, s, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    Head  = Common->Head ;	/* size n+1, 初始化全部为空 */
    Iwork = Common->Iwork ;
    Next  = Iwork ;		    /* size n (i/i/l) */
    Pstack = Iwork + n ;	/* size n (i/i/l) */

    /* ---------------------------------------------------------------------- */
    /* 为每个节点构造子链接列表 */
    /* ---------------------------------------------------------------------- */

    if (Weight == NULL)
    {
        /* 以相反的顺序，因此每个列表中的孩子都按升序排列 */
        for (j = n-1 ; j >= 0 ; j--)
        {
            p = Parent [j] ;
            if (p >= 0 && p < ((Int) n))
            {
            /* 将j添加到节点p的孩子列表中 */
            Next [j] = Head [p] ;
            Head [p] = j ;
            }
        }

        /* Head [p] = j 如果j是p中最小（编号最小）的孩子 */
        /* Next [j1] = j2 如果j2是j1的第二大兄弟 */
    }
    else
    {

	/* 首先，根据权重构造一组链接列表。
	 * Whead [w] = j 如果节点j是桶w中的第一个节点。
	 * Next [j1] = j2 如果节点j2在链接列表中跟随j1。
	 */

	Int *Whead = Pstack ;	    /* 使用Pstack作为Whead的工作区 */

	for (w = 0 ; w < ((Int) n) ; w++)
	{
	    Whead [w] = EMPTY ;
	}
	/* 以向前的顺序执行，因此联系的节点按节点索引排序 */
	for (j = 0 ; j < ((Int) n) ; j++)
	{
	    p = Parent [j] ;
	    if (p >= 0 && p < ((Int) n))
	    {
		w = Weight [j] ;
		w = MAX (0, w) ;
		w = MIN (w, ((Int) n) - 1) ;
		/* 将节点j放在权重w的链接列表的开头 */
		Next [j] = Whead [w] ;
		Whead [w] = j ;
	    }
	}

	/* 遍历权重桶，将每个节点放在其父级列表中 */
	for (w = n-1 ; w >= 0 ; w--)
	{
	    for (j = Whead [w] ; j != EMPTY ; j = nextj)
	    {
		nextj = Next [j] ;
		/* 将节点j放在其父节点的链接列表中 */
		p = Parent [j] ;
		Next [j] = Head [p] ;
		Head [p] = j ;
	    }
	}

	/* 不再需要Whead ] */
	/* Head [p] = j 如果j是p的最轻孩子*/
	/* Next [j1] = j2 如果j2是j1的下一个最重兄弟 */
    }

    /* ---------------------------------------------------------------------- */
    /* 在消去树的每个根节点上启动DFS */
    /* ---------------------------------------------------------------------- */

    k = 0 ;
    for (j = 0 ; j < ((Int) n) ; j++)
    {
        if (Parent [j] == EMPTY)
        {
            /* j是树的根； 在这里启动DFS */
            k = dfs (j, k, Post, Head, Next, Pstack) ;
        }
    }

    /* 这通常已经是EMPTY，除非Parent无效 */
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	    Head [j] = EMPTY ;
    }

    return (k) ;
}


/**
 * @brief   后序遍历消去树中第k个节点的初始化工作
 * 
 */
static Int initialize_node  
(
    Int k,		        /* 在算法的第k步（和第k个节点） */
    Int Post [ ],	    /* 后序遍历消去树中的第k个节点 */
    Int Parent [ ],	    /* Parent [i]是i在消去树中的父母 */
    Int ColCount [ ],	/* ColCount [c]是节点c的当前权重 */
    Int PrevNbr [ ]	    /* 如果u最后一次在步骤k被考虑 PrevNbr [u] = k */
)
{
    Int p, parent ;
    /* 确定p，即后序遍历消去树中的第k个节点 */
    p = Post [k] ;
    /* 如果p不是消去树的根，则调整权重 */
    parent = Parent [p] ;
    if (parent != EMPTY)
    {
	ColCount [parent]-- ;
    }
    /* 标记节点p以排除自身边（p，p）*/
    PrevNbr [p] = k ;
    return (p) ;
}

/**
 * @brief   边（p，u）正在处理。 p <u是其祖先u在消去树中的后代。
 *          节点p是后序遍历消去树中的第k个节点。
 * 
 */
static void process_edge
(
    Int p,		        /* 处理矩阵的边（p，u） */
    Int u,
    Int k,		        /* 我们在后序遍历消去树中的第k个节点 */
    Int First [ ],	    /* First [i] = k 如果节点i的第一个后代的后序遍历为k */
    Int PrevNbr [ ],	/* u最后一次在步骤k被考虑 k = PrevNbr [u] */
    Int ColCount [ ],	/* ColCount [c] 是节点c的当前权重 */
    Int PrevLeaf [ ],	/* s = PrevLeaf [u] 表示s是在以u为根的子树中看到的最后一片叶子。  */
    Int RowCount [ ],	/* RowCount [i] L的第i行（包括对角线）中的非零数。如果为NULL，则不计算。 */
    Int SetParent [ ],	/* FIND / UNION数据结构，该结构形成一组树。
			             * 根i为i = SetParent [i]。遵循从i到包含i的子树的根q的路径，意味着q是i的SetParent代表。
			             * 意味着q是i的SetParent代表。
			             * 树中的所有节点的SetParent都可以等于根q；
			             * 树表示法可以节省时间。
			             * 当从i到其根q跟踪路径时，将重新遍历该路径以将整个路径的SetParent设置为根q。
			             */
    Int Level [ ]	 /* Level [i] = 从节点i到根的路径长度 */
)
{
    Int prevleaf, q, s, sparent ;
    if (First [p] > PrevNbr [u])
    {
	/* p是u子树的叶子 */
	ColCount [p]++ ;
	prevleaf = PrevLeaf [u] ;
	if (prevleaf == EMPTY)
	{
	    /* p是u子树的第一片叶子；
	     * RowCount将以消去树中路径的长度从p递增到u。 */
	    q = u ;
	}
	else
	{
	    /* q = FIND (prevleaf): 找到包含prevleaf的SetParent树的根q */
	    for (q = prevleaf ; q != SetParent [q] ; q = SetParent [q])
	    {
		;
	    }
	    /* 根q已找到；重新遍历路径并执行路径压缩 */
	    s = prevleaf ;
	    for (s = prevleaf ; s != q ; s = sparent)
	    {
		sparent = SetParent [s] ;
		SetParent [s] = q ;
	    }
	    /* 调整Rowcount和ColCount；RowCount将增加从p到SetParent根q的路径长度，并将q的ColCount减1。*/
	    ColCount [q]-- ;
	}
	if (RowCount != NULL)
	{
	    /* 如果正在计算RowCount，则将其增加从p到q的路径长度 */
	    RowCount [u] += (Level [p] - Level [q]) ;
	}
	/* p是u的子树的叶子，因此将PrevLeaf [u]标记为p */
	PrevLeaf [u] = p ;
    }
    /* 标志u已在步骤k处理 */
    PrevNbr [u] = k ;
}

/**
 * @brief   计算UNION（p，Parent[p]）
 * 
 */
static void finalize_node    
(
    Int p,
    Int Parent [ ],	    /* Parent [p] 是消去树中p的父结点 */
    Int SetParent [ ]	/* 参见上面的process_edge */
)
{
    /* 现在，以p为根的SetParent树中的所有节点都
     * 将节点Parent [p]作为其最终根。这将计算UNION（p，Parent[p]）*/
    if (Parent [p] != EMPTY)
    {
	SetParent [p] = Parent [p] ;
    }
}

/**
 * @brief   计算矩阵A或A*A'的Cholesky因子L的行数和列数。
 *          必须已经计算了消去树及其存储(请参阅SparseCore_etree和
 *          SparseCore_postorder)，并将其作为这个例程的输入。
 * 
 */
int SparseChol_rowcolcounts
(
    /* ---- input ---- */
    sparse_csc *A,	/* 分析的矩阵 */
    Int *fset,		    /* 0:(A->ncol)-1的子集 */
    size_t fsize,	    /* fset的大小 */
    Int *Parent,	    /* size nrow.  Parent [i] = p 如果p是i的父结点 */
    Int *Post,		    /* size nrow.  Post [k] = i 如果i是后序遍历etree中的第k个节点 */
    /* ---- output --- */
    Int *RowCount,	    /* size nrow. RowCount [i] = L的第i行中的项，包括对角线 */
    Int *ColCount,	    /* size nrow. ColCount [i] = L的第i列中的项，包括对角线 */
    Int *First,		    /* size nrow.  First [i] = k 是i的任何后代中最少的后序遍历 */
    Int *Level,		    /* size nrow.  Level [i] 是从i到根的路径的长度，Level[root] = 0 */
    /* --------------- */
    sparse_common *Common
)
{
    double fl, ff ;
    Int *Ap, *Ai, *Anz, *PrevNbr, *SetParent, *Head, *PrevLeaf, *Anext, *Ipost,
	*Iwork ;
    Int i, j, r, k, len, s, p, pend, inew, stype, nf, anz, inode, parent,
	nrow, ncol, packed, use_fset, jj ;
    size_t w ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    // RETURN_IF_NULL_COMMON (FALSE) ;
    // RETURN_IF_NULL (A, FALSE) ;
    // RETURN_IF_NULL (Parent, FALSE) ;
    // RETURN_IF_NULL (Post, FALSE) ;
    // RETURN_IF_NULL (ColCount, FALSE) ;
    // RETURN_IF_NULL (First, FALSE) ;
    // RETURN_IF_NULL (Level, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    stype = A->stype ;
    if (stype > 0)
    {
	/* 不支持对称上三角部分 */
	ERROR (SPARSE_INVALID, "symmetric upper not supported") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作区 */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;	/* A的行数 */
    ncol = A->ncol ;	/* A的列数 */

    /* w = 2*nrow + (stype ? 0 : ncol) */
    w = SparseCore_mult_size_t (nrow, 2, &ok) ;
    w = SparseCore_add_size_t (w, (stype ? 0 : ncol), &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (nrow, w, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;	/* A的列指针 大小为ncol+1 */
    Ai = A->i ;	/* A的行索引 大小为nz=Ap[ncol+1] */
    Anz = A->nz ;
    packed = A->packed ;

    /* ---------------------------------------------------------------------- */
    /* 获取工作空间 */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    SetParent = Iwork ;		                    /* size nrow (i/i/l) */
    PrevNbr   = Iwork + nrow ;	                /* size nrow (i/i/l) */
    Anext     = Iwork + 2*((size_t) nrow) ;     /* size ncol (i/i/l) (unsym only) */
    PrevLeaf  = Common->Flag ;	                /* size nrow */
    Head      = Common->Head ;	                /* size nrow+1 (unsym only)*/

    /* ---------------------------------------------------------------------- */
    /* 查找树中每个节点的第一个子孙和层次 */
    /* ---------------------------------------------------------------------- */

    /* First [i] = k 如果节点i的第一个后代的后序为k */
    /* Level [i] = 从节点i到根的路径长度 (Level [root] = 0) */

    for (i = 0 ; i < nrow ; i++)
    {
	First [i] = EMPTY ;
    }

    /* 消去树的后序遍历 */
    for (k = 0 ; k < nrow ; k++)
    {
	/* 消去树的节点i是后序遍历消去树中的第k个节点 */
	i = Post [k] ;

	/* 如果First [i]仍然为空，则i是一片叶子 */
	/* 如果i是叶子，则ColCount [i]从1开始，否则为0 */
	ColCount [i] = (First [i] == EMPTY) ? 1 : 0 ;

	/* 遍历从节点i到根的路径，如果找到已经定义了First [r]的节点r，则停止 */
	len = 0 ;
	for (r = i ; (r != EMPTY) && (First [r] == EMPTY) ; r = Parent [r])
	{
	    First [r] = k ;
	    len++ ;
	}
	if (r == EMPTY)
	{
	    /* 我们命中了一个根节点，其层次为零 */
	    len-- ;
	}
	else
	{
	    /* 我们在已定义Level[r]的节点r处停止 */
	    len += Level [r] ;
	}
	/* 重新遍历从节点i到r的路径；设置每个节点的层次 */
	for (s = i ; s != r ; s = Parent [s])
	{
	    Level [s] = len-- ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* AA'情况：根据第一个后序遍历行下标对A的列进行排序 */
    /* ---------------------------------------------------------------------- */

    fl = 0.0 ;
    if (stype == 0)
    {
	/* [ 使用PrevNbr [0..nrow-1]作为Ipost的工作区 */
	Ipost = PrevNbr ;
	/* Ipost [i] = k 如果i是后序遍历消去树中的第k个节点 */
	for (k = 0 ; k < nrow ; k++)
	{
	    Ipost [Post [k]] = k ;
	}
	use_fset = (fset != NULL) ;
	if (use_fset)
	{
	    nf = fsize ;
	    /* 清除Anext以检查fset */
	    for (j = 0 ; j < ncol ; j++)
	    {
		Anext [j] = -2 ;
	    }
	    /* 在A（post，f）的每一列中找到第一个后序遍历行，并将该列放在相应的链接列表中 */
	    for (jj = 0 ; jj < nf ; jj++)
	    {
		j = fset [jj] ;
		if (j < 0 || j > ncol || Anext [j] != -2)
		{
		    /* fset中超出范围或重复的条目 */
		    ERROR (SPARSE_INVALID, "fset invalid") ;
		    return (FALSE) ;
		}
		/* 将第j列标记为已看到 */
		Anext [j] = EMPTY ;
	    }
	}
	else
	{
	    nf = ncol ;
	}
	for (jj = 0 ; jj < nf ; jj++)
	{
	    j = (use_fset) ? (fset [jj]) : jj ;
	    /* j列在fset中；找到最小的行（如果有 */
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    ff = (double) MAX (0, pend - p) ;
	    fl += ff*ff + ff ;
	    if (pend > p)
	    {
		k = Ipost [Ai [p]] ;
		for ( ; p < pend ; p++)
		{
		    inew = Ipost [Ai [p]] ;
		    k = MIN (k, inew) ;
		}
		/* 将j列放在链接列表k中 */
		Anext [j] = Head [k] ;
		Head [k] = j ;
	    }
	}
	/* 反向后序遍历不再需要Ipost Head [k]包含所有列的链接列表，
     * 其中第一个后序遍历行下标等于k，对于k = 0到nrow-1 */
    }

    /* ---------------------------------------------------------------------- */
    /* 计算行数和节点权重 */
    /* ---------------------------------------------------------------------- */

    if (RowCount != NULL)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    RowCount [i] = 1 ;
	}
    }
    for (i = 0 ; i < nrow ; i++)
    {
	PrevLeaf [i] = EMPTY ;
	PrevNbr [i] = EMPTY ;
	SetParent [i] = i ;	/* 每个节点本身都在自己的集合中 */
    }

    if (stype != 0)
    {

	/* ------------------------------------------------------------------ */
	/* 对称情况: LL' = A */
	/* ------------------------------------------------------------------ */

	/* 还要确定triu（A）中的条目数 */
	anz = nrow ;
	for (k = 0 ; k < nrow ; k++)
	{
	    /* j是后序遍历消去树中的第k个节点 */
	    j = initialize_node (k, Post, Parent, ColCount, PrevNbr) ;

	    /* 对于A的第j列中对角线以下的所有非零A（i，j） */
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i > j)
		{
		    /* j是消去树（A）中i的后代 */
		    anz++ ;
		    process_edge (j, i, k, First, PrevNbr, ColCount,
			    PrevLeaf, RowCount, SetParent, Level) ;
		}
	    }
	    /* 更新 SetParent: UNION (j, Parent [j]) */
	    finalize_node (j, Parent, SetParent) ;
	}
	Common->anz = anz ;
    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* 非对称情况: LL' = AA' */
	/* ------------------------------------------------------------------ */

	for (k = 0 ; k < nrow ; k++)
	{
	    /* inode是后序遍历消去树中的第k个节点 */
	    inode = initialize_node (k, Post, Parent, ColCount, PrevNbr) ;

	    /* 对于第一个后序遍历行为k的所有列j： */
	    for (j = Head [k] ; j != EMPTY ; j = Anext [j])
	    {
		/* k是A列j中的第一个后序遍历行 */
		/* 对于列j中的所有行i: */
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    /* i已经在此步骤k中被考虑过 */
		    if (PrevNbr [i] < k)
		    {
			/* inode是消去树（AA'）中i的后代 */
			/* 处理边（inode，i）并将PrevNbr [i]设置为k */
			process_edge (inode, i, k, First, PrevNbr, ColCount,
				PrevLeaf, RowCount, SetParent, Level) ;
		    }
		}
	    }
	    /* 清除链接列表k */
	    Head [k] = EMPTY ;
	    /* 更新 SetParent: UNION (inode, Parent [inode]) */
	    finalize_node (inode, Parent, SetParent) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 完成计算列数 */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j < nrow ; j++)
    {
	parent = Parent [j] ;
	if (parent != EMPTY)
	{
	    /* 将j的ColCount添加到其父结点 */
	    ColCount [parent] += ColCount [j] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 清除工作空间 */
    /* ---------------------------------------------------------------------- */

    Common->mark = EMPTY ;
    /* CORE(clear_flag) (Common) ; */
    SPARSE_CLEAR_FLAG (Common) ;

    /* ---------------------------------------------------------------------- */
    /* 浮点计数和nnz（L）用于随后的LL'数值分解 */
    /* ---------------------------------------------------------------------- */

    /* 使用double以避免整数溢出。lnz不能为NaN */
    Common->aatfl = fl ;
    Common->lnz = 0. ;
    fl = 0 ;
    for (j = 0 ; j < nrow ; j++)
    {
	ff = (double) (ColCount [j]) ;
	Common->lnz += ff ;
	fl += ff*ff ;
    }

    Common->fl = fl ;

    return (TRUE) ;
}


#define SUBTREE \
    for ( ; p < pend ; p++) \
    { \
	i = Ai [p] ; \
	if (i <= k) \
	{ \
	    /* 将A或AN * A'的列分散到Wx和Wz */ \
	    SCATTER ; \
	    /* 从节点i开始并且遍历子树，在节点k处停止 */ \
	    for (len = 0 ; i < k && i != EMPTY && Flag [i] < mark ; i = parent) \
	    { \
		/* L（k，i）不为零，并且是第一次看到 */ \
		Stack [len++] = i ;	    /* 将i放在栈上 */ \
		Flag [i] = mark ;	    /* 将i标记为已访问 */ \
		parent = PARENT (i) ;   /* 遍历消去树到父结点 */ \
	    } \
	    /* 将路径向下移动到栈的底部 */ \
	    while (len > 0) \
	    { \
		Stack [--top] = Stack [--len] ; \
	    } \
	} \
	else if (sorted) \
	{ \
	    break ; \
	} \
    }


/* ========================================================================== */
/* === 模板 ============================================================= */
/* ========================================================================== */

#define REAL
#include "SparseChol_t_rowfac.c"

#define MASK
#define REAL
#include "SparseChol_t_rowfac.c"
#undef MASK

/* Compute the nonzero pattern of the solution to the lower triangular system
 * L(0:k-1,0:k-1) * x = A (0:k-1,k) if A is symmetric, or
 * L(0:k-1,0:k-1) * x = A (0:k-1,:) * A (:,k)' if A is unsymmetric.
 * This gives the nonzero pattern of row k of L (excluding the diagonal).
 * The pattern is returned postordered.
 */

/**
 * @brief   如果A是对称的，则计算下三角系统的解的
 *          非零模式L（0：k-1,0：k-1）* x = A（0：k-1，k）
 *          或L（0：k -1,0：k-1）* x = A（0：k-1，:) * A（：，k）'（如果A不对称）。
 *          这给出了L的k行（不包括对角线）的非零模式。该模式后序遍历后返回。
 * 
 */
int SparseChol_row_subtree
(
    /* ---- input ---- */
    sparse_csc *A,	/* 分析的矩阵 */
    sparse_csc *F,	/* 仅用于A * A'情况。 F = A'或AN（：，f）' */
    size_t krow,	    /* L的k行 */
    Int *Parent,	    /* 消去树 */
    /* ---- output --- */
    sparse_csc *R,	/* pattern of L(k,:), 1-by-n with R->nzmax >= n */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Rp, *Stack, *Flag, *Ap, *Ai, *Anz, *Fp, *Fi, *Fnz ;
    Int p, pend, parent, t, stype, nrow, k, pf, pfend, Fpacked, packed,
	sorted, top, len, i, mark ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    // RETURN_IF_NULL_COMMON (FALSE) ;
    // RETURN_IF_NULL (A, FALSE) ;
    // RETURN_IF_NULL (R, FALSE) ;
    // RETURN_IF_NULL (Parent, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (R, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    stype = A->stype ;
    if (stype == 0)
    {
	RETURN_IF_NULL (F, FALSE) ;
	RETURN_IF_XTYPE_INVALID (F, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    }
    if (krow >= A->nrow)
    {
	ERROR (SPARSE_INVALID, "subtree: k invalid") ;
	return (FALSE) ;
    }
    if (R->ncol != 1 || A->nrow != R->nrow || A->nrow > R->nzmax)
    {
	ERROR (SPARSE_INVALID, "subtree: R invalid") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作区 */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    SparseCore_allocate_work (nrow, 0, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    if (stype > 0)
    {
	/* 对称上半部情况：不需要F。可能为NULL */
	Fp = NULL ;
	Fi = NULL ;
	Fnz = NULL ;
	Fpacked = TRUE ;
    }
    else if (stype == 0)
    {
	/* 非对称情况：需要F */
	Fp = F->p ;
	Fi = F->i ;
	Fnz = F->nz ;
	Fpacked = F->packed ;
    }
    else
    {
	/* 不支持对称的下三角形式 */
	ERROR (SPARSE_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;
    sorted = A->sorted ;

    k = krow ;
    Stack = R->i ;

    /* ---------------------------------------------------------------------- */
    /* 获取工作空间 */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;	/* size nrow, 必须保持Flag [i] < mark */
    /* mark = CORE(clear_flag) (Common) ; */
    SPARSE_CLEAR_FLAG (Common) ;
    mark = Common->mark ;

    /* ---------------------------------------------------------------------- */
    /* 计算L（k，:)的模式 */
    /* ---------------------------------------------------------------------- */

    top = nrow ;		        /* 栈为空 */
    Flag [k] = mark ;		    /* 在栈中不包括对角线项 */

#define SCATTER			        /* 不分发数值 */
#define PARENT(i) Parent [i]	/* 将父结点用于消去树 */

    if (stype != 0)
    {
	/* 分发triu（A）的第k个col，得到模式L（k，:) */
	p = Ap [k] ;
	pend = (packed) ? (Ap [k+1]) : (p + Anz [k]) ;
	SUBTREE ;
    }
    else
    {
	/* 分发triu（beta*I+AA'）的第k个col，得到模式L（k，:) */
	pf = Fp [k] ;
	pfend = (Fpacked) ? (Fp [k+1]) : (pf + Fnz [k]) ;
	for ( ; pf < pfend ; pf++)
	{
	    /* 获取非零项F（t，k） */
	    t = Fi [pf] ;
	    p = Ap [t] ;
	    pend = (packed) ? (Ap [t+1]) : (p + Anz [t]) ;
	    SUBTREE ;
	}
    }

#undef SCATTER
#undef PARENT

    /* 将栈向上移动到R的第一部分 */
    len = nrow - top ;
    for (i = 0 ; i < len ; i++)
    {
	Stack [i] = Stack [top + i] ;
    }

    Rp = R->p ;
    Rp [0] = 0 ;
    Rp [1] = len ;
    R->sorted = FALSE ;

    SparseCore_clear_flag (Common) ;
    return (TRUE) ;
}

/* Compute the nonzero pattern of Y=L\B.  L must be simplicial, and B
 * must be a single sparse column vector with B->stype = 0.  The values of
 * B are not used; it just specifies a nonzero pattern.  The pattern of
 * Y is not sorted, but is in topological order instead (suitable for a
 * sparse forward/backsolve).
 */

/**
 * @brief   计算Y=L\B的非零pattern。L必须是单纯的，B必须是一个单一的稀疏列向量，
 *          并且B->stype = 0。没有使用B的值;它只是指定了一个非零pattern。
 *          Y的pattern没有排序，而是按照拓扑顺序(适用于稀疏的前/后求解)。
 * 
 */
int SparseChol_lsolve_pattern
(
    /* ---- input ---- */
    sparse_csc *B,	    /* 稀疏right-hand-site(单个稀疏列) */
    sparse_factor *L,	    /* parent(i)的因子L */
    /* ---- output --- */
    sparse_csc *Yset,   /* Y=L\B的pattern, n*1且Y->nzmax >= n */
    /* --------------- */
    sparse_common *Common
)
{
    size_t krow ;
    RETURN_IF_NULL (B, FALSE) ;
    krow = B->nrow ;
    return (SparseChol_row_lsubtree (B, NULL, 0, krow, L, Yset, Common)) ;
}

/**
 * @brief   与SparseCore_row_subtree相同，除了消除树是从L本身获得的，
 *          它是每列中的第一个非对角线条目。
 *          L必须是简单的，而不是超节点的。
 * 
 */
int SparseChol_row_lsubtree
(
    /* ---- input ---- */
    sparse_csc *A,	    /* 分析的矩阵 */
    Int *Fi, size_t fnz,    /* A'的第k行的非零模式，对于对称情况不需要。无需排序 */
    size_t krow,	        /* row k of L */
    sparse_factor *L,	    /* parent(i)的因子L  */
    /* ---- output --- */
    sparse_csc *R,	    /* L(k,:)的pattern, n*1且R->nzmax >= n */
    /* --------------- */
    sparse_common *Common
)
{
    Int *Rp, *Stack, *Flag, *Ap, *Ai, *Anz, *Lp, *Li, *Lnz ;
    Int p, pend, parent, t, stype, nrow, k, pf, packed, sorted, top, len, i,
	mark, ka ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    // RETURN_IF_NULL_COMMON (FALSE) ;
    // RETURN_IF_NULL (A, FALSE) ;
    // RETURN_IF_NULL (R, FALSE) ;
    // RETURN_IF_NULL (L, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (R, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (L, SPARSE_REAL, SPARSE_REAL, FALSE) ;

    nrow = A->nrow ;
    stype = A->stype ;
    if (stype < 0)
    {
	/* 不支持对称的下三角形式 */
	ERROR (SPARSE_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }

    if (krow > nrow)
    {
        ERROR (SPARSE_INVALID, "lsubtree: krow invalid") ;
        return (FALSE) ;
    }
    else if (krow == nrow)
    {
        /* 找到x=L\b的pattern，其中b=A(:,0) */
        k = nrow ;      /* 计算所有结果；不要停在SUBTREE */
        ka = 0 ;        /* use column A(:,0) */
        if (stype != 0 || A->ncol != 1)
        {
            /* A必须是非对称的（它是单个稀疏列向量） */
            ERROR (SPARSE_INVALID, "lsubtree: A invalid") ;
            return (FALSE) ;
        }
    }
    else
    {
        /* 如果A不对称，则使用A（：，k）和Fi查找L（k，:)的pattern */
        k = krow ;      /* 要计算L的哪一行 */
        ka = k ;        /* 使用A的哪一列 */
        if (stype == 0)
        {
            RETURN_IF_NULL (Fi, FALSE) ;
        }
    }

    if (R->ncol != 1 || nrow != R->nrow || nrow > R->nzmax ||
        ((krow == nrow || stype != 0) && ka >= A->ncol))
    {
	ERROR (SPARSE_INVALID, "lsubtree: R invalid") ;
	return (FALSE) ;
    }
    if (L->is_super)
    {
	ERROR (SPARSE_INVALID, "lsubtree: L invalid (cannot be supernodal)") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作区 */
    /* ---------------------------------------------------------------------- */

    SparseCore_allocate_work (nrow, 0, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;
    sorted = A->sorted ;

    Stack = R->i ;

    Lp = L->p ;
    Li = L->i ;
    Lnz = L->nz ;

    /* ---------------------------------------------------------------------- */
    /* 获取工作空间 */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;	/* size nrow, 必须保持Flag [i] < mark */
    mark = SparseCore_clear_flag (Common) ;

    /* ---------------------------------------------------------------------- */
    /* 计算L（k，:)的pattern */
    /* ---------------------------------------------------------------------- */

    top = nrow ;		/* 栈空 */
    if (k < nrow)
    {
        Flag [k] = mark ;       /* 在栈中不包括对角线项 */
    }

#define SCATTER			/* 不要分发数值 */
#define PARENT(i) (Lnz [i] > 1) ? (Li [Lp [i] + 1]) : EMPTY

    if (krow == nrow || stype != 0)
    {
	/* 分发triu（A）的第k个col，得到模式L（k，:) */
	p = Ap [ka] ;
	pend = (packed) ? (Ap [ka+1]) : (p + Anz [ka]) ;
	SUBTREE ;
    }
    else
    {
	/* 分发triu（beta*I+AA'）的第k个col，得到模式L（k，:) */
	for (pf = 0 ; pf < (Int) fnz ; pf++)
	{
	    /* 获取非零项F（t，k） */
	    t = Fi [pf] ;
	    p = Ap [t] ;
	    pend = (packed) ? (Ap [t+1]) : (p + Anz [t]) ;
	    SUBTREE ;
	}
    }

#undef SCATTER
#undef PARENT

    /* 将栈向上移动到R的第一部分 */
    len = nrow - top ;
    for (i = 0 ; i < len ; i++)
    {
	Stack [i] = Stack [top + i] ;
    }

    Rp = R->p ;
    Rp [0] = 0 ;
    Rp [1] = len ;
    R->sorted = FALSE ;

    SparseCore_clear_flag (Common) ;
    return (TRUE) ;
}

/**
 * @brief   这是通用用途的增量分解。
 * 
 */
int SparseChol_rowfac
(
    /* ---- input ---- */
    sparse_csc *A,	/* 要分解的矩阵 */
    sparse_csc *F,	/* 仅用于A * A'情况。F = A'或AN（：，f）' */
    double beta [2],	/* 分解beta*I+A或者beta*I+AA' */
    size_t kstart,	    /* 分解的第一行 */
    size_t kend,	    /* 分解的最后一行为kend-1 */
    /* ---- in/out --- */
    sparse_factor *L,
    /* --------------- */
    sparse_common *Common
)
{
    return (SparseChol_rowfac_mask2 (A, F, beta, kstart, kend, NULL, 0, NULL, L,
	Common)) ;
}

/**
 * @brief   这仅适用于LPDASA
 * 
 */
int SparseChol_rowfac_mask
(
    /* ---- input ---- */
    sparse_csc *A,	/* 要分解的矩阵 */
    sparse_csc *F,	/* 仅用于A * A'情况。 F=A' or A(:,f)' */
    double beta [2],	/* 分解beta*I+A或者beta*I+AA' */
    size_t kstart,	    /* 分解的第一行 */
    size_t kend,	    /* 分解的最后一行kend-1 */
    Int *mask,		    /* size A->nrow. 若mask[i] >= 0，第i行设置为零 */
    Int *RLinkUp,	    /* size A->nrow. 链接要计算的行列表 */
    /* ---- in/out --- */
    sparse_factor *L,
    /* --------------- */
    sparse_common *Common
)
{
    Int maskmark = 0 ;
    return (SparseChol_rowfac_mask2 (A, F, beta, kstart, kend, mask, maskmark,
        RLinkUp, L, Common)) ;
}

/**
 * @brief 这仅适用于LPDASA
 * 
 */
int SparseChol_rowfac_mask2
(
    /* ---- input ---- */
    sparse_csc *A,	/* 要分解的矩阵 */
    sparse_csc *F,	/* 仅用于A * A'情况。 F=A' or A(:,f)' */
    double beta [2],	/* factorize beta*I+A or beta*I+AA' */
    size_t kstart,	    /* 分解的第一行 */
    size_t kend,	    /* 分解的最后一行kend-1 */
    Int *mask,		    /* size A->nrow. 若mask[i] >= maskmark 第i行设置为零 */
    Int maskmark,       /* 用于mask[i]测试 */
    Int *RLinkUp,	    /* size A->nrow. 链接要计算的行 列表 */
    /* ---- in/out --- */
    sparse_factor *L,
    /* --------------- */
    sparse_common *Common
)
{
    Int n ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    // RETURN_IF_NULL_COMMON (FALSE) ;
    // RETURN_IF_NULL (A, FALSE) ;
    // RETURN_IF_NULL (L, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (A, SPARSE_REAL, SPARSE_REAL, FALSE) ;
    // RETURN_IF_XTYPE_INVALID (L, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    if (L->xtype != SPARSE_PATTERN && A->xtype != L->xtype)
    {
	ERROR (SPARSE_INVALID, "xtype of A and L do not match") ;
	return (FALSE) ;
    }
    if (L->is_super)
    {
	ERROR (SPARSE_INVALID, "can only do simplicial factorization");
	return (FALSE) ;
    }
    if (A->stype == 0)
    {
	RETURN_IF_NULL (F, FALSE) ;
	if (A->xtype != F->xtype)
	{
	    ERROR (SPARSE_INVALID, "xtype of A and F do not match") ;
	    return (FALSE) ;
	}
    }
    if (A->stype < 0)
    {
	/* 不支持对称的下三角形式 */
	ERROR (SPARSE_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }
    if (kend > L->n)
    {
	ERROR (SPARSE_INVALID, "kend invalid") ;
	return (FALSE) ;
    }
    if (A->nrow != L->n)
    {
	ERROR (SPARSE_INVALID, "dimensions of A and L do not match") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;
    Common->rowfacfl = 0 ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作区 */
    /* ---------------------------------------------------------------------- */

    /* Xwork的实际大小为n */
    n = L->n  ;

    /* s = ((A->xtype != SPARSE_REAL) ? 2:1)*n */
    s = SparseCore_mult_size_t (n, ((A->xtype != SPARSE_REAL) ? 2:1), &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (n, n, s, Common) ;
    if (Common->status < SPARSE_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 使用模板例程分解矩阵 */
    /* ---------------------------------------------------------------------- */

    if (RLinkUp == NULL)
    {

	switch (A->xtype)
	{
	    case SPARSE_REAL:
		ok = r_SparseChol_rowfac (A, F, beta, kstart, kend, L, Common) ;
		break ;
	}

    }
    else
    {

	switch (A->xtype)
	{
	    case SPARSE_REAL:
		ok = r_SparseChol_rowfac_mask (A, F, beta, kstart, kend,
		    mask, maskmark, RLinkUp, L, Common) ;
		break ;
	}
    }

    return (ok) ;
}

/* ========================================================================== */
/* === LMINMAX ============================================================== */
/* ========================================================================== */

/* 更新一个条目L（j，j）的lmin和lmax */

#define FIRST_LMINMAX(Ljj,lmin,lmax) \
{ \
    double ljj = Ljj ; \
    if (IS_NAN (ljj)) \
    { \
	return (0) ; \
    } \
    lmin = ljj ; \
    lmax = ljj ; \
}

#define LMINMAX(Ljj,lmin,lmax) \
{ \
    double ljj = Ljj ; \
    if (IS_NAN (ljj)) \
    { \
	return (0) ; \
    } \
    if (ljj < lmin) \
    { \
	lmin = ljj ; \
    } \
    else if (ljj > lmax) \
    { \
	lmax = ljj ; \
    } \
}

/**
 * @brief   
 * 
 * @return double 返回 min(diag(L)) / max(diag(L))
 */
double SparseCore_rcond 
(
    /* ---- input ---- */
    sparse_factor *L,
    /* --------------- */
    sparse_common *Common
)
{
    double lmin, lmax, rcond ;
    double *Lx ;
    Int *Lpi, *Lpx, *Super, *Lp ;
    Int n, e, nsuper, s, k1, k2, psi, psend, psx, nsrow, nscol, jj, j ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    // RETURN_IF_NULL_COMMON (EMPTY) ;
    // RETURN_IF_NULL (L, EMPTY) ;
    // RETURN_IF_XTYPE_INVALID (L, SPARSE_REAL, SPARSE_REAL, EMPTY) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    n = L->n ;
    if (n == 0)
    {
	return (1) ;
    }
    if (L->minor < L->n)
    {
	return (0) ;
    }

    e = 1 ;

    if (L->is_super)
    {
	/* L是超节点 */
	nsuper = L->nsuper ;	/* L中的超节点数 */
	Lpi = L->pi ;		    /* 整型pattern的列指针 */
	Lpx = L->px ;		    /* 数值的列指针 */
	Super = L->super ;	    /* 超节点大小 */
	Lx = L->x ;		        /* 数值 */
	FIRST_LMINMAX (Lx [0], lmin, lmax) ;	/* L的第一个对角线入口 */
	for (s = 0 ; s < nsuper ; s++)
	{
	    k1 = Super [s] ;		/* 超节点s的第一列 */
	    k2 = Super [s+1] ;		/* 超节点的最后一列是k2-1 */
	    psi = Lpi [s] ;		    /* 第一行下标为L-> s [psi] */
	    psend = Lpi [s+1] ;		/* 最后一行下标是L-> s [psend-1] */
	    psx = Lpx [s] ;		    /* 第一个数字输入为Lx [psx] */
	    nsrow = psend - psi ;	/* 超级节点是nsrow-by-nscol */
	    nscol = k2 - k1 ;
	    for (jj = 0 ; jj < nscol ; jj++)
	    {
		LMINMAX (Lx [e * (psx + jj + jj*nsrow)], lmin, lmax) ;
	    }
	}
    }
    else
    {
	/* L简单的 */
	Lp = L->p ;
	Lx = L->x ;
	if (L->is_ll)
	{
	    /* LL' */
	    FIRST_LMINMAX (Lx [Lp [0]], lmin, lmax) ;
	    for (j = 1 ; j < n ; j++)
	    {
		LMINMAX (Lx [e * Lp [j]], lmin, lmax) ;
	    }
	}
	else
	{
	    /* LDL' 分解, 对角线可能为负 */
	    FIRST_LMINMAX (fabs (Lx [Lp [0]]), lmin, lmax) ;
	    for (j = 1 ; j < n ; j++)
	    {
		LMINMAX (fabs (Lx [e * Lp [j]]), lmin, lmax) ;
	    }
	}
    }
    rcond = lmin / lmax ;
    if (L->is_ll)
    {
	rcond = rcond*rcond ;
    }
    return (rcond) ;
}
#endif
