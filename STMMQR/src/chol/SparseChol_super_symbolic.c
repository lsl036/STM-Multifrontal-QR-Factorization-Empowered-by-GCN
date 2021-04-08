/**
 * @file SparseChol_super_symbolic.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-22
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#include "Sparse_internal.h"
#include "SparseChol.h"
#include "tpsm_base.h"

/* 在对称情况下，从A的非零中遍历第k行子树(0:k1-1,k)，并将找到的新条目添加到l的第k行
 * 模式中。当前的超节点s包含对角块k1: k1-1，因此可以跳过它。
 *
 * 在非对称情况下，A*F的非零模式一次只计算一列(因此，该函数中花费的总时间以乘以A*F
 * 所花费的时间为界，如果A又高又瘦，这个时间可能很高)。第k列是A*F(:，k)，或者是
 * F(j,k)非零的所有列A(:，j)的并集。该例程对每个条目j调用一次，只需要访问上三角
 * 部分，因此只访问A (0:k1-1,j)，其中k1: k1- 2为当前超节点s的列(k在k1到k1-1的范围内)。
 */

static void subtree
(
    Int j,		
    Int k,
    Int Ap [ ],
    Int Ai [ ],
    Int Anz [ ],
    Int SuperMap [ ],
    Int Sparent [ ],
    Int mark,
    Int sorted,         
    Int k1,             
    Int Flag [ ],
    Int Ls [ ],
    Int Lpi2 [ ]
)
{
    Int p, pend, i, si ;
    p = Ap [j] ;
    pend = (Anz == NULL) ? (Ap [j+1]) : (p + Anz [j]) ;

    for ( ; p < pend ; p++)
    {
	i = Ai [p] ;
	if (i < k1)
	{
	    for (si = SuperMap [i] ; Flag [si] < mark ; si = Sparent [si])
	    {
		Ls [Lpi2 [si]++] = k ;
		Flag [si] = mark ;
	    }
	}
        else if (sorted)
        {
            break ;
        }
    }
}


#define FREE_WORKSPACE \
{ \
    /* SparseCore_clear_flag (Common) ; */ \
    SPARSE_CLEAR_FLAG (Common) ; \
    for (k = 0 ; k <= nfsuper ; k++) \
    { \
	Head [k] = EMPTY ; \
    } \
} \

// 确定每层超节点所在深度
void SuperTreeDepth(chol_supernode *Snode, int node_id, int depth)
{
    int i, j;
    Snode[node_id].Sdepth = depth;
    if ( Snode[node_id].Snchild != 0)
    {
        for (i = 0; i < Snode[node_id].Snchild; ++i)
        {
            j = Snode[node_id].Schild[i];
            SuperTreeDepth( Snode, j, depth+1);
        }
    }
}


int SparseChol_super_symbolic2
(
    int for_whom,       /* (0): SparseQR  (1): SparseChol */
    sparse_csc *A,	/* 将被分析的矩阵 */
    sparse_csc *F,	/* F = A' or A(:,f)' */
    Int *Parent,	/* 消去树 */
    sparse_factor *L,	/* 简易的符号分解因子作为输入，超节点符号分解因子作为输出 */
    sparse_common *Common
)
{
    double zrelax0, zrelax1, zrelax2, xxsize ;
    Int *Wi, *Wj, *Super, *Snz, *Ap, *Ai, *Flag, *Head, *Ls, *Lpi, *Lpx, *Fnz,
	*Sparent, *Anz, *SuperMap, *Merged, *Nscol, *Zeros, *Fp, *Fj,
	*ColCount, *Lpi2, *Lsuper, *Iwork ;
    Int nsuper, d, n, j, k, s, mark, parent, p, pend, k1, k2, packed, nscol,
	nsrow, ndrow1, ndrow2, stype, ssize, xsize, sparent, plast, slast,
	csize, maxcsize, ss, nscol0, nscol1, ns, nfsuper, newzeros, totzeros,
	merge, snext, esize, maxesize, nrelax0, nrelax1, nrelax2, Asorted ;
    size_t w ;
    int ok = TRUE, find_xsize;
    const char* env_max_bytes;
    size_t max_bytes;
    const char* env_max_fraction;
    double max_fraction;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, SPARSE_PATTERN, SPARSE_PATTERN, FALSE) ;
    stype = A->stype ;
    if (stype < 0)
    {
	ERROR (SPARSE_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }
    if (stype == 0)
    {
	RETURN_IF_NULL (F, FALSE) ;
    }
    if (L->is_super)
    {
	ERROR (SPARSE_INVALID, "L must be symbolic on input") ;
	return (FALSE) ;
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;

    /* w = 5*n */
    w = SparseCore_mult_size_t (n, 5, &ok) ;
    if (!ok)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    SparseCore_allocate_work (n, w, 0, Common) ;
    if (Common->status < SPARSE_OK)
    {
	/* 内存溢出 */
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 获取输入 */
    /* ---------------------------------------------------------------------- */

    /* 对于对称情况，A现在为A或triu(A(p,p))。对于非对称情况，它要么是A，要么是
    * A(p,f)，都是列形式。它可以是压缩的，也可以是未压缩的，也可以是排序的或未排序的。
    * 如果A是对称的，则下三角形部分中的项可能存在，但这些被忽略。*/

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;

    if (stype != 0)
    {
	Fp = NULL ;
	Fj = NULL ;
	Fnz = NULL ;
	packed = TRUE ;
    }
    else
    {
	Fp = F->p ;
	Fj = F->i ;
	Fnz = F->nz ;
	packed = F->packed ;
    }

    ColCount = L->ColCount ;

    nrelax0 = Common->nrelax [0] ; // 4
    nrelax1 = Common->nrelax [1] ; // 16
    nrelax2 = Common->nrelax [2] ; // 48

    zrelax0 = Common->zrelax [0] ; // 0.8
    zrelax1 = Common->zrelax [1] ; // 0.1
    zrelax2 = Common->zrelax [2] ; // 0.05

    zrelax0 = IS_NAN (zrelax0) ? 0 : zrelax0 ;
    zrelax1 = IS_NAN (zrelax1) ? 0 : zrelax1 ;
    zrelax2 = IS_NAN (zrelax2) ? 0 : zrelax2 ;

    /* ---------------------------------------------------------------------- */
    /* 获取工作空间 */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wi      = Iwork ;	    /* size n (i/l/l).  Lpi2 is i/l/l */
    Wj      = Iwork + n ;   /* size n (i/l/l).  Zeros is i/l/l */
    Sparent = Iwork + 2*((size_t) n) ; /* size nfsuper <= n [ */
    Snz     = Iwork + 3*((size_t) n) ; /* size nfsuper <= n [ */
    Merged  = Iwork + 4*((size_t) n) ; /* size nfsuper <= n [ */

    Flag = Common->Flag ;   /* size n */
    Head = Common->Head ;   /* size n+1 */

    /* ---------------------------------------------------------------------- */
    /* 找出基本的超节点 */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j < n ; j++)
    {
	    Wi [j] = 0 ;
    }
    for (j = 0 ; j < n ; j++)
    {
        parent = Parent [j] ;
        if (parent != EMPTY)
        {
            Wi [parent]++ ;     // Wi[j] 表示j节点的孩子数
        }
    }

    Super = Head ;  

    nfsuper = (n == 0) ? 0 : 1 ;	/* 基本超节点的数量 */
    Super [0] = 0 ;

    for (j = 1 ; j < n ; j++)
    {
	if (Parent [j-1] != j	   
	    || (ColCount [j-1] != ColCount [j] + 1) 
	    || Wi [j] > 1	    
	    )
	{
	    Super [nfsuper++] = j ;
	}
    }
    Super [nfsuper] = n ;

    Nscol = Wi ; 

    /* ---------------------------------------------------------------------- */
    /* 找出基本节点到超节点的映射 */
    /* ---------------------------------------------------------------------- */
    SuperMap = Wj ;

    for (s = 0 ; s < nfsuper ; s++)
    {
	for (k = Super [s] ; k < Super [s+1] ; k++)
	{
	    SuperMap [k] = s ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 重构基本超节点消去树 */
    /* ---------------------------------------------------------------------- */

    for (s = 0 ; s < nfsuper ; s++)
    {
	j = Super [s+1] - 1 ;
	parent = Parent [j] ;	
	Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
    }

    Zeros = Wj ;   
    /* ---------------------------------------------------------------------- */
    /* 松弛合并超节点，将少数零元作为非零元处理 */
    /* ---------------------------------------------------------------------- */
    // 初始化
    for (s = 0 ; s < nfsuper ; s++)
    {
        Merged [s] = EMPTY ;			
        Nscol [s] = Super [s+1] - Super [s] ;	
        Zeros [s] = 0 ;				
        Snz [s] = ColCount [Super [s]] ;  
    }

    for (s = nfsuper-2 ; s >= 0 ; s--)
    {
        double lnz1 ;

        ss = Sparent [s] ;
        if (ss == EMPTY)
        {
            continue ;
        }

	    for (ss = Sparent [s] ; Merged [ss] != EMPTY ; ss = Merged [ss]) ;
	    sparent = ss ;

	    /* ss是s的当前父节点 */
        for (ss = Sparent [s] ; Merged [ss] != EMPTY ; ss = snext)
        {
            snext = Merged [ss] ;
            Merged [ss] = sparent ;
        }

        /* 如果s+1不是s的当前父节点，不合并 */
        if (sparent != s+1)
        {
            continue ;
        }

        nscol0 = Nscol [s] ;	
        nscol1 = Nscol [s+1] ;
        ns = nscol0 + nscol1 ;

        totzeros = Zeros [s+1] ;	
        lnz1 = (double) (Snz [s+1]) ;	

	    /* 确定超节点s和s+1是否应该合并 */
        if (ns <= nrelax0)
        {
            merge = TRUE ;
        }
        else 
        {
            double lnz0 = Snz [s] ;	
            double xnewzeros = nscol0 * (lnz1 + nscol0 - lnz0) ; // 额外的零元

            newzeros = nscol0 * (Snz [s+1] + nscol0 - Snz [s]) ;
            if (xnewzeros == 0)// 没有新增的非零项
            {
            merge = TRUE ;
            }
            else
            {
            double xtotzeros = ((double) totzeros) + xnewzeros ;

            double xns = (double) ns ;
            double xtotsize  = (xns * (xns+1) / 2) + xns * (lnz1 - nscol1) ;
            double z = xtotzeros / xtotsize ; // 合并后的总零元除以总size，即非零元的占比

            Int totsize ;
            totsize  = (ns * (ns+1) / 2) + ns * (Snz [s+1] - nscol1) ;
            totzeros += newzeros ;

            merge = ((ns <= nrelax1 && z < zrelax0) ||
                (ns <= nrelax2 && z < zrelax1) ||
                        (z < zrelax2)) &&
                (xtotsize < Int_max / sizeof (double)) ;

            }
	    }

        if (merge)// 合并超节点
        {
            Zeros [s] = totzeros ;
            Merged [s+1] = s ;
            Snz [s] = nscol0 + Snz [s+1] ;
            Nscol [s] += Nscol [s+1] ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* 重构松弛超节点列表 */
    /* ---------------------------------------------------------------------- */

    nsuper = 0 ;
    for (s = 0 ; s < nfsuper ; s++)
    {
	if (Merged [s] == EMPTY)
	{
	    Super [nsuper] = Super [s] ;
	    Snz [nsuper] = Snz [s] ;
	    nsuper++ ;
	}
    }
    Super [nsuper] = n ;

    /* ---------------------------------------------------------------------- */
    /* 找出松弛节点到超节点的映射 */
    /* ---------------------------------------------------------------------- */
    /* SuperMap [k] = s 如果列k被包含在超节点s中 */
    for (s = 0 ; s < nsuper ; s++)
    {
        for (k = Super [s] ; k < Super [s+1] ; k++)
        {
            SuperMap [k] = s ;
        }
    }

    if (for_whom == 1) // Cholesky
    {
        chol_supernode *Snode;
        int tmp_i, tmp_j, max_depth = 0, *Snode_num, **Snode_distribution;
        // 初始化超节点结构
        Snode = TPSM_Malloc_Align (nsuper * sizeof(chol_supernode));
		TPSM_assert(Snode,1);
        for (s = 0 ; s < nsuper ; ++s)
        {
            Snode[s].Snchild = 0;
            Snode[s].Sdepth = 0;
            Snode[s].Schild = NULL;
        }
        /* ---------------------------------------------------------------------- */
        /* 重构松弛超节点消去树 */
        /* ---------------------------------------------------------------------- */
        int Sroot_num = 0;
        int Sroot[nsuper];
        for (s = 0 ; s < nsuper ; ++s)
        {
            j = Super [s+1] - 1 ;	
            parent = Parent [j] ;	
            Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
            if (Sparent [s] != EMPTY){
                ++(Snode[Sparent [s]].Snchild);
            }
            else
            {
                Sroot[Sroot_num] = s;
                Sroot_num++;
            }
            
        }
        // 初始化每个超节点的孩子节点
        for (s = 0 ; s < nsuper ; ++s)
        {
            int child_num = Snode[s].Snchild;
            Snode[s].Schild = (int*)TPSM_Malloc_Align (child_num * sizeof(int));
			TPSM_assert(Snode[s].Schild,1);
            for(tmp_i = 0; tmp_i < child_num; ++tmp_i){
                Snode[s].Schild[tmp_i] = -1;
            }
        }

        // 确定每个超节点的孩子节点
        for (s = 0 ; s < nsuper ; ++s)
        {
            if (Sparent [s] != EMPTY){
                int parent = Sparent [s];
                int child_num = Snode[parent].Snchild;
                for(tmp_i = 0; tmp_i < child_num; ++tmp_i){
                    if(Snode[parent].Schild[tmp_i] == -1){
                        Snode[parent].Schild[tmp_i] = s;
                        break;
                    }
                }
            }
        }
        L->Snode = Snode;

        // 从根节点开始递归，确定每层节点深度
        for(tmp_i = 0; tmp_i < Sroot_num; ++tmp_i)
        {
            // printf("root Snode is %d\n", Sroot[tmp_i]);
            SuperTreeDepth(Snode, Sroot[tmp_i], 0);
        }
        
        // 确定超节点消去树深度
        for (s = 0 ; s < nsuper ; ++s)
        {
            if (Snode[s].Sdepth > max_depth){
                max_depth = Snode[s].Sdepth;
            }
        }  

        // 分配空间并初始化
        Snode_num = (int*)TPSM_Malloc_Align((max_depth+1)*sizeof(int));
		TPSM_assert(Snode_num,1);
        for(tmp_i = 0; tmp_i < max_depth+1; ++tmp_i){
            Snode_num[tmp_i] = 0;
        }

        // 统计每层超节点数量
        for (s = 0 ; s < nsuper ; ++s)
        {
            ++Snode_num[Snode[s].Sdepth]; 
        }  
        // 分配空间并初始化
        Snode_distribution = (int**)TPSM_Malloc_Align((max_depth+1)*sizeof(int*));
		TPSM_assert(Snode_distribution,1);
        for(tmp_i = 0; tmp_i < max_depth+1; ++tmp_i)
        {
            Snode_distribution[tmp_i] = (int*)TPSM_Malloc_Align(Snode_num[tmp_i]*sizeof(int));
			TPSM_assert(Snode_distribution[tmp_i],1);
        }
        for(tmp_i = 0; tmp_i < max_depth+1; ++tmp_i)
        {
            for(tmp_j = 0; tmp_j < Snode_num[tmp_i]; ++tmp_j)
            {
                Snode_distribution[tmp_i][tmp_j] = -1;
            }
        }

        // 为每个超节点在每层深度分配位置
        for (s = 0 ; s < nsuper ; ++s)
        {
            for(tmp_i = 0; tmp_i < Snode_num[Snode[s].Sdepth]; ++tmp_i)
            {
                if(Snode_distribution[Snode[s].Sdepth][tmp_i] == -1)
                {
                    Snode_distribution[Snode[s].Sdepth][tmp_i] = s;
                    break;
                }
            }
        }  
        // int i;
        // for(i = 0 ; i <= max_depth; ++i)
        // {
        //     printf("\n\ndepth %d\n\n", i);
        //     for(j = 0; j < Snode_num[i]; ++j)
        //     {
        //         printf("%d ", Snode_distribution[i][j]);
        //     }
        // }

        L->max_depth = max_depth;
        L->Snode_num = Snode_num;
        L->Snode_distribution = Snode_distribution;
    }
    else // QR
    {
        /* ---------------------------------------------------------------------- */
        /* 重构松弛超节点消去树 */
        /* ---------------------------------------------------------------------- */
        for (s = 0 ; s < nsuper ; ++s)
        {
            j = Super [s+1] - 1 ;	/* 超节点s中的最后节点 */
            parent = Parent [j] ;	/* 最后节点的父节点 */
            Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
        }
    }



    /* ---------------------------------------------------------------------- */
    /* 确定L->s和L->x的大小 */
    /* ---------------------------------------------------------------------- */

    ssize = 0 ;
    xsize = 0 ;
    xxsize = 0 ;
    find_xsize = for_whom == SPARSE_ANALYZE_FOR_CHOLESKY ;
    for (s = 0 ; s < nsuper ; s++)
    {
	nscol = Super [s+1] - Super [s] ;
	nsrow = Snz [s] ;
	ssize += nsrow ;
    if (find_xsize)
    {
        xsize += nscol * nsrow ;
        xxsize += ((double) nscol) * ((double) nsrow) ;
    }
	if (ssize < 0 ||(find_xsize && xxsize > Int_max))
	{
	    ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	    FREE_WORKSPACE ;
	    return (FALSE) ;
	}
    }
    xsize = MAX (1, xsize) ;
    ssize = MAX (1, ssize) ;

    L->ssize = ssize ;
    L->xsize = xsize ;
    L->nsuper = nsuper ;

    SparseCore_change_factor (SPARSE_PATTERN, TRUE, TRUE, TRUE, TRUE, L, Common);

    if (Common->status < SPARSE_OK)
    {
	FREE_WORKSPACE ;
	return (FALSE) ;
    }

    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Ls [0] = 0 ;  
    Lsuper = L->super ;

    /* 复制松弛超节点列表到L中最终列表 */
    for (s = 0 ; s <= nsuper ; s++)
    {
	Lsuper [s] = Super [s] ;
    }

    Super = Lsuper ;	    /* Super现在是松弛超节点的清单 */

    /* ---------------------------------------------------------------------- */
    /*构造松弛超节点模式的列指针 (L->pi) */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (s = 0 ; s < nsuper ; s++)
    {
	Lpi [s] = p ;
	p += Snz [s] ;
    }
    Lpi [nsuper] = p ;

    /* ---------------------------------------------------------------------- */
    /* 构造超节点的值的指针 (L->px) */
    /* ---------------------------------------------------------------------- */

    if (for_whom == SPARSE_ANALYZE_FOR_CHOLESKY)
    {
        Lpx [0] = 0 ;
        p = 0 ;
        for (s = 0 ; s < nsuper ; s++)
        {
            nscol = Super [s+1] - Super [s] ;   /* s中的列数 */
            nsrow = Snz [s] ;           /* 非零行数 */
            Lpx [s] = p ;               /* s数值部分的指针 */
            p += nscol * nsrow ;
        }
        Lpx [s] = p ;
    }
    else   //QR
    {
        Lpx [0] = 123456 ;
    }

    /* ---------------------------------------------------------------------- */
    /* 构造松弛超节点模式的符号分析 (L->s) */
    /* ---------------------------------------------------------------------- */

    Lpi2 = Wi ;	   
    for (s = 0 ; s < nsuper ; s++)
    {
	Lpi2 [s] = Lpi [s] ;
    }

    Asorted = A->sorted ;

    for (s = 0 ; s < nsuper ; s++)
    {
	/* 超节点s是第k1到k2-1列。
	 * 计算 L (k1:k2-1,:)的非零模式 */
	k1 = Super [s] ;
	k2 = Super [s+1] ;
	for (k = k1 ; k < k2 ; k++)
	{
	    Ls [Lpi2 [s]++] = k ;
	}

	/* 计算每行k1到k2-1的非零模式 */
	for (k = k1 ; k < k2 ; k++)
	{
	    /* 计算L的第k行。在对称情况下，L(k，:)的模式是A(0:k,k)的非零模式
        中任意第i行在超节点etree上可达的节点集。在非对称情况下，A*A'的
        第k列的模式是每个非零F(j,k)的所有列A(0:k,j)的并集。*/

	    SPARSE_CLEAR_FLAG (Common) ;
	    mark = Common->mark ;
	    Flag [s] = mark ;

	    if (stype != 0)
	    {
            subtree (k, k, Ap, Ai, Anz, SuperMap, Sparent, mark,
                            Asorted, k1, Flag, Ls, Lpi2) ;
	    }
	    else
	    {
            p = Fp [k] ;
            pend = (packed) ? (Fp [k+1]) : (p + Fnz [k]) ;
            for ( ; p < pend ; p++)
            {
                subtree (Fj [p], k, Ap, Ai, Anz, SuperMap, Sparent, mark,
                    Asorted, k1, Flag, Ls, Lpi2) ;
            }
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 确定最大的更新矩阵 (L->maxcsize) */
    /* ---------------------------------------------------------------------- */

    /* 在分配和定义L->s之前可以确定maxcsize，这意味着符号和数值分解的所有内存需求
     * 都可以使用O(nnz(A)+O(n))空间计算。然而，这将需要大量的额外工作。上面的分析阶
     * 段需要进行复制，但不保留Ls;相反,对于每个超节点d,该算法将跟踪当前的s和slast,并
     * 在一个新行索引出现在超节点d时更新它们。另一种是计算只有L→分配失败,在这种情
     * 况下,下面的代码将被忽略。
     *
     * 超节点的csize是它对后续祖先超节点的最大贡献的大小。例如，假设下面图中的#的行对
     * 应于后续超节点的列，而点就是该祖先中的条目。
     * 
     *	    c
     *	    c c
     *	    c c c
     *	    x x x
     *	    x x x
     *	    # # #   .
     *	    # # #   . .
     *	    * * *   . .
     *	    * * *   . .
     *	    * * *   . .
     *	            . .
     *
     * 那么对于这个更新，csize是3乘2，或者6，因为有3行*，这是更新中的行数，有2行#，
     * 这是更新中的列数。超节点的csize是所有祖先超节点的最大贡献。对于整个矩阵，
     * maxcsize具有任意超节点的最大大小的粗略上界。这个绑定是松散的，因为贡献必须
     * 小于它正在更新的祖先超节点的大小。具有一个超节点的完全稠密矩阵的maxcsize为零。
     *
     * maxesize是该求解所需的工作空间E的列维。E的大小为nrhs*maxesize，其中nrhs
     * 是右手边的列数。maxesize是任何超节点的最大esize。超节点的esize是它包含的行
     * 索引数，不包括超节点本身的列索引。下面的例子中，esize为4:
     *	    c
     *	    c c
     *	    c c c
     *	    x x x
     *	    x x x
     *	    x x x
     *	    x x x
     *
     * maxesize不会比n大.
     */

    maxcsize = 1 ;
    maxesize = 1 ;

    if (for_whom == SPARSE_ANALYZE_FOR_CHOLESKY)
    {
        for (d = 0 ; d < nsuper ; d++)
        {
            nscol = Super [d+1] - Super [d] ;
            p = Lpi [d] + nscol ;
            plast = p ;
            pend = Lpi [d+1] ;
            esize = pend - p ;
            maxesize = MAX (maxesize, esize) ;
            slast = (p == pend) ? (EMPTY) : (SuperMap [Ls [p]]) ;
            for ( ; p <= pend ; p++)
            {
                s = (p == pend) ? (EMPTY) : (SuperMap [Ls [p]]) ;
                if (s != slast)
                {
                    /* 行i是一个新的超节点的开始 */
                    ndrow1 = p - plast ;
                    ndrow2 = pend - plast ;
                    csize = ndrow2 * ndrow1 ;
                    maxcsize = MAX (maxcsize, csize) ;
                    plast = p ;
                    slast = s ;
                }
            }
        }
    }

    L->maxcsize = maxcsize ;
    L->maxesize = maxesize ;
    L->is_super = TRUE ;

    FREE_WORKSPACE ;
    return (TRUE) ;
}

/* 分析A, AA'，或A(:，f)*A(:，f)'，以准备超节点数值分解。 */
int SparseChol_super_symbolic
(
    sparse_csc *A,
    sparse_csc *F,	
    Int *Parent,	
    sparse_factor *L,	
    sparse_common *Common
)
{
    return (SparseChol_super_symbolic2 (SPARSE_ANALYZE_FOR_CHOLESKY,
        A, F, Parent, L, Common)) ;
}

