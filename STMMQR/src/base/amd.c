#include "amd_internal.h"

/**
 * 函数名: AMD_1	
 * 参数1: n					input ,	n > 0
 * 参数2: Ap				input , size: n+1, 不修改
 * 参数3: Ai				input , size: nz = Ap[n], 不修改
 * 参数4: P 				output, size: n, 置换数组
 * 参数5: Pinv 				output, size: n, 反转置换数组
 * 参数6: Len 				input , size: n
 * 参数7: slen 						slen >= sum (Len [0..n-1]) + 7n, ideally slen = 1.2 * sum (Len) + 8n
 * 参数8: S 						size: slen, 工作空间
 * 参数9: Control 			input , size: AMD_CONTROL
 * 参数10: Info 			output, size: AMD_INFO
 * 
 * 功能: 对于稀疏矩阵A, 构造A+A'. 并执行AMD重排序
 * 
 * 详情: n×n稀疏矩阵A可能是不对称的。它以matlab样式的压缩列形式存储，
 * 		每个列中都有已排序的行索引，并且没有重复条目。可能会出现对角线条目，
 * 		但它们会被忽略。A的j列的行索引存储在Ai [Ap [j] ... Ap [j+1]-1]。
 * 		Ap[0]必须为零，nz = Ap [n]为A中的项目数，矩阵的大小n ≥ 0 。
 * 
 * 注意事项: 这个例程之前必须调用AMD_aat，它计算A+A'中每行/列中的条目数， 不包括对角线。
 * 			Len [j]，在输入上，是A+A'的第j行/列的项目数。这个例程构造矩阵A+A'，
 * 			然后调用AMD_2。没有执行错误检查(这是在AMD_valid中完成的)。
 **/

GLOBAL void AMD_1
(
    Int n,		
    const Int Ap [ ],	
    const Int Ai [ ],	
    Int P [ ],		
    Int Pinv [ ],	
    Int Len [ ],	
    Int slen,		
    Int S [ ],		
    double Control [ ],	
    double Info [ ]	
)
{
    Int i, j, k, p, pfree, iwlen, pj, p1, p2, pj2, *Iw, *Pe, *Nv, *Head,
	*Elen, *Degree, *s, *W, *Sp, *Tp ;

    /* ----------------- */
    /*  为AMD 2构造矩阵  */
    /* ----------------- */

    ASSERT (n > 0) ;

    iwlen = slen - 6*n ;
    s = S ;
    Pe = s ;	    s += n ;
    Nv = s ;	    s += n ;
    Head = s ;	    s += n ;
    Elen = s ;	    s += n ;
    Degree = s ;    s += n ;
    W = s ;	    s += n ;
    Iw = s ;	    s += iwlen ;

    ASSERT (AMD_valid (n, n, Ap, Ai) == AMD_OK) ;

    /* 构造A + A' 的指针 */
    Sp = Nv ;			/* 使用 Nv 和 W 作为 Sp 和 Tp[ 的工作区 */
    Tp = W ;
    pfree = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	Pe [j] = pfree ;
	Sp [j] = pfree ;
	pfree += Len [j] ;
    }

    /* 注意，对 iwlen 的这个限制比AMD_2中的要求的要稍微严格一些。
     * AMD_2可以在完全没有足够空间的情况下运行，但速度会非常慢。
	 * 为了获得更好的性能，至少需要大小为n的足够空间。 */
    ASSERT (iwlen >= pfree + n) ;

    for (k = 0 ; k < n ; k++)
    {
	AMD_DEBUG1 (("Construct row/column k= "ID" of A+A'\n", k))  ;
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;

	/* 构造 A+A' */
	for (p = p1 ; p < p2 ; )
	{
	    /* 读取A的上三角部分 */
	    j = Ai [p] ;
	    ASSERT (j >= 0 && j < n) ;
	    if (j < k)
	    {
		/* 元素 A (j,k) 属于 A 的严格上三角部分 */
		ASSERT (Sp [j] < (j == n-1 ? pfree : Pe [j+1])) ;
		ASSERT (Sp [k] < (k == n-1 ? pfree : Pe [k+1])) ;
		Iw [Sp [j]++] = k ;
		Iw [Sp [k]++] = j ;
		p++ ;
	    }
	    else if (j == k)
	    {
		/* 跳过对角线元素 */
		p++ ;
		break ;
	    }
	    else /* j > k */
	    {
		/* 第一个低于对角线的元素 */
		break ;
	    }
	    /* 读取A的下三角部分, 在第 j 列中读到第k行停止。
	     * 从上次扫描结束的地方开始 */
	    ASSERT (Ap [j] <= Tp [j] && Tp [j] <= Ap [j+1]) ;
	    pj2 = Ap [j+1] ;
	    for (pj = Tp [j] ; pj < pj2 ; )
	    {
		i = Ai [pj] ;
		ASSERT (i >= 0 && i < n) ;
		if (i < k)
		{
		    /* A (i,j) 只是下三角部分的元素 */
		    ASSERT (Sp [i] < (i == n-1 ? pfree : Pe [i+1])) ;
		    ASSERT (Sp [j] < (j == n-1 ? pfree : Pe [j+1])) ;
		    Iw [Sp [i]++] = j ;
		    Iw [Sp [j]++] = i ;
		    pj++ ;
		}
		else if (i == k)
		{
		    /* 元素 A (k,j) 在下三角部分 and A (j,k) 在上三角部分 */
		    pj++ ;
		    break ;
		}
		else /* i > k */
		{
		    /* 当 k 前进到 i 时，之后再考虑*/
		    break ;
		}
	    }
	    Tp [j] = pj ;
	}
	Tp [k] = p ;
    }

    /* 清除剩余的不匹配条目 */
    for (j = 0 ; j < n ; j++)
    {
	for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
	{
	    i = Ai [pj] ;
	    ASSERT (i >= 0 && i < n) ;
	    /* A (i,j) 只是下三角部分的元素 */
	    ASSERT (Sp [i] < (i == n-1 ? pfree : Pe [i+1])) ;
	    ASSERT (Sp [j] < (j == n-1 ? pfree : Pe [j+1])) ;
	    Iw [Sp [i]++] = j ;
	    Iw [Sp [j]++] = i ;
	}
    }

    /* Tp 和 Sp 不再需要 ] */

    /* --------------------------------------------------------------------- */
    /* 重序矩阵 */
    /* --------------------------------------------------------------------- */

    AMD_2 (n, Pe, Iw, Len, iwlen, pfree,
	Nv, Pinv, P, Head, Elen, Degree, W, Control, Info) ;
}


/* ========================================================================= */
/* === clear_flag ========================================================== */
/* ========================================================================= */

static Int clear_flag (Int wflg, Int wbig, Int W [ ], Int n)
{
    Int x ;
    if (wflg < 2 || wflg >= wbig)
    {
	for (x = 0 ; x < n ; x++)
	{
	    if (W [x] != 0) W [x] = 1 ;
	}
	wflg = 2 ;
    }
    /*  此时, W [0..n-1] < wflg 保持 */
    return (wflg) ;
}


/* ========================================================================= */
/* === AMD_2 =============================================================== */
/* ========================================================================= */
/* 函数名: AMD_2  
 * 参数1: n 				Input: 矩阵A是 n * n 的,且 n > 0
 * 参数2: Pe[]				Input/Output: size n 的一个整型数组. 输入时表示第i行在 Iw 中的索引
 * 										  如果第i行没有非对角条目，则忽略Pe [i]。
 * 										  因此，对于非空行，Pe [i]必须在0到pfree-1的范围内。
 *  	
 * 		在执行时，它用来保存超级变量 和 基本元素。
 * 	     - 超级变量i ： 对超级变量i描述的索引Iw。 一个超级变量表示矩阵的一行或多行
 * 						具有相近的非零元分布模式，此时 Pe [i] >= 0.
 * 
 * 		 - 非超级变量i：如果i被另一个超变量j吸收，则Pe [i] = FLIP (j);
 * 						其中FLIP (j)定义为(-(j)-2)。第j行和第i行有相同的模式。
 * 						注意，j以后可能被吸收到另一个超变量j2中，此时Pe [i]仍然是FLIP (j);
 * 						Pe [j] = FLIP (j2) < EMPTY ;在amd_internal.h中，EMPTY被定义为(-1)
 * 
 * 参数3: Iw[]				Input: 大小为iwlen的工作空间，Iw[0 .. pfree-1] 保存输入矩阵
 * 参数4: Len[]				Input: Len[0..n-1]: 输入的 行/列 i 的长度
 * 								   在输入上，Len [i]保存矩阵第i行中元素的数量，
 * 								   不包括对角元素。Len的内容在输出中未定义。
 * 参数5: iwlen				Input: Iw的长度， iwlen >= pfree + n,否则，将会发生过度的压缩
 * 参数6: pfree				Input: 在数组的尾部Iw[pfree .. iwlen-1] 输入时为空.
 * 								   矩阵存储在Iw[0 .. pfree-1], 附加数据被放置在Iw中,
 * 								   pfree被修改.所以Iw[iwren-1]一直是 Iw 中未使用的部分。
 *   	
 * 			注意:输入矩阵会在计算过程中被覆盖。Iw 的内容在输出中没有定义。
 *  
 * 7个大小为n的工作区，未在输入中定义:
 * 参数7: Nv[]				Output: size n 输出每个超节点的大小
 * 									在执行过程中，abs(Nv [i])等于主超变量i所代表的主元数。
 * 									输出时，Nv [i]表示原始矩阵的超级 row/column i所表示的数
 * 									对于非主行/列，Nv [i] = 0。
 * 
 * 参数8: Next[]			Output: size n 输出逆置换
 * 参数9: Last[]			Output: size n 输出重排列
 * 参数10: Head []			Output: size n Head [deg]是degree表中的第一个超变量
 * 参数11: Elen []			Output: size n 执行时，Elen [i]为超变量i列表中的元素个数，
 * 									当e成为元素时，设置Elen [e] = FLIP (esize)，其中
 * 									esize 是element 的大小(主元素+非主元素)
 * 参数12: Degree[]			Output: size n 一个大小为n的整数数组。如果i是一个超变量，
 * 									那么Degree[i]保持当前第i行外部度的近似值(一个上界)。
 * 									外部度是第i行中非零的数量，减去ABS (Nv [i])，即对角线部分。
 * 参数13: W []				Output: size n 标志数组W决定元素和变量的状态，以及元素的外部度。
 * 
 * 控制参数和输出统计:	
 * 参数14: Control[]		Input: 大小为AMD_CONTROL的double数组,其中包含影响排序计算方式的输入参数。
 * 							       如果为空，则使用默认设置。
 * 参数15: Info[]			Ouput: AMD_INFO 数组
 * 
 * 功能: 在对称稀疏矩阵 A 上执行AMD排序，然后使用AMD_postorder例程（基于深度优先）
 * 		  对组合树进行后序排序。
 * 
 */
GLOBAL void AMD_2
(
    Int n,		
    Int Pe [ ],		
    Int Iw [ ],		
    Int Len [ ],	
    Int iwlen,		
    Int pfree,		
    Int Nv [ ],		
    Int Next [ ],	
    Int Last [ ],	
    Int Head [ ],
    Int Elen [ ],	
    Int Degree [ ],
    Int W [ ],
    double Control [ ],	
    double Info [ ]
)
{

/*
 * 给出了对称矩阵 A 的非零模式的一种表示（排除对角线）执行近似最小度排序
 * 来计算主元顺序，这样Cholesky因子a = LL'中的非零(填充)的引入保持在低水平。  
 * 在每一步中，选择的枢轴是在外部度上具有最小UMFAPACK/ ma38样式上限的枢轴。
 * 这个例程可以选择性地执行聚合吸收(就像MC47B在Harwell子例程库中所做的那样)。
 *
 * 这里实现的近似度算法是对MA38和UMFPACK (Davis和Duff的不对称模式多额包)
 * 中的度更新算法的对称模拟。该程序基于Iain Duff和John Reid提出的MA27最小度
 * 排序算法。
 *

 * ----------------------------------------------------------------------------
 * LOCAL INTEGERS:
 * ----------------------------------------------------------------------------
 */

    Int deg, degme, dext, lemax, e, elenme, eln, i, ilast, inext, j,
	jlast, jnext, k, knt1, knt2, knt3, lenj, ln, me, mindeg, nel, nleft,
	nvi, nvj, nvpiv, slenme, wbig, we, wflg, wnvi, ok, ndense, ncmpa,
	dense, aggressive ;

    unsigned Int hash ;	    /* unsigned, so that hash % n is well defined.*/

/*
 * deg:		the degree of a variable or element
 * degme:	size, |Lme|, of the current element, me (= Degree [me])
 * dext:	external degree, |Le \ Lme|, of some element e
 * lemax:	largest |Le| seen so far (called dmax in Fortran version)
 * e:		an element
 * elenme:	the length, Elen [me], of element list of pivotal variable
 * eln:		the length, Elen [...], of an element list
 * hash:	the computed value of the hash function
 * i:		a supervariable
 * ilast:	the entry in a link list preceding i
 * inext:	the entry in a link list following i
 * j:		a supervariable
 * jlast:	the entry in a link list preceding j
 * jnext:	the entry in a link list, or path, following j
 * k:		the pivot order of an element or variable
 * knt1:	loop counter used during element construction
 * knt2:	loop counter used during element construction
 * knt3:	loop counter used during compression
 * lenj:	Len [j]
 * ln:		length of a supervariable list
 * me:		current supervariable being eliminated, and the current
 *		    element created by eliminating that supervariable
 * mindeg:	current minimum degree
 * nel:		number of pivots selected so far
 * nleft:	n - nel, the number of nonpivotal rows/columns remaining
 * nvi:		the number of variables in a supervariable i (= Nv [i])
 * nvj:		the number of variables in a supervariable j (= Nv [j])
 * nvpiv:	number of pivots in current element
 * slenme:	number of variables in variable list of pivotal variable
 * wbig:	= (INT_MAX - n) for the int version, (Sparse_long_max - n)
 *                  for the Sparse_long version.  wflg is not allowed to
 *                  be >= wbig.
 * we:		W [e]
 * wflg:	used for flagging the W array.  See description of Iw.
 * wnvi:	wflg - Nv [i]
 * x:		either a supervariable or an element
 *
 * ok:		true if supervariable j can be absorbed into i
 * ndense:	number of "dense" rows/columns
 * dense:	rows/columns with initial degree > dense are considered "dense"
 * aggressive:	true if aggressive absorption is being performed
 * ncmpa:	number of garbage collections

 * ----------------------------------------------------------------------------
 * LOCAL DOUBLES, used for statistical output only (except for alpha):
 * ----------------------------------------------------------------------------
 */

    double f, r, ndiv, s, nms_lu, nms_ldl, dmax, alpha, lnz, lnzme ;

/*
 * f:		nvpiv
 * r:		degme + nvpiv
 * ndiv:	number of divisions for LU or LDL' factorizations
 * s:		number of multiply-subtract pairs for LU factorization, for the
 *		    current element me
 * nms_lu	number of multiply-subtract pairs for LU factorization
 * nms_ldl	number of multiply-subtract pairs for LDL' factorization
 * dmax:	the largest number of entries in any column of L, including the
 *		    diagonal
 * alpha:	"dense" degree ratio
 * lnz:		the number of nonzeros in L (excluding the diagonal)
 * lnzme:	the number of nonzeros in L (excl. the diagonal) for the
 *		    current element me

 * ----------------------------------------------------------------------------
 * LOCAL "POINTERS" (indices into the Iw array)
 * ----------------------------------------------------------------------------
*/

    Int p, p1, p2, p3, p4, pdst, pend, pj, pme, pme1, pme2, pn, psrc ;

/*
 * Any parameter (Pe [...] or pfree) or local variable starting with "p" (for
 * Pointer) is an index into Iw, and all indices into Iw use variables starting
 * with "p."  The only exception to this rule is the iwlen input argument.
 *
 * p:           pointer into lots of things
 * p1:          Pe [i] for some variable i (start of element list)
 * p2:          Pe [i] + Elen [i] -  1 for some variable i
 * p3:          index of first supervariable in clean list
 * p4:		
 * pdst:        destination pointer, for compression
 * pend:        end of memory to compress
 * pj:          pointer into an element or variable
 * pme:         pointer into the current element (pme1...pme2)
 * pme1:        the current element, me, is stored in Iw [pme1...pme2]
 * pme2:        the end of the current element
 * pn:          pointer into a "clean" variable, also used to compress
 * psrc:        source pointer, for compression
*/

/* ========================================================================= */
/*  INITIALIZATIONS */
/* ========================================================================= */

    /* Note that this restriction on iwlen is slightly more restrictive than
     * what is actually required in AMD_2.  AMD_2 can operate with no elbow
     * room at all, but it will be slow.  For better performance, at least
     * size-n elbow room is enforced. */
    ASSERT (iwlen >= pfree + n) ;
    ASSERT (n > 0) ;

    /* initialize output statistics */
    lnz = 0 ;
    ndiv = 0 ;
    nms_lu = 0 ;
    nms_ldl = 0 ;
    dmax = 1 ;
    me = EMPTY ;

    mindeg = 0 ;
    ncmpa = 0 ;
    nel = 0 ;
    lemax = 0 ;

    /* get control parameters */
    if (Control != (double *) NULL)
    {
	alpha = Control [AMD_DENSE] ;
	aggressive = (Control [AMD_AGGRESSIVE] != 0) ;
    }
    else
    {
	alpha = AMD_DEFAULT_DENSE ;
	aggressive = AMD_DEFAULT_AGGRESSIVE ;
    }
    /* Note: if alpha is NaN, this is undefined: */
    if (alpha < 0)
    {
	/* only remove completely dense rows/columns */
	dense = n-2 ;
    }
    else
    {
	dense = alpha * sqrt ((double) n) ;
    }
    dense = MAX (16, dense) ;
    dense = MIN (n,  dense) ;
    AMD_DEBUG1 (("\n\nAMD (debug), alpha %g, aggr. "ID"\n",
	alpha, aggressive)) ;

    for (i = 0 ; i < n ; i++)
    {
	Last [i] = EMPTY ;
	Head [i] = EMPTY ;
	Next [i] = EMPTY ;
	/* if separate Hhead array is used for hash buckets: *
	Hhead [i] = EMPTY ;
	*/
	Nv [i] = 1 ;
	W [i] = 1 ;
	Elen [i] = 0 ;
	Degree [i] = Len [i] ;
    }

    /* initialize wflg */
    wbig = Int_MAX - n ;
    wflg = clear_flag (0, wbig, W, n) ;

    /* --------------------------------------------------------------------- */
    /* initialize degree lists and eliminate dense and empty rows */
    /* --------------------------------------------------------------------- */

    ndense = 0 ;

    for (i = 0 ; i < n ; i++)
    {
	deg = Degree [i] ;
	ASSERT (deg >= 0 && deg < n) ;
	if (deg == 0)
	{

	    /* -------------------------------------------------------------
	     * we have a variable that can be eliminated at once because
	     * there is no off-diagonal non-zero in its row.  Note that
	     * Nv [i] = 1 for an empty variable i.  It is treated just
	     * the same as an eliminated element i.
	     * ------------------------------------------------------------- */

	    Elen [i] = FLIP (1) ;
	    nel++ ;
	    Pe [i] = EMPTY ;
	    W [i] = 0 ;

	}
	else if (deg > dense)
	{

	    /* -------------------------------------------------------------
	     * Dense variables are not treated as elements, but as unordered,
	     * non-principal variables that have no parent.  They do not take
	     * part in the postorder, since Nv [i] = 0.  Note that the Fortran
	     * version does not have this option.
	     * ------------------------------------------------------------- */

	    AMD_DEBUG1 (("Dense node "ID" degree "ID"\n", i, deg)) ;
	    ndense++ ;
	    Nv [i] = 0 ;		/* do not postorder this node */
	    Elen [i] = EMPTY ;
	    nel++ ;
	    Pe [i] = EMPTY ;

	}
	else
	{

	    /* -------------------------------------------------------------
	     * place i in the degree list corresponding to its degree
	     * ------------------------------------------------------------- */

	    inext = Head [deg] ;
	    ASSERT (inext >= EMPTY && inext < n) ;
	    if (inext != EMPTY) Last [inext] = i ;
	    Next [i] = inext ;
	    Head [deg] = i ;

	}
    }

/* ========================================================================= */
/* WHILE (selecting pivots) DO */
/* ========================================================================= */

    while (nel < n)
    {

/* ========================================================================= */
/* GET PIVOT OF MINIMUM DEGREE */
/* ========================================================================= */

	/* ----------------------------------------------------------------- */
	/* find next supervariable for elimination */
	/* ----------------------------------------------------------------- */

	ASSERT (mindeg >= 0 && mindeg < n) ;
	for (deg = mindeg ; deg < n ; deg++)
	{
	    me = Head [deg] ;
	    if (me != EMPTY) break ;
	}
	mindeg = deg ;
	ASSERT (me >= 0 && me < n) ;
	AMD_DEBUG1 (("=================me: "ID"\n", me)) ;

	/* ----------------------------------------------------------------- */
	/* remove chosen variable from link list */
	/* ----------------------------------------------------------------- */

	inext = Next [me] ;
	ASSERT (inext >= EMPTY && inext < n) ;
	if (inext != EMPTY) Last [inext] = EMPTY ;
	Head [deg] = inext ;

	/* ----------------------------------------------------------------- */
	/* me represents the elimination of pivots nel to nel+Nv[me]-1. */
	/* place me itself as the first in this set. */
	/* ----------------------------------------------------------------- */

	elenme = Elen [me] ;
	nvpiv = Nv [me] ;
	ASSERT (nvpiv > 0) ;
	nel += nvpiv ;

/* ========================================================================= */
/* CONSTRUCT NEW ELEMENT */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * At this point, me is the pivotal supervariable.  It will be
	 * converted into the current element.  Scan list of the pivotal
	 * supervariable, me, setting tree pointers and constructing new list
	 * of supervariables for the new element, me.  p is a pointer to the
	 * current position in the old list.
	 * ----------------------------------------------------------------- */

	/* flag the variable "me" as being in Lme by negating Nv [me] */
	Nv [me] = -nvpiv ;
	degme = 0 ;
	ASSERT (Pe [me] >= 0 && Pe [me] < iwlen) ;

	if (elenme == 0)
	{

	    /* ------------------------------------------------------------- */
	    /* construct the new element in place */
	    /* ------------------------------------------------------------- */

	    pme1 = Pe [me] ;
	    pme2 = pme1 - 1 ;

	    for (p = pme1 ; p <= pme1 + Len [me] - 1 ; p++)
	    {
		i = Iw [p] ;
		ASSERT (i >= 0 && i < n && Nv [i] >= 0) ;
		nvi = Nv [i] ;
		if (nvi > 0)
		{

		    /* ----------------------------------------------------- */
		    /* i is a principal variable not yet placed in Lme. */
		    /* store i in new list */
		    /* ----------------------------------------------------- */

		    /* flag i as being in Lme by negating Nv [i] */
		    degme += nvi ;
		    Nv [i] = -nvi ;
		    Iw [++pme2] = i ;

		    /* ----------------------------------------------------- */
		    /* remove variable i from degree list. */
		    /* ----------------------------------------------------- */

		    ilast = Last [i] ;
		    inext = Next [i] ;
		    ASSERT (ilast >= EMPTY && ilast < n) ;
		    ASSERT (inext >= EMPTY && inext < n) ;
		    if (inext != EMPTY) Last [inext] = ilast ;
		    if (ilast != EMPTY)
		    {
			Next [ilast] = inext ;
		    }
		    else
		    {
			/* i is at the head of the degree list */
			ASSERT (Degree [i] >= 0 && Degree [i] < n) ;
			Head [Degree [i]] = inext ;
		    }
		}
	    }
	}
	else
	{

	    /* ------------------------------------------------------------- */
	    /* construct the new element in empty space, Iw [pfree ...] */
	    /* ------------------------------------------------------------- */

	    p = Pe [me] ;
	    pme1 = pfree ;
	    slenme = Len [me] - elenme ;

	    for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++)
	    {

		if (knt1 > elenme)
		{
		    /* search the supervariables in me. */
		    e = me ;
		    pj = p ;
		    ln = slenme ;
		    AMD_DEBUG2 (("Search sv: "ID" "ID" "ID"\n", me,pj,ln)) ;
		}
		else
		{
		    /* search the elements in me. */
		    e = Iw [p++] ;
		    ASSERT (e >= 0 && e < n) ;
		    pj = Pe [e] ;
		    ln = Len [e] ;
		    AMD_DEBUG2 (("Search element e "ID" in me "ID"\n", e,me)) ;
		    ASSERT (Elen [e] < EMPTY && W [e] > 0 && pj >= 0) ;
		}
		ASSERT (ln >= 0 && (ln == 0 || (pj >= 0 && pj < iwlen))) ;

		/* ---------------------------------------------------------
		 * search for different supervariables and add them to the
		 * new list, compressing when necessary. this loop is
		 * executed once for each element in the list and once for
		 * all the supervariables in the list.
		 * --------------------------------------------------------- */

		for (knt2 = 1 ; knt2 <= ln ; knt2++)
		{
		    i = Iw [pj++] ;
		    ASSERT (i >= 0 && i < n && (i == me || Elen [i] >= EMPTY));
		    nvi = Nv [i] ;
		    AMD_DEBUG2 ((": "ID" "ID" "ID" "ID"\n",
				i, Elen [i], Nv [i], wflg)) ;

		    if (nvi > 0)
		    {

			/* ------------------------------------------------- */
			/* compress Iw, if necessary */
			/* ------------------------------------------------- */

			if (pfree >= iwlen)
			{

			    AMD_DEBUG1 (("GARBAGE COLLECTION\n")) ;

			    /* prepare for compressing Iw by adjusting pointers
			     * and lengths so that the lists being searched in
			     * the inner and outer loops contain only the
			     * remaining entries. */

			    Pe [me] = p ;
			    Len [me] -= knt1 ;
			    /* check if nothing left of supervariable me */
			    if (Len [me] == 0) Pe [me] = EMPTY ;
			    Pe [e] = pj ;
			    Len [e] = ln - knt2 ;
			    /* nothing left of element e */
			    if (Len [e] == 0) Pe [e] = EMPTY ;

			    ncmpa++ ;	/* one more garbage collection */

			    /* store first entry of each object in Pe */
			    /* FLIP the first entry in each object */
			    for (j = 0 ; j < n ; j++)
			    {
				pn = Pe [j] ;
				if (pn >= 0)
				{
				    ASSERT (pn >= 0 && pn < iwlen) ;
				    Pe [j] = Iw [pn] ;
				    Iw [pn] = FLIP (j) ;
				}
			    }

			    /* psrc/pdst point to source/destination */
			    psrc = 0 ;
			    pdst = 0 ;
			    pend = pme1 - 1 ;

			    while (psrc <= pend)
			    {
				/* search for next FLIP'd entry */
				j = FLIP (Iw [psrc++]) ;
				if (j >= 0)
				{
				    AMD_DEBUG2 (("Got object j: "ID"\n", j)) ;
				    Iw [pdst] = Pe [j] ;
				    Pe [j] = pdst++ ;
				    lenj = Len [j] ;
				    /* copy from source to destination */
				    for (knt3 = 0 ; knt3 <= lenj - 2 ; knt3++)
				    {
					Iw [pdst++] = Iw [psrc++] ;
				    }
				}
			    }

			    /* move the new partially-constructed element */
			    p1 = pdst ;
			    for (psrc = pme1 ; psrc <= pfree-1 ; psrc++)
			    {
				Iw [pdst++] = Iw [psrc] ;
			    }
			    pme1 = p1 ;
			    pfree = pdst ;
			    pj = Pe [e] ;
			    p = Pe [me] ;

			}

			/* ------------------------------------------------- */
			/* i is a principal variable not yet placed in Lme */
			/* store i in new list */
			/* ------------------------------------------------- */

			/* flag i as being in Lme by negating Nv [i] */
			degme += nvi ;
			Nv [i] = -nvi ;
			Iw [pfree++] = i ;
			AMD_DEBUG2 (("     s: "ID"     nv "ID"\n", i, Nv [i]));

			/* ------------------------------------------------- */
			/* remove variable i from degree link list */
			/* ------------------------------------------------- */

			ilast = Last [i] ;
			inext = Next [i] ;
			ASSERT (ilast >= EMPTY && ilast < n) ;
			ASSERT (inext >= EMPTY && inext < n) ;
			if (inext != EMPTY) Last [inext] = ilast ;
			if (ilast != EMPTY)
			{
			    Next [ilast] = inext ;
			}
			else
			{
			    /* i is at the head of the degree list */
			    ASSERT (Degree [i] >= 0 && Degree [i] < n) ;
			    Head [Degree [i]] = inext ;
			}
		    }
		}

		if (e != me)
		{
		    /* set tree pointer and flag to indicate element e is
		     * absorbed into new element me (the parent of e is me) */
		    AMD_DEBUG1 ((" Element "ID" => "ID"\n", e, me)) ;
		    Pe [e] = FLIP (me) ;
		    W [e] = 0 ;
		}
	    }

	    pme2 = pfree - 1 ;
	}

	/* ----------------------------------------------------------------- */
	/* me has now been converted into an element in Iw [pme1..pme2] */
	/* ----------------------------------------------------------------- */

	/* degme holds the external degree of new element */
	Degree [me] = degme ;
	Pe [me] = pme1 ;
	Len [me] = pme2 - pme1 + 1 ;
	ASSERT (Pe [me] >= 0 && Pe [me] < iwlen) ;

	Elen [me] = FLIP (nvpiv + degme) ;
	/* FLIP (Elen (me)) is now the degree of pivot (including
	 * diagonal part). */

	/* ----------------------------------------------------------------- */
	/* make sure that wflg is not too large. */
	/* ----------------------------------------------------------------- */

	/* With the current value of wflg, wflg+n must not cause integer
	 * overflow */

	wflg = clear_flag (wflg, wbig, W, n) ;

/* ========================================================================= */
/* COMPUTE (W [e] - wflg) = |Le\Lme| FOR ALL ELEMENTS */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * Scan 1:  compute the external degrees of previous elements with
	 * respect to the current element.  That is:
	 *       (W [e] - wflg) = |Le \ Lme|
	 * for each element e that appears in any supervariable in Lme.  The
	 * notation Le refers to the pattern (list of supervariables) of a
	 * previous element e, where e is not yet absorbed, stored in
	 * Iw [Pe [e] + 1 ... Pe [e] + Len [e]].  The notation Lme
	 * refers to the pattern of the current element (stored in
	 * Iw [pme1..pme2]).   If aggressive absorption is enabled, and
	 * (W [e] - wflg) becomes zero, then the element e will be absorbed
	 * in Scan 2.
	 * ----------------------------------------------------------------- */

	AMD_DEBUG2 (("me: ")) ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    eln = Elen [i] ;
	    AMD_DEBUG3 ((""ID" Elen "ID": \n", i, eln)) ;
	    if (eln > 0)
	    {
		/* note that Nv [i] has been negated to denote i in Lme: */
		nvi = -Nv [i] ;
		ASSERT (nvi > 0 && Pe [i] >= 0 && Pe [i] < iwlen) ;
		wnvi = wflg - nvi ;
		for (p = Pe [i] ; p <= Pe [i] + eln - 1 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    AMD_DEBUG4 (("    e "ID" we "ID" ", e, we)) ;
		    if (we >= wflg)
		    {
			/* unabsorbed element e has been seen in this loop */
			AMD_DEBUG4 (("    unabsorbed, first time seen")) ;
			we -= nvi ;
		    }
		    else if (we != 0)
		    {
			/* e is an unabsorbed element */
			/* this is the first we have seen e in all of Scan 1 */
			AMD_DEBUG4 (("    unabsorbed")) ;
			we = Degree [e] + wnvi ;
		    }
		    AMD_DEBUG4 (("\n")) ;
		    W [e] = we ;
		}
	    }
	}
	AMD_DEBUG2 (("\n")) ;

/* ========================================================================= */
/* DEGREE UPDATE AND ELEMENT ABSORPTION */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * Scan 2:  for each i in Lme, sum up the degree of Lme (which is
	 * degme), plus the sum of the external degrees of each Le for the
	 * elements e appearing within i, plus the supervariables in i.
	 * Place i in hash list.
	 * ----------------------------------------------------------------- */

	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n && Nv [i] < 0 && Elen [i] >= 0) ;
	    AMD_DEBUG2 (("Updating: i "ID" "ID" "ID"\n", i, Elen[i], Len [i]));
	    p1 = Pe [i] ;
	    p2 = p1 + Elen [i] - 1 ;
	    pn = p1 ;
	    hash = 0 ;
	    deg = 0 ;
	    ASSERT (p1 >= 0 && p1 < iwlen && p2 >= -1 && p2 < iwlen) ;

	    /* ------------------------------------------------------------- */
	    /* scan the element list associated with supervariable i */
	    /* ------------------------------------------------------------- */

	    /* UMFPACK/MA38-style approximate degree: */
	    if (aggressive)
	    {
		for (p = p1 ; p <= p2 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    if (we != 0)
		    {
			/* e is an unabsorbed element */
			/* dext = | Le \ Lme | */
			dext = we - wflg ;
			if (dext > 0)
			{
			    deg += dext ;
			    Iw [pn++] = e ;
			    hash += e ;
			    AMD_DEBUG4 ((" e: "ID" hash = "ID"\n",e,hash)) ;
			}
			else
			{
			    /* external degree of e is zero, absorb e into me*/
			    AMD_DEBUG1 ((" Element "ID" =>"ID" (aggressive)\n",
				e, me)) ;
			    ASSERT (dext == 0) ;
			    Pe [e] = FLIP (me) ;
			    W [e] = 0 ;
			}
		    }
		}
	    }
	    else
	    {
		for (p = p1 ; p <= p2 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    if (we != 0)
		    {
			/* e is an unabsorbed element */
			dext = we - wflg ;
			ASSERT (dext >= 0) ;
			deg += dext ;
			Iw [pn++] = e ;
			hash += e ;
			AMD_DEBUG4 (("	e: "ID" hash = "ID"\n",e,hash)) ;
		    }
		}
	    }

	    /* count the number of elements in i (including me): */
	    Elen [i] = pn - p1 + 1 ;

	    /* ------------------------------------------------------------- */
	    /* scan the supervariables in the list associated with i */
	    /* ------------------------------------------------------------- */

	    /* The bulk of the AMD run time is typically spent in this loop,
	     * particularly if the matrix has many dense rows that are not
	     * removed prior to ordering. */
	    p3 = pn ;
	    p4 = p1 + Len [i] ;
	    for (p = p2 + 1 ; p < p4 ; p++)
	    {
		j = Iw [p] ;
		ASSERT (j >= 0 && j < n) ;
		nvj = Nv [j] ;
		if (nvj > 0)
		{
		    /* j is unabsorbed, and not in Lme. */
		    /* add to degree and add to new list */
		    deg += nvj ;
		    Iw [pn++] = j ;
		    hash += j ;
		    AMD_DEBUG4 (("  s: "ID" hash "ID" Nv[j]= "ID"\n",
				j, hash, nvj)) ;
		}
	    }

	    /* ------------------------------------------------------------- */
	    /* update the degree and check for mass elimination */
	    /* ------------------------------------------------------------- */

	    /* with aggressive absorption, deg==0 is identical to the
	     * Elen [i] == 1 && p3 == pn test, below. */
	    ASSERT (IMPLIES (aggressive, (deg==0) == (Elen[i]==1 && p3==pn))) ;

	    if (Elen [i] == 1 && p3 == pn)
	    {

		/* --------------------------------------------------------- */
		/* mass elimination */
		/* --------------------------------------------------------- */

		/* There is nothing left of this node except for an edge to
		 * the current pivot element.  Elen [i] is 1, and there are
		 * no variables adjacent to node i.  Absorb i into the
		 * current pivot element, me.  Note that if there are two or
		 * more mass eliminations, fillin due to mass elimination is
		 * possible within the nvpiv-by-nvpiv pivot block.  It is this
		 * step that causes AMD's analysis to be an upper bound.
		 *
		 * The reason is that the selected pivot has a lower
		 * approximate degree than the true degree of the two mass
		 * eliminated nodes.  There is no edge between the two mass
		 * eliminated nodes.  They are merged with the current pivot
		 * anyway.
		 *
		 * No fillin occurs in the Schur complement, in any case,
		 * and this effect does not decrease the quality of the
		 * ordering itself, just the quality of the nonzero and
		 * flop count analysis.  It also means that the post-ordering
		 * is not an exact elimination tree post-ordering. */

		AMD_DEBUG1 (("  MASS i "ID" => parent e "ID"\n", i, me)) ;
		Pe [i] = FLIP (me) ;
		nvi = -Nv [i] ;
		degme -= nvi ;
		nvpiv += nvi ;
		nel += nvi ;
		Nv [i] = 0 ;
		Elen [i] = EMPTY ;

	    }
	    else
	    {

		/* --------------------------------------------------------- */
		/* update the upper-bound degree of i */
		/* --------------------------------------------------------- */

		/* the following degree does not yet include the size
		 * of the current element, which is added later: */

		Degree [i] = MIN (Degree [i], deg) ;

		/* --------------------------------------------------------- */
		/* add me to the list for i */
		/* --------------------------------------------------------- */

		/* move first supervariable to end of list */
		Iw [pn] = Iw [p3] ;
		/* move first element to end of element part of list */
		Iw [p3] = Iw [p1] ;
		/* add new element, me, to front of list. */
		Iw [p1] = me ;
		/* store the new length of the list in Len [i] */
		Len [i] = pn - p1 + 1 ;

		/* --------------------------------------------------------- */
		/* place in hash bucket.  Save hash key of i in Last [i]. */
		/* --------------------------------------------------------- */

		/* NOTE: this can fail if hash is negative, because the ANSI C
		 * standard does not define a % b when a and/or b are negative.
		 * That's why hash is defined as an unsigned Int, to avoid this
		 * problem. */
		hash = hash % n ;
		ASSERT (((Int) hash) >= 0 && ((Int) hash) < n) ;

		/* if the Hhead array is not used: */
		j = Head [hash] ;
		if (j <= EMPTY)
		{
		    /* degree list is empty, hash head is FLIP (j) */
		    Next [i] = FLIP (j) ;
		    Head [hash] = FLIP (i) ;
		}
		else
		{
		    /* degree list is not empty, use Last [Head [hash]] as
		     * hash head. */
		    Next [i] = Last [j] ;
		    Last [j] = i ;
		}

		/* if a separate Hhead array is used: *
		Next [i] = Hhead [hash] ;
		Hhead [hash] = i ;
		*/

		Last [i] = hash ;
	    }
	}

	Degree [me] = degme ;

	/* ----------------------------------------------------------------- */
	/* Clear the counter array, W [...], by incrementing wflg. */
	/* ----------------------------------------------------------------- */

	/* make sure that wflg+n does not cause integer overflow */
	lemax =  MAX (lemax, degme) ;
	wflg += lemax ;
	wflg = clear_flag (wflg, wbig, W, n) ;
	/*  at this point, W [0..n-1] < wflg holds */

/* ========================================================================= */
/* SUPERVARIABLE DETECTION */
/* ========================================================================= */

	AMD_DEBUG1 (("Detecting supervariables:\n")) ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    AMD_DEBUG2 (("Consider i "ID" nv "ID"\n", i, Nv [i])) ;
	    if (Nv [i] < 0)
	    {
		/* i is a principal variable in Lme */

		/* ---------------------------------------------------------
		 * examine all hash buckets with 2 or more variables.  We do
		 * this by examing all unique hash keys for supervariables in
		 * the pattern Lme of the current element, me
		 * --------------------------------------------------------- */

		/* let i = head of hash bucket, and empty the hash bucket */
		ASSERT (Last [i] >= 0 && Last [i] < n) ;
		hash = Last [i] ;

		/* if Hhead array is not used: */
		j = Head [hash] ;
		if (j == EMPTY)
		{
		    /* hash bucket and degree list are both empty */
		    i = EMPTY ;
		}
		else if (j < EMPTY)
		{
		    /* degree list is empty */
		    i = FLIP (j) ;
		    Head [hash] = EMPTY ;
		}
		else
		{
		    /* degree list is not empty, restore Last [j] of head j */
		    i = Last [j] ;
		    Last [j] = EMPTY ;
		}

		/* if separate Hhead array is used: *
		i = Hhead [hash] ;
		Hhead [hash] = EMPTY ;
		*/

		ASSERT (i >= EMPTY && i < n) ;
		AMD_DEBUG2 (("----i "ID" hash "ID"\n", i, hash)) ;

		while (i != EMPTY && Next [i] != EMPTY)
		{

		    /* -----------------------------------------------------
		     * this bucket has one or more variables following i.
		     * scan all of them to see if i can absorb any entries
		     * that follow i in hash bucket.  Scatter i into w.
		     * ----------------------------------------------------- */

		    ln = Len [i] ;
		    eln = Elen [i] ;
		    ASSERT (ln >= 0 && eln >= 0) ;
		    ASSERT (Pe [i] >= 0 && Pe [i] < iwlen) ;
		    /* do not flag the first element in the list (me) */
		    for (p = Pe [i] + 1 ; p <= Pe [i] + ln - 1 ; p++)
		    {
			ASSERT (Iw [p] >= 0 && Iw [p] < n) ;
			W [Iw [p]] = wflg ;
		    }

		    /* ----------------------------------------------------- */
		    /* scan every other entry j following i in bucket */
		    /* ----------------------------------------------------- */

		    jlast = i ;
		    j = Next [i] ;
		    ASSERT (j >= EMPTY && j < n) ;

		    while (j != EMPTY)
		    {
			/* ------------------------------------------------- */
			/* check if j and i have identical nonzero pattern */
			/* ------------------------------------------------- */

			AMD_DEBUG3 (("compare i "ID" and j "ID"\n", i,j)) ;

			/* check if i and j have the same Len and Elen */
			ASSERT (Len [j] >= 0 && Elen [j] >= 0) ;
			ASSERT (Pe [j] >= 0 && Pe [j] < iwlen) ;
			ok = (Len [j] == ln) && (Elen [j] == eln) ;
			/* skip the first element in the list (me) */
			for (p = Pe [j] + 1 ; ok && p <= Pe [j] + ln - 1 ; p++)
			{
			    ASSERT (Iw [p] >= 0 && Iw [p] < n) ;
			    if (W [Iw [p]] != wflg) ok = 0 ;
			}
			if (ok)
			{
			    /* --------------------------------------------- */
			    /* found it!  j can be absorbed into i */
			    /* --------------------------------------------- */

			    AMD_DEBUG1 (("found it! j "ID" => i "ID"\n", j,i));
			    Pe [j] = FLIP (i) ;
			    /* both Nv [i] and Nv [j] are negated since they */
			    /* are in Lme, and the absolute values of each */
			    /* are the number of variables in i and j: */
			    Nv [i] += Nv [j] ;
			    Nv [j] = 0 ;
			    Elen [j] = EMPTY ;
			    /* delete j from hash bucket */
			    ASSERT (j != Next [j]) ;
			    j = Next [j] ;
			    Next [jlast] = j ;

			}
			else
			{
			    /* j cannot be absorbed into i */
			    jlast = j ;
			    ASSERT (j != Next [j]) ;
			    j = Next [j] ;
			}
			ASSERT (j >= EMPTY && j < n) ;
		    }

		    /* -----------------------------------------------------
		     * no more variables can be absorbed into i
		     * go to next i in bucket and clear flag array
		     * ----------------------------------------------------- */

		    wflg++ ;
		    i = Next [i] ;
		    ASSERT (i >= EMPTY && i < n) ;

		}
	    }
	}
	AMD_DEBUG2 (("detect done\n")) ;

/* ========================================================================= */
/* RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVARIABLES FROM ELEMENT */
/* ========================================================================= */

	p = pme1 ;
	nleft = n - nel ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    nvi = -Nv [i] ;
	    AMD_DEBUG3 (("Restore i "ID" "ID"\n", i, nvi)) ;
	    if (nvi > 0)
	    {
		/* i is a principal variable in Lme */
		/* restore Nv [i] to signify that i is principal */
		Nv [i] = nvi ;

		/* --------------------------------------------------------- */
		/* compute the external degree (add size of current element) */
		/* --------------------------------------------------------- */

		deg = Degree [i] + degme - nvi ;
		deg = MIN (deg, nleft - nvi) ;
		ASSERT (IMPLIES (aggressive, deg > 0) && deg >= 0 && deg < n) ;

		/* --------------------------------------------------------- */
		/* place the supervariable at the head of the degree list */
		/* --------------------------------------------------------- */

		inext = Head [deg] ;
		ASSERT (inext >= EMPTY && inext < n) ;
		if (inext != EMPTY) Last [inext] = i ;
		Next [i] = inext ;
		Last [i] = EMPTY ;
		Head [deg] = i ;

		/* --------------------------------------------------------- */
		/* save the new degree, and find the minimum degree */
		/* --------------------------------------------------------- */

		mindeg = MIN (mindeg, deg) ;
		Degree [i] = deg ;

		/* --------------------------------------------------------- */
		/* place the supervariable in the element pattern */
		/* --------------------------------------------------------- */

		Iw [p++] = i ;

	    }
	}
	AMD_DEBUG2 (("restore done\n")) ;

/* ========================================================================= */
/* FINALIZE THE NEW ELEMENT */
/* ========================================================================= */

	AMD_DEBUG2 (("ME = "ID" DONE\n", me)) ;
	Nv [me] = nvpiv ;
	/* save the length of the list for the new element me */
	Len [me] = p - pme1 ;
	if (Len [me] == 0)
	{
	    /* there is nothing left of the current pivot element */
	    /* it is a root of the assembly tree */
	    Pe [me] = EMPTY ;
	    W [me] = 0 ;
	}
	if (elenme != 0)
	{
	    /* element was not constructed in place: deallocate part of */
	    /* it since newly nonprincipal variables may have been removed */
	    pfree = p ;
	}

	/* The new element has nvpiv pivots and the size of the contribution
	 * block for a multifrontal method is degme-by-degme, not including
	 * the "dense" rows/columns.  If the "dense" rows/columns are included,
	 * the frontal matrix is no larger than
	 * (degme+ndense)-by-(degme+ndense).
	 */

	if (Info != (double *) NULL)
	{
	    f = nvpiv ;
	    r = degme + ndense ;
	    dmax = MAX (dmax, f + r) ;

	    /* number of nonzeros in L (excluding the diagonal) */
	    lnzme = f*r + (f-1)*f/2 ;
	    lnz += lnzme ;

	    /* number of divide operations for LDL' and for LU */
	    ndiv += lnzme ;

	    /* number of multiply-subtract pairs for LU */
	    s = f*r*r + r*(f-1)*f + (f-1)*f*(2*f-1)/6 ;
	    nms_lu += s ;

	    /* number of multiply-subtract pairs for LDL' */
	    nms_ldl += (s + lnzme)/2 ;
	}

    }

/* ========================================================================= */
/* DONE SELECTING PIVOTS */
/* ========================================================================= */

    if (Info != (double *) NULL)
    {

	/* count the work to factorize the ndense-by-ndense submatrix */
	f = ndense ;
	dmax = MAX (dmax, (double) ndense) ;

	/* number of nonzeros in L (excluding the diagonal) */
	lnzme = (f-1)*f/2 ;
	lnz += lnzme ;

	/* number of divide operations for LDL' and for LU */
	ndiv += lnzme ;

	/* number of multiply-subtract pairs for LU */
	s = (f-1)*f*(2*f-1)/6 ;
	nms_lu += s ;

	/* number of multiply-subtract pairs for LDL' */
	nms_ldl += (s + lnzme)/2 ;

	/* number of nz's in L (excl. diagonal) */
	Info [AMD_LNZ] = lnz ;

	/* number of divide ops for LU and LDL' */
	Info [AMD_NDIV] = ndiv ;

	/* number of multiply-subtract pairs for LDL' */
	Info [AMD_NMULTSUBS_LDL] = nms_ldl ;

	/* number of multiply-subtract pairs for LU */
	Info [AMD_NMULTSUBS_LU] = nms_lu ;

	/* number of "dense" rows/columns */
	Info [AMD_NDENSE] = ndense ;

	/* largest front is dmax-by-dmax */
	Info [AMD_DMAX] = dmax ;

	/* number of garbage collections in AMD */
	Info [AMD_NCMPA] = ncmpa ;

	/* successful ordering */
	Info [AMD_STATUS] = AMD_OK ;
    }

/* ========================================================================= */
/* POST-ORDERING */
/* ========================================================================= */

/* -------------------------------------------------------------------------
 * Variables at this point:
 *
 * Pe: holds the elimination tree.  The parent of j is FLIP (Pe [j]),
 *	or EMPTY if j is a root.  The tree holds both elements and
 *	non-principal (unordered) variables absorbed into them.
 *	Dense variables are non-principal and unordered.
 *
 * Elen: holds the size of each element, including the diagonal part.
 *	FLIP (Elen [e]) > 0 if e is an element.  For unordered
 *	variables i, Elen [i] is EMPTY.
 *
 * Nv: Nv [e] > 0 is the number of pivots represented by the element e.
 *	For unordered variables i, Nv [i] is zero.
 *
 * Contents no longer needed:
 *	W, Iw, Len, Degree, Head, Next, Last.
 *
 * The matrix itself has been destroyed.
 *
 * n: the size of the matrix.
 * No other scalars needed (pfree, iwlen, etc.)
 * ------------------------------------------------------------------------- */

    /* restore Pe */
    for (i = 0 ; i < n ; i++)
    {
	Pe [i] = FLIP (Pe [i]) ;
    }

    /* restore Elen, for output information, and for postordering */
    for (i = 0 ; i < n ; i++)
    {
	Elen [i] = FLIP (Elen [i]) ;
    }

/* Now the parent of j is Pe [j], or EMPTY if j is a root.  Elen [e] > 0
 * is the size of element e.  Elen [i] is EMPTY for unordered variable i. */

/* ========================================================================= */
/* compress the paths of the variables */
/* ========================================================================= */

    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] == 0)
	{

	    /* -------------------------------------------------------------
	     * i is an un-ordered row.  Traverse the tree from i until
	     * reaching an element, e.  The element, e, was the principal
	     * supervariable of i and all nodes in the path from i to when e
	     * was selected as pivot.
	     * ------------------------------------------------------------- */

	    AMD_DEBUG1 (("Path compression, i unordered: "ID"\n", i)) ;
	    j = Pe [i] ;
	    ASSERT (j >= EMPTY && j < n) ;
	    AMD_DEBUG3 (("	j: "ID"\n", j)) ;
	    if (j == EMPTY)
	    {
		/* Skip a dense variable.  It has no parent. */
		AMD_DEBUG3 (("      i is a dense variable\n")) ;
		continue ;
	    }

	    /* while (j is a variable) */
	    while (Nv [j] == 0)
	    {
		AMD_DEBUG3 (("		j : "ID"\n", j)) ;
		j = Pe [j] ;
		AMD_DEBUG3 (("		j:: "ID"\n", j)) ;
		ASSERT (j >= 0 && j < n) ;
	    }
	    /* got to an element e */
	    e = j ;
	    AMD_DEBUG3 (("got to e: "ID"\n", e)) ;

	    /* -------------------------------------------------------------
	     * traverse the path again from i to e, and compress the path
	     * (all nodes point to e).  Path compression allows this code to
	     * compute in O(n) time.
	     * ------------------------------------------------------------- */

	    j = i ;
	    /* while (j is a variable) */
	    while (Nv [j] == 0)
	    {
		jnext = Pe [j] ;
		AMD_DEBUG3 (("j "ID" jnext "ID"\n", j, jnext)) ;
		Pe [j] = e ;
		j = jnext ;
		ASSERT (j >= 0 && j < n) ;
	    }
	}
    }

/* ========================================================================= */
/* postorder the assembly tree */
/* ========================================================================= */

    AMD_postorder (n, Pe, Nv, Elen,
	W,			/* output order */
	Head, Next, Last) ;	/* workspace */

/* ========================================================================= */
/* compute output permutation and inverse permutation */
/* ========================================================================= */

    /* W [e] = k means that element e is the kth element in the new
     * order.  e is in the range 0 to n-1, and k is in the range 0 to
     * the number of elements.  Use Head for inverse order. */

    for (k = 0 ; k < n ; k++)
    {
	Head [k] = EMPTY ;
	Next [k] = EMPTY ;
    }
    for (e = 0 ; e < n ; e++)
    {
	k = W [e] ;
	ASSERT ((k == EMPTY) == (Nv [e] == 0)) ;
	if (k != EMPTY)
	{
	    ASSERT (k >= 0 && k < n) ;
	    Head [k] = e ;
	}
    }

    /* construct output inverse permutation in Next,
     * and permutation in Last */
    nel = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	e = Head [k] ;
	if (e == EMPTY) break ;
	ASSERT (e >= 0 && e < n && Nv [e] > 0) ;
	Next [e] = nel ;
	nel += Nv [e] ;
    }
    ASSERT (nel == n - ndense) ;

    /* order non-principal variables (dense, & those merged into supervar's) */
    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] == 0)
	{
	    e = Pe [i] ;
	    ASSERT (e >= EMPTY && e < n) ;
	    if (e != EMPTY)
	    {
		/* This is an unordered variable that was merged
		 * into element e via supernode detection or mass
		 * elimination of i when e became the pivot element.
		 * Place i in order just before e. */
		ASSERT (Next [i] == EMPTY && Nv [e] > 0) ;
		Next [i] = Next [e] ;
		Next [e]++ ;
	    }
	    else
	    {
		/* This is a dense unordered variable, with no parent.
		 * Place it last in the output order. */
		Next [i] = nel++ ;
	    }
	}
    }
    ASSERT (nel == n) ;

    AMD_DEBUG2 (("\n\nPerm:\n")) ;
    for (i = 0 ; i < n ; i++)
    {
	k = Next [i] ;
	ASSERT (k >= 0 && k < n) ;
	Last [k] = i ;
	AMD_DEBUG2 (("   perm ["ID"] = "ID"\n", k, i)) ;
    }
}

/* ========================================================================= */
/* === AMD_aat ============================================================= */
/* ========================================================================= */

/* AMD_aat:  compute the symmetry of the pattern of A, and count the number of
 * nonzeros each column of A+A' (excluding the diagonal).  Assumes the input
 * matrix has no errors, with sorted columns and no duplicates
 * (AMD_valid (n, n, Ap, Ai) must be AMD_OK, but this condition is not
 * checked).
 */

GLOBAL size_t AMD_aat	/* returns nz in A+A' */
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Len [ ],	/* Len [j]: length of column j of A+A', excl diagonal*/
    Int Tp [ ],		/* workspace of size n */
    double Info [ ]
)
{
    Int p1, p2, p, i, j, pj, pj2, k, nzdiag, nzboth, nz ;
    double sym ;
    size_t nzaat ;

    if (Info != (double *) NULL)
    {
	/* clear the Info array, if it exists */
	for (i = 0 ; i < AMD_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
	Info [AMD_STATUS] = AMD_OK ;
    }

    for (k = 0 ; k < n ; k++)
    {
	Len [k] = 0 ;
    }

    nzdiag = 0 ;
    nzboth = 0 ;
    nz = Ap [n] ;

    for (k = 0 ; k < n ; k++)
    {
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;
	AMD_DEBUG2 (("\nAAT Column: "ID" p1: "ID" p2: "ID"\n", k, p1, p2)) ;

	/* construct A+A' */
	for (p = p1 ; p < p2 ; )
	{
	    /* scan the upper triangular part of A */
	    j = Ai [p] ;
	    if (j < k)
	    {
		/* entry A (j,k) is in the strictly upper triangular part,
		 * add both A (j,k) and A (k,j) to the matrix A+A' */
		Len [j]++ ;
		Len [k]++ ;
		AMD_DEBUG3 (("    upper ("ID","ID") ("ID","ID")\n", j,k, k,j));
		p++ ;
	    }
	    else if (j == k)
	    {
		/* skip the diagonal */
		p++ ;
		nzdiag++ ;
		break ;
	    }
	    else /* j > k */
	    {
		/* first entry below the diagonal */
		break ;
	    }
	    /* scan lower triangular part of A, in column j until reaching
	     * row k.  Start where last scan left off. */
	    ASSERT (Tp [j] != EMPTY) ;
	    ASSERT (Ap [j] <= Tp [j] && Tp [j] <= Ap [j+1]) ;
	    pj2 = Ap [j+1] ;
	    for (pj = Tp [j] ; pj < pj2 ; )
	    {
		i = Ai [pj] ;
		if (i < k)
		{
		    /* A (i,j) is only in the lower part, not in upper.
		     * add both A (i,j) and A (j,i) to the matrix A+A' */
		    Len [i]++ ;
		    Len [j]++ ;
		    AMD_DEBUG3 (("    lower ("ID","ID") ("ID","ID")\n",
			i,j, j,i)) ;
		    pj++ ;
		}
		else if (i == k)
		{
		    /* entry A (k,j) in lower part and A (j,k) in upper */
		    pj++ ;
		    nzboth++ ;
		    break ;
		}
		else /* i > k */
		{
		    /* consider this entry later, when k advances to i */
		    break ;
		}
	    }
	    Tp [j] = pj ;
	}
	/* Tp [k] points to the entry just below the diagonal in column k */
	Tp [k] = p ;
    }

    /* clean up, for remaining mismatched entries */
    for (j = 0 ; j < n ; j++)
    {
	for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
	{
	    i = Ai [pj] ;
	    /* A (i,j) is only in the lower part, not in upper.
	     * add both A (i,j) and A (j,i) to the matrix A+A' */
	    Len [i]++ ;
	    Len [j]++ ;
	    AMD_DEBUG3 (("    lower cleanup ("ID","ID") ("ID","ID")\n",
		i,j, j,i)) ;
	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the symmetry of the nonzero pattern of A */
    /* --------------------------------------------------------------------- */

    /* Given a matrix A, the symmetry of A is:
     *	B = tril (spones (A), -1) + triu (spones (A), 1) ;
     *  sym = nnz (B & B') / nnz (B) ;
     *  or 1 if nnz (B) is zero.
     */

    if (nz == nzdiag)
    {
	sym = 1 ;
    }
    else
    {
	sym = (2 * (double) nzboth) / ((double) (nz - nzdiag)) ;
    }

    nzaat = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	nzaat += Len [k] ;
    }

    AMD_DEBUG1 (("AMD nz in A+A', excluding diagonal (nzaat) = %g\n",
	(double) nzaat)) ;
    AMD_DEBUG1 (("   nzboth: "ID" nz: "ID" nzdiag: "ID" symmetry: %g\n",
		nzboth, nz, nzdiag, sym)) ;

    if (Info != (double *) NULL)
    {
	Info [AMD_STATUS] = AMD_OK ;
	Info [AMD_N] = n ;
	Info [AMD_NZ] = nz ;
	Info [AMD_SYMMETRY] = sym ;	    /* symmetry of pattern of A */
	Info [AMD_NZDIAG] = nzdiag ;	    /* nonzeros on diagonal of A */
	Info [AMD_NZ_A_PLUS_AT] = nzaat ;   /* nonzeros in A+A' */
    }

    return (nzaat) ;
}

/* ========================================================================= */
/* === AMD_control ========================================================= */
/* ========================================================================= */

/* User-callable.  Prints the control parameters for AMD.  See amd.h
 * for details.  If the Control array is not present, the defaults are
 * printed instead.
 */

GLOBAL void AMD_control
(
    double Control [ ]
)
{
    double alpha ;
    Int aggressive ;

    if (Control != (double *) NULL)
    {
	alpha = Control [AMD_DENSE] ;
	aggressive = Control [AMD_AGGRESSIVE] != 0 ;
    }
    else
    {
	alpha = AMD_DEFAULT_DENSE ;
	aggressive = AMD_DEFAULT_AGGRESSIVE ;
    }

    SPARSE_PRINTF ((
        "\nAMD version %d.%d.%d, %s: approximate minimum degree ordering\n"
	"    dense row parameter: %g\n", AMD_MAIN_VERSION, AMD_SUB_VERSION,
	AMD_SUBSUB_VERSION, AMD_DATE, alpha)) ;

    if (alpha < 0)
    {
	SPARSE_PRINTF (("    no rows treated as dense\n")) ;
    }
    else
    {
	SPARSE_PRINTF ((
	"    (rows with more than max (%g * sqrt (n), 16) entries are\n"
	"    considered \"dense\", and placed last in output permutation)\n",
	alpha)) ;
    }

    if (aggressive)
    {
	SPARSE_PRINTF (("    aggressive absorption:  yes\n")) ;
    }
    else
    {
	SPARSE_PRINTF (("    aggressive absorption:  no\n")) ;
    }

    SPARSE_PRINTF (("    size of AMD integer: %d\n\n", sizeof (Int))) ;
}

/* ========================================================================= */
/* === AMD_defaults ======================================================== */
/* ========================================================================= */

/* User-callable.  Sets default control parameters for AMD.  See amd.h
 * for details.
 */

/* ========================================================================= */
/* === AMD defaults ======================================================== */
/* ========================================================================= */

GLOBAL void AMD_defaults
(
    double Control [ ]
)
{
    Int i ;

    if (Control != (double *) NULL)
    {
	for (i = 0 ; i < AMD_CONTROL ; i++)
	{
	    Control [i] = 0 ;
	}
	Control [AMD_DENSE] = AMD_DEFAULT_DENSE ;
	Control [AMD_AGGRESSIVE] = AMD_DEFAULT_AGGRESSIVE ;
    }
}

/* ========================================================================= */
/* === AMD_info ============================================================ */
/* ========================================================================= */

/* User-callable.  Prints the output statistics for AMD.  See amd.h
 * for details.  If the Info array is not present, nothing is printed.
 */

#define PRI(format,x) { if (x >= 0) { SPARSE_PRINTF ((format, x)) ; }}

GLOBAL void AMD_info
(
    double Info [ ]
)
{
    double n, ndiv, nmultsubs_ldl, nmultsubs_lu, lnz, lnzd ;

    SPARSE_PRINTF (("\nAMD version %d.%d.%d, %s, results:\n",
	AMD_MAIN_VERSION, AMD_SUB_VERSION, AMD_SUBSUB_VERSION, AMD_DATE)) ;

    if (!Info)
    {
	return ;
    }

    n = Info [AMD_N] ;
    ndiv = Info [AMD_NDIV] ;
    nmultsubs_ldl = Info [AMD_NMULTSUBS_LDL] ;
    nmultsubs_lu = Info [AMD_NMULTSUBS_LU] ;
    lnz = Info [AMD_LNZ] ;
    lnzd = (n >= 0 && lnz >= 0) ? (n + lnz) : (-1) ;

    /* AMD return status */
    SPARSE_PRINTF (("    status: ")) ;
    if (Info [AMD_STATUS] == AMD_OK)
    {
	SPARSE_PRINTF (("OK\n")) ;
    }
    else if (Info [AMD_STATUS] == AMD_OUT_OF_MEMORY)
    {
	SPARSE_PRINTF (("out of memory\n")) ;
    }
    else if (Info [AMD_STATUS] == AMD_INVALID)
    {
	SPARSE_PRINTF (("invalid matrix\n")) ;
    }
    else if (Info [AMD_STATUS] == AMD_OK_BUT_JUMBLED)
    {
	SPARSE_PRINTF (("OK, but jumbled\n")) ;
    }
    else
    {
	SPARSE_PRINTF (("unknown\n")) ;
    }

    /* statistics about the input matrix */
    PRI ("    n, dimension of A:                                  %.20g\n", n);
    PRI ("    nz, number of nonzeros in A:                        %.20g\n",
	Info [AMD_NZ]) ;
    PRI ("    symmetry of A:                                      %.4f\n",
	Info [AMD_SYMMETRY]) ;
    PRI ("    number of nonzeros on diagonal:                     %.20g\n",
	Info [AMD_NZDIAG]) ;
    PRI ("    nonzeros in pattern of A+A' (excl. diagonal):       %.20g\n",
	Info [AMD_NZ_A_PLUS_AT]) ;
    PRI ("    # dense rows/columns of A+A':                       %.20g\n",
	Info [AMD_NDENSE]) ;

    /* statistics about AMD's behavior  */
    PRI ("    memory used, in bytes:                              %.20g\n",
	Info [AMD_MEMORY]) ;
    PRI ("    # of memory compactions:                            %.20g\n",
	Info [AMD_NCMPA]) ;

    /* statistics about the ordering quality */
    SPARSE_PRINTF (("\n"
	"    The following approximate statistics are for a subsequent\n"
	"    factorization of A(P,P) + A(P,P)'.  They are slight upper\n"
	"    bounds if there are no dense rows/columns in A+A', and become\n"
	"    looser if dense rows/columns exist.\n\n")) ;

    PRI ("    nonzeros in L (excluding diagonal):                 %.20g\n",
	lnz) ;
    PRI ("    nonzeros in L (including diagonal):                 %.20g\n",
	lnzd) ;
    PRI ("    # divide operations for LDL' or LU:                 %.20g\n",
	ndiv) ;
    PRI ("    # multiply-subtract operations for LDL':            %.20g\n",
	nmultsubs_ldl) ;
    PRI ("    # multiply-subtract operations for LU:              %.20g\n",
	nmultsubs_lu) ;
    PRI ("    max nz. in any column of L (incl. diagonal):        %.20g\n",
	Info [AMD_DMAX]) ;

    /* total flop counts for various factorizations */

    if (n >= 0 && ndiv >= 0 && nmultsubs_ldl >= 0 && nmultsubs_lu >= 0)
    {
	SPARSE_PRINTF (("\n"
	"    chol flop count for real A, sqrt counted as 1 flop: %.20g\n"
	"    LDL' flop count for real A:                         %.20g\n"
	"    LDL' flop count for complex A:                      %.20g\n"
	"    LU flop count for real A (with no pivoting):        %.20g\n"
	"    LU flop count for complex A (with no pivoting):     %.20g\n\n",
	n + ndiv + 2*nmultsubs_ldl,
	    ndiv + 2*nmultsubs_ldl,
	  9*ndiv + 8*nmultsubs_ldl,
	    ndiv + 2*nmultsubs_lu,
	  9*ndiv + 8*nmultsubs_lu)) ;
    }
}

/* ========================================================================= */
/* === AMD_order =========================================================== */
/* ========================================================================= */

/* User-callable AMD minimum degree ordering routine.  See amd.h for
 * documentation.
 */

GLOBAL Int AMD_order
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int P [ ],
    double Control [ ],
    double Info [ ]
)
{
    Int *Len, *S, nz, i, *Pinv, info, status, *Rp, *Ri, *Cp, *Ci, ok ;
    size_t nzaat, slen ;
    double mem = 0 ;

    /* clear the Info array, if it exists */
    info = Info != (double *) NULL ;
    if (info)
    {
	for (i = 0 ; i < AMD_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
	Info [AMD_N] = n ;
	Info [AMD_STATUS] = AMD_OK ;
    }

    /* make sure inputs exist and n is >= 0 */
    if (Ai == (Int *) NULL || Ap == (Int *) NULL || P == (Int *) NULL || n < 0)
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;	    /* arguments are invalid */
    }

    if (n == 0)
    {
	return (AMD_OK) ;	    /* n is 0 so there's nothing to do */
    }

    nz = Ap [n] ;
    if (info)
    {
	Info [AMD_NZ] = nz ;
    }
    if (nz < 0)
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;
    }

    /* check if n or nz will cause size_t overflow */
    if (((size_t) n) >= SIZE_T_MAX / sizeof (Int)
     || ((size_t) nz) >= SIZE_T_MAX / sizeof (Int))
    {
	if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;	    /* problem too large */
    }

    /* check the input matrix:	AMD_OK, AMD_INVALID, or AMD_OK_BUT_JUMBLED */
    status = AMD_valid (n, n, Ap, Ai) ;

    if (status == AMD_INVALID)
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;	    /* matrix is invalid */
    }

    /* allocate two size-n integer workspaces */
    Len  = SparseBase_malloc (n, sizeof (Int)) ;
    Pinv = SparseBase_malloc (n, sizeof (Int)) ;
    mem += n ;
    mem += n ;
    if (!Len || !Pinv)
    {
	/* :: out of memory :: */
	SparseBase_free (Len) ;
	SparseBase_free (Pinv) ;
	if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;
    }

    if (status == AMD_OK_BUT_JUMBLED)
    {
	/* sort the input matrix and remove duplicate entries */
	AMD_DEBUG1 (("Matrix is jumbled\n")) ;
	Rp = SparseBase_malloc (n+1, sizeof (Int)) ;
	Ri = SparseBase_malloc (nz,  sizeof (Int)) ;
	mem += (n+1) ;
	mem += MAX (nz,1) ;
	if (!Rp || !Ri)
	{
	    /* :: out of memory :: */
	    SparseBase_free (Rp) ;
	    SparseBase_free (Ri) ;
	    SparseBase_free (Len) ;
	    SparseBase_free (Pinv) ;
	    if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	    return (AMD_OUT_OF_MEMORY) ;
	}
	/* use Len and Pinv as workspace to create R = A' */
	AMD_preprocess (n, Ap, Ai, Rp, Ri, Len, Pinv) ;
	Cp = Rp ;
	Ci = Ri ;
    }
    else
    {
	/* order the input matrix as-is.  No need to compute R = A' first */
	Rp = NULL ;
	Ri = NULL ;
	Cp = (Int *) Ap ;
	Ci = (Int *) Ai ;
    }

    /* --------------------------------------------------------------------- */
    /* determine the symmetry and count off-diagonal nonzeros in A+A' */
    /* --------------------------------------------------------------------- */

    nzaat = AMD_aat (n, Cp, Ci, Len, P, Info) ;
    AMD_DEBUG1 (("nzaat: %g\n", (double) nzaat)) ;
    ASSERT ((MAX (nz-n, 0) <= nzaat) && (nzaat <= 2 * (size_t) nz)) ;

    /* --------------------------------------------------------------------- */
    /* allocate workspace for matrix, elbow room, and 6 size-n vectors */
    /* --------------------------------------------------------------------- */

    S = NULL ;
    slen = nzaat ;			/* space for matrix */
    ok = ((slen + nzaat/5) >= slen) ;	/* check for size_t overflow */
    slen += nzaat/5 ;			/* add elbow room */
    for (i = 0 ; ok && i < 7 ; i++)
    {
	ok = ((slen + n) > slen) ;	/* check for size_t overflow */
	slen += n ;			/* size-n elbow room, 6 size-n work */
    }
    mem += slen ;
    ok = ok && (slen < SIZE_T_MAX / sizeof (Int)) ; /* check for overflow */
    ok = ok && (slen < Int_MAX) ;	/* S[i] for Int i must be OK */
    if (ok)
    {
	S = SparseBase_malloc (slen, sizeof (Int)) ;
    }
    AMD_DEBUG1 (("slen %g\n", (double) slen)) ;
    if (!S)
    {
	/* :: out of memory :: (or problem too large) */
	SparseBase_free (Rp) ;
	SparseBase_free (Ri) ;
	SparseBase_free (Len) ;
	SparseBase_free (Pinv) ;
	if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;
    }
    if (info)
    {
	/* memory usage, in bytes. */
	Info [AMD_MEMORY] = mem * sizeof (Int) ;
    }

    /* --------------------------------------------------------------------- */
    /* order the matrix */
    /* --------------------------------------------------------------------- */

    AMD_1 (n, Cp, Ci, P, Pinv, Len, slen, S, Control, Info) ;

    /* --------------------------------------------------------------------- */
    /* free the workspace */
    /* --------------------------------------------------------------------- */

    SparseBase_free (Rp) ;
    SparseBase_free (Ri) ;
    SparseBase_free (Len) ;
    SparseBase_free (Pinv) ;
    SparseBase_free (S) ;
    if (info) Info [AMD_STATUS] = status ;
    return (status) ;	    /* successful ordering */
}

/* ========================================================================= */
/* === AMD_post_tree ======================================================= */
/* ========================================================================= */

/* Post-ordering of a supernodal elimination tree.  */

GLOBAL Int AMD_post_tree
(
    Int root,			/* root of the tree */
    Int k,			/* start numbering at k */
    Int Child [ ],		/* input argument of size nn, undefined on
				 * output.  Child [i] is the head of a link
				 * list of all nodes that are children of node
				 * i in the tree. */
    const Int Sibling [ ],	/* input argument of size nn, not modified.
				 * If f is a node in the link list of the
				 * children of node i, then Sibling [f] is the
				 * next child of node i.
				 */
    Int Order [ ],		/* output order, of size nn.  Order [i] = k
				 * if node i is the kth node of the reordered
				 * tree. */
    Int Stack [ ]		/* workspace of size nn */
)
{
    Int f, head, h, i ;

#if 0
    /* --------------------------------------------------------------------- */
    /* recursive version (Stack [ ] is not used): */
    /* --------------------------------------------------------------------- */

    /* this is simple, but can caouse stack overflow if nn is large */
    i = root ;
    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
    {
	k = AMD_post_tree (f, k, Child, Sibling, Order, Stack, nn) ;
    }
    Order [i] = k++ ;
    return (k) ;
#endif

    /* --------------------------------------------------------------------- */
    /* non-recursive version, using an explicit stack */
    /* --------------------------------------------------------------------- */

    /* push root on the stack */
    head = 0 ;
    Stack [0] = root ;

    while (head >= 0)
    {
	/* get head of stack */
	ASSERT (head < nn) ;
	i = Stack [head] ;
	AMD_DEBUG1 (("head of stack "ID" \n", i)) ;
	ASSERT (i >= 0 && i < nn) ;

	if (Child [i] != EMPTY)
	{
	    /* the children of i are not yet ordered */
	    /* push each child onto the stack in reverse order */
	    /* so that small ones at the head of the list get popped first */
	    /* and the biggest one at the end of the list gets popped last */
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		head++ ;
		ASSERT (head < nn) ;
		ASSERT (f >= 0 && f < nn) ;
	    }
	    h = head ;
	    ASSERT (head < nn) ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (h > 0) ;
		Stack [h--] = f ;
		AMD_DEBUG1 (("push "ID" on stack\n", f)) ;
		ASSERT (f >= 0 && f < nn) ;
	    }
	    ASSERT (Stack [h] == i) ;

	    /* delete child list so that i gets ordered next time we see it */
	    Child [i] = EMPTY ;
	}
	else
	{
	    /* the children of i (if there were any) are already ordered */
	    /* remove i from the stack and order it.  Front i is kth front */
	    head-- ;
	    AMD_DEBUG1 (("pop "ID" order "ID"\n", i, k)) ;
	    Order [i] = k++ ;
	    ASSERT (k <= nn) ;
	}

    }
    return (k) ;
}

/* ========================================================================= */
/* === AMD_postorder ======================================================= */
/* ========================================================================= */

/* Perform a postordering (via depth-first search) of an assembly tree. */

GLOBAL void AMD_postorder
(
    /* inputs, not modified on output: */
    Int nn,		/* nodes are in the range 0..nn-1 */
    Int Parent [ ],	/* Parent [j] is the parent of j, or EMPTY if root */
    Int Nv [ ],		/* Nv [j] > 0 number of pivots represented by node j,
			 * or zero if j is not a node. */
    Int Fsize [ ],	/* Fsize [j]: size of node j */

    /* output, not defined on input: */
    Int Order [ ],	/* output post-order */

    /* workspaces of size nn: */
    Int Child [ ],
    Int Sibling [ ],
    Int Stack [ ]
)
{
    Int i, j, k, parent, frsize, f, fprev, maxfrsize, bigfprev, bigf, fnext ;

    for (j = 0 ; j < nn ; j++)
    {
	Child [j] = EMPTY ;
	Sibling [j] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* place the children in link lists - bigger elements tend to be last */
    /* --------------------------------------------------------------------- */

    for (j = nn-1 ; j >= 0 ; j--)
    {
	if (Nv [j] > 0)
	{
	    /* this is an element */
	    parent = Parent [j] ;
	    if (parent != EMPTY)
	    {
		/* place the element in link list of the children its parent */
		/* bigger elements will tend to be at the end of the list */
		Sibling [j] = Child [parent] ;
		Child [parent] = j ;
	    }
	}
    }

    /* --------------------------------------------------------------------- */
    /* place the largest child last in the list of children for each node */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < nn ; i++)
    {
	if (Nv [i] > 0 && Child [i] != EMPTY)
	{

	    /* find the biggest element in the child list */
	    fprev = EMPTY ;
	    maxfrsize = EMPTY ;
	    bigfprev = EMPTY ;
	    bigf = EMPTY ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (f >= 0 && f < nn) ;
		frsize = Fsize [f] ;
		if (frsize >= maxfrsize)
		{
		    /* this is the biggest seen so far */
		    maxfrsize = frsize ;
		    bigfprev = fprev ;
		    bigf = f ;
		}
		fprev = f ;
	    }
	    ASSERT (bigf != EMPTY) ;

	    fnext = Sibling [bigf] ;

	    AMD_DEBUG1 (("bigf "ID" maxfrsize "ID" bigfprev "ID" fnext "ID
		" fprev " ID"\n", bigf, maxfrsize, bigfprev, fnext, fprev)) ;

	    if (fnext != EMPTY)
	    {
		/* if fnext is EMPTY then bigf is already at the end of list */

		if (bigfprev == EMPTY)
		{
		    /* delete bigf from the element of the list */
		    Child [i] = fnext ;
		}
		else
		{
		    /* delete bigf from the middle of the list */
		    Sibling [bigfprev] = fnext ;
		}

		/* put bigf at the end of the list */
		Sibling [bigf] = EMPTY ;
		ASSERT (Child [i] != EMPTY) ;
		ASSERT (fprev != bigf) ;
		ASSERT (fprev != EMPTY) ;
		Sibling [fprev] = bigf ;
	    }

	}
    }

    /* --------------------------------------------------------------------- */
    /* postorder the assembly tree */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < nn ; i++)
    {
	Order [i] = EMPTY ;
    }

    k = 0 ;

    for (i = 0 ; i < nn ; i++)
    {
	if (Parent [i] == EMPTY && Nv [i] > 0)
	{
	    AMD_DEBUG1 (("Root of assembly tree "ID"\n", i)) ;
	    k = AMD_post_tree (i, k, Child, Sibling, Order, Stack
		) ;
	}
    }
}

/* ========================================================================= */
/* === AMD_preprocess ====================================================== */
/* ========================================================================= */

/* Sorts, removes duplicate entries, and transposes from the nonzero pattern of
 * a column-form matrix A, to obtain the matrix R.  The input matrix can have
 * duplicate entries and/or unsorted columns (AMD_valid (n,Ap,Ai) must not be
 * AMD_INVALID).
 *
 * This input condition is NOT checked.  This routine is not user-callable.
 */

/* AMD_preprocess does not check its input for errors or allocate workspace.
 * On input, the condition (AMD_valid (n,n,Ap,Ai) != AMD_INVALID) must hold.
 */

GLOBAL void AMD_preprocess
(
    Int n,		/* input matrix: A is n-by-n */
    const Int Ap [ ],	/* size n+1 */
    const Int Ai [ ],	/* size nz = Ap [n] */

    /* output matrix R: */
    Int Rp [ ],		/* size n+1 */
    Int Ri [ ],		/* size nz (or less, if duplicates present) */

    Int W [ ],		/* workspace of size n */
    Int Flag [ ]	/* workspace of size n */
)
{

    /* --------------------------------------------------------------------- */
    /* local variables */
    /* --------------------------------------------------------------------- */

    Int i, j, p, p2 ;

    ASSERT (AMD_valid (n, n, Ap, Ai) != AMD_INVALID) ;

    /* --------------------------------------------------------------------- */
    /* count the entries in each row of A (excluding duplicates) */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
	W [i] = 0 ;		/* # of nonzeros in row i (excl duplicates) */
	Flag [i] = EMPTY ;	/* Flag [i] = j if i appears in column j */
    }
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		W [i]++ ;	    /* one more entry in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the row pointers for R */
    /* --------------------------------------------------------------------- */

    Rp [0] = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	Rp [i+1] = Rp [i] + W [i] ;
    }
    for (i = 0 ; i < n ; i++)
    {
	W [i] = Rp [i] ;
	Flag [i] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* construct the row form matrix R */
    /* --------------------------------------------------------------------- */

    /* R = row form of pattern of A */
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		Ri [W [i]++] = j ;  /* put col j in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }
}

/* ========================================================================= */
/* === AMD_valid =========================================================== */
/* ========================================================================= */

/* Check if a column-form matrix is valid or not.  The matrix A is
 * n_row-by-n_col.  The row indices of entries in column j are in
 * Ai [Ap [j] ... Ap [j+1]-1].  Required conditions are:
 *
 *	n_row >= 0
 *	n_col >= 0
 *	nz = Ap [n_col] >= 0	    number of entries in the matrix
 *	Ap [0] == 0
 *	Ap [j] <= Ap [j+1] for all j in the range 0 to n_col.
 *      Ai [0 ... nz-1] must be in the range 0 to n_row-1.
 *
 * If any of the above conditions hold, AMD_INVALID is returned.  If the
 * following condition holds, AMD_OK_BUT_JUMBLED is returned (a warning,
 * not an error):
 *
 *	row indices in Ai [Ap [j] ... Ap [j+1]-1] are not sorted in ascending
 *	    order, and/or duplicate entries exist.
 *
 * Otherwise, AMD_OK is returned.
 *
 * In v1.2 and earlier, this function returned TRUE if the matrix was valid
 * (now returns AMD_OK), or FALSE otherwise (now returns AMD_INVALID or
 * AMD_OK_BUT_JUMBLED).
 */

GLOBAL Int AMD_valid
(
    /* inputs, not modified on output: */
    Int n_row,		/* A is n_row-by-n_col */
    Int n_col,
    const Int Ap [ ],	/* column pointers of A, of size n_col+1 */
    const Int Ai [ ]	/* row indices of A, of size nz = Ap [n_col] */
)
{
    Int nz, j, p1, p2, ilast, i, p, result = AMD_OK ;

    if (n_row < 0 || n_col < 0 || Ap == NULL || Ai == NULL)
    {
	return (AMD_INVALID) ;
    }
    nz = Ap [n_col] ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
	AMD_DEBUG0 (("column 0 pointer bad or nz < 0\n")) ;
	return (AMD_INVALID) ;
    }
    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	AMD_DEBUG2 (("\nColumn: "ID" p1: "ID" p2: "ID"\n", j, p1, p2)) ;
	if (p1 > p2)
	{
	    /* column pointers must be ascending */
	    AMD_DEBUG0 (("column "ID" pointer bad\n", j)) ;
	    return (AMD_INVALID) ;
	}
	ilast = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    AMD_DEBUG3 (("row: "ID"\n", i)) ;
	    if (i < 0 || i >= n_row)
	    {
		/* row index out of range */
		AMD_DEBUG0 (("index out of range, col "ID" row "ID"\n", j, i));
		return (AMD_INVALID) ;
	    }
	    if (i <= ilast)
	    {
		/* row index unsorted, or duplicate entry present */
		AMD_DEBUG1 (("index unsorted/dupl col "ID" row "ID"\n", j, i));
		result = AMD_OK_BUT_JUMBLED ;
	    }
	    ilast = i ;
	}
    }
    return (result) ;
}

