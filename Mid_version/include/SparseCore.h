#ifndef SPARSE_CORE_H
#define SPARSE_CORE_H

/* 
 * SparseCore假设C编译器是 ANSI C89 兼容的.  它没有使用 ANSI C99 的特性。
 */
#include "SparseCore.h"
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

/* ========================================================================== */
/* === SparseCore objects ====================================================== */
/* ========================================================================== */

/* 每个SparseCore对象对应一个编号 */

#define SPARSE_COMMON 0
#define SPARSE_CSC 1
#define SPARSE_FACTOR 2
#define SPARSE_DENSE 3
#define SPARSE_TRIPLET 4

/* ========================================================================== */
/* === SparseCore Common ======================================================= */
/* ========================================================================== */

/* itype 定义了 integer 使用的类型: */
#define SPARSE_INT 0		/* int 类型 */
#define SPARSE_LONG 2		/* Sparse_long 类型*/

/* 所有SparseCore例程的所有参数的itype必须匹配。*/

/* dtype定义了什么是数值类型(只支持双精度） */
#define SPARSE_DOUBLE 0	/* 所有的数值类型为 double */ 

/*  SparseCore 函数中 所有的数值参数类型必须和 dtype 匹配.
 * 标量的浮点值一直是作为双精度数组传递*/

/* xtype 定义所使用的数据的类型 */
#define SPARSE_PATTERN 0	/* 仅有分布模式，没有数值 */
#define SPARSE_REAL 1		/* 实数矩阵 */

/* sparse_common 的相关定义: */
#define SPARSE_MAXMETHODS 9	/* SparseCore_analyze 可以使用的最大方法数目 */

/* Common->status 的值.  0表示成功, 负数表示ERROR, 正数为WARNING */
#define SPARSE_OK 0			    /* success */
#define SPARSE_NOT_INSTALLED (-1)	/* failure: 没有安装方法 */
#define SPARSE_OUT_OF_MEMORY (-2)	/* failure: 内存越界 */
#define SPARSE_TOO_LARGE (-3)		/* failure: 整数溢出 */
#define SPARSE_INVALID (-4)		/* failure: 非法输入 */
#define SPARSE_NOT_POSDEF (1)		/* warning: 矩阵非正定 */
#define SPARSE_DSMALL (2)		/* warning: LDL'或diag(L)或LL'中 D 的绝对值过小 */

/* 排序方法 (L->ordering 也会使用) */
#define SPARSE_NATURAL 0	/* 自然排序 */
#define SPARSE_GIVEN 1		/* 使用给定的排列 */
#define SPARSE_AMD 2		/* 使用最小度排序(AMD) */
#define SPARSE_COLAMD 5	/* 对A使用最小度排序(AMD), A*A' 使用COLAMD 排序*/

/* POSTORDERED不是一种方法，而是自然排序后加上加权后序排序的结果。
 * 也用在 L->ordering, 而不是 method [ ].ordering. */
#define SPARSE_POSTORDERED 6	/* 自然排序, 后序次序。 */

/* 超节点策略(Common->supernodal) */
#define SPARSE_SIMPLICIAL 0	/* 恒为 simplicial 分解 */
#define SPARSE_AUTO 1		/* 基于矩阵选择 simpl/super */
#define SPARSE_SUPERNODAL 2	/* 恒为 supernodal 分解 */

typedef struct sparse_common_struct
{
    /* ---------------------------------- */
    /*       符号分解/数值分解的参数      */
    /* ----------------------------------*/

    double dbound ;	/* LDL'分解和update/downdate/rowadd/rowdel的
                     * 对角线项 D (或LL'分解的对角线L) 的最小绝对值。
	                 *  0 - dbound 范围内的元素值 被替换为 dbound.
	                 * -dbound — 0 范围内的元素值 被替换为 -dbound. 
                     * 如果dbound <= 0 , 则不做修改; 默认:0  */

    double grow0 ;	/* 对于一个简单分解, L->i 和 L->x 在必要时可以扩大.  grow0是它增长的因子。
    * 对于初始空间，L的大小是所需空间的 MAX(1,grow0) 倍。如果 L 用完了空间，那么L的新空间大小
	* 为所需空间的 MAX(1.2,grow0) 倍。如果您不打算在 Modify 模块中修改LDL的因子分解，
    * 请将grow0设置为0 (或将grow2设置为0，见下文)。默认值: 1.2 */

    double grow1 ;

    size_t grow2 ;	/* 对于一个简单分解, L中的每个列j初始化空间为 grow1* L->ColCount[j] + grow2
    * 如果 grow0 < 1, grow1 < 1, 或者 grow2 == 0, 那么空间大小精确等于 L->ColCount[j].
	* 如果列 j 耗尽了空间，那么它增长为 grow1*need + grow2 ，其中需要的是这一列的非零元总数#
	* 如果您不打算在Modify模块中修改分解, 设置 grow2 = 0， 默认:grow1 = 1.2, grow2 = 5. */

    size_t maxrank ;	/* 最大更新/降级的 rank。  有效值: 2，4 或8
			 *  值如果 < 2 则设定为2, > 8 则设定为 8.
	* 如果不是2的整次幂的话，就四舍五入到更高的2的整数次幂。
    * 为更新/降级 分配大小为 nrow-by-maxrank 的double类型工作区(Xwork)
	* 如果rank k被需求更新/降级 而 k > maxrank, 它是按照maxrank的步骤完成的。  
    * 默认: 8, 这是最快的情况. 可以通过将maxrank设置为2或4来减少内存使用。
	*/

    double supernodal_switch ;	/* supernodal vs simplicial 分解方法 */
    int supernodal ;		/* 如果 Common->supernodal <= SPARSE_SIMPLICIAL
	* (0) 则 SparseCore_analyze 执行simplicial 分析；
	* 如果 >= SPARSE_SUPERNODAL (2), 则执行supernodal 分析；
	* 如果 == SPARSE_AUTO (1) 且flop/nnz(L) < Common->supernodal_switch
	* 则执行 simplicial 分析.  否则执行supernodal分析.
	* 默认:  AUTO. supernodal_switch = 40 */

    int final_asis ;	/* 如果为 TRUE, 则忽略其他 final_* 参数(除final_pack)
			 * 该参数 在完成时保持原样。  默认: TRUE.*/

    int final_super ;	/* 如果为 TRUE, 当超节点分解完成时，留下一个超节点形式的因子。
             * 如果为 FALSE, 则分解完成时转换为一个 simplicial 因子。 Default: TRUE */

    int final_ll ;	/* 如果为 TRUE, 完成后将因子保留为LL'形式。
			 * 否则，保存为 LDL' 形式.  Default: FALSE */

    int final_pack ;	/* 如果为 TRUE, 完成后打包列.  如果 TRUE 且 SparseCore_factorize
			 *  由符号分析的 L 调用, L使用L->ColCount精确地分配所需的空间
	* 如果您计划修改分解，请将Common->final_pack 设置为FALSE，
	* 由于更新的原因，每一列都有一些额外的空闲空间用于将来的填充增长。 Default: TRUE */

    int final_monotonic ;   /* 如果为 TRUE, 确保列号在完成时是单调递增的. Default: TRUE */

    int final_resymbol ;/* 如果 SparseCore_factorize 执行了超节点分解, final_resymbol 为真 且
	* final_super 为假 ( 转换为simplicial 数值分解), 然后，从数字上去除松弛超节点合并产生的零项.
	* 这一做法不会删除由于精确数字取消而归零的项，
	* 因为这么做将会中断 update/downdate rowadd/rowdel 例程 Default: FALSE. */

    /* 松弛超节点合并参数 */
    double zrelax [3] ;
    size_t nrelax [3] ;

	/* 让 ns 作为为两个相邻超节点的列总数。 
	 * 让 z  为合并后的两个超节点中零项的部分(z包含以前合并后的零项).  
     * 两个超节点合并，如果满足:
	 *    (ns <= nrelax [0]) || (no new zero entries added) ||
	 *    (ns <= nrelax [1] && z < zrelax [0]) ||
	 *    (ns <= nrelax [2] && z < zrelax [1]) || (z < zrelax [2])
	 *
	 * 默认参数导致以下规则::
	 *    (ns <= 4) || (no new zero entries added) ||
	 *    (ns <= 16 && z < 0.8) || (ns <= 48 && z < 0.1) || (z < 0.05)
	 */

    int prefer_upper ;	    /* SparseCore_analyze 和 SparseCore_factorize work 在
			     * 当对称矩阵以上三角形式存储时，且使用填充减少排序时，工作速度最快。  
                 * 当矩阵按原样排列时，且当对称矩阵为下三角形式时，工作速度最快。
	             * 这个参数只影响 SparseCore_read 返回对称矩阵的方式.
	             * 如果为 TRUE（默认）, 对称矩阵始终返回 上三角 形式  */

    int quick_return_if_not_posdef ;	/* 如果为 TRUE, 如果矩阵非正定，
                     * 则超节点数值分解将快速返回  Default: FALSE. */

    int prefer_binary ;	    /* SparseCore_read_triplet 将一个为PATTERN的对称矩阵转换为实矩阵
	* 如果 prefer_binary 为 FALSE, 对角线项被设置为1 +行/列的度数, 而非对角线项被设置为-1 
    * (如果对角无零，则得到一个正定矩阵).  大多数对称图形是正定矩阵的图形。
	* 如果为 TRUE, 然后返回矩阵的对应每个元素都是1。 Default: FALSE.  */

    /* ---------------------------------------------------------------------- */
    /* 打印和错误处理选项 */
    /* ---------------------------------------------------------------------- */

    int print ;		/* Default: 3 */
    int precise ;	/* 如果为 TRUE,则打印16位数字。否则打印5位 */

    int try_catch ;	/* 如果为 TRUE, 则忽略errors; SparseCore 在try/catch 块中
	                 * 没有打印错误消息，也没有调用Common->error_handler函数 Default: FALSE.*/

    void (*error_handler) (int status, const char *file,
        int line, const char *message) ;

	/* Common->error_handler 是用户的错误处理例程. 如果不是NULL，这个例程将在SparseCore模块
	 * 中有 error 出现时调用.  status 可以是 SPARSE_OK (0), 负数为一个致命的error, 正数
	 * 对应一个 warning. file是一个字符串，其中包含发生错误的源代码文件的名称, line是该文件
     * 中的行号。message 是更详细地描述错误的字符串。 */

    /* ------------------------------ */
    /*            排序选项            */
    /* ------------------------- ---- */

    /* SparseCore_analyze 例程会尝试多种不同的排序方法然后搜索最好的那一个.
     * 它同样会使用不同的参数设置来尝试一种排序方法多次. 默认情况下使用两个排序:
     * 用户的排列(如果提供)或者 AMD, 这是最快的排序且一般提供良好的填充。
     *
     * 如果您知道最适合您的矩阵的排序方法，则可以将Common->nmethods设为1，
     * 并将Common->method[0]设为对应该方法的参数。
     *
     * 为了尝试其它的办法, 设置 Common->nmethods 为你想要尝试的排序参数
     * SparseCore_defaults例程中描述了一组默认方法及其参数，并在这里进行总结:
     *
     *	    Common->method [i]:
     *	    i = 0: 用户提供的排序 (SparseCore_analyze_p only)
     *	    i = 1: AMD (for both A and A*A')
     *	    i = 4: 自然排序
     *	    i = 8: AMD for A, COLAMD for A*A'
     *
     * 在调用完 SparseCore_start 或 SparseCore_default 之后，
     * 您可以通过修改Common.method[ ... ] 来修改希望尝试的一组方法.
     *
     * 例如，使用AMD，然后加权后序 排序:
     *
     *	    Common->nmethods = 1 ;
     *	    Common->method [0].ordering = SPARSE_AMD ;
     *	    Common->postorder = TRUE ;
     *
     * 使用自然排序 (不需要后序):
     *
     *	    Common->nmethods = 1 ;
     *	    Common->method [0].ordering = SPARSE_NATURAL ;
     *	    Common->postorder = FALSE ;
     *
     * 如果你要对上百个或更多具有相同非零模式的矩阵进行因式分解，
     * 你可能希望花费大量的时间来找到一个好的排序。
     * 在这种情况下, 尝试设定 Common->nmethods 为 9.
     * 花费在SparseCore_analysis上的时间会非常长，但是您只需要调用它一次。
     *
     * SparseCore_analyze 设置 Common->current 为 0-nmethods-1 中的一个值.
     * 每个排序方法都使用此参数定义的一组选项。
     */

    int nmethods ;	/* 要尝试的排序方法数目.  Default: 0.
			 * nmethods = 0 是特殊情况.  SparseCore_analyze 将尝试用户提供的排序(如果提供)和AMD
	* 设定 f1 和lnz 是AMD排序中 L 的浮点数统计以及非零元数目
	* 设定 anz 是对称矩阵 A 上/下三角部分 的非零元数目。 如果 fl/lnz < 500 or lnz/anz < 5
	* 那么这是一个好的排序. 使用找到的最佳排序，如果 nmethods > 0, 则使用的 methods由 method[ ]
	* 数组给出.  默认排序组合中的方法是 (1) 使用给定的排序（如果给出）.(2) 使用AMD
	* 允许的最大值是SparseCore MAXMETHODS。 */

    int current ;	/* 目前正在尝试的方法.  Default: 0.  范围为 [0, nemethods-1] */

    int selected ;	/* 找到的最好方法 */

    /* 一套排序方法和参数 */
    struct SparseCore_method_struct
    {
	/* 这一个方法的信息统计 */
	double lnz ;	    /* 对于“纯” L, nnz(L)不包括超节点合并中的零 */

	double fl ;	    /* 对一个 "pure", 实数简单 LL' 分解的浮点数统计，合并不会产生额外工作
			         *  减去 n 以得到 LDL' 浮点数统计. */

	/* 排序方法参数 */
	double prune_dense ;/* AMD, SYMAMD, CSYMAMD
			     * 稠密 行/列 控制. 对于一个 n-by-n dd阿对称矩阵， 行/列 中项数目大于
                 * MAX(16, prune_dense * sqrt (n)) 的行/列会在排序之前被删除。
                 * 它们出现在重排序的矩阵的末尾。
	    *
	    * 如果 prune_dense < 0,只有完全稠密的 行/列 被删除。
	    *
	    * 该参数也是COLAMD 和CCOLAMD的稠密列控制选项。 
        * 对于一个 m-by-n 的矩阵, 列中的项大于MAX (16, prune_dense * sqrt (MIN (m,n)))
	    * 会在排序之前被删除。它们出现在重排序的矩阵的末尾。
	    * SparseCore 分解 A*A', 所以它利用A'调用 COLAMD 和 CCOLAMD , 而不用 A.
	    * 因此, 这一参数影响 SparseCore中矩阵的 稠密“行”控制 
        * 以及 COLAMD 和 CCOLAMD 模块中的稠密“列”控制 
	    *
	    * 删除密集的行和列会提高排序方法的运行时间。  
        * 它对排序的质量有一定的影响(通常是很小的，有时好，有时坏)。
	    * Default: 10. */

	double prune_dense2 ;/* 对 COLAMD 和 CCOLAMD 的稠密行控制.
			    *  对于m-by-n 矩阵中的行，如果项多于 MAX (16, dense2 * sqrt (n))
	    * 则在排序前移除此行。  SparseCore的 矩阵在利用 COLAMD 或 CCOLAMD 排序前会转置，
	    * 所以实际上这个控制了 SparseCore 矩阵的稠密“列”, 以及COLAMD 或 CCOLAMD 矩阵的稠密"行"
	    *
	    * 如果 prune_dense2 < 0, 只有完全稠密的 行/列 被删除。
	    *
	    * Default: -1.  注意，这不是COLAMD和CCOLAMD的默认值。
        *  -1 最适合Cholesky。 10 则对LU最好。 */

	// double nd_oksep ;   /* in NESDIS, when a node separator is computed, it
	// 		     * discarded if nsep >= nd_oksep*n, where nsep is
	//     * the number of nodes in the separator, and n is the size of the
	//     * graph being cut.  Valid range is 0 to 1.  If 1 or greater, the
	//     * separator is discarded if it consists of the entire graph.
	//     * Default: 1 */

	// double other_1 [4] ; /* future expansion */

	// size_t nd_small ;    /* do not partition graphs with fewer nodes than
	// 		     * nd_small, in NESDIS.  Default: 200 (same as
	// 		     * METIS) */

	// size_t other_2 [4] ; /* future expansion */

	int aggressive ;    /* AMD、COLAMD、SYMAMD、CCOLAMD、CSYMAMD的聚合吸收 Default: TRUE */

	int order_for_lu ;  /* CCOLAMD可以被优化为LU或Cholesky 分解生成一个排序。
        * SparseCore 只做 Cholesky 分解。 然而，你可能想去用SparseCore作为CCOLAMD的接口，
        * 但是可以将它用于您自己的 LU 分解。在本例中，order_for_lu应该设置为FALSE
	    * 在SparseCore本身中分解时，您应该***永远不要***设置该参数为FALSE。Default:TRUE.*/

	// int nd_compress ;   /* If TRUE, compress the graph and subgraphs before
	// 		     * partitioning them in NESDIS.  Default: TRUE */

	// int nd_camd ;	    /* If 1, follow the nested dissection ordering
	// 		     * with a constrained minimum degree ordering that
	//     * respects the partitioning just found (using CAMD).  If 2, use
	//     * CSYMAMD instead.  If you set nd_small very small, you may not need
	//     * this ordering, and can save time by setting it to zero (no
	//     * constrained minimum degree ordering).  Default: 1. */

	// int nd_components ; /* The nested dissection ordering finds a node
	// 		     * separator that splits the graph into two parts,
	//     * which may be unconnected.  If nd_components is TRUE, each of
	//     * these connected components is split independently.  If FALSE,
	//     * each part is split as a whole, even if it consists of more than
	//     * one connected component.  Default: FALSE */

	/* 使用的 填充简化排序 */
	int ordering ;

	// size_t other_3 [4] ; /* future expansion */

    } method [SPARSE_MAXMETHODS + 1] ;

    int postorder ;	/* 如果为 TRUE, SparseCore_analyze 使用消去树的加权后序遍历来执行排序。
	* 提高supernode 合并。不影响基本的nnz(L)和浮点数计数.  Default: TRUE. */

    /* ---------------------------------------------------------------------- */
    /* 工作区 */
    /* ---------------------------------------------------------------------- */

    /* SparseCore 有几个例程花费的时间比他们需要的工作空间的大小要少。
     * 分配和初始化工作空间将主导运行时间，除非工作空间只分配和初始化一次。
     * SparseCore在需要时分配这个空间，并在调用SparseCore时将其保存在这里。
     * SparseCore_start 将这些指针设置为 NULL(这就是为什么它必须是SparseCore中调用的第一个例程)。
     * SparseCore_finish 释放掉这些工作区 (这就是为什么它必须是SparseCore中调用的最后一个例程)。
     */

    size_t nrow ;	/* Flag 和 Head 的规模 */
    Sparse_long mark ;	/* 标记数组里的标记值*/
    size_t iworksize ;	/* Iwork的大小. 上限: 6*nrow+ncol */
    size_t xworksize ;	/* Xwork的大小, 以 bytes 为单位.
			 * 针对 update/downdate 为 maxrank*nrow*sizeof(double) .
			 * 否则为 2*nrow*sizeof(double) */

    /* 初始化工作区: 调用 SparseCore 之间必要的内容 */
    void *Flag ;	/* 规模为 nrow, 一个整数矩阵. 在各种调用SparseCore 例程间清理
			 * (Flag [i] < mark) */

    void *Head ;	/* 规模为 nrow+1, 一个整数矩阵. 在各种调用SparseCore 例程间清理
			 * (Head [i] = EMPTY) */

    void *Xwork ; 	/* 一个double数组. size可变. 对大部分例程来说为 nrow 大小
			 * (SparseCore_rowfac, SparseCore_add, SparseCore_aat, SparseCore_norm, SparseCore_ssmult)
	* 对于 SparseCore_rowadd 和 SparseCore_rowdel它的规模为 2*nrow.  对于 SparseCore_updown,
	* 其大小为 maxrank*nrow 其中 maxrank 是 2, 4,或 8. 在各种调用SparseCore 例程间清理 (设置为0). */

    /* 未初始化的工作空间，在调用SparseCore之间不需要的内容  */
    void *Iwork ;	/* 大小为 iworksize, 对大部分例程来说是 2*nrow+ncol ,
			 * SparseCore_analyze的最大值为6*nrow+ncol */

    int itype ;		/* 如果为 SPARSE_LONG, Flag, Head, 和 Iwork 都是
                         * Sparse_long.  否则所有的这些参数都为 int. */

    int dtype ;		/* double */

	/* Common->itype 和 Common->dtype 用来确定所有稀疏矩阵、三元组、稠密矩阵和
     * 使用这个公共结构创建的因子的types. 所有SparseCore例程的所有参数的itypes
     * 和dtype必须匹配。  */

    int no_workspace_reallocate ;   /* 这是内置的标志, 用作SparseCore_analyze的预防措施。
	* 通常为 false.  如果为true，SparseCore_allocate_work 不允许重新分配任何工作空间;
	* 他们必须使用Common中现存的工作空间 (Iwork, Flag, Head,and Xwork)*/

    /* ---------------------------------------------------------------------- */
    /* 统计信息 */
    /* ---------------------------------------------------------------------- */

    /* fl 和 lnz 只在 Cholesky模块中的 SparseCore_analyze 和 SparseCore_rowcolcounts设置,
     * modfl只在 Modify模块中设定 */

    int status ;	    /* 错误码 */
    double fl ;		    /* LL' 最近一次分析的浮点数统计 */
    double lnz ;	    /* L 中的基本非零元 */
    double anz ;	    /* tril(A)中的非零元数目，如果 A 是对称/下三角;
			     * triu(A)中的非零元数目 如果 A 是对称/上三角,;或者 tril(A*A')中的非零元数目
			     * 如果A是非对称矩阵, 在最后一次调用SparseCore分析. */
    size_t malloc_count ;   /* 已经 malloc的对象数目减去 已经 free 掉的对象数目*/
    size_t memory_usage ;   /* 峰值内存使用量(以字节为单位) */
    size_t memory_inuse ;   /* 当前内存使用情况(单位为字节) */

    double nrealloc_col ;   /* 列重新分配的# */
    double nrealloc_factor ;/* 由于列重新分配 而造成的factor重新分配的# */
    double ndbounds_hit ;   /* 被dbound对角修改的时间 */

    double rowfacfl ;	    /* 最后一次调用 SparseCore_rowfac 的 flops 数 */
    double aatfl ;	    /* 计算 A(:,f)*A(:,f)' 的 flops 数 */

    

    int blas_ok ;           /* 如果BLAS的 int 溢出，则为FALSE ; 否则为 TRUE */

    /* ---------------------------------------------------------------------- */
    /* HnuSparseQR 控制参数 */
    /* ---------------------------------------------------------------------- */

    double SPQR_grain ;      /* 任务规模 >= max (total flops / grain) */
    double SPQR_small ;      /* 任务规模 >= small */
    int SPQR_shrink ;        /* 控件堆栈再分配的 方法 */

    /* ---------------------------------------------------------------------- */
    /* HnuSparseQR 统计信息 */
    /* ---------------------------------------------------------------------- */

    double SPQR_flopcount ;         /* SPQR的浮点计算次数 */

    double SPQR_flopcount_bound ;   /* 浮点计算量的上界 */
    double SPQR_tol_used ;          /* 使用的阈值 */
    double SPQR_norm_E_fro ;        /* 删除项的范数 */

    /* 范围是 SPQR_istat [0:9] */
    Sparse_long SPQR_istat [10] ;

    double thread_pool_time;   //线程池创建时间

} sparse_common ;


/* 关于超节点分析 */
#define SPARSE_ANALYZE_FOR_SPQR     0
#define SPARSE_ANALYZE_FOR_CHOLESKY 1

/* -------------------------------------------------------------------------- */
/* SparseCore_start:  使用 SparseCore 必须最优先调用的例程 */
/* -------------------------------------------------------------------------- */

int SparseCore_start (sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_finish:  使用 SparseCore 必须最后调用的例程 */
/* -------------------------------------------------------------------------- */

int SparseCore_finish (sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_defaults:  恢复默认设置 */
/* -------------------------------------------------------------------------- */

int SparseCore_defaults (sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_maxrank:  为 update/downdate 返回有效的最大rank */
/* -------------------------------------------------------------------------- */

size_t SparseCore_maxrank (size_t, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_allocate_work:  给Common结构体 分配空间 */
/* -------------------------------------------------------------------------- */

int SparseCore_allocate_work (size_t, size_t, size_t, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_free_work:  释放 Common 结构体中的空间   */
/* -------------------------------------------------------------------------- */

int SparseCore_free_work (sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_clear_flag:  清除 Common 结构体中的 Flag 工作区  */
/* -------------------------------------------------------------------------- */

/* use a macro for speed */
#define SPARSE_CLEAR_FLAG(Common) \
{ \
    Common->mark++ ; \
    if (Common->mark <= 0) \
    { \
	Common->mark = EMPTY ; \
	CORE(clear_flag) (Common) ; \
    } \
}

Sparse_long SparseCore_clear_flag (sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_error:  当 SparseCore 遇到错误时调用 */
/* -------------------------------------------------------------------------- */

int SparseCore_error (int, const char *, int, const char *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_dbound:  仅供SparseCore内部使用 */
/* -------------------------------------------------------------------------- */

double SparseCore_dbound (double, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_hypot:  准确计算√(x*x + y*y) */
/* -------------------------------------------------------------------------- */

double SparseCore_hypot (double, double) ;


/* ========================================================================== */
/* === Core/sparse_csc ================================================== */
/* ========================================================================== */

/* 以压缩列形式存储的稀疏矩阵. */

typedef struct sparse_csc_struct
{
    size_t nrow ;	/* 矩阵大小是nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;	/* 矩阵中元素的最大数目 */

    /* 指向 int 或 Sparse_long 的指针 */
    void *p ;		/* p [0..ncol], 列指针 */
    void *i ;		/* i [0..nzmax-1], 行索引 */

    /* 仅对非压缩矩阵生效 */
    void *nz ;		/* nz [0..ncol-1], 每一列的非零元 #. 在压缩格式下
			 *  列 j 的非零元素对应位置为 A->i [A->p [j] ... A->p [j+1]-1] 
	* 在非压缩格式下, 列 j 的非零元素为 A->i [A->p [j] ... A->p [j]+A->nz[j]-1]
	* 在这两种情况下, 数值(如果存在)位于数组 x 的相应位置 */

    /* double 类型 */
    void *x ;		/* 大小nzmax或2*nzmax，如果存在 */
    void *z ;		/* 大小nzmax，如果存在 */

    int stype ;		/* 描述考虑矩阵的哪些部分:
			 *
	* 0:  非对称矩阵: 使用上三角部分和下三角部分(矩阵实际上可能在模式和值上是对称的
    * ，但是这两个部分都是显式存储和使用的).  可以是正方形或长方形。
	*  
	* >0: 矩阵是正方形和对称的,使用上三角部分.
	*     下三角部分的元素忽略.
	* <0: 矩阵是正方形和对称的,使用下三角部分.
	*     上三角部分的元素忽略.
	*
	* 注意到 stype>0 和 stype<0 SparseCore_sparse和SparseCore_triplet是不同的。
    * 有关更多细节，请参阅SparseCore_triplet数据结构。
	*/

    int itype ;		/* SPARSE_INT:     p, i, nz 为 int.
			 * SPARSE_LONG:    p, i, nz 为 Sparse_long */

    int xtype ;		/* pattern, real */
    int dtype ;		/* x 和 z 为double */
    int sorted ;	/* 如果列已排序为TRUE, 否则为FALSE */
    int packed ;	/* 如果是压缩格式为TRUE( 不需要nz ), 否则为FALSE( 需要nz ) */

} sparse_csc ;

typedef struct SparseCore_descendant_score_t {
  double score;
  Sparse_long d;
} descendantScore;

/* 用于用qsort排序子代超节点 */

int SparseCore_score_comp (struct SparseCore_descendant_score_t *i,
			       struct SparseCore_descendant_score_t *j);

/* -------------------------------------------------------------------------- */
/* SparseCore_allocate_sparse:  分配稀疏矩阵空间 */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_allocate_sparse (size_t, size_t, size_t, int, int,
    int, int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_free_sparse:  释放稀疏矩阵空间 */
/* -------------------------------------------------------------------------- */

int SparseCore_free_sparse (sparse_csc **, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_reallocate_sparse:  改变稀疏矩阵的大小( # 条目) */
/* -------------------------------------------------------------------------- */

int SparseCore_reallocate_sparse ( size_t, sparse_csc *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_nnz:  返回稀疏矩阵中的非零元数目 */
/* -------------------------------------------------------------------------- */

Sparse_long SparseCore_nnz (sparse_csc *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_speye:  稀疏单位矩阵 */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_speye (size_t, size_t, int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_spzeros:  稀疏 零矩阵 */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_spzeros (size_t, size_t, size_t, int,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_transpose:  稀疏矩阵转置 */
/* -------------------------------------------------------------------------- */

/* 返回 A' “values”参数为0,1, 分别表示模式转置、数组转置(A.'). */

sparse_csc *SparseCore_transpose (sparse_csc *, int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_transpose_unsym:  一个非对称稀疏矩阵的转置 */
/* -------------------------------------------------------------------------- */

/* 计算F = A'， A (:， F)'，或A (p, F)'，其中A不对称且F已经分配。
 * 	参见SparseCore_transpose获得一个更简单的例程 */

int SparseCore_transpose_unsym (sparse_csc *, int, Sparse_long *,
    Sparse_long *, size_t, sparse_csc *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_transpose_sym:  一个对称稀疏矩阵的转置 */
/* -------------------------------------------------------------------------- */

/* 计算F = A'或A (p,p)'，其中A是对称的，F已经分配。
 * 参见SparseCore_transpose获得一个更简单的例程。 */

int SparseCore_transpose_sym (sparse_csc *, int, Sparse_long *,
    sparse_csc *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_ptranspose:  转置一个稀疏矩阵 */
/* -------------------------------------------------------------------------- */

/*如果A是对称的，返回A'或A(p,p)'。如果A不对称，返回A'， A(:，f)'，或A(p,f)' */

sparse_csc *SparseCore_ptranspose (sparse_csc *, int, Sparse_long *,
    Sparse_long *, size_t, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_sort: 对稀疏矩阵的每一列中的行索引排序 */
/* -------------------------------------------------------------------------- */

int SparseCore_sort (sparse_csc *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_band:  C = tril (triu (A,k1), k2) */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_band (sparse_csc *, Sparse_long,
    Sparse_long, int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_band_inplace:  A = tril (triu (A,k1), k2) */
/* -------------------------------------------------------------------------- */

int SparseCore_band_inplace (Sparse_long, Sparse_long, int,
    sparse_csc *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_aat:  C = A*A' or A(:,f)*A(:,f)' */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_aat (sparse_csc *, Sparse_long *, size_t,
    int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_copy_sparse:  C = A,创建一个稀疏矩阵的精确副本 */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_copy_sparse (sparse_csc *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_copy:  C = A, 可能会改变 stype */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_copy (sparse_csc *, int, int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_add: C = alpha*A + beta*B */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_add (sparse_csc *, sparse_csc *, double *,
    double *, int, int, sparse_common *) ;

// 超节点结构体，保存当前在消去树中的深度和孩子节点在后续遍历中的位置
typedef struct chol_supernode_struct
{
    int Snchild;
    int Sdepth;
    int *Schild;
} chol_supernode ;

/* ========================================================================== */
/* === Core/sparse_factor ================================================== */
/* ========================================================================== */

/* 符号分解或数值分解, 简单分解/ 超节点分解.
 * 在所有情况下，L的列中的行下标都是有序的. */
typedef struct sparse_factor_struct
{
    /* ---------------------------------------------------------------------- */
    /* 对于简单分解和超节点分解 */
    /* ---------------------------------------------------------------------- */

    size_t n ;		/* L 是 n*n 的 */

    size_t minor ;	/* 如果分解失败, L->minor 表示失败在第几列(0, n-1)
			 * 如果minor值为n 表示分解成功 或者矩阵还没有被分解。 */

    /* ---------------------------------------------------------------------- */
    /* 符号排序与分析 */
    /* ---------------------------------------------------------------------- */

    void *Perm ;	/* size n, 重排列使用 */
    void *ColCount ;	/* size n, 简单分解 L 的列数统计 */

    void *IPerm ;       /* size n, 逆排列。只有在使用Bset时才由cholmod_solve2创建. */
    /* ---------------------------------------------------------------------- */
    /* 简单分解  */
    /* ---------------------------------------------------------------------- */

    size_t nzmax ;	/*  i 和 x 的大小*/

    void *p ;		/* p [0..ncol], 列指针 */
    void *i ;		/* i [0..nzmax-1], 行索引 */
    void *x ;		/* x [0..nzmax-1], 值。 */
    void *z ;
    void *nz ;		/* nz [0..ncol-1], 每列中的非零项.
			 * i [p [j] ... p [j]+nz[j]-1] 包含行索引,
			 * 且 数值 也在 x 的对应位置上。
			 *   i [p [k]] 的值 恒为 k  */

    void *next ;	/* size ncol+2. next [j] 在 i/x 中表示下一列  */
    void *prev ;	/* size ncol+2. prev [j] 在 i/x 中表示上一列.
			 * 列表的头是 ncol+1 ，尾是 ncol . */

    /* ---------------------------------------------------------------------- */
    /* 超节点分解 */
    /* ---------------------------------------------------------------------- */

    /* 注意到 L->x 由 简单分解的数据结构共享.  L->x 对于简单分解 规模为L->nzmax,
     * 对于超节点分解 规模为L->xsize  */

    size_t nsuper ;	/* 超节点数目 */
    size_t ssize ;	/* size of s, 超节点的整数部分 */
    size_t xsize ;	/* size of x, 超节点的实数部分 */
    size_t maxcsize ;	/* 最大更新矩阵的size */
    size_t maxesize ;	/* 超节点行中的最大元素值 # 不包括三角形部分 */

    void *super ;	/* size nsuper+1, 每个supernode节点的第一列列号 */
    void *pi ;		/* size nsuper+1, 指向整数模式的指针 */
    void *px ;		/* size nsuper+1, 指向实数模式的指针 */
    void *s ;		/* size ssize, 超节点的整数部分 */

    int max_depth;     /* 松弛超节点消去树深度 */
    int *Snode_num;    /* 消去树每层深度包含节点数 */
    int **Snode_distribution;   /* 消去树每层深度包含节点位置 */
    chol_supernode *Snode;      /* 超节点消去树，从根节点出发 */

    /* ---------------------------------------------------------------------- */
    /* 分解类型 */
    /* ---------------------------------------------------------------------- */

    int ordering ;	/* 使用的排序方法 */

    int is_ll ;		/* 如果做 LL' 为TRUE , 否则是 FALSE */
    int is_super ;	/* 如果做超节点分解则为TRUE, 否则为 FALSE */
    int is_monotonic ;	/*  如果L的列是0...n-1顺序的，则为TRUE。
			 * 只适用于 简单 数值类型 */

    /* SparseCore_factor可以表示8种类型的因子对象(只使用了6种) :
     *
     * 数值类型(xtype不是SPARSE_PATTERN)
     * --------------------------------------------
     *
     * 简单LDL':  (is_ll FALSE, is_super FALSE).  
     * 使用上面的简单组件(nzmax, p, i, x, z, nz, next，和prev)
     * 以压缩列的形式存储。L的单位对角线不被存储，而D被存储在它所在的位置。
     * 没有超级节点。
     *
     * 简单LL': (is_ll TRUE, is_super FALSE).  
     * 使用与简单LDL'相同的存储方案，只是D没有出现。
     * L的每一列的第一项是L的那一列的对角项。
     *
     * 超节点LL': (is_ll TRUE, is_super TRUE).  
     * 一个超节点因子，使用上面描述的
     * 超节点组件(nsuper, ssize, xsize, maxcsize, maxesize, super, pi, px, s, x和z)。
     * 
     * 符号类型 (xtype为SPARSE_PATTERN)
     * -----------------------------------------
     *
     * 简单LDL': (is_ll FALSE, is_super FALSE).  
     * 除了Perm和ColCount，什么也没有.
     *
     * 简单LL': (is_ll TRUE, is_super FALSE).  
     * 除了is_ll标志，与单纯的LDL'相同
     *
     * 超节点LL': (is_ll TRUE, is_super TRUE).  
     * 一个超节点符号分解。除了数值(x和z)外，
     * 所有的超节点因数分解都存在简单的符号信息(Perm和ColCount)。
     */

    int itype ; /* 整数数组是Perm, ColCount, p, i, nz, next, prev, super, pi, px和s。
		 * SPARSE_INT：所有这些都是int数组.
		 * SPARSE_LONG:    所有整型都是Sparse_long. */
    int xtype ; /* 二进制或者实数 */
    int dtype ; /* x z为double */

} sparse_factor ;

/* -------------------------------------------------------------------------- */
/* SparseCore_allocate_factor: 指定一个分解(符号LL'或LDL') */
/* -------------------------------------------------------------------------- */

sparse_factor *SparseCore_allocate_factor (size_t, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_free_factor:  释放一个分解 */
/* -------------------------------------------------------------------------- */

int SparseCore_free_factor (sparse_factor **, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_reallocate_factor:  更改分解中的#项 */
/* -------------------------------------------------------------------------- */

int SparseCore_reallocate_factor (size_t, sparse_factor *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_change_factor:  改变因子的类型(例如，LDL'变为LL') */
/* -------------------------------------------------------------------------- */

int SparseCore_change_factor ( int, int, int, int, int, sparse_factor *,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_pack_factor:  压缩一个因子的列 */
/* -------------------------------------------------------------------------- */

/* 压缩一个简单因子的列。与SparseCore_change_factor不同，
 * 它可以压缩一个因子的列，即使它们不是按照它们的自然顺序存储的(非单调)。 */

int SparseCore_pack_factor (sparse_factor *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_reallocate_column:  调整因子的单个列的大小 */
/* -------------------------------------------------------------------------- */

int SparseCore_reallocate_column (size_t, size_t, sparse_factor *,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_factor_to_sparse:  创建一个因子的稀疏矩阵副本 */
/* -------------------------------------------------------------------------- */

/* 只对数值因子起作用，而不是符号因子 */

sparse_csc *SparseCore_factor_to_sparse (sparse_factor *,
	sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_copy_factor:  创造一个因子的副本 */
/* -------------------------------------------------------------------------- */

sparse_factor *SparseCore_copy_factor (sparse_factor *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_factor_xtype: 改变因子的xtype */
/* -------------------------------------------------------------------------- */

int SparseCore_factor_xtype (int, sparse_factor *, sparse_common *) ;


/* ========================================================================== */
/* === Core/dense_array =================================================== */
/* ========================================================================== */

/* 一种面向列的稠密矩阵。它没有itype，因为它不包含整数。第i行和第j列中的项位于x [i+j*d]中。
 */

typedef struct dense_array_struct
{
    size_t nrow ;	
    size_t ncol ;   /* 矩阵大小为nrow*ncol */

    size_t nzmax ;	/* 矩阵中元素的最大数目 */
    size_t d ;		/* 前导维度 (d >= nrow) */
    void *x ;		/* 如果存在，则大小为nzmax或者2*nzmax */
    void *z ;		/* 如果存在，则大小为nzmax */
    int xtype ;		/* pattern, real, complex, zomplex */
    int dtype ;		/* x z为double */

} dense_array ;

/* -------------------------------------------------------------------------- */
/* SparseCore_allocate_dense:  分配一个稠密矩阵(内容未初始化) */
/* -------------------------------------------------------------------------- */

dense_array *SparseCore_allocate_dense (size_t, size_t, size_t, int,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_zeros: 分配一个稠密矩阵并且置为0 */
/* -------------------------------------------------------------------------- */

dense_array *SparseCore_zeros (size_t, size_t, int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_ones: 分配一个稠密矩阵并且全部置为1 */
/* -------------------------------------------------------------------------- */

dense_array *SparseCore_ones (size_t, size_t, int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_eye: 分配一个稠密矩阵并且设置为单位矩阵 */
/* -------------------------------------------------------------------------- */

dense_array *SparseCore_eye (size_t, size_t, int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_free_dense:  释放一个稠密矩阵 */
/* -------------------------------------------------------------------------- */

int SparseCore_free_dense (dense_array **, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_ensure_dense:  确保一个稠密矩阵有给定的大小和数据类型 */
/* -------------------------------------------------------------------------- */

dense_array *SparseCore_ensure_dense (dense_array **, size_t, size_t, size_t,
    int, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_sparse_to_dense:  创建稀疏矩阵的稠密矩阵副本 */
/* -------------------------------------------------------------------------- */

dense_array *SparseCore_sparse_to_dense (sparse_csc *,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_dense_to_sparse:  创建一个稠密矩阵的稀疏矩阵副本 */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_dense_to_sparse (dense_array *, int,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_copy_dense:  创建一个稠密矩阵的副本 */
/* -------------------------------------------------------------------------- */

dense_array *SparseCore_copy_dense (dense_array *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_copy_dense2:  拷贝稠密矩阵（提前分配好空间） */
/* -------------------------------------------------------------------------- */

int SparseCore_copy_dense2 (dense_array *, dense_array *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_dense_xtype: 改变稠密矩阵的数据类型 */
/* -------------------------------------------------------------------------- */

int SparseCore_dense_xtype (int, dense_array *, sparse_common *) ;


/* ========================================================================== */
/* === Core/sparse_triplet ================================================= */
/* ========================================================================== */

/* 以三元组形式存储的稀疏矩阵 */

typedef struct sparse_triplet_struct
{
    size_t nrow ;	
    size_t ncol ;   /* 矩阵的大小为nrow*ncol */

    size_t nzmax ;	/* 矩阵中元素的最大数目 */
    size_t nnz ;	/* 矩阵中的非零元数目 */

    void *i ;		/* i [0..nzmax-1], 行索引 */
    void *j ;		/* j [0..nzmax-1], 列索引 */
    void *x ;		/* 如果存在，则大小为nzmax或者2*nzmax */
    void *z ;		/* 如果存在，则大小为nzmax */

    int stype ;		/* 描述考虑矩阵的哪些部分:
			 *
	* 0:  矩阵是“不对称的”:同时使用上三角部分和下三角部分(实际上，矩阵在模式和
    * 值上可能是对称的，但这两个部分都是显式存储和使用的)。可以是正方形或长方形。
    * 
	* >0: 对称方阵。当矩阵转换为SparseCore_sparse形式时，下三角部分中的项被转置并
    * 添加到上三角部分。
    * 
	* <0: 对称方阵。当矩阵转换为SparseCore_sparse形式时，上三角部分中的项被转置并
    * 添加到下三角部分。
	*
	* 注意，对于SparseCore_sparse和SparseCore_triplet, stype>0和stype<&lt;>0是不同的。
    * 原因很简单。您可以通过逆置换将一个对称三联体矩阵的行和列索引替换为新的行和列索引。
    * 假设P = L->Perm是你的排列，Pinv是一个大小为n的数组。假设一个对称矩阵a由一个三元
    * 组矩阵T表示，其元素只在上三角形部分。那么以下代码:
	*
	*	Ti = T->i ;
	*	Tj = T->j ;
	*	for (k = 0 ; k < n  ; k++) Pinv [P [k]] = k ;
	*	for (k = 0 ; k < nz ; k++) Ti [k] = Pinv [Ti [k]] ;
	*	for (k = 0 ; k < nz ; k++) Tj [k] = Pinv [Tj [k]] ;
	*
	* 创建C=P*A*P'的三重形式。但是，如果T最初只包含上三角项(T->stype = 1)，那么在置换
    * 之后，它同时包含上三角项和下三角项。在构造A的SparseCore_sparse形式时，应该调换这些项，
    * 这正是SparseCore_triplet_to_sparse所做的。因此:
	*
	*	C = SparseCore_triplet_to_sparse (T, 0, &Common) ;
	*
	* 将会返回矩阵C = P*A*P'.
	*
	* 由于三元组矩阵T的生成非常简单，所以在将T转换为SparseCore_sparse形式之前，删除不需要的
    * 项是非常容易的。因此，如果您将这些条目包含在T中，SparseCore假设一定有原因(如上面的原因)。
    * 因此，三元矩阵中的任何项都不会被忽略。
	*/

    int itype ; /* SPARSE_LONG: i j为Sparse_long.  其余为int */
    int xtype ; /* pattern, real, complex, zomplex */
    int dtype ; /* x z为double */

} sparse_triplet ;

/* -------------------------------------------------------------------------- */
/* SparseCore_allocate_triplet:  分配一个三元组数组 */
/* -------------------------------------------------------------------------- */

sparse_triplet *SparseCore_allocate_triplet (size_t, size_t, size_t, int, int,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_free_triplet:  释放一个三元组数组 */
/* -------------------------------------------------------------------------- */

int SparseCore_free_triplet (sparse_triplet **, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_reallocate_triplet:  更改三元组矩阵中的项# */
/* -------------------------------------------------------------------------- */

int SparseCore_reallocate_triplet (size_t, sparse_triplet *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_sparse_to_triplet:  创建一个稀疏矩阵的三元组矩阵副本 */
/* -------------------------------------------------------------------------- */

sparse_triplet *SparseCore_sparse_to_triplet (sparse_csc *,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_triplet_to_sparse:  创建一个三元组矩阵的稀疏矩阵的副本 */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_triplet_to_sparse (sparse_triplet *, size_t,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_copy_triplet:  创建一个三元组矩阵的副本 */
/* -------------------------------------------------------------------------- */

sparse_triplet *SparseCore_copy_triplet (sparse_triplet *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_triplet_xtype: 改变三元组矩阵的xtype */
/* -------------------------------------------------------------------------- */

int SparseCore_triplet_xtype (int, sparse_triplet *, sparse_common *) ;


/* ========================================================================== */
/* === Core/SparseCore_memory ================================================== */
/* ========================================================================== */

/* 用户可以使用它们，就像malloc和free一样。您甚至可以malloc对象并使用SparseCore_free
 * 安全地释放它，反之亦然(除了内存使用统计数据将被破坏)。这些例程确实与malloc和free不同。
 * 例如，如果给SparseCore_free一个空指针，它什么也不做(不像ANSI free)。如果给定一个非空指
 * 针和一个非零大小的指针，SparseCore_realloc不会返回NULL，即使它失败了(它返回原始指针并
 * 设置一个Common->status内的错误代替代)。
 *
 * SparseCore跟踪分配的内存量，因此SparseCore_free例程还会获取被释放对象的大小。这只用于统计。
 * 如果您(SparseCore的用户)传递了错误的大小，唯一的后果是内存使用统计信息将被破坏。
 */

void *SparseCore_malloc (size_t, size_t, sparse_common *) ;

void *SparseCore_calloc (size_t, size_t, sparse_common *) ;

void *SparseCore_free (size_t, size_t, void *, sparse_common *) ;

void *SparseCore_realloc (size_t, size_t, void *, size_t *, sparse_common *) ;

int SparseCore_realloc_multiple (size_t, int, int, void **, void **, void **,
    void **, size_t *, sparse_common *) ;

double SparseCore_norm_dense( dense_array *X, int norm, sparse_common *Common);

double SparseCore_norm_sparse( sparse_csc *A, int norm, sparse_common *Common);

int SparseCore_sdmult( sparse_csc *A, int transpose, double alpha [2], double beta [2], 
    dense_array *X, dense_array *Y, sparse_common *Common);

/* -------------------------------------------------------------------------- */
/* SparseCore_check_common:  检查common结构体 */
/* -------------------------------------------------------------------------- */

int SparseCore_check_common (sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_print_common:  输出Common对象 */
/* -------------------------------------------------------------------------- */

int SparseCore_print_common (const char *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_cpu_stats:  输出 cpu 信息 */
/* -------------------------------------------------------------------------- */

int SparseCore_cpu_stats (sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_check_sparse:  检查一个稀疏矩阵 */
/* -------------------------------------------------------------------------- */

int SparseCore_check_sparse (sparse_csc *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_print_sparse  输出 稀疏矩阵  */
/* -------------------------------------------------------------------------- */

int SparseCore_print_sparse (sparse_csc *, const char *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_check_dense:  检查一个稠密矩阵  */
/* -------------------------------------------------------------------------- */

int SparseCore_check_dense (dense_array *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_print_dense:  打印一个稠密矩阵 */
/* -------------------------------------------------------------------------- */

int SparseCore_print_dense (dense_array *, const char *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_check_factor:  检查 L */
/* -------------------------------------------------------------------------- */

int SparseCore_check_factor (sparse_factor *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_print_factor:  打印 分解因子L */
/* -------------------------------------------------------------------------- */

int SparseCore_print_factor (sparse_factor *, const char *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_check_triplet:  检查一个三元组格式的稀疏矩阵 */
/* -------------------------------------------------------------------------- */

int SparseCore_check_triplet (sparse_triplet *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_print_triplet:  打印一个三元组矩阵 */
/* -------------------------------------------------------------------------- */

int SparseCore_print_triplet (sparse_triplet *, const char *, sparse_common *);

/* -------------------------------------------------------------------------- */
/* SparseCore_check_subset:  检查一个子集 */
/* -------------------------------------------------------------------------- */

int SparseCore_check_subset (Sparse_long *, Sparse_long, size_t,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_print_subset:  打印子集 */
/* -------------------------------------------------------------------------- */

int SparseCore_print_subset (Sparse_long *, Sparse_long, size_t,
    const char *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_check_perm:  检查一个重排列 */
/* -------------------------------------------------------------------------- */

int SparseCore_check_perm (Sparse_long *, size_t, size_t, sparse_common *);

/* -------------------------------------------------------------------------- */
/* SparseCore_print_perm:  打印一个置换向量 */
/* -------------------------------------------------------------------------- */

int SparseCore_print_perm (Sparse_long *, size_t, size_t, const char *,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_check_parent:  检查消去树 */
/* -------------------------------------------------------------------------- */

int SparseCore_check_parent (Sparse_long *, size_t, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_print_parent 打印parent */
/* -------------------------------------------------------------------------- */

int SparseCore_print_parent (Sparse_long *, size_t, const char *,
    sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_read_sparse: 从文件中读取稀疏矩阵  */
/* -------------------------------------------------------------------------- */

sparse_csc *SparseCore_read_sparse (FILE *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_read_triplet: 从文件中读取三元组矩阵 */
/* -------------------------------------------------------------------------- */

sparse_triplet *SparseCore_read_triplet (FILE *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_read_dense: 从文件中读取稠密矩阵 */
/* -------------------------------------------------------------------------- */

dense_array *SparseCore_read_dense (FILE *, sparse_common *) ; 

/* -------------------------------------------------------------------------- */
/* SparseCore_read_matrix: 从文件中读取 稀疏 / 稠密 矩阵 */
/* -------------------------------------------------------------------------- */

void *SparseCore_read_matrix (FILE *, int, int *, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_write_sparse: 将稀疏矩阵写入文件 */
/* -------------------------------------------------------------------------- */

int SparseCore_write_sparse (FILE *, sparse_csc *, sparse_csc *,
    const char *c, sparse_common *) ;

/* -------------------------------------------------------------------------- */
/* SparseCore_write_dense: 将稠密矩阵写入文件 */
/* -------------------------------------------------------------------------- */

int SparseCore_write_dense (FILE *, dense_array *, const char *,
    sparse_common *) ;

/* ========================================================================== */
/* === symmetry types ======================================================= */
/* ========================================================================== */

#define SPARSE_MM_RECTANGULAR 1
#define SPARSE_MM_UNSYMMETRIC 2
#define SPARSE_MM_SYMMETRIC 3
#define SPARSE_MM_HERMITIAN 4
#define SPARSE_MM_SKEW_SYMMETRIC 5
#define SPARSE_MM_SYMMETRIC_POSDIAG 6
#define SPARSE_MM_HERMITIAN_POSDIAG 7

/* ========================================================================== */
/* === Numerical relop macros =============================================== */
/* ========================================================================== */

/* 这些宏正确地处理了NaN情况 */

#define SPARSE_IS_NAN(x)	((x) != (x))
#define SPARSE_IS_ZERO(x)	((x) == 0.)
#define SPARSE_IS_NONZERO(x)	((x) != 0.)
#define SPARSE_IS_LT_ZERO(x)	((x) < 0.)
#define SPARSE_IS_GT_ZERO(x)	((x) > 0.)
#define SPARSE_IS_LE_ZERO(x)	((x) <= 0.)



#endif
