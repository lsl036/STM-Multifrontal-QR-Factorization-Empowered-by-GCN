// =============================================================================
// === SparseQR_struct.h =======================================================
// =============================================================================


#ifndef SparseQR_STRUCT_H
#define SparseQR_STRUCT_H

#include "SparseQR_definitions.h"
#include "Sparse.h"


// =============================================================================
// === qr_symbolic ===========================================================
// =============================================================================

// 符号对象，分解期间不会改变。符号对象只依赖于输入矩阵的模式，而不依赖于它的值。 
// 结构体内容也不会因为用于秩检测的列旋转而改变。 
// 这使得并行性更容易管理，因为所有线程都可以在不同步的情况下访问这个对象。

// 其中 用户输入的矩阵A 是 m-by-n 规模的，非零元有 anz 个. nf <= MIN(m,n) 是
// 波前矩阵的数目， rnz <= nnz(R) 是表示R的超节点形式的列索引（在R的每个块首行
// 中给每个非关键列一个 Long 索引）


typedef struct qr_symbolic_struct
{

    // -------------------------------------------------------------------------
    // 输入矩阵的行格式及其排序
    // -------------------------------------------------------------------------
    /*
        符号分析时，构造S = A(P,Q)的非零模式，其中A是用户的输入矩阵。
        它的数值也被构造，但它们不作为符号对象的一部分。
        矩阵S以面向行的形式存储。S的行根据它们最左边的列索引(通过PLinv)进行排序。 
        S的每一行的列索引都严格按照升序排列，即使输入矩阵A不需要排序。
    */

    Sparse_long m, n, anz ; 

    Sparse_long *Sp ;       // S的行指针

    Sparse_long *Sj ;       // S的列索引 （ anz = Sp[m]）

    Sparse_long *Qfill ;    // 减少填充的列排序矩阵 （S 相比于 A 的置换矩阵）
                        // Qfill[k] = j 如果矩阵A的第k列是矩阵S的第j列

    Sparse_long *PLinv ;    // 将S = A(P,Q)按照最左列索引的递增顺序排列的逆行置换。
                        // PLinv [i] = k 如果矩阵A的第i行是矩阵S的第k行.

    Sparse_long *Sleft ;    // S最左列是第j列的行索引如下给出 Sleft [j] ... Sleft [j+1]-1.
            //  这也可以是空的(也就是说，Sleft [j]可以等于Sleft [j+1])。
            //   Sleft[n] 为S的非空行数，Sleft [n+1] == m。
            //  也就是说，Sleft [n] … Sleft [n+1]-1给出了S的空行，如果有的话。

    // -------------------------------------------------------------------------
    // 波前矩阵： 分布模式 和 树
    // -------------------------------------------------------------------------

    // 每一个波前矩阵 规模为 fm-by-fn ，主列共有 fnpiv 列。fn 列索引由一组大小为
    // fnpiv 的主列给出，由 Super 定义， 伴随分布模式 Rj [ Rp[f] ... Rp[f+1]-1 ].

    // 波前阵的行索引不保存， 如果Household向量不保存的话则行索引不保存。
    // 如果 Household 向量需要保存， 行索引则会在数值分解的过程中动态计算。

    Sparse_long nf ;        // 波前矩阵的数目
    Sparse_long maxfn ;     // 任意波前矩阵 的 列的最大 #

    // parent, child and childp 定义了行合并树
    Sparse_long *Parent ;   // size nf+1
    Sparse_long *Child ;    // size nf+1
    Sparse_long *Childp ;   // size nf+2

    // 波前矩阵 f 的父节点是 Parent [f]，或者 为空 如果f =nf
    // f 的孩子 可以通过数组获取： Child [Childp [f] ... Childp [f+1]-1].

    // 树中的节点 nf 是一个占位符，它不代表一个波前矩阵。 所有波前矩阵树（森林）根的
    // 父节点都是这个节点 nf。 所以树节点 0：nf 是真实的树，共有一个父节点 nf。

    Sparse_long *Super ;    // 大小为nf+1. Super[f]给出了波前阵F的第一个主列
        //  这个主列在S中，因此波前阵F的 主列数目一共有Super [f+1] - Super [f] 

    Sparse_long *Rp ;       // Rp [f+1] - Rp [f] 表示波前阵f的列数
    Sparse_long *Rj ;       // R的压缩超节点形式，大小为rjsize

    Sparse_long *Post ;     // 波前树的后序， f=Post[k] 给出了后序遍历树的第k个节点

    Sparse_long rjsize ;    

    Sparse_long do_rank_detection ; 

    // 剩下的参数取决于是否允许 秩检查
    Sparse_long maxstack  ; // 最大堆栈大小(顺序情况)
    Sparse_long hisize ;    // Hii 的大小 （定义在numeric 结构体中）

    Sparse_long keepH ;     // 为真则保存Household变换H

    Sparse_long *Hip ;      // 大小为nf+1.如果H保存，波前阵f的行索引
        //  在Hii [Hip [f] ... Hip [f] + Hm [f]] 中
        //  Hii 和Hm 保存在数值分解结构体中

        //    每一个波前阵都对应R的一个行块。 R 中fn条列的索引
        //   由 Rj [Rp [f] ... Rp [f+1]-1] 给出，其中第一个主列fp索引
        //  是 Super [f] ... Super [f+1]-1 。
        //   剩余列索引 Rj[...] 是非主列，范围是Super[f+1] ~ n.
        //     R 的行数目最大为fp，如果矩阵中出现死列则有可能更少。
        //   贡献块C 的列数 恒为 cn = fn - fp ,其中 fn = Rp [f+1] - Rp [f].

    Sparse_long ntasks ;    // 任务图中的任务数目
    Sparse_long ns ;        // 栈的数量

    // -------------------------------------------------------------------------
    // 剩下的 QR 符号对象仅在 ntasks > 1 的时候存在 ( 并行的时候起作用)
    // -------------------------------------------------------------------------

    // 任务树（节点 0：ntasks），包含占位节点
    Sparse_long *TaskChildp ;      
    Sparse_long *TaskChild ;        

    Sparse_long *TaskStack ;      

    // 每个任务的 波前阵 列表
    Sparse_long *TaskFront ;                  
    Sparse_long *TaskFrontp  ;      

    Sparse_long *On_stack  ;        //  波前阵f 在 On_stack[f] 中

    // 每个栈的大小
    Sparse_long *Stack_maxstack ;  

    // 每个波前阵的行数
    Sparse_long *Fm ;              

    // 每个波前阵中贡献块C 的行数
    Sparse_long *Cm ;              

}qr_symbolic ;


// =============================================================================
// === qr_numeric ============================================================
// =============================================================================

// 数值对象， 包含了三角/梯形 因子R，和可选的Householder 向量H
typedef struct qr_numeric_struct
{

    // -------------------------------------------------------------------------
    // Numeric R factor     数值因子R 
    // -------------------------------------------------------------------------

    double **Rblock ;    // 大小为nf。R[f] 是一个(double *)指针指向波前阵 f 中的R块
                        //  它是一个大小为 Rm(f)-by-Rn(f) 的上梯形，但仅
                        //     以列压缩格式保存上三角部分。

    double **Stacks ;   //  大小为ns。一个堆栈数组，保存 R 和 H 因子以及当前在顶部的波前矩阵F
                        //  接下来是空区域，然后底部是先前的波前阵的贡献块C。当分解完成后，只有
                        //  在栈顶的 R 和 H 部分保留。

    Sparse_long *Stack_size ;   // 大小为ns; Stack_size [s] 是 Stacks[s] 的大小

    Sparse_long hisize ;        //  Hii 的大小       

    Sparse_long n ;             // 矩阵A的规模 m-by-n
    Sparse_long m ;
    Sparse_long nf ;            //  波前矩阵 数目 nf
    Sparse_long ntasks ;        //  任务图中实际上使用的tasks
    Sparse_long ns ;            //  栈 数目
    Sparse_long maxstack ;      //  串行的 栈的大小，如果使用

    // -------------------------------------------------------------------------
    // 用于 秩检查 且 m < n 的情况
    // -------------------------------------------------------------------------

    char *Rdead ;       //   如果k是一个死主列，则Rdead[k] = 1，否则Rdead[k] = 0.
                        //   如果没有死列，这个指针为NULL。 如果 m < n. 则至少有n-m个死列


    Sparse_long rank ;      //   活的 主列数目 （rank）
    Sparse_long rank1 ;     //   A的第一个ntol列中 活的 主列数目

    Sparse_long maxfrank ;  //    任意 R 块中的最大行数

    double norm_E_fro ; //    w的2范数，死列2范数的向量

    // -------------------------------------------------------------------------
    // 保存Householder向量 的相关参数
    // -------------------------------------------------------------------------

    Sparse_long keepH ;     

    Sparse_long rjsize ;    //    HStair 和Htau的大小

    Sparse_long *HStair ;   //   大小为rjsize, staircase是每列中最后一个非零项的行索引
                        //     数组 Hstair[Rp [f] ...Rp [f+1]-1 ] 给出了波前阵 f 每一列的 staircase

    double *HTau ;       //  数组HTau [Rp [f] ... Rp [f+1]-1] 给出 波前矩阵 f 每一列的 Householder系数

    Sparse_long *Hii ;      //  大小为hisize(在符号分解中定义). H的行索引

    Sparse_long *HPinv ;    //  大小为m，一种行置换。如果A和H 的 行i 是 R的行k，则 HPinv[i] = k 
                        //    这个排序包含了 QRsym->PLinv，
                        //    以及在分解过程中通过行主元 排序 构造的置换。

    Sparse_long *Hm ;      // 大小为nf，Hm [f] 表示波前阵 f 的行数
    Sparse_long *Hr ;      // 大小为nf，Hr [f] 表示波前阵 f 中 R块的行数
    Sparse_long maxfm ;    // max (Hm [0:nf-1])

} qr_numeric;


// =============================================================================
// === SparseQR_factorization =============================================
// =============================================================================

// 矩阵 A 或[A B] 包含单例的符号分解和数值分解 的混合结构体

typedef struct SparseQR_factorization_struct
{
    // A 或 [A Binput] 在移除单例后的 QR 分解
    double tol ;        
    qr_symbolic *QRsym ;
    qr_numeric  *QRnum ;

    // 单例， 行压缩CSR格式； R 是 n1rows-by-n ,把这个格式转置就是L矩阵了
    Sparse_long *R1p ;      
    Sparse_long *R1j ;
    double *R1x ;
    Sparse_long r1nz ;      

    // 合并单例 ， 填充简化排序。
    Sparse_long *Q1fill ;
    Sparse_long *P1inv ;
    Sparse_long *HP1inv ;   //  如果n1cols == 0,则为NULL
                        //     在这种情况下，QRnum->HPinv服务于它的位置。

    //   如果满秩 QR->RANK == A->ncol ,则Rmap 和 RmapInv 是NULL
    Sparse_long *Rmap ;     //  Rmap大小为n。 如果R的第j列是第k个有效列，且k < QR->rank，
                        //         则Rmap[j]=k; 如果j是一个死列，那么 k >= QR->rank . 


    Sparse_long *RmapInv ;

    Sparse_long n1rows ;    //   [A B]的单例行数
    Sparse_long n1cols ;    //  [A B]的单例列数

    Sparse_long narows ;    //   A 的行数
    Sparse_long nacols ;    //    A 的列数
    Sparse_long rank ;      //   A的秩估计( n1rows + QRnum->rank1 )
                                
    double Ana_time;
    double Fac_time;

    int allow_tol ;    
} SparseQR_factorization;


#endif
