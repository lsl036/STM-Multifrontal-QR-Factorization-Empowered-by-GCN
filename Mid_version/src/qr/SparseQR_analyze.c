/******************************************************************************
 * VERSION: 1.1
 * DATE:    2020年9月27日
 * FILE:    SparseQR_analyze.c
 * BRIEF:   符号分解
 * FUNCTION: analyze函数，用来计算矩阵的填充简化排序，符号分解和预估数值分解计算量
 *****************************************************************************/

/*给定一个稀疏m×n矩阵a的非零模式，分析它以进行后续的数值分解。
    这个函数只作用于A的分布模式;它不需要实例化。*/

#include "SparseQR.h"
#include <stdbool.h>
// =============================================================================

#define FREE_WORK \
    SparseCore_free_factor (&Sc, cc) ; \
    SparseCore_free (2*(nf+1), sizeof (double), Flops,         cc) ; \
    SparseCore_free (ns+2,     sizeof (Long),    Stack_stack,   cc) ; \
    SparseCore_free (nf,       sizeof (Long),    Rh,            cc) ; \
    SparseCore_free (ntasks,   sizeof (Long),    TaskParent,    cc) ;

// =============================================================================
// === qr_analyze ============================================================
// =============================================================================

qr_symbolic *qr_analyze
(
    
    sparse_csc *A,
    int ordering,          
    Long *Quser,            
    int do_rank_detection,  
    sparse_common *cc
)
{
    qr_symbolic *QRsym ;
    Long *Parent, *Child, *Childp, *W, *Rj, *Rp, *Super, *Stair, *Fmap, *Sleft,
        *Post, *Ap, *Weight, *On_stack, *Task, *TaskParent,
        *TaskChildp, *TaskChild, *Fm, *Cm, *TaskFront, *TaskFrontp, *Rh,
        *Stack_stack, *Stack_maxstack, *Hip,
        *TaskStack, *InvPost ;
    Long nf, f, j, col1, col2, p, p1, p2, t, parent, anz, fp, csize_max,
        fmc, fnc, fpc, cm, cn, ci, fm, fn, cm_min, cm_max, csize_min, kf,
        rm, rn, col, c, pc, rsize, maxfn, csize, m, n, klast,
        stack, maxstack, rxsize, hisize,
        rhxsize, ctot, fsize, ns, ntasks, task ;
    sparse_csc *AT ;
    sparse_factor *Sc ;
    int ok = TRUE, do_parallel_analysis ;
    double total_flops = 0 ;
    double *Flops, *Flops_subtree ;
    Long *Sp, *Sj;

    // -------------------------------------------------------------------------
    // 获取输入
    // -------------------------------------------------------------------------

    m = A->nrow ;  // 行号
    n = A->ncol ;  // 列号
    Ap = (Long *) A->p ; //  Ap [0...n], 列偏移指针
    anz = Ap [n] ;       // Ap[n] 为非零元素总数

    // 并行情况。 如果 SPQR_grain <= 1 或者总计算量低于 SPQR_small 则用串行,
    // 合适的cc->SPQR_grain值大约是内核数量的2倍。
    
    do_parallel_analysis = (cc->SPQR_grain > 1) ; 

    // -------------------------------------------------------------------------
    // 给这两个函数 qr_analyze 和 SparseChol_postorder 分配内存空间
    // -------------------------------------------------------------------------

    nf = n ;    

    SparseCore_allocate_work (n+1, MAX (m, 2*(n+1) + 2*(nf+2)) + 1, 0, cc) ;

    // 工作空间之后分配
    Rh = NULL ;
    Flops = NULL ;
    Flops_subtree = NULL ;
    Sc = NULL ;
    Stack_stack = NULL ;
    TaskParent = NULL ;
    ns = 0 ;                    
    ntasks = 0 ;               

    if (cc->status < SPARSE_OK) // 内存检查
    {
        // 内存溢出
        FREE_WORK ;
        return (NULL) ;
    }

    // ----------------------------------------------------
    // A'A的超节点Cholesky排序和分析    1.2
    // ----------------------------------------------------

    AT = SparseCore_transpose (A, 0, cc) ;   // AT = 转置 (A') 

    // 保存当前设置
    Long save [6] ;
    save [0] = cc->supernodal ; // SPARSE_AUTO
    save [1] = cc->nmethods ;  // 0, 默认策略
    save [2] = cc->postorder ; // TRUE
    save [3] = cc->method [0].ordering ; // SPARSE_GIVEN
    save [4] = cc->method [1].ordering ; // SPARSE_AMD
    save [5] = cc->method [2].ordering ; // 默认为 SPARSE_AMD

    // 采用supernode 分析 来寻找波前矩阵
    cc->supernodal = SPARSE_SUPERNODAL ;

    /* -------------------------------------------------------------- */
    /*  确定排序方法，这里是GIVEN(无单例，对A) 或者 FIXED(有单例，对Y) */
    /* -------------------------------------------------------------- */

    if (ordering == QR_ORDERING_NATURAL || (ordering == QR_ORDERING_GIVEN && Quser == NULL))
    {
        ordering = QR_ORDERING_FIXED ;
    }

    // if (ordering == QR_ORDERING_BEST)
    // {
    //     ordering = QR_ORDERING_CHOL ;
    //     cc->postorder = TRUE ;
    //     cc->nmethods = 2 ;
    //     cc->method [0].ordering = SPARSE_COLAMD ;
    //     cc->method [1].ordering = SPARSE_AMD ;
    // }


    // if (ordering == QR_ORDERING_BESTAMD)
    // {
    //     ordering = QR_ORDERING_CHOL ;
    //     cc->postorder = TRUE ;
    //     cc->nmethods = 2 ;
    //     cc->method [0].ordering = SPARSE_COLAMD ;
    //     cc->method [1].ordering = SPARSE_AMD ;
    // }

    /* ------------------- */
    /*     设置排序方法    */
    /* ------------------ */

    if (ordering == QR_ORDERING_FIXED)
    {
        cc->nmethods = 1 ;
        cc->method [0].ordering = SPARSE_NATURAL ;
        cc->postorder = FALSE ;
        Quser = NULL ;
    }
    else if (ordering == QR_ORDERING_GIVEN) // 进到这里
    {

        cc->nmethods = 1 ;
        cc->method [0].ordering = SPARSE_GIVEN ;
        cc->postorder = FALSE ;
    }
    // else if (ordering == QR_ORDERING_CHOL)
    // {

    //     Quser = NULL ;
    // }
    // else if (ordering == QR_ORDERING_AMD)
    // {

    //     cc->nmethods = 1 ;
    //     cc->method [0].ordering = SPARSE_AMD ;
    //     cc->postorder = TRUE ;
    //     Quser = NULL ;
    // }
    // else 
    // {
    //     ordering = QR_ORDERING_COLAMD ;
    //     cc->nmethods = 1 ;
    //     cc->method [0].ordering = SPARSE_COLAMD ;
    //     cc->postorder = TRUE ;
    //     Quser = NULL ;
    // }


    // 多波前QR 排序 与分析
    Sc = SparseChol_analyze_p2 (SPARSE_ANALYZE_FOR_SPQR,
        AT, (Sparse_long *) Quser, NULL, 0, cc) ;

    // 记录实际使用的排序
    if (Sc != NULL)
    {
        switch (Sc->ordering)
        {
            case SPARSE_NATURAL: ordering = QR_ORDERING_NATURAL ; break ;
            case SPARSE_GIVEN:   ordering = QR_ORDERING_GIVEN   ; break ;
            case SPARSE_AMD:     ordering = QR_ORDERING_AMD     ; break ;
            case SPARSE_COLAMD:  ordering = QR_ORDERING_COLAMD  ; break ;
        }
    }

    cc->SPQR_istat [7] = ordering ;
    // printf("SparseChol_analyze_p2 ordering = %d\n", ordering);  // p2使用的排序方法

    // 恢复HNUCHOL设置
    cc->supernodal              = save [0] ;
    cc->nmethods                = save [1] ;
    cc->postorder               = save [2] ;
    cc->method [0].ordering     = save [3] ;
    cc->method [1].ordering     = save [4] ;
    cc->method [2].ordering     = save [5] ;

    SparseCore_free_sparse (&AT, cc) ;      

    if (cc->status < SPARSE_OK) // 检查内存
    {
        // 越界
        FREE_WORK ;
        return (NULL) ;
    }

    if (Sc == NULL || !(Sc->is_super) || !(Sc->is_ll)) //释放空间
    {
        SparseCore_free_factor (&Sc, cc) ;
        FREE_WORK ;

        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // extract the contents of HNUCHOL's supernodal factorization
    // 提取 HNUCHOL 超节点分解的内容 1.3.2
    // -------------------------------------------------------------------------

    QRsym = (qr_symbolic *) SparseCore_malloc (1, sizeof (qr_symbolic), cc) ;

    if (cc->status < SPARSE_OK)
    {
        // 内存越界
        FREE_WORK ;
        return (NULL) ;
    }

    QRsym->m = m ;
    QRsym->n = n ;
    QRsym->do_rank_detection = do_rank_detection ;
    QRsym->anz = anz ;
    QRsym->nf = nf = Sc->nsuper ;       // 超节点数目
    QRsym->rjsize = Sc->ssize ;         // 超节点R的int部分的大小
    QRsym->keepH = TRUE ;

    QRsym->Qfill = (Long *) Sc->Perm ;           
    Sc->Perm = NULL ;

    QRsym->Super = Super = (Long *) Sc->super ;  // 超级是大小nf+1
    Sc->super = NULL ;

    QRsym->Rp = Rp = (Long *) Sc->pi ;           
    Sc->pi = NULL ;

    QRsym->Rj = Rj = (Long *) Sc->s ;            
    Sc->s = NULL ;

    Sc->ColCount = NULL ;

    SparseCore_free_factor (&Sc, cc) ;

    // -------------------------------------------------------------------------
    // 其余属性赋值
    // -------------------------------------------------------------------------
    // 波前矩阵 f 的父节点是 Parent [f]，或者 为空 如果f =nf
    // f 的孩子 可以通过数组获取： Child [Childp [f] ... Childp [f+1]-1].
    QRsym->Parent = Parent = (Long *) SparseCore_malloc (nf+1, sizeof(Long), cc);
    QRsym->Childp = Childp = (Long *) SparseCore_calloc (nf+2, sizeof(Long), cc);
    QRsym->Child  = Child  = (Long *) SparseCore_calloc (nf+1, sizeof(Long), cc);
    
    // 波前树的后序， f=Post[k] 给出了后序遍历树的第k个节点
    QRsym->Post   = Post   = (Long *) SparseCore_malloc (nf+1, sizeof(Long), cc);
    // 一个行置换操作。 PLinv [i] = k 如果矩阵A的第i行是矩阵S的第k行.
    QRsym->PLinv           = (Long *) SparseCore_malloc (m,    sizeof(Long), cc);
    // Sleft[j] ... Sleft[j+1]-1 保存矩阵S中最左列为第j列的行索引
    QRsym->Sleft  = Sleft  = (Long *) SparseCore_malloc (n+2,  sizeof(Long), cc);
    // CSC 格式下。 Sp为非零元行号， Sj为列偏移。 非零元anz = Sp[m]
    QRsym->Sp     = Sp     = (Long *) SparseCore_malloc (m+1,  sizeof(Long), cc);
    QRsym->Sj     = Sj     = (Long *) SparseCore_malloc (anz,  sizeof(Long), cc);
    // 每个波前阵的行数
    QRsym->Fm              = (Long *) SparseCore_malloc (nf+1, sizeof(Long), cc);
    // 每个波前阵中贡献块C 的行数
    QRsym->Cm              = (Long *) SparseCore_malloc (nf+1, sizeof(Long), cc);

    QRsym->Hip = Hip = (Long *) SparseCore_malloc (nf+1, sizeof (Long), cc) ;

    // 稍后分配(如果没有并行性则跳过) 对栈和子节点赋值
    QRsym->TaskChild  = NULL ;
    QRsym->TaskChildp = NULL ;
    QRsym->TaskFront = NULL ;
    QRsym->TaskFrontp = NULL ;
    QRsym->TaskStack = NULL ;
    QRsym->On_stack = NULL ;
    QRsym->Stack_maxstack = NULL ;
    QRsym->ntasks = EMPTY ;         
    QRsym->ns = EMPTY ;

    // -------------------------------------------------------------------------
    // 分配并行分析所需的工作区
    // -------------------------------------------------------------------------

    if (do_parallel_analysis) 
    {
        // 计算浮点计算量的空间
        Flops = (double *) SparseCore_malloc (2*(nf+1), sizeof (double), cc) ;
        Flops_subtree = Flops + (nf+1) ;
        // 开辟保存 Household 的空间
        Rh = (Long *) SparseCore_malloc (nf, sizeof (Long), cc) ;

    }

    if (cc->status < SPARSE_OK)
    {

        qr_freesym (&QRsym, cc) ;
        FREE_WORK ;
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // 构造超节点树 及 其后序排序分布 1.3.3
    // -------------------------------------------------------------------------
    /*  在cholmod_super_symbolic中松弛超节点的Sparent与parent相同，但没有保留。
        在这里需要重新构建。未来工作:告诉CHOLMOD要保留Sparent;可以删除其中一些代码。 
        注意，cholmod_analyze 按照 R 行数递增增的顺序对子节点进行 后序排列,
        因此较大的子元素在后面，波前阵f的最大子元素是波前阵f-1。
    */

    W = (Long *) cc->Iwork ;

    // 记录每一列所属于的 超节点f
    for (f = 0 ; f < nf ; f++)
    {
        
        col1 = Super [f] ;
        col2 = Super [f+1] ;
        for (j = col1 ; j < col2 ; j++) 
        {
            W [j] = f ;  // 表示第 j 列属于超节点f 
        }
    }

    //  寻找父节点
    for (f = 0 ; f < nf ; f++)
    {
        col1 = Super [f] ;          //记录一个f 的起始列号
        col2 = Super [f+1] ;        //记录一个f 的最终列号
        p1 = Rp [f] ;                  
        p2 = Rp [f+1] ;
        fp = col2 - col1 ;          // 列数
        p = p1 + fp ;                   
        if (p < p2)                 // Rj细节 不清楚
        {
            j = Rj [p] ;               
            
            parent = W [j] ;            
        }
        else
        {
            parent = nf ;              
        }
        Parent [f] = parent ;    // 记录下父子关系
        Childp [parent]++ ;             
    }
    Parent [nf] = EMPTY ;              


    // -------------------------------------------------------------------------
    //  把波前树 进行后序排列。
    // -------------------------------------------------------------------------

    // 使用 W (2*(nf+1):3*(nf+1)-1) 为节点的权重 
    Weight = W + 2*(nf+1) ;

    for (f = 0 ; f < nf ; f++)
    {
        
        Weight [f] = (Rp [f+1] - Rp [f]) ;

    }
    Weight [nf] = 1 ;                   

    // 确保HNUCHOL不会重新分配工作空间
    cc->no_workspace_reallocate = TRUE ;

    SparseChol_postorder ((Sparse_long *) Parent, nf+1,
        (Sparse_long *) Weight, (Sparse_long *) Post, cc) ;

    cc->no_workspace_reallocate = FALSE ;

    // -------------------------------------------------------------------------
    // 构造每个波前阵的子节点列表
    // -------------------------------------------------------------------------

    // 用 cumsum 替代 Childp
    // 在输入时，Childp [0:n-1]包含计数。
    // 在输出时，Childp [k]被替换为输入计数Childp[0:k-1]的和。
    qr_cumsum (nf+1, Childp) ;

    for (kf = 0 ; kf < nf ; kf++)
    {
        c = Post [kf] ;
        
        parent = Parent [c] ;
        Child [Childp [parent]++] = c ;
    }

    qr_shift (nf+1, Childp) ;
    

    // -------------------------------------------------------------------------
    // S = A (P,Q) in row form, where P = leftmost column sort of A(:,Q)
    // 在行形式上 S = A (P,Q) , 其中 P =  A(:,Q) 基于最左列的排序 1.3.4
    // -------------------------------------------------------------------------

    //使用 W [0:m-1] 工作区 :
    qr_stranspose1 (A, QRsym->Qfill, QRsym->Sp, QRsym->Sj, QRsym->PLinv,
        Sleft, W) ;

    // -------------------------------------------------------------------------
    // 确定浮点计算数， 波前阵规模 和 连续内存使用情况   1.3.5
    // -------------------------------------------------------------------------

    Fm = QRsym->Fm ;
    Cm = QRsym->Cm ;

    Stair = W ;            
    Fmap = Stair + (n+1) ; 

    stack = 0 ;            
    maxstack = 0 ;         

    rxsize = 0 ;           
    rhxsize = 0 ;           
    maxfn = 0 ;            

    

    for (kf = 0 ; ok && kf < nf ; kf++)  // 计算波前阵的规模
    {
        f = Post [kf] ;
        
        // ---------------------------------------------------------------------
        //  分析波前阵F的分解
        // ---------------------------------------------------------------------

        col1 = Super [f] ;          
        col2 = Super [f+1] ;
        p1 = Rp [f] ;              
        p2 = Rp [f+1] ;
        fp = col2 - col1 ;         
        fn = p2 - p1 ;              
        maxfn = MAX (maxfn, fn) ;
        
        // ---------------------------------------------------------------------
        // 为波前阵F创建Fmap
        // ---------------------------------------------------------------------

        for (p = p1, j = 0 ; p < p2 ; p++, j++)
        {
            col = Rj [p] ;             
            Fmap [col] = j ;
            
        }

        // 如果do_rank_detection是假的，那么没有可能的主列故障;在这种情况下，分析是精确的。
        // 如果do_rank_detection为真，则考虑最坏情况的主列故障;在这种情况下，分析给出了一个上限。
        // ---------------------------------------------------------------------
        //  初始化 波前阵F 的staircase
        // ---------------------------------------------------------------------

        //  用S的原始行初始化staircase
        for (j = 0 ; j < fp ; j++)
        {
            col = j + col1 ;
            Stair [j] = Sleft [col+1] - Sleft [col] ;

        }

        // 来自子节点的贡献块将在这里添加
        for ( ; j < fn ; j++)
        {
            Stair [j] = 0 ;
        }

        // ---------------------------------------------------------------------
        // 组装每一个子节点
        // ---------------------------------------------------------------------

        ctot = 0 ;
        for (p = Childp [f] ; p < Childp [f+1] ; p++)
        {
            c = Child [p] ;                
            
            pc = Rp [c] ;
            fnc = Rp [c+1] - pc ;           
            fpc = Super [c+1] - Super [c] ; 
            cn = fnc - fpc ;                
            fmc = Fm [c] ;                 

            if (do_rank_detection)
            {
                cm = MIN (fmc, cn) ;        
            }
            else
            {

                Long rc = MIN (fmc, fpc) ;  
                cm = MAX (fmc - rc, 0) ;
                cm = MIN (cm, cn) ;        
            }
            

            for (ci = 0 ; ci < cm ; ci++)
            {
                col = Rj [pc + fpc + ci] ;  
                j = Fmap [col] ;            

                Stair [j]++ ;              

            }

            /*
                跟踪所有孩子的C块的总大小。这个子元素的C块最多有cm行，并且总是有(fnc-fpc)列。
                不能发生长溢出，因为子元素的csize < fsize，并且已经检查了fsize。
            */
            csize = cm*(cm+1)/2 + cm*(cn-cm) ;
            ctot += csize ;

        }

        // ---------------------------------------------------------------------
        // 将Stair替换为cumsum (Stair)，并查找F的# rows
        // ---------------------------------------------------------------------

        fm = 0 ;
        for (j = 0 ; j < fn ; j++)
        {
            fm += Stair [j] ;
            Stair [j] = fm ;
            
        }
        

        // ---------------------------------------------------------------------
        // 确定大小F的上界
        // ---------------------------------------------------------------------

        // 波前阵 F 最多 fm-by-fn, 有fp个主列

        fsize = qr_mult (fm, fn, &ok) ;
        
        if (!ok)
        {
            // 问题太大
            break ;
        }

        Fm [f] = fm ;       

        // ---------------------------------------------------------------------
        // 如果不保留H，确定R规模的上界
        // ---------------------------------------------------------------------

        rm = MIN (fm, fp) ;     
        rn = fn ;              
        
        rsize = rm*(rm+1)/2 + rm*(rn-rm) ;  
        
        rxsize += rsize ;

        // ---------------------------------------------------------------------
        // 确定 C 规模的上界
        // ---------------------------------------------------------------------

        cn = fn - fp ;         

        cm_max = fm ;
        cm_max = MIN (cm_max, cn) ;        

        // 没有 关键 故障
        cm_min = MAX (fm - rm, 0) ;
        cm_min = MIN (cm_min, cn) ;         

        // 不能出现长溢出:
        csize_max = cm_max*(cm_max+1)/2 + cm_max*(cn-cm_max) ;
        csize_min = cm_min*(cm_min+1)/2 + cm_min*(cn-cm_min) ;
        csize = do_rank_detection ? csize_max : csize_min ;

        Cm [f] = do_rank_detection ? cm_max : cm_min ;

        

        // ---------------------------------------------------------------------
        // 如果保留H，则确定浮点数计算和 R+H 规模的上限
        // ---------------------------------------------------------------------

        double fflops = 0 ;                 
        Long rhsize = 0 ;                   
        for (j = 0 ; j < fn ; j++)
        {
            t = MAX (j+1, Stair [j]) ;      
            t = MIN (t, fm) ;              
            
            rhsize += t ;                  
            if (t > j)
            {
                double h = (t-j) ;          
                fflops += 3*h               
                    + 4*h*(fn-j-1) ;        
            }
        }

        rhsize -= csize_min ;               

        if ( do_parallel_analysis)
        {
            Rh [f] = rhsize ;
        }

        rhxsize += rhsize ;

        

        // ---------------------------------------------------------------------
        // 找到节点 f 和所有他后代的 浮点计算量
        // ---------------------------------------------------------------------

        total_flops += fflops ; // 从孩子节点算起， total_flops保存从子节点累加上来的浮点数计算量
        if (do_parallel_analysis)
        {
            Flops [f] = fflops ;            // 保存本front的浮点数计算量
            for (p = Childp [f] ; p < Childp [f+1] ; p++)
            {
                fflops += Flops_subtree [Child [p]] ;
            }
            Flops_subtree [f] = fflops ;   // 消去树上的计算量保存为自己和孩子的累计和
        }

        // ---------------------------------------------------------------------
        // 预估连续情况下的 栈使用
        // ---------------------------------------------------------------------

        stack = qr_add (stack, fsize, &ok) ;
        maxstack = MAX (maxstack, stack) ;

        stack -= ctot ;

        stack = qr_add (stack, csize, &ok) ;
        maxstack = MAX (maxstack, stack) ;

        stack -= fsize ;

        stack += rhsize ;

    }

    // -------------------------------------------------------------------------
    // 计算Hip
    // -------------------------------------------------------------------------

    hisize = 0 ;           

        for (f = 0 ; f < nf ; f++)
        {

            fm = Fm [f] ;          
            Hip [f] = hisize ;
            hisize = qr_add (hisize, fm, &ok) ;
        }
        Hip [nf] = hisize ;


    // -------------------------------------------------------------------------
    // 检查越界
    // -------------------------------------------------------------------------

    if (!ok)
    {
        qr_freesym (&QRsym, cc) ;
        FREE_WORK ;
        return (NULL) ;
    }
    
    // -------------------------------------------------------------------------
    // 最终确定浮点运算次数 ， 如果问题太小关闭并行
    // -------------------------------------------------------------------------
    //printf("flag1:do_parallel_analysis= %d \n" ,do_parallel_analysis);
    if (do_parallel_analysis)
    {
        Flops [nf] = 0 ;
        Flops_subtree [nf] = total_flops ; // 根节点，总浮点计算数
        //printf("whole tree total_flops= %lf \n" , total_flops);
        if (total_flops < cc->SPQR_small)         
        {
            do_parallel_analysis = FALSE ;
        }
    }

    // -------------------------------------------------------------------------
    // 串行
    // -------------------------------------------------------------------------


    QRsym->maxstack = maxstack ;    
    QRsym->hisize = hisize ;       
    QRsym->maxfn = maxfn ;         

    cc->SPQR_flopcount_bound = total_flops ;   

    cc->SPQR_istat [0] = rxsize ;              
    cc->SPQR_istat [1] = rhxsize - rxsize ;    
    cc->SPQR_istat [2] = nf ;                 
    cc->SPQR_istat [3] = 1 ;                   

    // -------------------------------------------------------------------------
    // 如果没有并行分析，快速返回
    // -------------------------------------------------------------------------
    if (!do_parallel_analysis)
    {
        QRsym->ntasks = 1 ;
        QRsym->ns = 1 ;
        FREE_WORK ;
        return (QRsym) ; // 单任务串行这里就返回了
    }

    // =========================================================================
    // ===  并行分析    ======================================
    // =========================================================================

    // 如果 f 的子树 计算量 > big_flops, 则波前阵f 是 大 的

    double big_flops = total_flops / cc->SPQR_grain ;

    big_flops = MAX (big_flops, cc->SPQR_small) ;

    double small_flops = big_flops / 32 ;
    // -------------------------------------------------------------------------
    //  为每个节点初始化任务分配
    // -------------------------------------------------------------------------

    //  Stair, Fmap 不再需要
    //  还需要Fm、Cm、Rh、Flops和Flops_subtree。

    // 由于不再需要Stair和Fmap，我们可以使用该工作区来执行任务和InvPost。Fm和Cm在W[0:…]仍然需要。
    Task = W ;                     
    InvPost = Task + (nf+1) ;       

    for (Long k = 0 ; k <= nf ; k++)
    {
        f = Post [k] ;
        InvPost [f] = k ;
    }



#define PENDING (-2)
#define TASK_IS_PENDING(f) (Task [f] <= PENDING)
#define SMALL_TASK(f) (TASK_IS_PENDING (f) ? ((- Task [f]) - 3) : EMPTY) 
#define ASSIGN_SMALL_TASK(f,task) { Task [f] = -(task) - 3 ; }

    for (f = 0 ; f < nf ; f++)
    {

        // 一个大节点被标记为挂起，并且不会在第一次分配。
        // Task [f] = PENDING意味着以f为根的子树中有大于big_flops的工作。
        Task [f] = (Flops_subtree [f] > big_flops) ? PENDING : EMPTY ;

    }
    Task [nf] = PENDING ; // 占位符节点必定是大节点


    // -------------------------------------------------------------------------
    // 为任务分配小的子树
    // -------------------------------------------------------------------------

    ntasks = 0 ;

    for (kf = 0 ; kf < nf ; kf++)
    {
        Long fstart = Post [kf] ;
        if (Task [fstart] != EMPTY)
        {
            continue ;
        }


        // ---------------------------------------------------------------------
        // 在根节点遍历f，在一个挂起的(大的)节点停止
        // ---------------------------------------------------------------------

        for (f = fstart ; Task [Parent [f]] == EMPTY ; f = Parent [f])
        {
            
            ;
        }

        Long flast = f ;
        parent = Parent [flast] ;


        // ---------------------------------------------------------------------
        // 启动一个新任务，或与之前的同级任务合并
        // ---------------------------------------------------------------------

        double flops = Flops_subtree [flast] ;
        if (SMALL_TASK (parent) == EMPTY)
        {
            // 悬挂的父节点 flast 没有子节点小任务，分配一个新task
            task = ntasks++ ;
            if (flops <= small_flops)
            {
                // 这个节点本身是小的，分配它
                ASSIGN_SMALL_TASK (parent, task) ;
                Flops_subtree [parent] = flops ;
            }
        }
        else
        {
            
            task = SMALL_TASK (parent) ;
            Flops_subtree [parent] += flops ;
            
            if (flops > small_flops)
            {
                
                ASSIGN_SMALL_TASK (parent, task) ;
                Flops_subtree [parent] = 0 ;
                
            }
        }

        // ---------------------------------------------------------------------
        // 将以flast为根的子树中的所有节点放置到当前任务中
        // ---------------------------------------------------------------------

        klast = InvPost [flast] ;
        for (Long k = kf ; k <= klast ; k++)
        {
            f = Post [k] ;
            Task [f] = task ;
        }
    }

    // -------------------------------------------------------------------------
    // 完成对大节点的任务分配
    // -------------------------------------------------------------------------

    for (f = 0 ; f < nf ; f++)    // 任务划分 
    {
        if (TASK_IS_PENDING (f))
        {
            // f 没有子任务，或者如果它的子节点分配给不同任务，
            // 则创建一个 f 分配的新任务。
            if ((Childp [f+1] - Childp [f]) == 0)
            {
                task = ntasks++ ;
                
            }
            else
            {
                task = EMPTY ;
                for (p = Childp [f] ; p < Childp [f+1] ; p++)
                {
                    c = Child [p] ;
                    
                    if (task == EMPTY)
                    {
                        // 第一个孩子的任务 Task[c]
                        task = Task [c] ;
                    }
                    else if (Task [c] != task)   //孩子所在的任务不一样
                    {

                        task = ntasks++ ;
                        
                        break ;
                    }
                }
            }
            
            Task [f] = task ;
        }
    }

    // -------------------------------------------------------------------------
    // 如果只有一个真正的任务找到 则 快速返回
    // -------------------------------------------------------------------------
    if (ntasks <= 1)
    {
        // 没有并行 found
        QRsym->ntasks = 1 ;
        QRsym->ns = 1 ;
        FREE_WORK ;
        return (QRsym) ;
    }

    // 接下来是并行的分析
    // -------------------------------------------------------------------------
    //  打包任务调度
    // -------------------------------------------------------------------------

    Task [nf] = ntasks++ ;      // 赋予一个统一的父节点

    // 浮点统计不再需要
    SparseCore_free (2*(nf+1), sizeof (double), Flops, cc) ;
    Flops = NULL ;
    Flops_subtree = NULL ;

    // 在波前树中，每一个波前阵对应一个节点(0~nf-1),加上一个额外的占位节点nf作为树的根
    // 所以一共有nf+1个节点。 注意一棵树可以有nf=0(矩阵含有行m=0)。每个节点最多一个
    // 任务，包括占位符
    QRsym->ntasks = ntasks ;
    cc->SPQR_istat [3] = ntasks ;

    // -------------------------------------------------------------------------
    // 为并行分析分配空间
    // -------------------------------------------------------------------------

    // 在分解过程中，只需获得任务中的第一个f的On_stack [f]。

    //TaskParent是临时工作区:
    TaskParent = (Long *) SparseCore_malloc (ntasks,   sizeof (Long), cc) ;

    TaskChildp = (Long *) SparseCore_calloc (ntasks+2, sizeof (Long), cc) ;
    TaskChild  = (Long *) SparseCore_calloc (ntasks+1, sizeof (Long), cc) ;
    TaskFront  = (Long *) SparseCore_malloc (nf+1,     sizeof (Long), cc) ;
    TaskFrontp = (Long *) SparseCore_calloc (ntasks+2, sizeof (Long), cc) ;
    TaskStack  = (Long *) SparseCore_malloc (ntasks+1, sizeof (Long), cc) ;
    On_stack   = (Long *) SparseCore_malloc (nf+1,     sizeof (Long), cc) ;

    // 每个任务的波前矩阵列表
    QRsym->TaskFront  = TaskFront ;  
    QRsym->TaskFrontp = TaskFrontp ;
    // 任务树
    QRsym->TaskChild  = TaskChild ;
    QRsym->TaskChildp = TaskChildp ;
    QRsym->TaskStack  = TaskStack ;
    // On_stack[f] 保存front f
    QRsym->On_stack   = On_stack ;

    if (cc->status < SPARSE_OK) // 空间检查
    {
        qr_freesym (&QRsym, cc) ;
        FREE_WORK ;
        return (NULL) ;
    }
    
    // -------------------------------------------------------------------------
    //  构建任务树
    // -------------------------------------------------------------------------

    for (task = 0 ; task < ntasks ; task++)
    {
        TaskParent [task] = EMPTY ;
    }

    for (f = 0 ; f < nf ; f++)
    {
        Long my_task = Task [f] ;
        Long parent_task = Task [Parent [f]] ;
        
        if (my_task != parent_task)
        {
            TaskParent [my_task] = parent_task ;
        }
    }


    
    // 统计每一个task的子 task 数目
    for (task = 0 ; task < ntasks ; task++)
    {
        parent = TaskParent [task] ;
        if (parent != EMPTY)
        {
            TaskChildp [parent]++ ;
        }
    }

    //构建 TaskChildp 指针
    qr_cumsum (ntasks, TaskChildp) ;

    // 创建孩子清单
    for (Long child_task = 0 ; child_task < ntasks ; child_task++)
    {
        parent = TaskParent [child_task] ;
        if (parent != EMPTY)
        {
            TaskChild [TaskChildp [parent]++] =  child_task ;
        }
    }

    qr_shift (ntasks, TaskChildp) ;

    // -------------------------------------------------------------------------
    // 为每项任务列出波前阵列表
    // -------------------------------------------------------------------------

    // 计算每个任务的fronts数目
    for (f = 0 ; f < nf ; f++)
    {
        TaskFrontp [Task [f]]++ ;
    }

    // 构建 TaskFrontp 指针
    qr_cumsum (ntasks+1, TaskFrontp) ;

    for (kf = 0 ; kf < nf ; kf++)
    {
        f = Post [kf] ;
        task = Task [f] ;
        TaskFront [TaskFrontp [task]++] = f ;
    }

    qr_shift (ntasks+1, TaskFrontp) ;


    // -------------------------------------------------------------------------
    //  每个task分一个stack
    // -------------------------------------------------------------------------


    for (task = 0 ; task < ntasks ; task++)
    {
        TaskStack [task] = EMPTY ;
    }

    for (Long task_start = 0 ; task_start < ntasks ; task_start++)
    {
        if (TaskStack [task_start] == EMPTY)
        {
            Long s = ns++ ;
            for (task = task_start ;
                task != EMPTY && TaskStack [task] == EMPTY ;
                task = TaskParent [task])
            {
                TaskStack [task] = s ;
            }
        }
        
    }

    QRsym->ns = ns ;

    for (task = 0 ; task < ntasks ; task++)
    {
        Long s = TaskStack [task] ;

        for (p = TaskFrontp [task] ; p < TaskFrontp [task+1] ; p++)
        {
            f = TaskFront [p] ;

            
            On_stack [f] = s ;
        }
    }

    // -------------------------------------------------------------------------
    // 分配QRsym的其余部分
    // -------------------------------------------------------------------------

    // 临时工作空间:
    // Stack_stack (s): 当前栈使用
    Stack_stack = (Long *) SparseCore_calloc (ns+2, sizeof (Long), cc) ;

    // 永久部分:
    // Stack_maxstack (s): 最大堆栈使用，如果H没有保持
    Stack_maxstack = (Long *) SparseCore_calloc (ns+2, sizeof (Long), cc) ;

    QRsym->Stack_maxstack = Stack_maxstack ;

    if (cc->status < SPARSE_OK)
    {
        qr_freesym (&QRsym, cc) ;
        FREE_WORK ;
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // 找到栈的大小
    // -------------------------------------------------------------------------

    for (kf = 0 ; kf < nf ; kf++)
    {

        // ---------------------------------------------------------------------
        // 分析堆栈 s 上波前矩阵F的 分解
        // ---------------------------------------------------------------------

        f = Post [kf] ;
        Long s = On_stack [f] ;

        fp = Super [f+1] - Super[f] ;
        fn = Rp [f+1] - Rp [f] ;    

        

        // ---------------------------------------------------------------------
        // 组装每个子堆栈(忽略其他堆栈中的子堆栈)
        // ---------------------------------------------------------------------

        ctot = 0 ;
        for (p = Childp [f] ; p < Childp [f+1] ; p++)
        {
            c = Child [p] ;                 
            if (On_stack [c] != s)
            {
                
                continue ;
            }
            fnc = Rp [c+1] - Rp [c] ;       
            fpc = Super [c+1] - Super [c] ; 
            cn = fnc - fpc ;            
            cm = Cm [c] ;                 

            csize = cm*(cm+1)/2 + cm*(cn-cm) ;
            ctot += csize ;
            
        }

        // ---------------------------------------------------------------------
        // 确定F大小的上界
        // ---------------------------------------------------------------------

        fm = Fm [f] ;                      
        fsize = fm * fn ;

        // ---------------------------------------------------------------------
        //  确定R大小的上界
        // ---------------------------------------------------------------------

        rm = MIN (fm, fp) ;     
        rn = fn ;               
        
        rsize = rm*(rm+1)/2 + rm*(rn-rm) ;

        // ---------------------------------------------------------------------
        // 确定C大小的上界
        // ---------------------------------------------------------------------

        cn = fn - fp ;         

        // 有秩检测和可能的主列故障
        cm_max = fm ;
        cm_max = MIN (cm_max, cn) ;        

        // 没有主列故障
        cm_min = MAX (fm - rm, 0) ;
        cm_min = MIN (cm_min, cn) ;         

        // 不能出现长溢出:
        csize_max = cm_max*(cm_max+1)/2 + cm_max*(cn-cm_max) ;
        csize_min = cm_min*(cm_min+1)/2 + cm_min*(cn-cm_min) ;
        csize = do_rank_detection ? csize_max : csize_min ;

        // ---------------------------------------------------------------------
        // 预估栈的使用情况
        // ---------------------------------------------------------------------

        Long ss = Stack_stack [s] ;  
        Long sm = Stack_maxstack [s] ; 

        // 在栈s上分配波前阵F
        ss += fsize ;
        sm = MAX (sm, ss) ;

        // 组装和删除s上的子元素，然后分解F
        ss -= ctot ;
        

        // 为波前阵F创建C 
        ss += csize ;
        sm = MAX (sm, ss) ;

        // 删除F
        ss -= fsize ;

        ss += Rh [f] ;

        Stack_stack [s] = ss ;      
        Stack_maxstack [s] = sm ;   
    }


    // -------------------------------------------------------------------------
    // 释放空间， 返回结果
    // -------------------------------------------------------------------------

    FREE_WORK ;
    return (QRsym) ;
} // end of analyze

void qr_stranspose1( sparse_csc *A, Long *Qfill, Long *Sp, Long *Sj, 
                        Long *PLinv, Long *Sleft, Long *W )
{
    Long i, j, p, pend, t, k, row, col, kstart, s, m, n, *Ap, *Ai ;

    // -------------------------------------------------------------------------
    // 获取输入
    // -------------------------------------------------------------------------

    m = A->nrow ;
    n = A->ncol ;
    Ap = (Long *) A->p ;
    Ai = (Long *) A->i ;

    // -------------------------------------------------------------------------
    // 清除逆置换
    // -------------------------------------------------------------------------

    for (i = 0 ; i < m ; i++)
    {
        PLinv [i] = EMPTY ;            
        
    }

    // -------------------------------------------------------------------------
    // 计算A的每一行中的条目并找到PLinv
    // -------------------------------------------------------------------------

    k = 0 ;
    for (col = 0 ; col < n ; col++)    
    {
        j = Qfill ? Qfill [col] : col ;         
        
        kstart = k ;
        pend = Ap [j+1] ;
        for (p = Ap [j] ; p < pend ; p++)
        {
            i = Ai [p] ;                
            row = PLinv [i] ;           
            if (row == EMPTY)
            {
                row = k++ ;
                PLinv [i] = row ;
                W [row] = 1 ;           
            }
            else
            {
                W [row]++ ;             
            }
        }
        Sleft [col] = k - kstart ;  
    }

    // -------------------------------------------------------------------------
    // 计算 Sleft = cumsum ([0 Sleft])
    // -------------------------------------------------------------------------

    s = 0 ;
    for (col = 0 ; col < n ; col++)     
    {
        t = s ;
        s += Sleft [col] ;
        Sleft [col] = t ;
    }
    Sleft [n] = k ;
    Sleft [n+1] = m ;

    // -------------------------------------------------------------------------
    // 完成行置换，以防A的行是空的
    // -------------------------------------------------------------------------

    if (k < m)
    {
        for (i = 0 ; i < m ; i++)
        {
            if (PLinv [i] == EMPTY)     
            {
                row = k++ ;             
                PLinv [i] = row ;
                W [row] = 0 ;          
            }
        }
    }


    // -------------------------------------------------------------------------
    // 计算S的行指针(和W中的一个副本)
    // -------------------------------------------------------------------------

    p = 0 ;
    for (row = 0 ; row < m ; row++)
    {
        t = p ;
        p += W [row] ;
        W [row] = t ;                   
        Sp [row] = t ;                 
    }
    Sp [m] = p ;

    // -------------------------------------------------------------------------
    // create S=A(p,q) '，或者S=A(p,q)，如果S被认为是行形式
    // -------------------------------------------------------------------------

    for (col = 0 ; col < n ; col++)     
    {
        j = Qfill ? Qfill [col] : col ; 
        pend = Ap [j+1] ;
        for (p = Ap [j] ; p < pend ; p++)
        {
            i = Ai [p] ;                
            row = PLinv [i] ;           
            s = W [row]++ ;             
            Sj [s] = col ;
        }
    }
}
