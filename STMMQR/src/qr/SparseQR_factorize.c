/******************************************************************************
 * SPQR, Copyright 2008-2016 by Timothy A. Davis.
 * All Rights Reserved.
 * SPQR is available under alternate licenses, contact T. Davis for details.
 * Availability:

    http://www.suitesparse.com

 * Used by permission.
 
 * ChangeLog: 1. Merged a part of QR numerical factorization functions, 
              and implemented data affinity in get_Work.
              2. Modified the parallel mode of the nodes in task tree 
              (in qr_factorize).
 *****************************************************************************/
#ifdef PRINT_TIME
#include <sys/time.h>
#endif
#include "SparseQR.h"
#include "tpsm.h"

// #define FCHUNK 64       // 32 Householder 块规模
// #define SMALL 5000
// #define MINCHUNK 8
// #define MINCHUNK_RATIO 4  // 4


#define FREE_WORK \
{ \
    free_Work (Work, ns, n, maxfn, wtsize, cc) ; \
    if (freeA) SparseCore_free_sparse (Ahandle, cc) ; \
    SparseCore_free (anz, sizeof (double), Sx, cc) ; \
    Sx = NULL ; \
    SparseCore_free (ns, sizeof (qr_work), Work, cc) ; \
    Work = NULL ; \
    SparseCore_free (nf+1, sizeof (double *), Cblock, cc) ; \
    Cblock = NULL ; \
}
// =============================
// ======  获取分块参数  =======
// =============================
/**
* @brief 读环境变量值
* @param env 环境变量名称.
* @return 环境变量值(int类型)
*/
int CHUNK_Read_Env(char *env)
{
  char *p;
  //getenv函数是gcc自带的，再stdlib.h中
  if (( p = getenv(env)))
  	return (atoi(p));
  else
	return(0);
}

int chunk_getSettings
(
    size_t const _FCHUNK_size,
    size_t const _SMALL_size,
    size_t const _MINCHUNK_size,
    size_t const _MINCHUNK_RATIO
)
{
    int edit = 0;
    // 接收是否需要调整分块参数
    // printf("Want to edit chunk settings? 1 for yes, 0 for no: ");
    // scanf("%d", &edit);

    if (edit) // 需要读环境变量输入参数
    {
        if (CHUNK_Read_Env("FCHUNK"))
            FCHUNK = CHUNK_Read_Env("FCHUNK") ;
        else
            FCHUNK = _FCHUNK_size;
        
        if (CHUNK_Read_Env("SMALL"))
            SMALL = CHUNK_Read_Env("SMALL") ;
        else
            SMALL = _SMALL_size;
        
        if (CHUNK_Read_Env("MINCHUNK"))
            MINCHUNK = CHUNK_Read_Env("MINCHUNK") ;
        else
            MINCHUNK = _MINCHUNK_size;

        if (CHUNK_Read_Env("MINCHUNK_RATIO"))
            MINCHUNK_RATIO = CHUNK_Read_Env("MINCHUNK_RATIO") ;
        else
            MINCHUNK_RATIO = _MINCHUNK_RATIO;
        // size_t t1,t2,t3,t4;
        // printf("Please input FCHUNK, SMALL, MINCHUNK and RATIO:");
        // scanf("%d %d %d %d", &t1, &t2, &t3, &t3);
        // FCHUNK = t1;
        // SMALL = t2;
        // MINCHUNK = t3;
        // MINCHUNK_RATIO = t4;
        /* code */
    } 
    else // 采用默认参数
    {
        FCHUNK = _FCHUNK_size;
        SMALL = _SMALL_size;
        MINCHUNK = _MINCHUNK_size;
        MINCHUNK_RATIO = _MINCHUNK_RATIO;
    }

    return 0;
}

// =============================================================================
// === get_Work ================================================================
// =============================================================================

// 开辟 factorize 和 kernel 中的每个堆栈都需要的工作空间
// 元素 Work[s]包含堆栈s的工作空间，对于每个堆栈0到ns-1。

qr_work *get_Work 
(
    Long ns,            // 栈的数目
    Long n,             // A的列数
    Long maxfn,         // 所有波前阵的最大列数
    Long keepH,         // 是否保留H
    Long fchunk,
    Long *p_wtsize,     // 给每个 WTwork 的 空间大小
    sparse_common *cc
)
{
    int ok = TRUE ;
    qr_work *Work ;
    Long wtsize ;
    *p_wtsize = 0 ;

    wtsize = qr_mult (fchunk + (keepH ? 0:1), maxfn, &ok) ;

    // 开辟 ns 个 qr_work 空间
    Work = (qr_work *)    
        SparseCore_malloc (ns, sizeof (qr_work), cc) ;

    //检查
    if (!ok || cc->status < SPARSE_OK)
    {
        SparseCore_free (ns, sizeof (qr_work), Work, cc) ;
        
        return (NULL) ;
    }
    #ifdef NUMA_ALLOC
        size_t numa_node_rank = 0;
    #endif
    
    for (Long stack = 0 ; stack < ns ; stack++)
    {
        #ifdef NUMA_ALLOC
            numa_node_rank = stack % 4;
            Work [stack].Fmap = (Long *) SparseCore_malloc_onnode (n, sizeof (Long), numa_node_rank, cc) ;
            Work [stack].Cmap = (Long *) SparseCore_malloc_onnode (maxfn, sizeof(Long), numa_node_rank, cc);
            Work [stack].WTwork =
            (double *) SparseCore_malloc_onnode (wtsize, sizeof (double), numa_node_rank, cc) ;
        #else
            Work [stack].Fmap = (Long *) SparseCore_malloc (n, sizeof (Long), cc) ;
            Work [stack].Cmap = (Long *) SparseCore_malloc (maxfn, sizeof(Long), cc);
            Work [stack].WTwork =
            (double *) SparseCore_malloc (wtsize, sizeof (double), cc) ;
        #endif
        if (keepH)
        {
            // Staircase 是 H的 永久部分，走的这里
            Work [stack].Stair1 = NULL ;
        }
        else
        {
            // 给 H 开辟 Staircase 空间
            Work [stack].Stair1 =
                (Long *) SparseCore_malloc (maxfn, sizeof (Long), cc) ;
        }
        
        
        Work [stack].sumfrank = 0 ;
        Work [stack].maxfrank = 0 ;

        Work [stack].wscale = 0 ;
        Work [stack].wssq   = 0 ;
    }

    *p_wtsize = wtsize ;
    return (Work) ;
}


// =============================================================================
// === free_Work ===============================================================
// =============================================================================

// 释放work的内容，但不释放工作数组本身
void free_Work
(
    qr_work *Work,
    Long ns,            
    Long n,             
    Long maxfn,         
    Long wtsize,        
    sparse_common *cc
)
{
    if (Work != NULL)
    {
        for (Long stack = 0 ; stack < ns ; stack++)
        {
            #ifdef NUMA_ALLOC
            SparseCore_numa_free(n,      sizeof (Long),   Work [stack].Fmap,   cc) ;
            SparseCore_numa_free(maxfn,  sizeof (Long),   Work [stack].Cmap,   cc) ;
            SparseCore_numa_free(wtsize, sizeof (double), Work [stack].WTwork, cc) ;
            #else
            SparseCore_free (n,      sizeof (Long),   Work [stack].Fmap,   cc) ;
            SparseCore_free (maxfn,  sizeof (Long),   Work [stack].Cmap,   cc) ;
            SparseCore_free (wtsize, sizeof (double), Work [stack].WTwork, cc) ;
            #endif
            SparseCore_free (maxfn,  sizeof (Long),   Work [stack].Stair1, cc) ;
            
            Work [stack].Fmap = NULL ;
            Work [stack].Cmap = NULL ;
            Work [stack].Stair1 = NULL ;
            Work [stack].WTwork = NULL ;
        }
    }
}


// =============================================================================
// === qr_factorize ==========================================================
// =============================================================================
// 请注意，代码的第一部分分配工作空间和产生的QR因子。
// 然后它在这个空间内完成所有的工作，而不重新分配。
// 分配的内存块的大小是预先知道的，从符号分析，不依赖于特定的数值。完成后，它释放其工作区。
qr_numeric *qr_factorize
(
    sparse_csc **Ahandle, Long freeA, 
    double tol, Long ntol, qr_symbolic *QRsym, sparse_common *cc )
{
    Long *Wi, *Qfill, *PLinv, *Cm, *Sp, *Stack_size,
        *TaskFront, *TaskFrontp, *TaskStack, *Stack_maxstack ;
    double *Sx, **Rblock, **Cblock, **Stacks ;
    qr_numeric *QRnum ;
    Long nf, m, n, anz, fchunk, maxfn, rank, maxfrank, rjsize, rank1,
        maxstack,j, wtsize, stack, ns, ntasks, keepH, hisize ;
    char *Rdead ;
    sparse_csc *A ;
    qr_work *Work ;
    #ifdef PRINT_TIME
    // double timeStart, timeEnd;
    // struct timeval tv;
    #endif

    // -------------------------------------------------------------------------
    // 获取符号对象的输入和内容
    // -------------------------------------------------------------------------
    
    // chunk_getSettings(32, 5000, 4, 4);

    if (QRsym == NULL)
    {
        if (freeA)
        {
            SparseCore_free_sparse (Ahandle, cc) ;
        }
        return (NULL) ;
    }

    A = *Ahandle ;

    nf = QRsym->nf ;             
    m = QRsym->m ;              
    n = QRsym->n ;
    anz = QRsym->anz ;          

    keepH = QRsym->keepH ;

    rjsize = QRsym->rjsize ;

    Sp = QRsym->Sp ;               
    Qfill = QRsym->Qfill ;         
    PLinv = QRsym->PLinv ;          

    ns = QRsym->ns ;                
    ntasks = QRsym->ntasks ;      

    maxfn  = QRsym->maxfn ;        
    hisize = QRsym->hisize ;       

    // 并行分解中的 任务Task包含的波前阵信息与栈信息·
    TaskFrontp = QRsym->TaskFrontp ;
    TaskFront  = QRsym->TaskFront ;
    TaskStack  = QRsym->TaskStack ;

    maxstack = QRsym->maxstack ;
    Stack_maxstack = QRsym->Stack_maxstack ;

    if (!(QRsym->do_rank_detection))
    {
        // 如果分析中没有考虑到秩检测，则禁用秩检
        tol = -1 ;
    }

    // -------------------------------------------------------------------------
    // 分配空间
    // -------------------------------------------------------------------------

    SparseCore_allocate_work (0, MAX (m,nf), 0, cc) ;

    Wi = (Long *) cc->Iwork ;  
    Cm = Wi ;                  // 规模为 nf

    // Cblock 是所有线程共享的工作区， 应该是 nf 个波前阵的贡献块C
    Cblock = (double **) SparseCore_malloc (nf+1, sizeof (double *), cc) ;

    Work = NULL ;              
    fchunk = MIN (m, FCHUNK) ;
    wtsize = 0 ;

    // -------------------------------------------------------------------------
    // 创建S = AE
    // -------------------------------------------------------------------------
    // Sx 的空间在这里开辟， S = A(p,q) 
    Sx = (double *) SparseCore_malloc (anz, sizeof (double), cc) ;

    if (cc->status == SPARSE_OK)
    {
        // 这个函数使用了 Wi 作为工作空间 ( Iwork(0:m-1))
        qr_stranspose2 (A, Qfill, Sp, PLinv, Sx, Wi) ;
        // 之后 Wi 不再需要
    }

    // -------------------------------------------------------------------------
    // A 不再需要， 释放空间
    // -------------------------------------------------------------------------

    if (freeA)
    {
        SparseCore_free_sparse (Ahandle, cc) ;
    }

    if (cc->status < SPARSE_OK)
    {
        FREE_WORK ;
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // 分配数值对象, 给 QRnum 分配空间
    // -------------------------------------------------------------------------

    QRnum = (qr_numeric *)
        SparseCore_malloc (1, sizeof (qr_numeric), cc) ;

    if (cc->status < SPARSE_OK)
    {
        FREE_WORK ;
        return (NULL) ;
    }

    Rblock     = (double **) SparseCore_malloc (nf, sizeof (double *), cc) ;
    Rdead      = (char *)   SparseCore_calloc (n,  sizeof (char),    cc) ;

    Stacks     = (double **) SparseCore_calloc (ns, sizeof (double *), cc) ;
    Stack_size = (Long *)   SparseCore_calloc (ns, sizeof (Long),    cc) ;

    QRnum->Rblock     = Rblock ;     // 253 行
    QRnum->Rdead      = Rdead ;      // 254 行
    QRnum->Stacks     = Stacks ;     // 310 行  
    QRnum->Stack_size = Stack_size ;  // 308 行

    // 给每个波前阵的 Stair, Tau, Hii 分配永久空间， 就在这里分配的空间
    QRnum->HStair= (Long *)  SparseCore_malloc (rjsize, sizeof (Long),  cc) ;
    QRnum->HTau  = (double *) SparseCore_malloc (rjsize, sizeof (double), cc) ;
    QRnum->Hii   = (Long *)  SparseCore_malloc (hisize, sizeof (Long),  cc) ;
    QRnum->Hm    = (Long *)  SparseCore_malloc (nf,     sizeof (Long),  cc) ;
    QRnum->Hr    = (Long *)  SparseCore_malloc (nf,     sizeof (Long),  cc) ;
    QRnum->HPinv = (Long *)  SparseCore_malloc (m,      sizeof (Long),  cc) ;

    QRnum->n = n ;
    QRnum->m = m ;
    QRnum->nf = nf ;
    QRnum->rjsize = rjsize ;
    QRnum->hisize = hisize ;
    QRnum->keepH = keepH ;
    QRnum->maxstack = maxstack ;
    QRnum->ns = ns ;
    QRnum->ntasks = ntasks ;
    QRnum->maxfm = EMPTY ;     

    if (cc->status < SPARSE_OK)
    {
        qr_freenum (&QRnum, cc) ;
        FREE_WORK ;
        return (NULL) ;
    }


    // -------------------------------------------------------------------------
    // 分配工作空间
    // -------------------------------------------------------------------------

    Work = get_Work (ns, n, maxfn, keepH, fchunk, &wtsize, cc) ;

    // -------------------------------------------------------------------------
    // 分配初始化每个栈
    // -------------------------------------------------------------------------
#ifdef NUMA_ALLOC
    size_t numa_node_rank = 0;
#endif
    // double timeStart, timeEnd;
    // struct timeval tv;

    if (cc->status == SPARSE_OK)
    {
        // gettimeofday(&tv, NULL);
        // timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
        for (stack = 0 ; stack < ns ; stack++)
        {
            double *Stack ;
            size_t stacksize = (ntasks == 1) ?
                maxstack : Stack_maxstack [stack] ;
            Stack_size [stack] = stacksize ;
            // 给 ns 个栈开辟空间， 这里也可以分NUMA 开辟空间
            #ifdef NUMA_ALLOC
                numa_node_rank = stack % 4;
                Stack = (double *) SparseCore_malloc_onnode (stacksize, sizeof (double), numa_node_rank, cc) ;
                Work [stack].numa_node_rank = numa_node_rank;
            #else
                Stack = (double *) SparseCore_malloc (stacksize, sizeof (double), cc) ;
            #endif
            Stacks [stack] = Stack ;
            Work [stack].Stack_head = Stack ;
            Work [stack].Stack_top  = Stack + stacksize ;
        }
        // timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        // printf("Stack allocate time: %f\n\n", (timeEnd-timeStart));
    }
    // printf("ns = %d, ntasks = %d \n", ns, ntasks);

    // -------------------------------------------------------------------------
    // 如果内存不足，则使用顺序 case和fchunk = 1
    // -------------------------------------------------------------------------

    if (cc->status < SPARSE_OK)
    {
        printf(" 内存不足，串行使用！\n");
        if (Stacks != NULL)
        {
            for (stack = 0 ; stack < ns ; stack++)
            {
                size_t stacksize = (ntasks == 1) ?
                    maxstack : Stack_maxstack [stack] ;
                SparseCore_free (stacksize, sizeof (double), Stacks [stack], cc) ;
            }
        }
        SparseCore_free (ns, sizeof (double *), Stacks,     cc) ;
        SparseCore_free (ns, sizeof (Long),    Stack_size, cc) ;

        free_Work (Work, ns, n, maxfn, wtsize, cc) ;
        SparseCore_free (ns, sizeof (qr_work), Work, cc) ;

        ns = 1 ;
        ntasks = 1 ;
        fchunk = 1 ;
        cc->status = SPARSE_OK ;
        Work = get_Work (ns, n, maxfn, keepH, fchunk, &wtsize, cc) ;
        Stacks     = (double **) SparseCore_calloc (ns, sizeof (double *), cc) ;
        Stack_size = (Long *)   SparseCore_calloc (ns, sizeof (Long),    cc) ;
        QRnum->Stacks     = Stacks ;
        QRnum->Stack_size = Stack_size ;
        if (cc->status == SPARSE_OK)
        {
            double *Stack ;
            Stack_size [0] = maxstack ;
            Stack = (double *) SparseCore_malloc (maxstack, sizeof (double), cc) ;
            Stacks [0] = Stack ;
            Work [0].Stack_head = Stack ;
            Work [0].Stack_top  = Stack + maxstack ;
        }
    }

    QRnum->ns = ns ;
    QRnum->ntasks = ntasks ;

    // -------------------------------------------------------------------------
    // 检查
    // -------------------------------------------------------------------------

    if (cc->status < SPARSE_OK)
    {
        qr_freenum (&QRnum, cc) ;
        FREE_WORK ;
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // 创建 Blob , kernel 分解需要
    // -------------------------------------------------------------------------

    qr_blob Blob ;
    Blob.QRsym = QRsym ;
    Blob.QRnum = QRnum ;
    Blob.tol = tol ;
    Blob.Work = Work ;  //栈在这里
    Blob.Cm = Cm ;      // cc 里的空间，已经开完了
    Blob.Cblock = Cblock ;   // 206 行
    Blob.Sx = Sx ;
    Blob.ntol = ntol ;
    Blob.fchunk = fchunk ;
    Blob.cc = cc ;

    // -------------------------------------------------------------------------
    // 初始化 "pure" 浮点计数
    // -------------------------------------------------------------------------

    cc->SPQR_flopcount = 0 ;

    // -------------------------------------------------------------------------
    // 数值QR分解
    // -------------------------------------------------------------------------

    #ifdef PRINT_TIME
    // gettimeofday(&tv, NULL);
    // timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
    #endif
    if (ntasks == 1)
    {
        
        qr_kernel (0, &Blob) ;     
    }
    else
    {
// #ifdef MULTI
        qr_multithreads(ntasks, &Blob);
// #else

//     for (Long id = 0 ; id < ntasks-1 ; id++)
//     {
//         qr_kernel (id, &Blob) ;
//     }
// #endif
    }
    #ifdef PRINT_TIME
    // gettimeofday(&tv, NULL);
    // timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
    // printf ("REAL kernel_sum time: %lf\n", timeEnd - timeStart);
    #endif
    // -------------------------------------------------------------------------
    // 检查是否越界
    // -------------------------------------------------------------------------

    if (cc->status < SPARSE_OK)
    {
        qr_freenum (&QRnum, cc) ;
        FREE_WORK ;
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // 最终的秩
    // -------------------------------------------------------------------------
    #ifdef PRINT_TIME
    // gettimeofday(&tv, NULL);
    // timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
    #endif
    rank = 0 ;
    maxfrank = 1 ;
    for (stack = 0 ; stack < ns ; stack++)
    {
        rank += Work [stack].sumfrank ;
        maxfrank = MAX (maxfrank, Work [stack].maxfrank) ;
    }
    QRnum->rank = rank ;         // 稀疏矩阵A的 rank         
    QRnum->maxfrank = maxfrank ; // 所有波前阵的最大rank

    // -------------------------------------------------------------------------
    // 最终死列的二范数
    // -------------------------------------------------------------------------

    // double wscale = 0 ;
    // double wssq = 1 ;
    // for (stack = 0 ; stack < ns ; stack++)
    // {
    //     double ws = Work [stack].wscale ;
    //     double wq = Work [stack].wssq ;
    //     if (wq != 0)
    //     {
    //         double wk = ws * sqrt (wq) ;
    //         if (wscale < wk)
    //         {
    //             double rr = wscale / wk ;
    //             wssq = 1 + wssq * rr * rr ;
    //             wscale = wk ;
    //         }
    //         else
    //         {
    //             double rr = wk / wscale ;
    //             wssq += rr * rr ;
    //         }
    //     }
    // }
    // QRnum->norm_E_fro = wscale * sqrt (wssq) ;
    // cc->SPQR_norm_E_fro = QRnum->norm_E_fro ;

    // -------------------------------------------------------------------------
    // 释放空间
    // -------------------------------------------------------------------------

    free_Work (Work, ns, n, maxfn, wtsize, cc) ; 
    if (freeA) SparseCore_free_sparse (Ahandle, cc) ; 
    SparseCore_free (anz, sizeof (double), Sx, cc) ; 
    Sx = NULL ;

    // -------------------------------------------------------------------------
    // 压缩栈
    // -------------------------------------------------------------------------
    Long any_moved = FALSE ;

    int shrink = cc->SPQR_shrink ; // 默认为 1 
    //int shrink = 0;
    if (shrink > 0)
    {
        for (stack = 0 ; stack < ns ; stack++)
        {
            size_t stacksize = Stack_size [stack] ;
            double *Stack = Stacks [stack] ;
            
            size_t newstacksize = Work [stack].Stack_head - Stack ;
            
            if (shrink > 1)
            {
                Cblock [stack] = (double *) SparseCore_malloc (newstacksize,
                    sizeof (double), cc) ;
                if (Cblock [stack] == NULL)
                {
                    cc->status = SPARSE_OK ;
                    Cblock [stack] = Stack ;
		    cc->memory_inuse +=
                        ((newstacksize-stacksize) * sizeof (double)) ;
                }
                else
                {
                    memcpy (Cblock [stack], Stack, newstacksize*sizeof(double)) ;
                    #ifdef NUMA_ALLOC
                    SparseCore_numa_free(stacksize, sizeof (double), Stack, cc) ;
                    #else
                    SparseCore_free (stacksize, sizeof (double), Stack, cc) ;
                    #endif
                }
                stacksize = newstacksize ;
            }
            else
            {
                #ifdef NUMA_ALLOC
                Cblock [stack] =  
                    (double *) SparseCore_numa_realloc (
                    newstacksize,   
                    sizeof (double), 
                    Stack,          
                    &stacksize,     
                    cc) ;
                #else
                Cblock [stack] =  
                    (double *) SparseCore_realloc (
                    newstacksize,   
                    sizeof (double), 
                    Stack,          
                    &stacksize,     
                    cc) ;
                #endif
            }
            Stack_size [stack] = stacksize ;
            any_moved = any_moved || (Cblock [stack] != Stack) ;
        }
    }

    // -------------------------------------------------------------------------
    // 调整R 块的指针
    // -------------------------------------------------------------------------

    if (any_moved)
    {
        for (Long task = 0 ; task < ntasks ; task++)
        {
            Long kfirst, klast ;
            if (ntasks == 1)
            {
                kfirst = 0 ;
                klast = nf ;
                stack = 0 ;
            }
            else
            {
                kfirst = TaskFrontp [task] ;
                klast  = TaskFrontp [task+1] ;
                stack  = TaskStack [task] ;
            }
            double *Old_Stack = Stacks [stack] ;
            double *New_Stack = Cblock [stack] ;
            if (New_Stack != Old_Stack)
            {
                for (Long kf = kfirst ; kf < klast ; kf++)
                {
                    Long f = (ntasks == 1) ? kf : TaskFront [kf] ;
                    Rblock [f] = New_Stack + (Rblock [f] - Old_Stack) ;
                }
            }
        }
        for (stack = 0 ; stack < ns ; stack++)
        {
            Stacks [stack] = Cblock [stack] ;
        }
    }
    
    // -------------------------------------------------------------------------
    // 释放剩余空间
    // -------------------------------------------------------------------------
    

    SparseCore_free (ns, sizeof (qr_work), Work, cc) ; 
    Work = NULL ; 
    SparseCore_free (nf+1, sizeof (double *), Cblock, cc) ; 
    Cblock = NULL ;
    // -------------------------------------------------------------------------
    // 提取H的隐式行置换
    // -------------------------------------------------------------------------

    // 这一步必须串行 
    if (keepH)
    {
        qr_hpinv (QRsym, QRnum, Wi) ;
        // 完成 H 分布模式 中隐含的行置换
    }

    // -------------------------------------------------------------------------
    // 找到秩，返回结果
    // -------------------------------------------------------------------------

    if (ntol >= n)
    {
        rank1 = rank ;
    }
    else
    {
        rank1 = 0 ;
        for (j = 0 ; j < ntol ; j++)
        {
            if (!Rdead [j])
            {
                rank1++ ;
            }
        }
    }
    QRnum->rank1 = rank1 ;
    #ifdef PRINT_TIME
    // gettimeofday(&tv, NULL);
    // timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
    // printf ("final operation time: %lf\n", timeEnd - timeStart);
    #endif
    return (QRnum) ;
} // qr_factorize

// =============================================================================
// === qr_stranspose2 ========================================================
// =============================================================================
// 用行压缩形式构造S = A (p,q)的数值
void qr_stranspose2(
    sparse_csc *A, Long *Qfill, Long *Sp, Long *PLinv, double *Sx, Long *W 
)
{
    Long i, j, p, pend, row, col, s, m, n, *Ap, *Ai ;
    double *Ax ;

    m = A->nrow ;
    n = A->ncol ;
    Ap = (Long *) A->p ;
    Ai = (Long *) A->i ;
    Ax = (double *) A->x ;

    for (row = 0 ; row < m ; row++)
    {
        W [row] = Sp [row] ;
    }

    for (col = 0 ; col < n ; col++)     
    {
        j = Qfill ? Qfill [col] : col ; 
        pend = Ap [j+1] ;
        for (p = Ap [j] ; p < pend ; p++)
        {
            i = Ai [p] ;             // 用列索引找到行号  i
            row = PLinv [i] ;        // 找到置换后的实际行 row
            s = W [row]++ ;          // 利用符号分解的Sp,找到S 的行指针s 
            Sx [s] = Ax [p] ;        // 给 Sx 赋值
        }
    }
} 

// =============================================================================
// === qr_kernel =============================================================
// =============================================================================
// 在单个任务中分解所有波前阵
void qr_kernel( 
    Long task, qr_blob *Blob
)
{
    #ifdef PRINT_PACK_TIME
    double timeStart, timeEnd;
    struct timeval tv;
    #endif
    // ---------------------
    // 获取 Blob
    // ---------------------
    qr_symbolic *          QRsym = Blob->QRsym ;
    qr_numeric *   QRnum = Blob->QRnum ;
    double                   tol = Blob->tol ;
    Long                     ntol = Blob->ntol ;
    Long                     fchunk = Blob->fchunk ;
    qr_work *      Work = Blob->Work ;
    Long *                   Cm = Blob->Cm ;
    double **                 Cblock = Blob->Cblock ;
    double *                  Sx = Blob->Sx ;
    sparse_common *         cc = Blob->cc ;

    // ---------------------
    // 获取 QRsym 中的对象
    // ---------------------
    Long *  Super = QRsym->Super ;  // 波前阵的 主列索引 Super[f+1] - Super[f] 为f的主列数目
    Long *  Rp = QRsym->Rp ;        // 波前阵 的列索引 Rp [f+1] - Rp [f] 表示波前阵 f 的列数
    Long *  Rj = QRsym->Rj ;        // 波前阵 的非主列 索引
    Long *  Sleft = QRsym->Sleft ;  // S的最左列是第j列的行索引  Sleft[j] ... Sleft[j+1]-1
    Long *  Sp = QRsym->Sp ;        // S的行指针    
    Long *  Sj = QRsym->Sj ;        // S的列索引   
    Long *  Child = QRsym->Child ;  // 表示列消去树的 孩子关系 
    Long *  Childp = QRsym->Childp ;   // child的索引 
    Long    maxfn  = QRsym->maxfn ;    // 波前矩阵中 列数目的最大值 
    Long    nf = QRsym->nf ;           // 波前矩阵的数目

    Long *  Hip = QRsym->Hip ;       //H保存时，波前阵f的行 数据如下保存：Hii [Hip [f] ... Hip [f] + Hm [f]]

    Long *  TaskFront = QRsym->TaskFront ;    // 保存每个任务中的波前阵 数组，size:nf+1
    Long *  TaskFrontp = QRsym->TaskFrontp ;   // TaskFront的索引 
    Long *  TaskStack = QRsym->TaskStack ;     // 任务栈 TaskStack
    Long *  On_stack = QRsym->On_stack ;       // On_stack[f] 给波前阵f 记录栈号  

    Long *  Post = QRsym->Post ;            // 给串行使用的，大小为nf    

    // ---------------------
    // 获取 QRnum 中的对象
    // ---------------------
    double ** Rblock = QRnum->Rblock ; // 大小为nf。R[f] 是一个(double *)指针指向波前阵 f 中的R块
    char *   Rdead = QRnum->Rdead ;    // 如果k是一个死主列，则Rdead[k] = 1，否则Rdead[k] = 0.
                                       // 如果没有死列，这个指针为NULL。 如果 m < n. 则至少有n-m个死列
    Long *   HStair = QRnum->HStair ;  // 数组 Hstair[Rp [f] ...Rp [f+1]-1 ] 给出了波前阵 f 每一列的 staircase
    double *  HTau = QRnum->HTau ;   //  数组HTau [Rp [f] ... Rp [f+1]-1] 给出 波前矩阵 f 每一列的 Householder系数
    Long *   Hii = QRnum->Hii ;      // 大小为Hisize， 保存波前阵 f 包含的 行号
    Long *   Hm = QRnum->Hm ;        // 大小为nf， Hm[f] 表示波前阵 f 的行数
    Long *   Hr = QRnum->Hr ;        // 大小为nf， Hr[f] 表示波前阵 f 中 R 块的行数
    Long     keepH = QRnum->keepH ;  // 是否保留 H，这里保留
    Long     ntasks = QRnum->ntasks ;// 符号分析划分的 tasks 数目    

    Long stack, kfirst, klast ;

    if (ntasks == 1)
    {
        kfirst = 0 ;
        klast  = nf ;
        stack  = 0 ;
    }
    else
    {
        kfirst = TaskFrontp [task] ;   // 取task的第一个波前阵
        klast  = TaskFrontp [task+1] ; // 取task的最后一个波前阵
        stack  = TaskStack [task] ;    // 取task的 栈号
    }

    double * Stack_head = Work [stack].Stack_head ; // 栈底
    double * Stack_top = Work [stack].Stack_top ;   // 栈顶

    //如果H保存，则Tau和Stair将指向QRnum中的永久空间
    double * Tau = keepH ? NULL : Work [stack].WTwork ;
    Long *  Stair = keepH ? NULL : Work [stack].Stair1 ;
    double * W = Work [stack].WTwork + (keepH ? 0 : maxfn) ;

    Long *  Fmap = Work [stack].Fmap ;
    Long *  Cmap = Work [stack].Cmap ;

    Long    sumfrank = Work [stack].sumfrank ; //用于保存这个task里波前阵的秩和
    Long    maxfrank = Work [stack].maxfrank ; //用于保存这个task里的最大rank

    double wscale = Work [stack].wscale ;
    double wssq   = Work [stack].wssq   ;

    #ifdef PRINT_PACK_TIME
    
    double assemble_time = 0, pack_time = 0;
    #endif
    for (Long kf = kfirst ; kf < klast ; kf++) // 遍历这个task中的波前阵
    {
        Long f = (ntasks == 1) ? Post [kf] : TaskFront [kf] ; // 取波前阵号 f

        if (keepH)
        {
            Stair = HStair + Rp [f] ; // 靠索引Rp[f]找到列号，然后找到Staircase起点
            Tau = HTau + Rp [f] ;     // 同上靠索引Rp[f]， 找到对应的 HTau 起点
        }
        // 计算波前阵F 的行数，初始化staircase 和 Fmap
        Long fm = qr_fsize (f, Super, Rp, Rj, Sleft, Child, Childp, Cm,
            Fmap, Stair) ;
        Long fn = Rp [f+1] - Rp [f] ;     // 波前阵 f 的列数
        Long col1 = Super [f] ;           // 主列起点   
        Long fp = Super [f+1] - col1 ;    // 波前阵 f 的主列数   
        Long fsize = fm * fn ;
        if (keepH)
        {
            Hm [f] = fm ;
        }

        double *F = Stack_head ;
        Rblock [f] = F ;                     
        Stack_head += fsize ;  // 给Rblock 留空间
        #ifdef PRINT_PACK_TIME
        gettimeofday(&tv, NULL);
        timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
        #endif
        // 把子节点组装进来
        qr_assemble (f, fm, keepH,
            Super, Rp, Rj, Sp, Sj, Sleft, Child, Childp,
            Sx, Fmap, Cm, Cblock,
            Hr, Stair, Hii, Hip, F, Cmap) ;  
        #ifdef PRINT_PACK_TIME
        gettimeofday(&tv, NULL);
        timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        assemble_time += (timeEnd - timeStart);
        #endif

        for (Long p = Childp [f] ; p < Childp [f+1] ; p++)
        {
            Long c = Child [p] ;
            if (ntasks == 1 || On_stack [c] == stack)
            {
                Long ccsize = qr_csize (c, Rp, Cm, Super) ;
                Stack_top = MAX (Stack_top, Cblock [c] + ccsize) ;
            }
        }
        #ifdef PRINT_TIME
        // gettimeofday(&tv, NULL);
        // timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
        #endif
        Long frank = qr_front (fm, fn, fp, tol, ntol - col1,
            fchunk, F, Stair, Rdead + col1, Tau, W,
            &wscale, &wssq, cc) ;
        #ifdef PRINT_TIME
        // gettimeofday(&tv, NULL);
        // timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        // front_time += (timeEnd - timeStart);
        #endif
        sumfrank += frank ;
        maxfrank = MAX (maxfrank, frank) ;

        Long csize = qr_fcsize (fm, fn, fp, frank) ;
        Stack_top -= csize ;

        #ifdef PRINT_PACK_TIME
        gettimeofday(&tv, NULL);
        timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
        #endif
        Cblock [f] = Stack_top ;
        Cm [f] = qr_cpack (fm, fn, fp, frank, F, Cblock [f]) ;

        Long rm ;
        Long rsize = qr_rhpack (keepH, fm, fn, fp, Stair, F, F, &rm) ;
        #ifdef PRINT_PACK_TIME
        gettimeofday(&tv, NULL);
        timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        pack_time += (timeEnd - timeStart);
        #endif
        if (keepH)
        {
            Hr [f] = rm ;
        }

        Stack_head -= fsize ;             
        Stack_head += rsize ;             

    }
    #ifdef PRINT_PACK_TIME
    printf("stack = %d, assemble time = %lf, pack time = %lf\n", stack, assemble_time, pack_time);
    #endif
    Work [stack].Stack_head = Stack_head  ;
    Work [stack].Stack_top = Stack_top ;
    Work [stack].sumfrank = sumfrank ;  
    Work [stack].maxfrank = maxfrank ;  

    Work [stack].wscale = wscale ;
    Work [stack].wssq   = wssq   ;
} 

// =============================================================================
// === qr_hpinv ==============================================================
// =============================================================================
// 完成h模式中隐含的行置换。这必须在所有线程完成分解矩阵并找到所有死列之后按顺序完成。也决定了QRnum→maxfm。
void qr_hpinv
(
    qr_symbolic *QRsym,
    qr_numeric *QRnum,
    Long *W           
)
{
    Long *Hi, *Hii, *Hip, *HPinv, *Hr, *Super, *Rp, *Hm, *Sleft, *PLinv ;
    Long nf, m, n, f, rm, i, row1, row2, fm, fn, fp, cm, cn, maxfm ;

    nf = QRsym->nf ;
    m = QRsym->m ;
    n = QRsym->n ;
    Hr = QRnum->Hr ;
    Hm = QRnum->Hm ;
    Hii = QRnum->Hii ;
    Hip = QRsym->Hip ;
    HPinv = QRnum->HPinv ;
    Super = QRsym->Super ;
    Rp = QRsym->Rp ;
    Sleft = QRsym->Sleft ;
    PLinv = QRsym->PLinv ;
    maxfm = 0 ;

    row1 = 0 ;                              
    row2 = m ;                              

    for (i = Sleft [n] ; i < m ; i++)
    {
        W [i] = (--row2) ;
    }

    for (f = 0 ; f < nf ; f++)              
    {
        Hi = &Hii [Hip [f]] ;              
        rm = Hr [f] ;                      
        for (i = 0 ; i < rm ; i++)
        {
            W [Hi [i]] = row1++ ;
        }

        fp = Super [f+1] - Super [f] ;
        fn = Rp [f+1] - Rp [f] ;
        fm = Hm [f] ;
        maxfm = MAX (maxfm, fm) ;
        cn = fn - fp ;
        cm = MIN (fm - rm, cn) ;
        for (i = fm-1 ; i >= rm + cm ; i--)
        {
            W [Hi [i]] = (--row2) ;
        }
    }
    QRnum->maxfm = maxfm ;

    for (i = 0 ; i < m ; i++)
    {
        HPinv [i] = W [PLinv [i]] ;
    }

    for (f = 0 ; f < nf ; f++)
    {
        Hi = &Hii [Hip [f]] ;                   
        fm = Hm [f] ;
        for (i = 0 ; i < fm ; i++)
        {
            Hi [i] = W [Hi [i]] ;
        }

    }
} 

// =============================================================================
// === qr_fsize ==============================================================
// =============================================================================
// 计算波前阵F的行数, 初始化他的阶梯staircase，创建Fmap
Long qr_fsize   (  
    Long f, Long *Super, Long *Rp, Long *Rj, 
    Long *Sleft, Long *Child,  Long *Childp, 
    Long *Cm, Long *Fmap, Long *Stair 
)
{
    Long col1, col2, p1, p2, fp, fn, fm, col, p, j, c, pc, cm, ci, t, fpc ;

    // -------------------------------------------------------------------------
    // 获取矩阵F
    // -------------------------------------------------------------------------

    col1 = Super [f] ;      
    col2 = Super [f+1] ;
    p1 = Rp [f] ;          
    p2 = Rp [f+1] ;
    fp = col2 - col1 ;      
    fn = p2 - p1 ;          

    // -------------------------------------------------------------------------
    // 创建波前阵 F 的 Fmap
    // -------------------------------------------------------------------------

    for (p = p1, j = 0 ; p < p2 ; p++, j++)
    {
        col = Rj [p] ;              
        Fmap [col] = j ;
    }

    // -------------------------------------------------------------------------
    // 初始化波前阵F的阶梯staircase
    // -------------------------------------------------------------------------

    for (j = 0 ; j < fp ; j++)
    {
        col = j + col1 ;
        Stair [j] = Sleft [col+1] - Sleft [col] ;
    }

    for ( ; j < fn ; j++)
    {
        Stair [j] = 0 ;
    }

    // -------------------------------------------------------------------------
    // 为每个孩子建造阶梯
    // -------------------------------------------------------------------------

    for (p = Childp [f] ; p < Childp [f+1] ; p++)
    {
        c = Child [p] ;                 
        pc = Rp [c] ;                   
        cm = Cm [c] ;                  

        fpc = Super [c+1] - Super [c] ; 
    
        pc += fpc ;                     
        for (ci = 0 ; ci < cm ; ci++)
        {
            col = Rj [pc + ci] ;        
            j = Fmap [col] ;            
            //Stair [j]++ ; 
            __sync_fetch_and_add(&Stair[j],1);              
        }
    }

    // -------------------------------------------------------------------------
    // 将Stair替换为cumsum ([0 Stair])，并查找波前阵F的对应行
    // -------------------------------------------------------------------------

    fm = 0 ;
    for (j = 0 ; j < fn ; j++)
    {
        t = fm ;
        fm += Stair [j] ;
        Stair [j] = t ;
    }

    return (fm) ;
} //end of qr_fsize

// ==========================================
// === qr_assemble ==========================
// ==========================================
// 组装函数，把 C 组装到 波前阵F
void qr_assemble(
    /* inputs, not modified */
    Long f, Long fm, int keepH, Long *Super,
    Long *Rp, Long *Rj, Long *Sp, Long *Sj, Long *Sleft,
    Long *Child, Long *Childp, double *Sx, Long *Fmap,
    Long *Cm, double **Cblock,
    Long *Hr, Long *Stair, Long *Hii,  Long *Hip, double *F, Long *Cmap
)
{
    double *Fi, *Fj, *C ;
    Long k, fsize, fn, col1, col2, p, p1, p2, fp, j, leftcol, row, col, i,
        cm, cn, ci, cj, c, pc, fnc, fpc, rmc ;
    Long *Hi = NULL, *Hichild = NULL ;

    /* ---------------------------------------------------------------------- */
    /* 获取波前阵F */
    /* ---------------------------------------------------------------------- */

    col1 = Super [f] ;      
    col2 = Super [f+1] ;
    p1 = Rp [f] ;           
    p2 = Rp [f+1] ;
    fp = col2 - col1 ;      
    fn = p2 - p1 ;          

    fsize = fm * fn ;
    for (k = 0 ; k < fsize ; k++)
    {
        F [k] = 0 ;
    }

    Hi = &Hii [Hip [f]] ;   

    /* ---------------------------------------------------------------------- */
    /* 组装 S 中最左列是主列的行到 F 中 */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < fp ; k++)
    {
        leftcol = k + col1 ;

        for (row = Sleft [leftcol] ; row < Sleft [leftcol+1] ; row++)
        {
            i = Stair [k]++ ;
            Fi = F + i ;       

            for (p = Sp [row] ; p < Sp [row+1] ; p++)
            {
                col = Sj [p] ;          
                j = Fmap [col] ;       
                Fi [j*fm] = Sx [p] ;   
            }

            
            Hi [i] = row ;        
            
        }
    }

    /* ---------------------------------------------------------------------- */
    /* 组装每一个孩子 */
    /* ---------------------------------------------------------------------- */

    for (p = Childp [f] ; p < Childp [f+1] ; p++)
    {

        /* ------------------------------------------------------------------ */
        /* 获取孩子节点索引*/
        /* ------------------------------------------------------------------ */

        c = Child [p] ;                
        pc = Rp [c] ;                   
        cm = Cm [c] ;                   
        fnc = Rp [c+1] - pc ;          
        fpc = Super [c+1] - Super [c] ;
        cn = fnc - fpc ;               
        pc += fpc ;                   
        C = Cblock [c] ;               

        rmc = Hr [c] ;

        Hichild = &Hii [Hip [c] + rmc] ;
        

        /* ------------------------------------------------------------------ */
        /* 构造 Cmap */
        /* ------------------------------------------------------------------ */

        for (ci = 0 ; ci < cm ; ci++)
        {
            col = Rj [pc + ci] ;        
            j = Fmap [col] ;            
            i = Stair [j]++ ;           
            Cmap [ci] = i ;            
            
            Hi [i] = Hichild [ci] ;
            
        }

        /* ------------------------------------------------------------------ */
        /* 将C的三角形部分复制到F (cm乘cm，上三角形) */
        /* ------------------------------------------------------------------ */

        for (cj = 0 ; cj < cm ; cj++)
        {
            
            col = Rj [pc + cj] ;            
            j = Fmap [col] ;                
            Fj = F + fm * j ;               
            for (ci = 0 ; ci <= cj ; ci++)
            {
                i = Cmap [ci] ;             
                Fj [i] = *(C++) ;           
            }
        }

        /* ------------------------------------------------------------------ */
        /* 将C的矩形部分复制到F (cm-by-(cn-cm)的大小) */
        /* ------------------------------------------------------------------ */

        for ( ; cj < cn ; cj++)
        {
            col = Rj [pc + cj] ;           
            j = Fmap [col] ;                
            Fj = F + fm * j ;              
            for (ci = 0 ; ci < cm ; ci++)
            {
                i = Cmap [ci] ;             
                Fj [i] = *(C++) ;           
            }
        }

    }

} //End of qr_assemble

// =============================================================================
// === qr_csize ==============================================================
// =============================================================================
// 返回C中子元素的项#
Long qr_csize (    
    Long c, Long *Rp, 
    Long *Cm, Long *Super )
{
    Long pc, cm, fnc, fpc, cn, csize ;

    pc = Rp [c] ;                   
    cm = Cm [c] ;                   
    fnc = Rp [c+1] - pc ;           
    fpc = Super [c+1] - Super [c] ; 
    cn = fnc - fpc ;                
    csize = (cm * (cm+1)) / 2 + cm * (cn - cm) ;
    return (csize) ;
} // End of qr_csize

// =============================================================================
// === qr_front ==============================================================
// =============================================================================
inline double qr_private_larfg (Long n, double *X, sparse_common *cc)
{
    double tau = 0 ;
    BLAS_INT N = n, one = 1 ;
    if (CHECK_BLAS_INT && !EQ (N,n))
    {
        cc->blas_ok = FALSE ;
    }
    if (!CHECK_BLAS_INT || cc->blas_ok)
    {
        #ifndef USEATL
            LAPACK_DLARFG (&N, X, X + 1, &one, &tau) ;
        #else
            HNU_larfg(N, X, X + 1, one, &tau) ;
        #endif
    }
    return (tau) ;
}

// 返回 tau
double qr_private_house  
(
    Long n,
    double *X,       
    sparse_common *cc
)
{
    return (qr_private_larfg (n, X, cc)) ;
}

inline void qr_private_larf (Long m, Long n, double *V, double tau,
    double *C, Long ldc, double *W, sparse_common *cc)
{
    BLAS_INT M = m, N = n, LDC = ldc, one = 1 ;
    char left = 'L' ;
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDC,ldc)))
    {
        cc->blas_ok = FALSE ;
        
    }
    if (!CHECK_BLAS_INT || cc->blas_ok)
    {
        #ifndef USEATL
            LAPACK_DLARF (&left, &M, &N, V, &one, &tau, C, &LDC, W) ;
        #else
            HNU_larf(CblasLeft, M, N, V, one, tau, C, LDC, W);    
        #endif
    }
}

void qr_private_apply1
(
    Long m,             
    Long n,
    Long ldc,           
    double *V,           
    double tau,          

    double *C,           
    double *W,           
    sparse_common *cc
)
{
    double vsave ;
    if (m <= 0 || n <= 0)
    {
        return ;        
    }
    vsave = V [0] ;     
    V [0] = 1 ;
    qr_private_larf (m, n, V, tau, C, ldc, W, cc) ;
    V [0] = vsave ;     
}

Long qr_front
(
    Long m,             
    Long n,
    Long npiv,          
    double tol,         
    Long ntol,          
    Long fchunk,        

    // input/output
    double *F,           
    Long *Stair,        
    char *Rdead,       

    // output, not defined on input
    double *Tau,         

    // workspace, undefined on input and output
    double *W,           

    // input/output
    double *wscale,
    double *wssq,

    sparse_common *cc
)
{
    double tau ;
    double wk ;
    double *V ;
    Long k, t, g, g1, nv, k1, k2, i, t0, vzeros, mleft, nleft, vsize, minchunk,
        rank ;

    npiv = MAX (0, npiv) ;  
    npiv = MIN (n, npiv) ;

    g1 = 0 ;                
    k1 = 0 ;                
    k2 = 0 ;
    V = F ;                 
    g = 0 ;                
    nv = 0 ;                
    vzeros = 0 ;            
    t = 0 ;                
    fchunk = MAX (fchunk, 1) ;
    minchunk = MAX (MINCHUNK, fchunk/MINCHUNK_RATIO) ; // max(4,32/4)
    rank = MIN (m,npiv) ;  


    ntol = MIN (ntol, npiv) ;   

    for (k = 0 ; k < n ; k++)
    {

        // ---------------------------------------------------------------------
        // 减少F的第k列，除去除“对角”F (g,k)之外的所有F (g,k)
        // ---------------------------------------------------------------------

        t0 = t ;           
        t = Stair [k] ;    

        if (g >= m)
        {
            for ( ; k < npiv ; k++)
            {
                Rdead [k] = 1 ;
                Stair [k] = 0 ;         
                Tau [k] = 0 ;
            }
            for ( ; k < n ; k++)
            {
                Stair [k] = m ;        
                Tau [k] = 0 ;
            }
            return (rank) ;
        }

        t = MAX (g+1,t) ;
        Stair [k] = t ;

        // ---------------------------------------------------------------------
        // 如果t比上一个t增长很多，现在对所有F应用H
        // ---------------------------------------------------------------------

        vzeros += nv * (t - t0) ;
        if (nv >= minchunk)
        {
            vsize = (nv*(nv+1))/2 + nv*(t-g1-nv) ;
            if (vzeros > MAX (16, vsize/2))
            {
                qr_larftb (
                    0,                          
                    t0-g1, n-k2, nv, m, m,
                    V,                          
                    &Tau [k1],                  
                    &F [INDEX (g1,k2,m)],      
                    W, cc) ;                   
                nv = 0 ;        
                vzeros = 0 ;
            }
        }

        // ---------------------------------------------------------------------
        // 找到一个减少k列的household反射
        // ---------------------------------------------------------------------

        tau = qr_private_house (t-g, &F [INDEX (g,k,m)], cc) ;

        // ---------------------------------------------------------------------
        // 检查第k列是否正确
        // ---------------------------------------------------------------------

        if (k < ntol && (wk = qr_abs (F [INDEX (g,k,m)])) <= tol)
        {

            // -----------------------------------------------------------------
            // 范数(F (g:t-1, k))太小;则第k个主列标记为死列
            // -----------------------------------------------------------------

            if (wk != 0)
            {
                if ((*wscale) == 0)
                {
                    (*wssq) = 1 ;
                }
                if ((*wscale) < wk)
                {
                    double rr = (*wscale) / wk ;
                    (*wssq) = 1 + (*wssq) * rr * rr ;
                    (*wscale) = wk ;
                }
                else
                {
                    double rr = wk / (*wscale) ;
                    (*wssq) += rr * rr ;
                }
            }

            for (i = g ; i < m ; i++)
            {

                F [INDEX (i,k,m)] = 0 ;
            }
            Stair [k] = 0 ;
            Tau [k] = 0 ;
            Rdead [k] = 1 ;

            if (nv > 0)
            {
                // 应用待定块的H 变换
                qr_larftb (
                    0,                          
                    t0-g1, n-k2, nv, m, m,
                    V,                          
                    &Tau [k1],                  
                    &F [INDEX (g1,k2,m)],       
                    W, cc) ;                   
                nv = 0 ;        
                vzeros = 0 ;
            }

        }
        else
        {

            // -----------------------------------------------------------------
            // 又发现了一个好的主列
            // -----------------------------------------------------------------

            Tau [k] = tau ;             
            if (nv == 0)
            {

                g1 = g ;                        
                k1 = k ;                        
                k2 = MIN (n, k+fchunk) ;      // FCHUNK 主要影响在这里
                V = &F [INDEX (g1,k1,m)] ;     

                mleft = m-g1 ;                  
                nleft = n-k1 ;                  
                if (mleft * (nleft-(fchunk+4)) < SMALL || mleft <= fchunk/2
                    || fchunk <= 1)
                {
                    k2 = n ;
                }
            }
            nv++ ;  

            FLOP_COUNT ((t-g) * (3 + 4 * (n-k-1))) ;

            // -----------------------------------------------------------------
            // 对当前面板应用第k个 Household 变换
            // -----------------------------------------------------------------

            qr_private_apply1 (t-g, k2-k-1, m, &F [INDEX (g,k,m)], tau,
                &F [INDEX (g,k+1,m)], W, cc) ;

            g++ ;   

            // -----------------------------------------------------------------
            // 如果达到panel最后，应用Household反射
            // -----------------------------------------------------------------

            if (k == k2-1 || g == m)            
            {
                qr_larftb (
                    0,                          
                    t-g1, n-k2, nv, m, m,
                    V,                         
                    &Tau [k1],                  
                    &F [INDEX (g1,k2,m)],       
                    W, cc) ;                    
                nv = 0 ;        
                vzeros = 0 ;
            }
        }

        // ---------------------------------------------------------------------
        // 确定主列的秩
        // ---------------------------------------------------------------------

        if (k == npiv-1)
        {
            
            rank = g ;
        }
    }

    if (CHECK_BLAS_INT && !cc->blas_ok)
    {
        
        return (0) ;
    }

    return (rank) ;
} //end of qr_front

// =============================================================================
// === qr_fcsize =============================================================
// =============================================================================
Long qr_fcsize    (
    Long m, Long n, 
    Long npiv, Long rank 
)
{
    Long cm, cn, csize ;
    cn = n - npiv ;                         
    cm = MIN (m-rank, cn) ;                 

    csize = (cm * (cm+1)) / 2 + cm * (cn - cm) ;
    return (csize) ;                        
} //end of qr_fcsize

// =============================================================================
// === qr_cpack ==============================================================
// =============================================================================
Long qr_cpack   (  
    Long m, Long n, Long npiv, Long rank,
    double *F, double *C 
)
{
    Long i, k, cm, cn ;

    // -------------------------------------------------------------------------
    // 获取输入
    // -------------------------------------------------------------------------

    cn = n - npiv ;                     
    cm = MIN (m-rank, cn) ;            
    if (cm <= 0 || cn <= 0)
    {
        return (0) ;                    
    }

    F += INDEX (rank,npiv,m) ;         

    // -------------------------------------------------------------------------
    // 压缩C的上三角形部分
    // -------------------------------------------------------------------------

    for (k = 0 ; k < cm ; k++)
    {
        for (i = 0 ; i <= k ; i++)
        {
            *(C++) = F [i] ;
        }
        F += m ;                       
    }

    // -------------------------------------------------------------------------
    // 压缩C的矩形部分
    // -------------------------------------------------------------------------

    for ( ; k < cn ; k++)
    {
        for (i = 0 ; i < cm ; i++)
        {
            *(C++) = F [i] ;
        }
        F += m ;                       
    }
    return (cm) ;                      
} //end of qr_cpack

// =============================================================================
// === qr_rhpack =============================================================
// =============================================================================
// 返回完整的 R+H 项
Long qr_rhpack   
(
    // input, not modified
    int keepH,              
    Long m,                
    Long n,                
    Long npiv,             
    Long *Stair,            

    // input, not modified (unless the pack occurs in-place)
    double *F,              

    // output, contents not defined on input
    double *R,               
    Long *p_rm             
)
{
    double *R0 = R ;
    Long i, k, h, t, rm ;

    // -------------------------------------------------------------------------
    // 获取输入
    // -------------------------------------------------------------------------

    if (m <= 0 || n <= 0)
    {
        *p_rm = 0 ;                     
        return (0) ;                   
    }


    // -------------------------------------------------------------------------
    // 压缩R的上三角形部分
    // -------------------------------------------------------------------------

    rm = 0 ;                            
    for (k = 0 ; k < npiv ; k++)
    {
        t = Stair [k] ;                
        if (t == 0)
        {
            t = rm ;                    
        }
        else if (rm < m)
        {
            rm++ ;                      
        }
        if (keepH)
        {
            
            for (i = 0 ; i < t ; i++)
            {
                *(R++) = F [i] ;
            }
        }
        else
        {
            for (i = 0 ; i < rm ; i++)
            {
                *(R++) = F [i] ;
            }
        }
        F += m ;                        
    }

    // -------------------------------------------------------------------------
    // 把R的矩形部分压缩
    // -------------------------------------------------------------------------

    h = rm ;                        
    for ( ; k < n ; k++)
    {

        for (i = 0 ; i < rm ; i++)
        {
            *(R++) = F [i] ;
        }

        if (keepH)
        {
            t = Stair [k] ;            
            h = MIN (h+1, m) ;          
            for (i = h ; i < t ; i++)
            {
                *(R++) = F [i] ;
            }
        }

        F += m ;                   
    }

    *p_rm = rm ;                       
    return (R-R0) ;                    
} // end of qr_rhpack

// =============================================================================
// === qr_larftb =============================================================
// =============================================================================
// 对矩阵应用一组Household反射。给定向量V和系数Tau，构造矩阵T，然后应用更新:
inline void qr_private_larft (char direct, char storev, Long n, Long k,
    double *V, Long ldv, double *Tau, double *T, Long ldt, sparse_common *cc)
{
    BLAS_INT N = n, K = k, LDV = ldv, LDT = ldt ;
    if (CHECK_BLAS_INT &&
        !(EQ (N,n) && EQ (K,k) && EQ (LDV,ldv) && EQ (LDT,ldt)))
    {
        cc->blas_ok = FALSE ;
    }
    if (!CHECK_BLAS_INT || cc->blas_ok)
    {
        #ifndef USEATL
            LAPACK_DLARFT (&direct, &storev, &N, &K, V, &LDV, Tau, T, &LDT) ;
        #else
            HNU_larftFC(LAForward, LAColumnStore, N, K, V, LDV, Tau, T, LDT);
        #endif
    }
}

inline void qr_private_larfb (char side, char trans, char direct, char storev,
    Long m, Long n, Long k, double *V, Long ldv, double *T, Long ldt, double *C,
    Long ldc, double *Work, Long ldwork, sparse_common *cc)
{
    BLAS_INT M = m, N = n, K = k, LDV = ldv, LDT = ldt, LDC = ldc,
        LDWORK = ldwork ;
    // 格式转换
#ifdef USEATL
    enum CBLAS_SIDE SIDE;
    enum CBLAS_TRANSPOSE TRANS;
    enum ATL_LADIRECT DIRECT;
    enum ATL_LASTOREV STOREV;
    if (side == 'L') SIDE = CblasLeft;
    else SIDE = CblasRight;

    if (trans =='T') TRANS = CblasTrans;
    else TRANS = CblasNoTrans;

    if (direct == 'F') DIRECT = LAForward;
    else DIRECT = LABackward;

    if (storev == 'C') STOREV = LAColumnStore;
    else STOREV = LARowStore;
#endif
    if (CHECK_BLAS_INT &&
        !(EQ (M,m) && EQ (N,n) && EQ (K,k) && EQ (LDV,ldv) &&
          EQ (LDT,ldt) && EQ (LDV,ldv) && EQ (LDWORK,ldwork)))
    {
        cc->blas_ok = FALSE ;
    }
    if (!CHECK_BLAS_INT || cc->blas_ok)
    {
        #ifndef USEATL
            LAPACK_DLARFB (&side, &trans, &direct, &storev, &M, &N, &K, V, &LDV,
                    T, &LDT, C, &LDC, Work, &LDWORK);
        #else
            HNU_larfb( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
                    T, LDT, C, LDC, Work, LDWORK);
        #endif
    }
}
// ====================================== //
void qr_larftb
(
    int method, Long m, Long n, Long k, 
    Long ldc, Long ldv, double *V, double *Tau, 
    double *C, double *W, sparse_common *cc
)
{
    double *T, *Work ;

    // -------------------------------------------------------------------------
    // 检查输入并分割工作区
    // -------------------------------------------------------------------------

    if (m <= 0 || n <= 0 || k <= 0)
    {
        return ; 
    }

    T = W ;             
    Work = W + k*k ;    

    // -------------------------------------------------------------------------
    // 构造并应用k×k上三角矩阵T
    // -------------------------------------------------------------------------

    if (method == QR_QTX)
    {
        qr_private_larft ('F', 'C', m, k, V, ldv, Tau, T, k, cc) ;

        qr_private_larfb ('L', 'T', 'F', 'C', m, n, k, V, ldv, T, k, C, ldc,
            Work, n, cc) ;
    }
    else if (method == QR_QX)
    {
        qr_private_larft ('F', 'C', m, k, V, ldv, Tau, T, k, cc) ;

        qr_private_larfb ('L', 'N', 'F', 'C', m, n, k, V, ldv, T, k, C, ldc,
            Work, n, cc) ;
    }
    else if (method == QR_XQT)
    {
        qr_private_larft ('F', 'C', n, k, V, ldv, Tau, T, k, cc) ;

        qr_private_larfb ('R', 'T', 'F', 'C', m, n, k, V, ldv, T, k, C, ldc,
            Work, m, cc) ;
    }
    else if (method == QR_XQ)
    {
        qr_private_larft ('F', 'C', n, k, V, ldv, Tau, T, k, cc) ;

        qr_private_larfb ('R', 'N', 'F', 'C', m, n, k, V, ldv, T, k, C, ldc,
            Work, m, cc) ;
    }
} //end of qr_larftb