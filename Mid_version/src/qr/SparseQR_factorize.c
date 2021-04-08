/******************************************************************************
 * VERSION: 1.1
 * DATE:    2020年9月27日
 * FILE:    SparseQR_factorize.c
 * BRIEF:   数值分解
 * FUNCTION: factorize函数，用来计算QR分解的数值分解部分，可以通过线程池任务级并行
 *****************************************************************************/

#include "SparseQR.h"
#include "tpsm.h"
#define FCHUNK 32       
#define SMALL 5000
#define MINCHUNK 4
#define MINCHUNK_RATIO 4

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

// =============================================================================
// === get_Work ================================================================
// =============================================================================

// 分配工作对象，它是一个结构体数组。
// 元素 Work[s]包含堆栈s的工作空间，对于每个堆栈0到ns-1。

qr_work *get_Work
(
    Long ns,            
    Long n,             
    Long maxfn,        
    Long keepH,         
    Long fchunk,
    Long *p_wtsize,     
    sparse_common *cc
)
{
    int ok = TRUE ;
    qr_work *Work ;
    Long wtsize ;
    *p_wtsize = 0 ;

    wtsize = qr_mult (fchunk + (keepH ? 0:1), maxfn, &ok) ;

    Work = (qr_work *)    
        SparseCore_malloc (ns, sizeof (qr_work), cc) ;

    if (!ok || cc->status < SPARSE_OK)
    {
        SparseCore_free (ns, sizeof (qr_work), Work, cc) ;
        
        return (NULL) ;
    }

    for (Long stack = 0 ; stack < ns ; stack++)
    {
        Work [stack].Fmap = (Long *) SparseCore_malloc (n, sizeof (Long), cc) ;
        Work [stack].Cmap = (Long *) SparseCore_malloc (maxfn, sizeof(Long), cc);
        if (keepH)
        {
            Work [stack].Stair1 = NULL ;
        }
        else
        {
            Work [stack].Stair1 =
                (Long *) SparseCore_malloc (maxfn, sizeof (Long), cc) ;
        }
        Work [stack].WTwork =
            (double *) SparseCore_malloc (wtsize, sizeof (double), cc) ;
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
            SparseCore_free (n,      sizeof (Long),   Work [stack].Fmap,   cc) ;
            SparseCore_free (maxfn,  sizeof (Long),   Work [stack].Cmap,   cc) ;
            SparseCore_free (maxfn,  sizeof (Long),   Work [stack].Stair1, cc) ;
            SparseCore_free (wtsize, sizeof (double), Work [stack].WTwork, cc) ;
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
    TPSM_t *tpool, sparse_csc **Ahandle, Long freeA, 
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

    // -------------------------------------------------------------------------
    // 获取符号对象的输入和内容
    // -------------------------------------------------------------------------

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
    Cm = Wi ;                  

    Cblock = (double **) SparseCore_malloc (nf+1, sizeof (double *), cc) ;

    Work = NULL ;              
    fchunk = MIN (m, FCHUNK) ;
    wtsize = 0 ;

    // -------------------------------------------------------------------------
    // 创建S = AE
    // -------------------------------------------------------------------------

    Sx = (double *) SparseCore_malloc (anz, sizeof (double), cc) ;

    if (cc->status == SPARSE_OK)
    {
        qr_stranspose2 (A, Qfill, Sp, PLinv, Sx, Wi) ;
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
    // 分配数值对象
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

    QRnum->Rblock     = Rblock ;
    QRnum->Rdead      = Rdead ;
    QRnum->Stacks     = Stacks ;
    QRnum->Stack_size = Stack_size ;


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
    // 串行，只有一个栈
    if (cc->status == SPARSE_OK)
    {
        for (stack = 0 ; stack < ns ; stack++)
        {
            double *Stack ;
            size_t stacksize = (ntasks == 1) ?
                maxstack : Stack_maxstack [stack] ;
            Stack_size [stack] = stacksize ;
            Stack = (double *) SparseCore_malloc (stacksize, sizeof (double), cc) ;
            Stacks [stack] = Stack ;
            Work [stack].Stack_head = Stack ;
            Work [stack].Stack_top  = Stack + stacksize ;
        }
    }


    // -------------------------------------------------------------------------
    // 如果内存不足，则使用顺序 case和fchunk = 1
    // -------------------------------------------------------------------------

    if (cc->status < SPARSE_OK)
    {
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
    Blob.Cm = Cm ;
    Blob.Cblock = Cblock ;
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


    if (ntasks == 1)
    {
        
        qr_kernel (0, &Blob) ;     
    }
    else
    {
// #ifdef MULTI
        qr_multithreads(ntasks, &Blob, tpool);
// #else

//     for (Long id = 0 ; id < ntasks-1 ; id++)
//     {
//         qr_kernel (id, &Blob) ;
//     }
// #endif
    }

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

    rank = 0 ;
    maxfrank = 1 ;
    for (stack = 0 ; stack < ns ; stack++)
    {
        rank += Work [stack].sumfrank ;
        maxfrank = MAX (maxfrank, Work [stack].maxfrank) ;
    }
    QRnum->rank = rank ;                  
    QRnum->maxfrank = maxfrank ;

    // -------------------------------------------------------------------------
    // 最终死列的二范数
    // -------------------------------------------------------------------------

    double wscale = 0 ;
    double wssq = 1 ;
    for (stack = 0 ; stack < ns ; stack++)
    {
        double ws = Work [stack].wscale ;
        double wq = Work [stack].wssq ;
        if (wq != 0)
        {
            double wk = ws * sqrt (wq) ;
            if (wscale < wk)
            {
                double rr = wscale / wk ;
                wssq = 1 + wssq * rr * rr ;
                wscale = wk ;
            }
            else
            {
                double rr = wk / wscale ;
                wssq += rr * rr ;
            }
        }
    }
    QRnum->norm_E_fro = wscale * sqrt (wssq) ;
    cc->SPQR_norm_E_fro = QRnum->norm_E_fro ;

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

    int shrink = cc->SPQR_shrink ;

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
                    SparseCore_free (stacksize, sizeof (double), Stack, cc) ;
                }
                stacksize = newstacksize ;
            }
            else
            {
                Cblock [stack] =  
                    (double *) SparseCore_realloc (
                    newstacksize,   
                    sizeof (double), 
                    Stack,          
                    &stacksize,     
                    cc) ;
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
    return (QRnum) ;
}

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
            i = Ai [p] ;                
            row = PLinv [i] ;           
            s = W [row]++ ;             
            Sx [s] = Ax [p] ;
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

    Long *  Super = QRsym->Super ;      
    Long *  Rp = QRsym->Rp ;           
    Long *  Rj = QRsym->Rj ;            
    Long *  Sleft = QRsym->Sleft ;      
    Long *  Sp = QRsym->Sp ;            
    Long *  Sj = QRsym->Sj ;           
    Long *  Child = QRsym->Child ;      
    Long *  Childp = QRsym->Childp ;    
    Long    maxfn  = QRsym->maxfn ;     
    Long    nf = QRsym->nf ;            

    Long *  Hip = QRsym->Hip ;          

    Long *  TaskFront = QRsym->TaskFront ;      
    Long *  TaskFrontp = QRsym->TaskFrontp ;    
    Long *  TaskStack = QRsym->TaskStack ;      
    Long *  On_stack = QRsym->On_stack ;        

    Long *  Post = QRsym->Post ;                

    double ** Rblock = QRnum->Rblock ;
    char *   Rdead = QRnum->Rdead ;
    Long *   HStair = QRnum->HStair ;
    double *  HTau = QRnum->HTau ;
    Long *   Hii = QRnum->Hii ;          
    Long *   Hm = QRnum->Hm ;
    Long *   Hr = QRnum->Hr ;
    Long     keepH = QRnum->keepH ;
    Long     ntasks = QRnum->ntasks ;    

    Long stack, kfirst, klast ;

    if (ntasks == 1)
    {
        kfirst = 0 ;
        klast  = nf ;
        stack  = 0 ;
    }
    else
    {
        kfirst = TaskFrontp [task] ;
        klast  = TaskFrontp [task+1] ;
        stack  = TaskStack [task] ;
    }

    double * Stack_head = Work [stack].Stack_head ;
    double * Stack_top = Work [stack].Stack_top ;

    double * Tau = keepH ? NULL : Work [stack].WTwork ;
    Long *  Stair = keepH ? NULL : Work [stack].Stair1 ;
    double * W = Work [stack].WTwork + (keepH ? 0 : maxfn) ;

    Long *  Fmap = Work [stack].Fmap ;
    Long *  Cmap = Work [stack].Cmap ;

    Long    sumfrank = Work [stack].sumfrank ;
    Long    maxfrank = Work [stack].maxfrank ;

    double wscale = Work [stack].wscale ;
    double wssq   = Work [stack].wssq   ;

    for (Long kf = kfirst ; kf < klast ; kf++)
    {
        Long f = (ntasks == 1) ? Post [kf] : TaskFront [kf] ;

        if (keepH)
        {
            Stair = HStair + Rp [f] ;
            Tau = HTau + Rp [f] ;
        }

        Long fm = qr_fsize (f, Super, Rp, Rj, Sleft, Child, Childp, Cm,
            Fmap, Stair) ;
        Long fn = Rp [f+1] - Rp [f] ;       
        Long col1 = Super [f] ;              
        Long fp = Super [f+1] - col1 ;       
        Long fsize = fm * fn ;
        if (keepH)
        {
            Hm [f] = fm ;
        }

        double *F = Stack_head ;
        Rblock [f] = F ;
        Stack_head += fsize ;

        qr_assemble (f, fm, keepH,
            Super, Rp, Rj, Sp, Sj, Sleft, Child, Childp,
            Sx, Fmap, Cm, Cblock,
            Hr, Stair, Hii, Hip, F, Cmap) ;

        for (Long p = Childp [f] ; p < Childp [f+1] ; p++)
        {
            Long c = Child [p] ;
            if (ntasks == 1 || On_stack [c] == stack)
            {
                Long ccsize = qr_csize (c, Rp, Cm, Super) ;
                Stack_top = MAX (Stack_top, Cblock [c] + ccsize) ;
            }
        }

        Long frank = qr_front (fm, fn, fp, tol, ntol - col1,
            fchunk, F, Stair, Rdead + col1, Tau, W,
            &wscale, &wssq, cc) ;

        sumfrank += frank ;
        maxfrank = MAX (maxfrank, frank) ;

        Long csize = qr_fcsize (fm, fn, fp, frank) ;
        Stack_top -= csize ;

        Cblock [f] = Stack_top ;
        Cm [f] = qr_cpack (fm, fn, fp, frank, F, Cblock [f]) ;

        Long rm ;
        Long rsize = qr_rhpack (keepH, fm, fn, fp, Stair, F, F, &rm) ;
        if (keepH)
        {
            Hr [f] = rm ;
        }

        Stack_head -= fsize ;             
        Stack_head += rsize ;             

    }

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
            Stair [j]++ ;               
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
    minchunk = MAX (MINCHUNK, fchunk/MINCHUNK_RATIO) ;
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
                k2 = MIN (n, k+fchunk) ;        
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