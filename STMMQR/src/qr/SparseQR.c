/******************************************************************************
 * VERSION: 1.1
 * DATE:    2020年9月27日
 * FILE:    SparseQR.c
 * BRIEF:   QR 分解函数
 * FUNCTION:SparseQR函数，用来计算矩阵的QR分解，并保留Household变换
 *****************************************************************************/
/*******************************
 *         INCLUDE
 ******************************/ 
#ifdef PRINT_TIME
#include <sys/time.h>
#endif
#include "SparseQR.h"
#include <float.h>
#include "tpsm.h"

// =============================================================================
// === SparseQR ============================================================
// =============================================================================
/*  计算稀疏矩阵的QR分解，包括符号矩阵和数字矩阵。这个函数利用列单例，因此不能分割为符号和数字阶段。
    function [C,R,E,X,H] = sparseQR (A,B)
    E = A 的填充简化排序, 先排除单例
    S = A*E

    对Y进行QR分解，得到因子R，Household 变换 H，
        并且将H应用到 B2 到 C2.
        X = E*(R\C)
    仅分解 时 B 为 NULL
    
    首先，找到稀疏矩阵A的列单例。如果一个矩阵的列单例及其对应的行重排列到左侧和顶部，
    且该矩阵不存在秩亏，则结果如下:
        c x x x x x x x
        . c x x x x x x
        . . c x x x x x
        . . . c x x x x
        . . . . a a a a
        . . . . a a a a
        . . . . a a a a
        . . . . a a a a
        . . . . a a a a
    
    c 表示一个列单例j 和它的对应行i。x表示单例行中的一个条目，a表示剩余矩阵.
    在 S2 子矩阵("a"项)中没有一列只含一个项，除非这项的值低于tol.

    如果一个单例元素c的大小低于tol，那么则不接受它。 如果 tol为负数，矩阵A中
    任意单例元素都是可接受的，甚至是显示存储的 0 项。

    如果矩阵在结构上秩亏，那么排序后的矩阵可能形如下，其中第三列称为"死列":
        c x x x x x x x
        . c x x x x x x
        . . . c x x x x
        . . . . a a a a
        . . . . a a a a
        . . . . a a a a
        . . . . a a a a
        . . . . a a a a
        . . . . a a a a

    在找到列单例后， 剩余的S2矩阵(连同R12矩阵)做一个填充简化排序，得到排序E。
    B 的 列不排序。

    矩阵 [AE B] 分为了2 x 3 的块形式:
        R11 R12 B1
        0   S2  B2
    
    R1 = [R11 R12] 是 A的行单例，[R11 ; 0] 是列单例。R1以行的形式存储(带有已排序的行)
    Y = [S2 B2] 以列的形式存储(如果输入矩阵A的列是有序的则它的列也是有序的).
    S2不包含空列，而且对于S2中任何只有一个项的列，其范数都小于tol。
    注意，Y的行是根据单例行进行排列的。
*/
SparseQR_factorization *SparseQR(
    int ordering, double tol,  
    sparse_csc *A, sparse_common *cc, char *result_name )
{
    #ifdef PRINT_TIME
    double timeStart, timeEnd;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
    #endif
    qr_symbolic *QRsym ;
    qr_numeric *QRnum ;
    SparseQR_factorization *QR ;
    Long *Yp, *Yi, *Q1fill, *R1p, *R1j, *P1inv, *Ap, *Ai ;
    double *Yx, *R1x, *Ax ;
    Long noY, anz, a2nz, r1nz, ynz, i, j, k, p, p2, bnz, py, n1rows,
        n1cols, n2, d, iold, inew, m, n ;
    sparse_csc *Y = NULL ;
    // -------------------------------------------------------------------------
    // 获取输入值并分配结果
    // -------------------------------------------------------------------------

    m = A->nrow ;
    n = A->ncol ;
    Ap = (Long *) A->p ;
    Ai = (Long *) A->i ;
    Ax = (double *) A->x ;

    QR = (SparseQR_factorization *)
        SparseCore_malloc (1, sizeof (SparseQR_factorization), cc) ;

    if (cc->status < SPARSE_OK)  
    {
        // 内存越界
        return (NULL) ;
    }
    QR->narows = m ;
    QR->nacols = n ;
    QR->n1rows = 0 ;
    QR->n1cols = 0 ;
    QR->r1nz = 0 ;
    r1nz = 0 ;
    QR->QRsym = NULL ;
    QR->QRnum = NULL ;

    QR->R1p = NULL ;
    QR->R1j = NULL ;
    QR->R1x = NULL ;
    QR->P1inv = NULL ;
    QR->Q1fill = NULL ;
    QR->Rmap = NULL ;
    QR->RmapInv = NULL ;
    QR->HP1inv = NULL ;

    QR->Ana_time = 0;
    QR->Fac_time = 0;

    // -------------------------------------------------------------------------
    // 检查阈值 tol
    // -------------------------------------------------------------------------

    if (tol <= -2)  
    {
        tol = qr_tol (A, cc) ;  
    }
    if (tol < 0)
    {
        QR->allow_tol = FALSE ;
        tol = EMPTY ;
    }
    else
    {
        QR->allow_tol = TRUE ;
    }
    QR->tol = tol ;
    
    // -------------------------------------------------------------------------
    // 寻找单例 ， 构造Y 的在矩阵A部分的列指针   1.1
    // -------------------------------------------------------------------------
    // 寻找列单例并排列矩阵 A     
    qr_1colamd (ordering, tol, A, &Q1fill,
        &R1p, &P1inv, &Y, &n1cols, &n1rows, cc, result_name) ; //1.2
    ordering = cc->SPQR_istat [7]  ;

    if (cc->status < SPARSE_OK)  // 再次检查内存
    {
        qr_freefac (&QR, cc) ;
        return (NULL) ;
    }

    // 通过之前寻找单例后 ，保存下来
    QR->R1p = R1p ;
    QR->P1inv = P1inv ;
    QR->Q1fill = Q1fill ; 
    QR->n1rows = n1rows ;
    QR->n1cols = n1cols ;

    noY = (Y == NULL) ;                         
    
    Yp = noY ? NULL : (Long *) Y->p ;
    anz = Ap [n] ;                              // A中的非零元素
    a2nz = noY ? anz : Yp [n-n1cols] ;          // S2中的非零元素
    n2 = n - n1cols ;                           // S2中的列数
    ynz = 0 ;
    // -------------------------------------------------------------------------
    // 构建 删除列单例后 的下部分Y的指针
    // -------------------------------------------------------------------------
    if (noY)
    {
        
    }
    else if (n1cols == 0)
    {
        ynz = a2nz ;
        Yp [(n-n1cols)] = ynz ;
        
    }
    else
    {
        ynz = a2nz ;
        Yp [(n-n1cols)] = ynz ;
    }
    
    // -------------------------------------------------------------------------
    //  给 Y 部分分配非零元
    // -------------------------------------------------------------------------
    if (noY) 
    {
        // 没有找到单例
        Yi = NULL ;
        Yx = NULL ;
    }
    else   
    {
        SparseCore_reallocate_sparse (ynz, Y, cc) ;
        Yi = (Long  *) Y->i ;
        Yx = (double *) Y->x ;
    }
    
    if (cc->status < SPARSE_OK) // 再次检查状态内存
    {
        // 内存越界
        qr_freefac (&QR, cc) ;
        SparseCore_free_sparse (&Y, cc) ;
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // 创建Y和R1的模式和值
    // -------------------------------------------------------------------------
    if (noY) 
    {

        // ---------------------------------------------------------------------
        // R1不存在
        // ---------------------------------------------------------------------

        R1j = NULL ;
        R1x = NULL ;

    }
    else if (n1cols == 0)
    {

        // ---------------------------------------------------------------------
        //  R1不存在
        // ---------------------------------------------------------------------

        R1j = NULL ;
        R1x = NULL ;

        // ---------------------------------------------------------------------
        // 构造 A 的这一部分: Y = [S B]
        // ---------------------------------------------------------------------

        py = 0 ;
        for (k = 0 ; k < n ; k++)
        {
            j = Q1fill ? Q1fill [k] : k ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                Yi [py] = Ai [p] ;
                Yx [py] = Ax [p] ;
                py++ ;
            }
        }
    }
    else    // 有列单例划分
    {
        r1nz = qr_cumsum (n1rows, R1p) ;      // 不能发生 长溢出

        // ---------------------------------------------------------------------
        // 分配R1空间
        // ---------------------------------------------------------------------

        R1j = (Long  *) SparseCore_malloc (r1nz, sizeof (Long ), cc) ;
        R1x = (double *) SparseCore_malloc (r1nz, sizeof (double), cc) ;
        QR->R1j = R1j ;
        QR->R1x = R1x ;
        QR->r1nz = r1nz ;

        if (cc->status < SPARSE_OK)
        {
            // 内存溢出
            qr_freefac (&QR, cc) ;
            SparseCore_free_sparse (&Y, cc) ;
            return (NULL) ;
        }

        // ---------------------------------------------------------------------
        // 扫描 A ,构造 R11
        // ---------------------------------------------------------------------

        for (k = 0 ; k < n1cols ; k++)
        {
            j = Q1fill ? Q1fill [k] : k ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                // A的第i行是在单列置换之后的第inew 行
                i = Ai [p] ;
                inew = P1inv [i] ;
                // A (i,j) 是单例所在的行.  转换成了 R1 (inew,k)
                p2 = R1p [inew]++ ;
                R1j [p2] = k ;
                R1x [p2] = Ax [p] ;
            }
        }

        // ---------------------------------------------------------------------
        // 扫描A，构造R12, Y的S2部分= [S2 B2]
        // ---------------------------------------------------------------------

        py = 0 ;
        for ( ; k < n ; k++)
        {
            j = Q1fill ? Q1fill [k] : k ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                // A的第i行是在单列置换之后的第inew 行
                i = Ai [p] ;
                inew = P1inv [i] ;
                if (inew < n1rows)
                {
                    // A (i,j) 是单例所在的行.  转换成了 R1 (inew,k)
                    p2 = R1p [inew]++ ;
                    R1j [p2] = k ;
                    R1x [p2] = Ax [p] ;
                }
                else
                {
                    // A (i,j) 不在单例中.  将他放置在Y (inew-n1rows, k-n1cols)
                    Yi [py] = inew - n1rows ;
                    Yx [py] = Ax [p] ;
                    py++ ;
                }
            }
        }

        // ---------------------------------------------------------------------
        //  恢复R1的行指针
        // ---------------------------------------------------------------------

        qr_shift (n1rows, R1p) ;
    }

    // -------------------------------------------------------------------------
    //  对A或Y进行QR分解
    // -------------------------------------------------------------------------
    
    if (noY) 
    {
        // 将A因式分解，在Q1fill中已经给出了填充减少排序
        QRsym = qr_analyze (A, QR_ORDERING_GIVEN, Q1fill,tol >= 0, cc, result_name) ; //符号分析
        #ifdef PRINT_TIME
        gettimeofday(&tv, NULL);
        timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        printf ("\nAnalyze time: %lf\n", timeEnd - timeStart);

        QR->Ana_time = timeEnd - timeStart;

        gettimeofday(&tv, NULL);
        timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
        #endif
        QRnum = qr_factorize ( &A, FALSE, tol, n, QRsym, cc) ; //数值分解
        #ifdef PRINT_TIME
        gettimeofday(&tv, NULL);
        timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        printf ("Factorize time: %lf\n", timeEnd - timeStart);

        QR->Fac_time = timeEnd - timeStart;
        #endif
    }
    else
    {
        QRsym = qr_analyze (Y, QR_ORDERING_FIXED, NULL, tol >= 0, cc, result_name) ; //符号分析
        #ifdef PRINT_TIME
        gettimeofday(&tv, NULL);
        timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        printf ("\nAnalyze time: %lf\n", timeEnd - timeStart);

        QR->Ana_time = timeEnd - timeStart;

        gettimeofday(&tv, NULL);
        timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
        #endif
        QRnum = qr_factorize ( &Y, TRUE, tol, n2, QRsym, cc) ;
        #ifdef PRINT_TIME
        gettimeofday(&tv, NULL);
        timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        printf ("Factorize time: %lf\n", timeEnd - timeStart);

        QR->Fac_time = timeEnd - timeStart;
        #endif
    }

    cc->SPQR_istat [7] = ordering ;

    QR->QRsym = QRsym ;
    QR->QRnum = QRnum ;

    if (cc->status < SPARSE_OK) //检查内存
    {
        // 内存溢出
        qr_freefac (&QR, cc) ;
        return (NULL) ;
    }

    cc->SPQR_istat [0] += r1nz ; 
    QR->rank = n1rows + QRnum->rank1 ;

    // -------------------------------------------------------------------------
    // 如果保存 H 且单例存在，则构造一个全局行置换
    // -------------------------------------------------------------------------
    if ( n1cols > 0)
    {

        Long kk ;
        Long *HP1inv, *HPinv ;
        QR->HP1inv = HP1inv = (Long *) SparseCore_malloc (m, sizeof (Long), cc) ;
        HPinv = QRnum->HPinv ;

        if (cc->status < SPARSE_OK)
        {
            // 内存溢出
            qr_freefac (&QR, cc) ;
            return (NULL) ;
        }

        for (i = 0 ; i < m ; i++)
        {
            // i 是A的一行，k是一个经过行单例重排序后的行索引。则kk是R的行索引
            k = P1inv ? P1inv [i] : i ;
            if (k < n1rows)
            {
                kk = k ;
            }
            else
            {
                kk = HPinv [k - n1rows] + n1rows ;
            }
            HP1inv [i] = kk ;
        }
    }

    // -------------------------------------------------------------------------
    // 如果 A 是秩亏的， 则找到压缩 R 的映射
    // -------------------------------------------------------------------------

    if (QR->rank < n && !qr_rmap (QR, cc))
    {
        qr_freefac (&QR, cc) ;
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // 保存输出数据
    // -------------------------------------------------------------------------

    cc->SPQR_istat [4] = QR->rank ;         //  A 的预估秩
    cc->SPQR_istat [5] = n1cols ;           // 列单例数目
    cc->SPQR_istat [6] = n1rows ;           // 行单例数目
    cc->SPQR_tol_used = tol ;               // 使用的tol

    return (QR) ;
} // end of SparseQR

// ===================================
// === SparseQR_free =================
// ===================================
// 释放 QR 结构体对象
int SparseQR_free(
    SparseQR_factorization **QR,sparse_common *cc
)
{
    qr_freefac (QR, cc) ;
    return (TRUE) ;
} //end of SparseQR_free

// ===================================
// =========== qr_freefac ============
// ===================================
void qr_freefac(
    SparseQR_factorization **QR_handle,

    // 工作区域
    sparse_common *cc
)
{
    SparseQR_factorization *QR ;
    Long n, m, n1rows, r1nz ;

    if (QR_handle == NULL || *QR_handle == NULL)
    {
        return ;
    }
    QR = *QR_handle ;

    n      = QR->nacols ;
    m      = QR->narows ;
    n1rows = QR->n1rows ;
    r1nz   = QR->r1nz ;

    qr_freenum (& (QR->QRnum), cc) ;
    qr_freesym (& (QR->QRsym), cc) ;

    SparseCore_free (n, sizeof (Long),  QR->Q1fill,  cc) ; 
    SparseCore_free (m,        sizeof (Long),  QR->P1inv,   cc) ;
    SparseCore_free (m,        sizeof (Long),  QR->HP1inv,  cc) ;
    SparseCore_free (n1rows+1, sizeof (Long),  QR->R1p,     cc) ;
    SparseCore_free (r1nz,     sizeof (Long),  QR->R1j,     cc) ;
    SparseCore_free (r1nz,     sizeof (double), QR->R1x,     cc) ;
    SparseCore_free (n,        sizeof (Long),  QR->Rmap,    cc) ;
    SparseCore_free (n,        sizeof (Long),  QR->RmapInv, cc) ;

    SparseCore_free (1, sizeof (SparseQR_factorization), QR, cc) ;
    *QR_handle = NULL ;
} // end of qr_freefac

// ======================================
// === qr_1colamd =======================
// ======================================
/*
    查找列单例，允许列排列。在列单例找到之后(先使用Q1fill排序)
    其余的列可以可选地通过COLAMD或CHOLMOD的内部排序方法进行排列。
    这个函数处理natural、COLAMD和CHOLMOD排序选项(而不是固定排序选项)。

    返回一个稀疏矩阵Y，其中包含分配和初始化的列指针，但没有值。
    Y有n-n1cols列，和m-n1rows行。B是空的，没有找到单例，Y是NULL。
*/
int qr_1colamd  // 如果OK为真，否则为假
(
    int ordering, double tol, sparse_csc *A, 
    Long **p_Q1fill, Long **p_R1p, Long **p_P1inv, 
    sparse_csc **p_Y, Long *p_n1cols, Long *p_n1rows, 
    sparse_common *cc, char *result_name
)
{
    Long *Q1fill, *Degree, *Qrows, *W, *Winv, *ATp, *ATj, *R1p, *P1inv, *Yp,
        *Ap, *Ai, *Work ;
    double *Ax ;
    Long p, d, j, i, k, n1cols, n1rows, row, pend, n2rows, n2cols = EMPTY,
        nz2, kk, p2, col2, ynz, fill_reducing_ordering, m, n, xtype, worksize ;
    sparse_csc *AT, *Y ;

    // -------------------------------------------------------------------------
    // 获取输入的值
    // -------------------------------------------------------------------------

    xtype = 1;

    m = A->nrow ;
    n = A->ncol ;
    Ap = (Long *) A->p ;
    Ai = (Long *) A->i ;
    Ax = (double *) A->x ;

    // 将输出设置为NULL
    *p_Q1fill = NULL ;
    *p_R1p = NULL ;
    *p_P1inv = NULL ;
    *p_Y = NULL ;
    *p_n1cols = EMPTY ;
    *p_n1rows = EMPTY ;

    // -------------------------------------------------------------------------
    // 分配结果Q1fill (Y, R1p, P1inv 在之后分配)
    // -------------------------------------------------------------------------

    Q1fill = (Long *) SparseCore_malloc (n, sizeof (Long), cc) ;

    // -------------------------------------------------------------------------
    // 分配空间
    // -------------------------------------------------------------------------

    fill_reducing_ordering = !
        ((ordering == QR_ORDERING_FIXED) ||
         (ordering == QR_ORDERING_GIVEN) ||
         (ordering == QR_ORDERING_NATURAL)) ;

    worksize = ((fill_reducing_ordering) ? 3:2) * n ;

    Work = (Long *) SparseCore_malloc (worksize, sizeof (Long), cc) ;
    Degree = Work ;        
    Qrows  = Work + n ;    
    Winv   = Qrows ;        
    W      = Qrows + n ;    

    if (cc->status < SPARSE_OK)
    {
        // 内存溢出; 返回
        SparseCore_free (worksize, sizeof (Long), Work, cc) ;
        SparseCore_free (n, sizeof (Long), Q1fill, cc) ;
        return (FALSE) ;
    }

    // -------------------------------------------------------------------------
    // 初始化具有空列的队列，以及只有一个条目的列
    // -------------------------------------------------------------------------

    n1cols = 0 ;
    n1rows = 0 ;

    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;
        d = Ap [j+1] - p ;
        if (d == 0)
        {
            // j是一个空列
            Q1fill [n1cols] = j ;
            Qrows [n1cols] = EMPTY ; //空列的行号为-1
            n1cols++ ;
            Degree [j] = EMPTY ;
        }
        else if (d == 1 && qr_abs (Ax [p]) > tol)
        {
            // j是一个列单例
            Q1fill [n1cols] = j ;   // 记录下列号
            Qrows [n1cols] = Ai [p] ;       // 这是一个副本 //记录单例行号
            n1cols++ ;   // 单例列加一
            Degree [j] = EMPTY ;
        }
        else
        {
            // j的度> 1，它还不是一个单例
            // 如果j在单例队列中则 Degree [j] = EMPTY, 否则 Degree [j]>1 是列的j 的度
            Degree [j] = d ;
        }
    }

    // -------------------------------------------------------------------------
    // 创建 AT = A'
    // -------------------------------------------------------------------------

    AT = SparseCore_transpose (A, 0, cc) ;       // 

    if (cc->status < SPARSE_OK)
    {
        // 内存溢出; 返回
        SparseCore_free (worksize, sizeof (Long), Work, cc) ;
        SparseCore_free (n, sizeof (Long), Q1fill, cc) ;
        return (FALSE) ;
    }

    ATp = (Long *) AT->p ;
    ATj = (Long *) AT->i ;

    // -------------------------------------------------------------------------
    // 通过宽度优先搜索删除列单例
    // -------------------------------------------------------------------------

    for (k = 0 ; k < n1cols ; k++)
    {

        // ---------------------------------------------------------------------
        // 从队列中获取新的单例
        // ---------------------------------------------------------------------
        #define col (Q1fill [k])

        row = Qrows [k] ;

        if (row == EMPTY || ATp [row] < 0)
        {

            // -----------------------------------------------------------------
            // col 是死的列单例; 移除行索引的副本
            // -----------------------------------------------------------------
            Qrows [k] = EMPTY ;
            row = EMPTY ;
        }
        else
        {

            // -----------------------------------------------------------------
            // col 是活的列单例; 从矩阵中删除它的行
            // -----------------------------------------------------------------

            n1rows++ ;
            p = ATp [row] ;
            ATp [row] = FLIP (p) ;          // 标记单例行
            pend = UNFLIP (ATp [row+1]) ;
            for ( ; p < pend ; p++)
            {
                // 在行被删除后寻找新的列单例
                j = ATj [p] ;
                d = Degree [j] ;
                if (d == EMPTY)
                {
                    // j 已经在单例列中
                    continue ;
                }
                d-- ;
                Degree [j] = d ;
                if (d == 0)
                {
                    // 一个新的死列单例
                    Q1fill [n1cols] = j ;
                    Qrows [n1cols] = EMPTY ;
                    n1cols++ ;
                    Degree [j] = EMPTY ;
                }
                else if (d == 1)
                {
                    // 一个活的列单例，找到其单独的活行
                    for (p2 = Ap [j] ; p2 < Ap [j+1] ; p2++)
                    {
                        i = Ai [p2] ;
                        if (ATp [i] >= 0 && qr_abs (Ax [p2]) > tol)
                        {
                            // i 可能出现在 Qrows [k+1:n1cols-1]
                            Q1fill [n1cols] = j ;
                            Qrows [n1cols] = i ;
                            n1cols++ ;
                            Degree [j] = EMPTY ;
                            break ;
                        }
                    }
                }
            }
        }
        #undef col
    }

    // -------------------------------------------------------------------------
    // 查找行 序 
    // -------------------------------------------------------------------------

    if (n1cols == 0)
    {

        // ---------------------------------------------------------------------
        // 矩阵中没有单例; 没有R1矩阵，没有P1inv 重排序
        // ---------------------------------------------------------------------

        R1p = NULL ;
        P1inv = NULL ;

    }
    else
    {

        // ---------------------------------------------------------------------
        // 构造行单例 排列
        // ---------------------------------------------------------------------

        // 分配结果数组R1p和P1inv的空间
        R1p   = (Long *) SparseCore_malloc (n1rows+1, sizeof (Long), cc) ;
        P1inv = (Long *) SparseCore_malloc (m,        sizeof (Long), cc) ;

        if (cc->status < SPARSE_OK)
        {
            // 越界，返回
            SparseCore_free_sparse (&AT, cc) ;
            SparseCore_free (worksize, sizeof (Long), Work, cc) ;
            SparseCore_free (n, sizeof (Long), Q1fill, cc) ;
            SparseCore_free (n1rows+1, sizeof (Long), R1p, cc) ;
            SparseCore_free (m,        sizeof (Long), P1inv, cc) ;
            return (FALSE) ;
        }

        kk = 0 ;
        for (k = 0 ; k < n1cols ; k++)
        {
            i = Qrows [k] ;
            if (i != EMPTY)
            {
                // 第i行是第k个单例行
                P1inv [i] = kk ;
                R1p [kk] = UNFLIP (ATp [i+1]) - UNFLIP (ATp [i]) ;
                kk++ ;
            }
        }
        for (i = 0 ; i < m ; i++)
        {
            if (ATp [i] >= 0)
            {
                P1inv [i] = kk ;
                kk++ ;
            }
        }

    }

    // -------------------------------------------------------------------------
    // 完成列排序
    // -------------------------------------------------------------------------

    if (!fill_reducing_ordering)
    {

        // ---------------------------------------------------------------------
        // 自然排序
        // ---------------------------------------------------------------------

        if (n1cols == 0)
        {

            // 没有单例，所以现在的自然顺序是0:n-1
            for (k = 0 ; k < n ; k++)
            {
                Q1fill [k] = k ;
            }

        }
        else
        {

            //首先出现单例列，然后出现非列单例列
            k = n1cols ;
            for (j = 0 ; j < n ; j++)
            {
                if (Degree [j] > 0)
                {
                    // 列j不是一个单列
                    Q1fill [k++] = j ;
                }
            }
        }

    }
    else
    {

        // ---------------------------------------------------------------------
        // 修剪子矩阵的填充减少排序
        // ---------------------------------------------------------------------

        if (n1cols == 0)
        {

            // -----------------------------------------------------------------
            // 没有单例;在整个矩阵上做填充减少排序
            // -----------------------------------------------------------------

            n2cols = n ;
            n2rows = m ;

        }
        else
        {

            // -----------------------------------------------------------------
            // 通过删除单例来创建精简矩阵以减少填充
            // -----------------------------------------------------------------

            // 查找原始列到经过修剪的列的映射
            n2cols = 0 ;
            for (j = 0 ; j < n ; j++)
            {
                if (Degree [j] > 0)
                {
                    // 列j不是一个单列
                    W [j] = n2cols++ ;
                }
                else
                {
                    // 列j是一个单列
                    W [j] = EMPTY ;
                }
            }

            // -----------------------------------------------------------------
            // 从A'中删除行和列单例
            // -----------------------------------------------------------------

            nz2 = 0 ;
            n2rows = 0 ;
            for (i = 0 ; i < m ; i++)
            {
                p = ATp [i] ;
                if (p >= 0)
                {
                    ATp [n2rows++] = nz2 ;
                    pend = UNFLIP (ATp [i+1]) ;
                    for (p = ATp [i] ; p < pend ; p++)
                    {
                        j = ATj [p] ;
                        ATj [nz2++] = W [j] ;
                    }
                }
            }
            ATp [n2rows] = nz2 ;
        }

        // ---------------------------------------------------------------------
        // 经过修剪的A'矩阵的 转置的填充减少排序
        // ---------------------------------------------------------------------

        AT->nrow = n2cols ;
        AT->ncol = n2rows ;
        // 保存当前的设置
        Long save [6] ;
        save [0] = cc->supernodal ;
        save [1] = cc->nmethods ;
        save [2] = cc->postorder ;
        save [3] = cc->method [0].ordering ;
        save [4] = cc->method [1].ordering ;
        save [5] = cc->method [2].ordering ;

        // 按照后续顺序排列列etree
        cc->postorder = TRUE ;

        /* --------------------------------- */
        // 选QR重排序方法：AMD\COLAMD\METIS  //
        /* --------------------------------- */
        if (ordering == QR_ORDERING_BEST)
        {
            ordering = QR_ORDERING_CHOL ;
            cc->nmethods = 4 ;
            cc->method [0].ordering = SPARSE_COLAMD ;
            cc->method [1].ordering = SPARSE_AMD ;
            cc->method [2].ordering = SPARSE_METIS ;
            cc->method [3].ordering = SPARSE_NESDIS ;
        }

        if (ordering == QR_ORDERING_DEFAULT)
        {
            
            // QR 默认 COLAMD
            ordering = QR_ORDERING_COLAMD ; 
            
            #ifdef all_methods_time
            // 测试所有方法的总时间
            ordering = QR_ORDERING_CHOL ;
            #endif

            #ifdef ONLY_METIS
            ordering = QR_ORDERING_ONLYMETIS;
            #endif
            #ifdef ONLY_AMD
            ordering = QR_ORDERING_AMD;
            #endif
            #ifdef ONLY_COLAMD
            ordering = QR_ORDERING_COLAMD;
            #endif
            #ifdef ONLY_NESDIS
            ordering = QR_ORDERING_NESDIS;
            #endif
        }
        double timeStart, timeEnd;
        struct timeval tv;
        /* --------------------------------- */
        // 根据ordering的值来调用重排序方法  //
        /* --------------------------------- */
        if (ordering == QR_ORDERING_AMD)
        {
            printf("QR use AMD \n");
            // 使用CHOL的接口到AMD 排序 A'*A
            SparseChol_amd (AT, NULL, 0, (Long *) (Q1fill + n1cols), cc) ;
        }
        
        else if (ordering == QR_ORDERING_COLAMD)
        {
            printf("QR use COLAMD \n");
            ordering = QR_ORDERING_COLAMD ; 
            SparseChol_colamd (AT, NULL, 0, TRUE,
                (Long *) (Q1fill + n1cols), cc) ; 
        }
        // #ifndef NMETIS
        else if (ordering == QR_ORDERING_ONLYMETIS)
        {
            printf("QR use METIS \n");
            
            SparseCore_metis (AT, NULL, 0, TRUE,
                (Long *) (Q1fill + n1cols), cc) ;
        }
        // #endif
        else if (ordering == QR_ORDERING_NESDIS)
        {
            printf("QR use NESDIS \n");
            Int *Cmember, *CParent;
            Cmember = (Int *)malloc( m * sizeof(Int)) ;
            CParent = (Int *)malloc( m * sizeof(Int)) ;
            SparseCore_nested_dissection (AT, NULL, 0, (Long *) (Q1fill + n1cols),
            CParent, Cmember, cc) ;
        }
        else if (ordering == QR_ORDERING_CHOL)
        {
            gettimeofday(&tv, NULL);
            timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
            //使用CHOL的内部排序(由cc定义)在
            // cc->supernodal = SPARSE_AUTO ;
            cc->postorder = TRUE ;
            sparse_factor *Sc ;
            Sc = SparseChol_analyze_p2 (1, AT, NULL, NULL, 0, cc, result_name) ;
            if (Sc != NULL)
            {
                // 复制perm(Sc->Perm [0:n2cols-1]) 到 Q1fill (n1cols:n)
                Long *Sc_perm = (Long *) Sc->Perm ;
                for (k = 0 ; k < n2cols ; k++)
                {
                    Q1fill [k + n1cols] = Sc_perm [k] ;
                }
                // HNUCHOL选择了一个排序;确定使用的排序
                switch (Sc->ordering)
                {
                    case SPARSE_AMD:    ordering = QR_ORDERING_AMD    ;break;//5
                    case SPARSE_COLAMD: ordering = QR_ORDERING_COLAMD ;break;//2
                    case SPARSE_METIS: ordering = QR_ORDERING_METIS ;break;//10
                    default: ordering = QR_ORDERING_NESDIS; break;//6
                }
            }
            gettimeofday(&tv, NULL);
            timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
            switch (ordering)
            {
            case QR_ORDERING_COLAMD:
                printf("Selected method is COLAMD\n");
                break;
            case QR_ORDERING_AMD:
                printf("Selected method is AMD\n");
                break;
            case QR_ORDERING_NESDIS:
                printf("Selected method is NESDIS\n");
                break;
            case QR_ORDERING_METIS:
                printf("Selected method is METIS\n");
                break;
            default:
                break;
            }
            

            printf("Brute-force method time = %lf\n",timeEnd - timeStart);
            FILE* fresult;
            fresult = fopen("./Results/Brute_force_time.txt","a+");
            fprintf (fresult,"%lf\n",timeEnd - timeStart);
            fclose(fresult);
            cc->SPQR_istat [7] = ordering ;
            SparseCore_free_factor (&Sc, cc) ;
        }

        cc->SPQR_istat [7] = ordering ; 
        // 保存值
        cc->supernodal              = save [0] ;
        cc->nmethods                = save [1] ;
        cc->postorder               = save [2] ;
        cc->method [0].ordering     = save [3] ;
        cc->method [1].ordering     = save [4] ;
        cc->method [2].ordering     = save [5] ;

        AT->nrow = n ;
        AT->ncol = m ;
    }

    // -------------------------------------------------------------------------
    // 释放 AT
    // -------------------------------------------------------------------------

    SparseCore_free_sparse (&AT, cc) ;   // ]

    // -------------------------------------------------------------------------
    //  检查排序方法是否成功
    // -------------------------------------------------------------------------

    if (cc->status < SPARSE_OK)
    {
        // 内存越界
        SparseCore_free (worksize, sizeof (Long), Work, cc) ;
        SparseCore_free (n, sizeof (Long), Q1fill, cc) ;
        SparseCore_free (n1rows+1, sizeof (Long), R1p, cc) ;
        SparseCore_free (m,        sizeof (Long), P1inv, cc) ;
        return (FALSE) ;
    }

    // -------------------------------------------------------------------------
    // 将填充减少顺序映射回A
    // -------------------------------------------------------------------------

    if (n1cols > 0 && fill_reducing_ordering)
    {
        for (j = 0 ; j < n ; j++)
        {
            // j为A的列。col2 = W [j]要么为空，要么为剪枝矩阵对应的列
            col2 = W [j] ;
            if (col2 != EMPTY)
            {
                Winv [col2] = j ;
            }
        }

        for (k = n1cols ; k < n ; k++)
        {
            //col2是修剪后的矩阵中的一列
            col2 = Q1fill [k] ;
            // j是A矩阵的相关列
            j = Winv [col2] ;
            Q1fill [k] = j ;
        }
    }

    // -------------------------------------------------------------------------
    // 查找Y = [A2 B2]的列指针;A2的列
    // -------------------------------------------------------------------------

    if (n1cols == 0 )
    {
        // A将被因式分解
        Y = NULL ;
    }
    else
    {
        // Y还没有分量;这里只把p置零和nzmax重新分配。 nnz(Y)稍后确定
        Y = SparseCore_allocate_sparse (m-n1rows, n-n1cols, 0,
            FALSE, TRUE, 0, xtype, cc) ;

        if (cc->status < SPARSE_OK)
        {
            // 内存越界
            SparseCore_free (worksize, sizeof (Long), Work, cc) ;
            SparseCore_free (n, sizeof (Long), Q1fill, cc) ;
            SparseCore_free (n1rows+1, sizeof (Long), R1p, cc) ;
            SparseCore_free (m,        sizeof (Long), P1inv, cc) ;
            return (FALSE) ;
        }

        Yp = (Long *) Y->p ; 

        ynz = 0 ;
        for (k = n1cols ; k < n ; k++)
        {
            j = Q1fill [k] ;
            d = Degree [j] ;
            Yp [k-n1cols] = ynz ;
            ynz += d ;
        }
        Yp [n-n1cols] = ynz ;
    }

    // -------------------------------------------------------------------------
    // 释放工作区并返回结果
    // -------------------------------------------------------------------------

    SparseCore_free (worksize, sizeof (Long), Work, cc) ;

    *p_Q1fill = Q1fill ;
    *p_R1p    = R1p ;
    *p_P1inv  = P1inv ;
    *p_Y      = Y ;
    *p_n1cols = n1cols ;
    *p_n1rows = n1rows ;
    return (TRUE) ;
} // end of qr_1colamd

// =============================================================================
// === qr_tol ================================================================
// =============================================================================
// 返回默认的列2-范数公差；返回默认的tol(如果错误，返回-1)
double qr_tol
(
    sparse_csc *A,
    sparse_common *cc
)
{
    double tol = (20 * ((double) A->nrow + (double) A->ncol) * DBL_EPSILON *
                  qr_maxcolnorm (A, cc));
    tol = MIN(tol, DBL_MAX);
    return (tol) ;
} // end of qr_tol

// =============================================================================
// === qr_shift ==============================================================
// =============================================================================
// 插入一个0作为向量的第一个元素，将所有其他元素向下移动一个位置。
// 如果X为空，则不执行任何操作。
void qr_shift(Long n,Long *X)
{
    Long k ;
    if (X != NULL)
    {
        for (k = n ; k >= 1 ; k--)
        {
            X [k] = X [k-1] ;
        }
        X [0] = 0 ;
    }
} // end of qr_shift

// =============================================================================
// === qr_freesym ============================================================
// =============================================================================
// 释放 QR 符号结构体
void qr_freesym( qr_symbolic **QRsym_handle, sparse_common *cc)
{
    qr_symbolic *QRsym ;
    Long m, n, anz, nf, rjsize, ns, ntasks ;

    if (QRsym_handle == NULL || *QRsym_handle == NULL)
    {
        // 无事可做;调用可能耗尽了内存
        return ;
    }
    QRsym = *QRsym_handle  ;

    m = QRsym->m ;
    n = QRsym->n ;
    nf = QRsym->nf ;
    anz = QRsym->anz ;
    rjsize = QRsym->rjsize ;

    SparseCore_free (n,      sizeof (Long), QRsym->Qfill, cc) ;
    SparseCore_free (nf+1,   sizeof (Long), QRsym->Super, cc) ;
    SparseCore_free (nf+1,   sizeof (Long), QRsym->Rp, cc) ;
    SparseCore_free (rjsize, sizeof (Long), QRsym->Rj, cc) ;
    SparseCore_free (nf+1,   sizeof (Long), QRsym->Parent, cc) ;
    SparseCore_free (nf+2,   sizeof (Long), QRsym->Childp, cc) ;
    SparseCore_free (nf+1,   sizeof (Long), QRsym->Child, cc) ;
    SparseCore_free (nf+1,   sizeof (Long), QRsym->Post, cc) ;
    SparseCore_free (m,      sizeof (Long), QRsym->PLinv, cc) ;
    SparseCore_free (n+2,    sizeof (Long), QRsym->Sleft, cc) ;
    SparseCore_free (m+1,    sizeof (Long), QRsym->Sp, cc) ;
    SparseCore_free (anz,    sizeof (Long), QRsym->Sj, cc) ;

    SparseCore_free (nf+1,   sizeof (Long), QRsym->Hip, cc) ;

    SparseCore_free (nf+1,   sizeof (Long), QRsym->Fm, cc) ;
    SparseCore_free (nf+1,   sizeof (Long), QRsym->Cm, cc) ;

    // 并行分析
    ntasks = QRsym->ntasks ;
    SparseCore_free (ntasks+2, sizeof (Long), QRsym->TaskChildp, cc) ;
    SparseCore_free (ntasks+1, sizeof (Long), QRsym->TaskChild, cc) ;
    SparseCore_free (nf+1,     sizeof (Long), QRsym->TaskFront, cc) ;
    SparseCore_free (ntasks+2, sizeof (Long), QRsym->TaskFrontp, cc) ;
    SparseCore_free (ntasks+1, sizeof (Long), QRsym->TaskStack, cc) ;
    SparseCore_free (nf+1,     sizeof (Long), QRsym->On_stack, cc) ;

    ns = QRsym->ns ;
    SparseCore_free (ns+2,     sizeof (Long), QRsym->Stack_maxstack, cc) ;

    SparseCore_free (1, sizeof (qr_symbolic), QRsym, cc) ;

    *QRsym_handle = NULL ;
} // end of qr_freesym

// =============================================================================
// === qr_freenum ============================================================
// =============================================================================
// 释放QR数值分解处理的结构体
void qr_freenum( 
            qr_numeric **QRnum_handle,sparse_common *cc )
{
    qr_numeric *QRnum ;
    Long nf, n, m, rjsize, hisize, ns, stack, maxstack ;

    if (QRnum_handle == NULL || *QRnum_handle == NULL)
    {
        return ;
    }
    QRnum = *QRnum_handle ;

    n  = QRnum->n ;
    m  = QRnum->m ;
    nf = QRnum->nf ;
    rjsize = QRnum->rjsize ;
    hisize = QRnum->hisize ;
    ns = QRnum->ns ;
    maxstack = QRnum->maxstack ;

    SparseCore_free (nf, sizeof (double *), QRnum->Rblock, cc) ;
    SparseCore_free (n,  sizeof (char),    QRnum->Rdead,  cc) ;


    // QRnum->H* 只有当H被保留时，项目才会存在
    SparseCore_free (rjsize, sizeof (Long),  QRnum->HStair,  cc) ;
    SparseCore_free (rjsize, sizeof (double), QRnum->HTau,    cc) ;
    SparseCore_free (nf,     sizeof (Long),  QRnum->Hm,      cc) ;
    SparseCore_free (nf,     sizeof (Long),  QRnum->Hr,      cc) ;
    SparseCore_free (hisize, sizeof (Long),  QRnum->Hii,     cc) ;
    SparseCore_free (m,      sizeof (Long),  QRnum->HPinv,   cc) ;

    // 释放每一个栈
    if (QRnum->Stacks != NULL)
    {
        Long *Stack_size = QRnum->Stack_size ;
        for (stack = 0 ; stack < ns ; stack++)
        {
            size_t s = Stack_size ? (Stack_size [stack]) : maxstack ;
            SparseCore_free (s, sizeof (double), QRnum->Stacks [stack], cc) ;
        }
    }
    SparseCore_free (ns, sizeof (double *), QRnum->Stacks, cc) ;
    SparseCore_free (ns, sizeof (Long), QRnum->Stack_size, cc) ;

    SparseCore_free (1, sizeof (qr_numeric), QRnum, cc) ;
    *QRnum_handle = NULL ;
} //end of qr_freenum

// =============================================================================
// === qr_cumsum =============================================================
// =============================================================================
Long qr_cumsum (Long n,Long *X)
{
    Long itot, t, x, k ;

    // -------------------------------------------------------------------------
    // X = cumsum ([0 X])
    // -------------------------------------------------------------------------

    itot = 0 ;
    if (X != NULL)
    {
        for (k = 0 ; k < n ; k++)
        {
            t = itot ;    
            x = X [k] ;
            itot += x ;
            X [k] = t ;
        }
        X [n] = itot ;
    }

    return (itot) ;
} //end of qr_cumsum

// =============================================================================
// === qr_rmap ===============================================================
// =============================================================================
int qr_rmap
(
    SparseQR_factorization *QR,
    sparse_common *cc
)
{
    Long n, j, i, p, n1rows, n1cols ;
    Long *Rmap, *RmapInv, *R1p, *R1j ;

    n = QR->nacols ;
    Rmap = QR->Rmap ;
    RmapInv = QR->RmapInv ;

    if (Rmap == NULL)
    {
        QR->Rmap    = Rmap    = (Long *) SparseCore_malloc (n, sizeof(Long), cc);
        QR->RmapInv = RmapInv = (Long *) SparseCore_malloc (n, sizeof(Long), cc);
        if (cc->status < SPARSE_OK)
        {
            // 内存越界
            return (FALSE) ;
        }
    }

    for (j = 0 ; j < n ; j++)
    {
        Rmap [j] = EMPTY ;
    }

    R1p = QR->R1p ;
    R1j = QR->R1j ;
    n1rows = QR->n1rows ;
    n1cols = QR->n1cols ;

    // 查找单例行的映射
    for (i = 0 ; i < n1rows ; i++)
    {
        // R的第i行是单例行;找到它对应的关键列。
        p = R1p [i] ;
        j = R1j [p] ;
        Rmap [j] = i ;
    }

    // 找到多波前R的关键行的映射
    char *Rdead = QR->QRnum->Rdead ;
    for (j = n1cols ; j < n ; j++)
    {
        if (!Rdead [j-n1cols])
        {
            Rmap [j] = i++ ;
        }
    }

    // 用R的死列完成映射，包括R的单元素部分和R的多波前部分
    for (j = 0 ; j < n ; j++)
    {
        if (Rmap [j] == EMPTY)
        {
            Rmap [j] = i++ ;
        }
    }

    //构造Rmap的逆
    for (j = 0 ; j < n ; j++)
    {
        i = Rmap [j] ;
        RmapInv [i] = j ;
    }
    return (TRUE) ;
}

// =============================================================================
// === qr_maxcolnorm =========================================================
// =============================================================================
inline double qr_private_nrm2 (Long n, double *X, sparse_common *cc)
{
    double norm = 0 ;
    BLAS_INT N = n, one = 1 ;
    if (CHECK_BLAS_INT && !EQ (N,n))
    {
        cc->blas_ok = FALSE ;
    }
    if (!CHECK_BLAS_INT || cc->blas_ok)
    {
        #ifndef USEATL
            norm = BLAS_DNRM2 (&N, X, &one) ;
        #else
            norm = hnu_nrm2(N, X, one);
        #endif
    }
    return (norm) ;
}


// =============================================================================
// === qr_maxcolnorm =========================================================
// =============================================================================

double qr_maxcolnorm
(
    sparse_csc *A,
    sparse_common *cc
)
{
    double norm, maxnorm ;
    Long j, p, len, n, *Ap ;
    double *Ax ;

    cc->blas_ok = TRUE ;
    n = A->ncol ;
    Ap = (Long *) A->p ;
    Ax = (double *) A->x ;

    maxnorm = 0 ;
    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;
        len = Ap [j+1] - p ;
        norm = qr_private_nrm2 (len, Ax + p, cc) ;
        maxnorm = MAX (maxnorm, norm) ;
    }

    if (CHECK_BLAS_INT && !cc->blas_ok)
    {
        return (EMPTY) ;
    }

    return (maxnorm) ;
}

/**************************************
 *            求解器需要  
 **************************************/
#define FREE_WORK \
{ \
    SparseCore_free_dense (&Zdense, cc) ; \
    SparseCore_free_dense (&Vdense, cc) ; \
    SparseCore_free_dense (&Wdense, cc) ; \
    SparseCore_free_dense (&Cdense, cc) ; \
    SparseCore_free (maxfn, sizeof (double), H_Tau,   cc) ; \
    SparseCore_free (maxfn, sizeof (Long),  H_start, cc) ; \
    SparseCore_free (maxfn, sizeof (Long),  H_end,   cc) ; \
}
#define HCHUNK 32

// =============================================================================
// === qr_private_get_H_vectors ==============================================
// =============================================================================

// 获取指向单个波前阵中的Householder向量的指针。
// 返回 F中的 Householder向量 #.
Long qr_private_get_H_vectors
(
    // inputs
    Long f,                 
    SparseQR_factorization *QR,

    // outputs
    double *H_Tau,           
    Long *H_start,          
    Long *H_end,            
    sparse_common *cc
)
{
    qr_symbolic *QRsym ;
    qr_numeric  *QRnum ;
    double *Tau ;
    Long *Rj, *Stair ;
    Long col1, col2, fp, pr, fn, fm, h, nh, p, rm, k, j, t, n1cols, n ;

    // -------------------------------------------------------------------------
    // 得到波前阵F的R块
    // -------------------------------------------------------------------------

    QRsym = QR->QRsym ;
    QRnum = QR->QRnum ;
    n1cols = QR->n1cols ;
    Rj = QRsym->Rj ;
    n = QR->nacols ;

    col1 = QRsym->Super [f] ;          
    col2 = QRsym->Super [f+1] ;         
    fp = col2 - col1 ;                  
    pr = QRsym->Rp [f] ;               
    fn = QRsym->Rp [f+1] - pr ;         

    Stair = QRnum->HStair + pr ;        
    Tau = QRnum->HTau + pr ;           

    fm = QRnum->Hm [f] ;                
    h = 0 ;                             
    nh = 0 ;                            
    p = 0 ;                            

    // -------------------------------------------------------------------------
    // 求R+H块中所有的Household向量
    // -------------------------------------------------------------------------

    rm = 0 ;                            // R块中的行数
    for (k = 0 ; k < fn && nh < fm ; k++)
    {
        if (k < fp)
        {
            // 波前阵F的关键列
            j = col1 + k ;
            
            t = Stair [k] ;                 // R+H 向量的长度
            
            if (t == 0)
            {
                p += rm ;                   // 没有H向量，跳过R
                continue ;
            }
            else if (rm < fm)
            {
                rm++ ;                      // k列没有死
            }
            h = rm ;                        //H向量从H行开始
            
        }
        else
        {
            // F前面的非关键的一列
            j = Rj [pr + k] ;
            
            t = Stair [k] ;                 // R+H 向量的程度
            
            h = MIN (h+1,fm) ;              
            
        }
        if (j+n1cols >= n) break ;          // 如果[A Binput]被分解
        p += rm ;                          
        H_Tau [nh] = Tau [k] ;              // 又找到了一个H向量
        H_start [nh] = p ;                  // H向量从这里开始
        p += (t-h) ;                        // 跳过H向量
        H_end [nh] = p ;                    // H向量在这里结束
        nh++ ;
        if (h == fm) break ;                // 这是最后一个H向量
    }

    return (nh) ;
}

// =============================================================================
// === qr_private_load_H_vectors =============================================
// =============================================================================

// 将Householder向量h1:h2-1加载到panel V中，返回V中的行元素#。

Long qr_private_load_H_vectors
(
    // input
    Long h1,            
    Long h2,
    Long *H_start,      
    Long *H_end,        
    double *R,           
    // output
    double *V,           
    sparse_common *cc
)
{
    // v 为最后一个H向量的长度
    Long v = H_end [h2-1] - H_start [h2-1] + (h2-h1) ;
    double *V1 = V ;
    for (Long h = h1 ; h < h2 ; h++)
    {
        Long i ;
        // V的此部分不被访问，仅用于测试:
        i = h-h1 ;
        V1 [i++] = 1 ;
        for (Long p = H_start [h] ; p < H_end [h] ; p++)
        {
            V1 [i++] = R [p] ;
        }
        for ( ; i < v ; i++)
        {
            V1 [i] = 0 ;
        }
        V1 += v ;
    }
    return (v) ;
}

// ===========================
// === qr_panel ===========
// ===========================
void qr_panel
(
    // input
    int method,         
    Long m,
    Long n,
    Long v,             
    Long h,             
    Long *Vi,           
    double *V,           
    double *Tau,         
    Long ldx,

    // input/output
    double *X,          

    // workspace
    double *C,           
    double *W,           

    sparse_common *cc
)
{
    double *C1, *X1 ;
    Long k, p, i ;

    // -------------------------------------------------------------------------
    // 将X收集到工作空间C中
    // -------------------------------------------------------------------------

    if (method == QR_QTX || method == QR_QX)
    {

        C1 = C ;
        X1 = X ;
        for (k = 0 ; k < n ; k++)
        {
            for (p = 0 ; p < v ; p++)
            {
                i = Vi [p] ;
                C1 [p] = X1 [i] ;
            }
            C1 += v ;
            X1 += ldx ;
        }
    }
    else 
    {

        C1 = C ;
        for (p = 0 ; p < v ; p++)
        {
            i = Vi [p] ;
            X1 = X + i*ldx ;
            for (k = 0 ; k < m ; k++)
            {
                C1 [k] = X1 [k] ;
            }
            C1 += m ;
        }
    }

    // -------------------------------------------------------------------------
    // 将household panel面板应用到C
    // -------------------------------------------------------------------------

    if (method == QR_QTX || method == QR_QX)
    {
        qr_larftb (method, v, n, h, v, v, V, Tau, C, W, cc) ;
    }
    else 
    {
        qr_larftb (method, m, v, h, m, v, V, Tau, C, W, cc) ;
    }

    // -------------------------------------------------------------------------
    // 把C 映射回X
    // -------------------------------------------------------------------------

    if (method == QR_QTX || method == QR_QX)
    {
        C1 = C ;
        X1 = X ;
        for (k = 0 ; k < n ; k++)
        {
            for (p = 0 ; p < v ; p++)
            {
                i = Vi [p] ;
                X1 [i] = C1 [p] ;
            }
            C1 += v ;
            X1 += ldx ;
        }
    }
    else
    {
        C1 = C ;
        for (p = 0 ; p < v ; p++)
        {
            i = Vi [p] ;
            X1 = X + i*ldx ;
            for (k = 0 ; k < m ; k++)
            {
                X1 [k] = C1 [k] ;
            }
            C1 += m ;
        }
    }
}

// ============================
// === qr_private_Happly ====
// ============================
// 输入SparseQR得到的QR, 将Housholder向量应用到稠密矩阵X上

void qr_private_Happly
(
    // inputs
    int method,            
    SparseQR_factorization *QR,
    Long hchunk,            

    // input/output
    Long m,
    Long n,
    double *X,               

    //工作区，没有在输入或输出上定义
    double *H_Tau,           
    Long *H_start,          
    Long *H_end,            
    double *V,               
    double *C,               
    double *W,               
    sparse_common *cc
)
{

    qr_symbolic *QRsym ;
    qr_numeric  *QRnum ;
    double **Rblock, *R, *X2 ;
    Long *Hii, *Hip, *Hi ;
    Long nf, f, nh, h1, h2, v, n1rows, m2, n2 ;

    // -------------------------------------------------------------------------
    // 获取QR分解的内容
    // -------------------------------------------------------------------------

    QRsym = QR->QRsym ;
    QRnum = QR->QRnum ;
    
    nf = QRsym->nf ;
    Rblock = QRnum->Rblock ;
    Hii = QRnum->Hii ;
    Hip = QRsym->Hip ;
    
    n1rows = QR->n1rows ;

    // -------------------------------------------------------------------------
    //  X (n1rows:m-1,:) or X (:,n1rows:n-1) 上的操作
    // -------------------------------------------------------------------------

    // Household向量必须向下移动n1row，超过单例行，这些单例行不属于QR分解的多波前部分。

    if (method == QR_QTX || method == QR_QX)
    {
        // Q*X or Q'*X
        X2 = X + n1rows ;
        n2 = n ;
        m2 = m - n1rows ;
    }
    else
    {
        X2 = X + n1rows * m ;
        n2 = n - n1rows ; 
        m2 = m ;
    }

    // -------------------------------------------------------------------------
    // 应用Household向量
    // -------------------------------------------------------------------------

    if (method == QR_QTX || method == QR_XQ)
    {

        // ---------------------------------------------------------------------
        // 正向应用Household变量
        // ---------------------------------------------------------------------

        for (f = 0 ; f < nf ; f++)
        {
            // 得到波前F的Household向量
            nh = qr_private_get_H_vectors (f, QR, H_Tau, H_start, H_end, cc) ;
            R = Rblock [f] ;
            Hi = &Hii [Hip [f]] ;               // H 的行下标列表

            // 应用Household向量，一次一个panel
            for (h1 = 0 ; h1 < nh ; h1 = h2)
            {
                // 从R到面板V的加载向量h1:h2-1，并应用他们
                h2 = MIN (h1 + hchunk, nh) ;
                v = qr_private_load_H_vectors (h1, h2, H_start, H_end, R, V,
                    cc) ;
                qr_panel (method, m2, n2, v, h2-h1, Hi+h1, V, H_Tau+h1,
                    m, X2, C, W, cc) ;
            }
        }
    }
    else
    {

        // ---------------------------------------------------------------------
        // 反向应用
        // ---------------------------------------------------------------------

        for (f = nf-1 ; f >= 0 ; f--)
        {

            // 得到波前F的Household向量
            nh = qr_private_get_H_vectors (f, QR, H_Tau, H_start, H_end, cc) ;
            R = Rblock [f] ;
            Hi = &Hii [Hip [f]] ;               // H的行下标列表

            // 应用户主向量，一次一个panel
            for (h2 = nh ; h2 > 0 ; h2 = h1)
            {
                // 从R到panel V的加载向量h1:h2-1，并应用
                h1 = MAX (h2 - hchunk, 0) ;
                v = qr_private_load_H_vectors (h1, h2, H_start, H_end, R, V,
                    cc) ;
                
                qr_panel (method, m2, n2, v, h2-h1, Hi+h1, V, H_Tau+h1,
                    m, X2, C, W, cc) ;
            }
        }
    }
}

/**************************************
 *           QR_qmult
 **************************************/
// 对一个稠密矩阵X 应用Householder变换
//  method QR_QTX (0): Y = Q'*X
//  method QR_QX  (1): Y = Q*X
//  method QR_XQT (2): Y = X*Q'
//  method QR_XQ  (3): Y = X*Q

dense_array *QR_qmult
(
    int method,             
    SparseQR_factorization *QR,
    dense_array *Xdense,  // 大小m×n，前导维ldx
    sparse_common *cc
)
{
    dense_array *Ydense, *Cdense, *Vdense, *Wdense, *Zdense ;
    double *X, *Y, *X1, *Y1, *Z1, *C, *V, *Z, *W, *H_Tau ;
    Long *HPinv, *H_start, *H_end ;
    Long i, k, mh, v, hchunk, ldx, m, n, maxfn, ok ;

    // -------------------------------------------------------------------------
    // 获取输入
    // -------------------------------------------------------------------------

    Long xtype = 1 ;

    cc->status = SPARSE_OK ;

    HPinv = (QR->n1cols > 0) ? QR->HP1inv : QR->QRnum->HPinv ;

    v = QR->QRnum->maxfm ;
    mh = QR->narows ;
    maxfn = QR->QRsym->maxfn ;

    X = (double *) Xdense->x ;
    m = Xdense->nrow ;
    n = Xdense->ncol ;
    ldx = Xdense->d ;

    if (method == QR_QTX || method == QR_QX)
    {
        // H和X的行必须是相同的
        if (mh != m)
        {
            return (NULL) ;
        }
    }
    else if (method == QR_XQT || method == QR_XQ)
    {
        // H的行和X的列必须是相同的
        if (mh != n)
        {
            return (NULL) ;
        }
    }
    else
    {
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // 分配结果 Y
    // -------------------------------------------------------------------------

    Ydense = SparseCore_allocate_dense (m, n, m, xtype, cc) ;
    if (cc->status < SPARSE_OK)
    {
        // 内存越界
        return (NULL) ;
    }
    Y = (double *) Ydense->x ;

    if (m == 0 || n == 0)
    {
        return (Ydense) ;    
    }

    // -------------------------------------------------------------------------
    // 分配空间
    // -------------------------------------------------------------------------

    Z = NULL ;
    Zdense = NULL ;
    ok = TRUE ;
    if (method == QR_QX || method == QR_XQT)
    {
        
        Zdense = SparseCore_allocate_dense (m, n, m, xtype, cc) ;
        ok = (Zdense != NULL) ;
    }

    hchunk = HCHUNK ;
    
    Cdense = SparseCore_allocate_dense (v, (method <= QR_QX) ? n : m,
        v, xtype, cc) ;
    Vdense = NULL ;
    Wdense = NULL ;

    H_Tau   = (double *) SparseCore_malloc (maxfn, sizeof (double), cc) ;
    H_start = (Long *)  SparseCore_malloc (maxfn, sizeof (Long),  cc) ;
    H_end   = (Long *)  SparseCore_malloc (maxfn, sizeof (Long),  cc) ;

    if (!ok || Cdense == NULL || cc->status < SPARSE_OK)
    {

        SparseCore_free_dense (&Ydense, cc) ;
        FREE_WORK ;
        return (NULL) ;
    }

    // -------------------------------------------------------------------------
    // X 复制到 Z
    // -------------------------------------------------------------------------

    if (method == QR_QX || method == QR_XQT)
    {
        Z = (double *) Zdense->x ;
        Z1 = Z ;
        X1 = X ;
        for (k = 0 ; k < n ; k++)
        {
            for (i = 0 ; i < m ; i++)
            {
                Z1 [i] = X1 [i] ;
            }
            X1 += ldx ;
            Z1 += m ;
        }
    }

    // -------------------------------------------------------------------------
    // 分配O (hchunk) 大小的工作区
    // -------------------------------------------------------------------------

    Vdense = SparseCore_allocate_dense (v, hchunk, v, xtype, cc) ;

    Wdense = SparseCore_allocate_dense (hchunk,
        hchunk + ((method <= QR_QX) ? n : m), hchunk, xtype, cc) ;

    // -------------------------------------------------------------------------
    // 如果内存不足，弃掉
    // -------------------------------------------------------------------------

    if (Vdense == NULL || Wdense == NULL)
    {
        // PUNT: 内存越界; 使用hchunk = 1再试一次
        cc->status = SPARSE_OK ;
        hchunk = 1 ;

        SparseCore_free_dense (&Vdense, cc) ;
        SparseCore_free_dense (&Wdense, cc) ;

        Vdense = SparseCore_allocate_dense (v, hchunk, v, xtype, cc) ;
        Wdense = SparseCore_allocate_dense (hchunk,
            hchunk + ((method <= QR_QX) ? n : m), hchunk, xtype, cc) ;

        if (Vdense == NULL || Wdense == NULL)
        {

            SparseCore_free_dense (&Ydense, cc) ;
            FREE_WORK ;
            return (NULL) ;
        }
    }

    V = (double *) Vdense->x ;
    C = (double *) Cdense->x ;
    W = (double *) Wdense->x ;

    // -------------------------------------------------------------------------
    // Y = Q'*X, Q*X, X*Q, or X*Q'
    // -------------------------------------------------------------------------

    if (method == QR_QTX)
    {

        // ---------------------------------------------------------------------
        // Y = Q'*X
        // ---------------------------------------------------------------------

        X1 = X ;
        Y1 = Y ;
        for (k = 0 ; k < n ; k++)
        {
            for (i = 0 ; i < m ; i++)
            {
                Y1 [HPinv [i]] = X1 [i] ;
            }
            X1 += ldx ;
            Y1 += m ;
        }

        // 将H作用于Y
        qr_private_Happly (method, QR, hchunk, m, n, Y, H_Tau, H_start, H_end,
            V, C, W, cc) ;

    }
    else if (method == QR_QX)
    {

        // ---------------------------------------------------------------------
        // Y = Q*X
        // ---------------------------------------------------------------------

        // 将H作用于Z
        qr_private_Happly (method, QR, hchunk, m, n, Z, H_Tau, H_start, H_end,
            V, C, W, cc) ;

        Z1 = Z ;
        Y1 = Y ;
        for (k = 0 ; k < n ; k++)
        {
            for (i = 0 ; i < m ; i++)
            {
                Y1 [i] = Z1 [HPinv [i]] ;
            }
            Z1 += m ;
            Y1 += m ;
        }

    }
    else if (method == QR_XQT)
    {

        // ---------------------------------------------------------------------
        // Y = X*Q'
        // ---------------------------------------------------------------------

        // 将H作用于Y
        qr_private_Happly (method, QR, hchunk, m, n, Z, H_Tau, H_start, H_end,
            V, C, W, cc) ;

        Y1 = Y ;
        for (k = 0 ; k < n ; k++)
        {
            Z1 = Z + HPinv [k] * m ;   
            for (i = 0 ; i < m ; i++)
            {
                Y1 [i] = Z1 [i] ;
            }
            Y1 += m ;
        }

    }
    else if (method == QR_XQ)
    {

        // ---------------------------------------------------------------------
        // Y = X*Q
        // ---------------------------------------------------------------------

        X1 = X ;
        for (k = 0 ; k < n ; k++)
        {
            Y1 = Y + HPinv [k] * m ;    
            for (i = 0 ; i < m ; i++)
            {
                Y1 [i] = X1 [i] ;
            }
            X1 += ldx ;
        }

        // apply H to Y
        qr_private_Happly (method, QR, hchunk, m, n, Y, H_Tau, H_start, H_end,
            V, C, W, cc) ;

    }

    // -------------------------------------------------------------------------
    // 释放空间 返回Y
    // -------------------------------------------------------------------------

    FREE_WORK ;

    if (CHECK_BLAS_INT && !cc->blas_ok)
    {

        SparseCore_free_dense (&Ydense, cc) ;
        return (NULL) ;
    }

    return (Ydense) ;
}
// ==================================================
// ================== QR_solve ======================
// ==================================================
// 调用 SparseQR_factorize 返回的QR对象的R因子来求解上下三角系统。
dense_array *QR_solve    
(
    // inputs, not modified:
    int system,                 
    SparseQR_factorization  *QR,    
    dense_array *B,           
    // workspace and parameters
    sparse_common *cc
)
{
    double *Bx ;
    dense_array *W, *X ;
    Long m, n, nrhs, ldb, ok ;

    // -------------------------------------------------------------------------
    // 获取输入参数
    // -------------------------------------------------------------------------

    Long xtype = 1 ;

    if (system < QR_RX_EQUALS_B || system > QR_RTX_EQUALS_ETB)
    {
        
        return (NULL) ;
    }
    m = QR->narows ;
    n = QR->nacols ;
    if ((Long) B->nrow != ((system <= QR_RETX_EQUALS_B) ? m : n))
    {
       
        return (NULL) ;
    }

    cc->status = SPARSE_OK ;

    nrhs = B->ncol ;
    Bx = (double *) B->x ;
    ldb = B->d ;

    if (system == QR_RX_EQUALS_B || system == QR_RETX_EQUALS_B)
    {

        // ---------------------------------------------------------------------
        // X = E*(R\B) or X=R\B
        // ---------------------------------------------------------------------

        Long *Rlive ;
        double **Rcolp ;
        X = SparseCore_allocate_dense (n, nrhs, n, xtype, cc) ;
        Long maxfrank = QR->QRnum->maxfrank  ;
        W = SparseCore_allocate_dense (maxfrank, nrhs, maxfrank, xtype, cc) ;
        Rlive = (Long *)   SparseCore_malloc (maxfrank, sizeof (Long),    cc) ;
        Rcolp = (double **) SparseCore_malloc (maxfrank, sizeof (double *), cc) ;
        ok = (X != NULL) && (W != NULL) && (cc->status == SPARSE_OK) ;
        if (ok)
        {
            qr_rsolve (QR, system == QR_RETX_EQUALS_B, nrhs, ldb, Bx,
                (double *) X->x, Rcolp, Rlive, (double *) W->x, cc) ;
        }
        SparseCore_free (maxfrank, sizeof (Long),    Rlive, cc) ;
        SparseCore_free (maxfrank, sizeof (double *), Rcolp, cc) ;
        SparseCore_free_dense (&W, cc) ;

    }
    else
    {

        // ---------------------------------------------------------------------
        // X = R'\(E'*B) or R'\B
        // ---------------------------------------------------------------------

        X = SparseCore_allocate_dense (m, nrhs, m, xtype, cc) ;
        ok = (X != NULL) ;
        if (ok)
        {
            qr_private_rtsolve (QR, system == QR_RTX_EQUALS_ETB, nrhs, ldb,
                Bx, (double *) X->x, cc) ;
        }
    }

    if (!ok)
    {
        
        SparseCore_free_dense (&X, cc) ;
        return (NULL) ;
    }

    return (X) ;
}

inline double qr_divide (double a, double b, sparse_common *cc)  // cc unused
{
    return (a/b) ;
}

// =============================================================================
// === qr_rsolve =============================================================
// =============================================================================

// 使用 SparseQR 得到的 QR 对象，求解 X = E*(R\B) 或 X=R\B  
void qr_rsolve
(
    // inputs
    SparseQR_factorization *QR,
    int use_Q1fill,         

    Long nrhs,              
    Long ldb,               
    double *B,               

    // output
    double *X,               

    // workspace
    double **Rcolp,         
    Long *Rlive,           
    double *W,              

    sparse_common *cc
)
{
    qr_symbolic *QRsym ;
    qr_numeric  *QRnum ;
    Long n1rows, n1cols, n ;
    Long *Q1fill, *R1p, *R1j ;
    double *R1x ;

    double xi ;
    double **Rblock, *R, *W1, *B1, *X1 ;
    Long *Rp, *Rj, *Super, *HStair, *Hm, *Stair ;
    char *Rdead ;
    Long nf, // m,
        rank, j, f, col1, col2, fp, pr, fn, rm, k, i, row1, row2, ii,
        keepH, fm, h, t, live, kk ;

    // -------------------------------------------------------------------------
    // 获取 QR 对象里面的值
    // -------------------------------------------------------------------------

    QRsym = QR->QRsym ;
    QRnum = QR->QRnum ;
    n1rows = QR->n1rows ;
    n1cols = QR->n1cols ;
    n = QR->nacols ;
    Q1fill = use_Q1fill ? QR->Q1fill : NULL ;
    R1p = QR->R1p ;
    R1j = QR->R1j ;
    R1x = QR->R1x ;

    keepH = QRnum->keepH ;

    nf = QRsym->nf ;

    Rblock = QRnum->Rblock ;
    Rp = QRsym->Rp ;
    Rj = QRsym->Rj ;
    Super = QRsym->Super ;
    Rdead = QRnum->Rdead ;
    rank = QR->rank ;   
                        
    HStair = QRnum->HStair ;
    Hm = QRnum->Hm ;

    // -------------------------------------------------------------------------
    // X = 0
    // -------------------------------------------------------------------------

    X1 = X ;
    for (kk = 0 ; kk < nrhs ; kk++)
    {
        for (i = 0 ; i < n ; i++)
        {
            X1 [i] = 0 ;
        }
        X1 += n ;
    }

    // ====================================================
    // === 用R的多波前行求解 ===============================
    // ====================================================

    Stair = NULL ;
    fm = 0 ;
    h = 0 ;
    t = 0 ;

    // 从row2 = QR-num->rank + n1rows开始，这是[A Binput]的组合R因子的最后一行

    row2 = QRnum->rank + n1rows ;
    for (f = nf-1 ; f >= 0 ; f--)
    {

        // ---------------------------------------------------------------------
        // 得到F前面的R块
        // ---------------------------------------------------------------------

        R = Rblock [f] ;
        col1 = Super [f] ;                  
        col2 = Super [f+1] ;                
        fp = col2 - col1 ;                  
        pr = Rp [f] ;                       
        fn = Rp [f+1] - pr ;                

        if (keepH)
        {
            Stair = HStair + pr ;           
            fm = Hm [f] ;                   
            h = 0 ;                         
        }

        // ---------------------------------------------------------------------
        // 在R或RH块中找到活的主列
        // ---------------------------------------------------------------------

        rm = 0 ;                            
        for (k = 0 ; k < fp ; k++)
        {
            j = col1 + k ;
            if (keepH)
            {
                t = Stair [k] ;             

                if (t == 0)
                {
                    live = FALSE ;          
                    t = rm ;                
                    h = rm ;
                }
                else
                {
                    live = (rm < fm) ;      
                    h = rm + 1 ;           
                }

            }
            else
            {
                live = (!Rdead [j])  ;
            }

            if (live)
            {
                // R (rm,k) 是对角线; rm 和 k 是本地数值
                // 跟踪指向第一个条目R(0,k)的指针
                Rcolp [rm] = R ;
                Rlive [rm] = j ;
                rm++ ;
            }
            else
            {
                // 计算基础解; 死列为0
                ii = Q1fill ? Q1fill [j+n1cols] : j+n1cols ;
                if (ii < n)
                {
                    for (kk = 0 ; kk < nrhs ; kk++)
                    {
                        X [INDEX (ii,kk,n)] = 0 ;
                    }
                }
            }

            R += rm + (keepH ? (t-h) : 0) ;
        }


        row1 = row2 - rm ;

        // ---------------------------------------------------------------------
        // 得到这些rm方程的右边
        // ---------------------------------------------------------------------

        
        W1 = W ;
        B1 = B ;
        for (kk = 0 ; kk < nrhs ; kk++)
        {
            for (i = 0 ; i < rm ; i++)
            {
                ii = row1 + i ;

                W1 [i] = (ii < rank) ? B1 [ii] : 0 ;
            }
            W1 += rm ;
            B1 += ldb ;
        }

        // ---------------------------------------------------------------------
        // 用R (W = W - R2*x2)的矩形部分求解
        // ---------------------------------------------------------------------

        for ( ; k < fn ; k++)
        {
            j = Rj [pr + k] ;
            
            ii = Q1fill ? Q1fill [j+n1cols] : j+n1cols ;
            
            if (ii >= n) break ;

            if (!Rdead [j])
            {
                W1 = W ;
                for (kk = 0 ; kk < nrhs ; kk++)
                {
                    xi = X [INDEX (ii,kk,n)] ;       
                    if (xi != (double) 0)
                    {
                        FLOP_COUNT (2*rm) ;
                        for (i = 0 ; i < rm ; i++)
                        {
                            W1 [i] -= R [i] * xi ;
                        }
                    }
                    W1 += rm ;
                }
            }

            // 转到R的下一列
            R += rm ;
            if (keepH)
            {
                t = Stair [k] ;             
                
                h = MIN (h+1, fm) ;         
                
                R += (t-h) ;
            }
        }

        // ---------------------------------------------------------------------
        // 用R的压缩上三角形部分求解
        // ---------------------------------------------------------------------

        for (k = rm-1 ; k >= 0 ; k--)
        {
            R = Rcolp [k] ;               
            j = Rlive [k] ;                
            ii = Q1fill ? Q1fill [j+n1cols] : j+n1cols ;
            
            if (ii < n)
            {
                W1 = W ;
                for (kk = 0 ; kk < nrhs ; kk++)
                {

                    xi = qr_divide (W1 [k], R [k], cc) ;
                    FLOP_COUNT (1) ;
                    X [INDEX(ii,kk,n)] = xi ;
                    if (xi != (double) 0)
                    {
                        FLOP_COUNT (2*k) ;
                        for (i = 0 ; i < k ; i++)
                        {
                            W1 [i] -= R [i] * xi ;
                        }
                    }
                    W1 += rm ;
                }
            }
        }

        // ---------------------------------------------------------------------
        // 由波前阵f-1 准备 R的块
        // ---------------------------------------------------------------------

        row2 = row1 ;
    }

    // ========================================================
    // === 用R的单例行来求解 ==================================
    // ========================================================

    FLOP_COUNT ((n1rows <= 0) ? 0 :
        nrhs * (n1rows + (2 * (R1p [n1rows] - n1rows)))) ;

    for (kk = 0 ; kk < nrhs ; kk++)
    {
        for (i = n1rows-1 ; i >= 0 ; i--)
        {
            // 得到右边的第 i 单例元素行
            double x = B [i] ;
            
            for (Long p = R1p [i] + 1 ; p < R1p [i+1] ; p++)
            {
                Long jnew = R1j [p] ;
                
                Long jold = Q1fill ? Q1fill [jnew] : jnew ;
                
                x -= R1x [p] * X [jold] ;
            }
            // 除以“对角线”(单元素条目本身)
            Long p = R1p [i] ;
            Long jnew = R1j [p] ;
            Long jold = Q1fill ? Q1fill [jnew] : jnew ;

            X [jold] = qr_divide (x, R1x [p], cc) ;
        }
        B += ldb ;
        X += n ;
    }
}

// ====================================
// === qr_private_rtsolve ===========
// ====================================
void qr_private_rtsolve
(
    // inputs
    SparseQR_factorization  *QR,
    int use_Q1fill,

    Long nrhs,              
    Long ldb,               
    double *B,               

    double *X,               

    sparse_common *cc
)
{
    double xi ;
    qr_symbolic *QRsym ;
    qr_numeric  *QRnum ;
    Long *R1p, *R1j, *Rmap, *Rp, *Rj, *Super, *HStair, *Hm, *Stair, *Q1fill,
        *RmapInv ;
    double *R1x, **Rblock, *R, *X1, *X2 ;
    char *Rdead ;
    Long i, j, k, m, n, p, kk, n1rows, n1cols, rank, nf, f, col1, col2, fp, pr,
        fn, rm, row1, keepH, fm, h, t, live, jj ;

    // -------------------------------------------------------------------------
    // 获取QR分解对象的内容
    // -------------------------------------------------------------------------

    QRsym = QR->QRsym ;
    QRnum = QR->QRnum ;
    n1rows = QR->n1rows ;
    n1cols = QR->n1cols ;
    n = QR->nacols ;
    m = QR->narows ;
    Q1fill = use_Q1fill ? QR->Q1fill : NULL ;
    R1p = QR->R1p ;
    R1j = QR->R1j ;
    R1x = QR->R1x ;
    Rmap = QR->Rmap ;
    RmapInv = QR->RmapInv ;
    rank = QR->rank ;
    

    keepH = QRnum->keepH ;

    nf = QRsym->nf ;
    Rblock = QRnum->Rblock ;
    Rp = QRsym->Rp ;
    Rj = QRsym->Rj ;
    Super = QRsym->Super ;
    Rdead = QRnum->Rdead ;

    HStair = QRnum->HStair ;
    Hm = QRnum->Hm ;

    // -------------------------------------------------------------------------
    // X = E'*B or X = B
    // -------------------------------------------------------------------------

    X1 = X ;
    if (rank == n)
    {
        
        for (kk = 0 ; kk < nrhs ; kk++)
        {
            for (k = 0 ; k < n ; k++)
            {
                X1 [k] = B [Q1fill ? Q1fill [k] : k] ;
            }
            for ( ; k < m ; k++)
            {
                X1 [k] = 0 ;
            }
            X1 += m ;
            B += ldb ;
        }
    }
    else
    {
        // R是压缩形式;使用到梯形的映射
        for (kk = 0 ; kk < nrhs ; kk++)
        {
            for (i = 0 ; i < rank ; i++)
            {
                k = RmapInv [i] ;
                
                Long knew = Q1fill ? Q1fill [k] : k ;
                
                X1 [i] = B [knew] ;
            }
            for ( ; i < m ; i++)
            {
                X1 [i] = 0 ;
            }
            X1 += m ;
            B += ldb ;
        }
    }

    // =======================================================
    // === 用R的单例行来求解 ==================================
    // =======================================================

    X1 = X ;
    for (kk = 0 ; kk < nrhs ; kk++)
    {
        for (i = 0 ; i < n1rows ; i++)
        {
            // 除以对角线(单元素条目本身)
            p = R1p [i] ;
            xi = qr_divide (X1 [i], R1x [p], cc) ;
            X1 [i] = xi ;

            // 用非对角元素求解
            for (p++ ; p < R1p [i+1] ; p++)
            {
                j = R1j [p] ;
                
                k = Rmap ? Rmap [j] : j ;
                
                if (k < rank)
                {

                    X1 [k] -= R1x [p] * xi ;
                }
            }
        }
        X1 += m ;
    }

    // ====================================================
    // === 用R的多波前行求解 ===============================
    // ====================================================
    Stair = NULL ;
    fm = 0 ;
    h = 0 ;
    t = 0 ;


    row1 = n1rows ;
    for (f = 0 ; f < nf ; f++)
    {

        // ---------------------------------------------------------------------
        // 从波前阵F 中获取 R 块
        // ---------------------------------------------------------------------

        R = Rblock [f] ;
        col1 = Super [f] ;                  
        col2 = Super [f+1] ;                
        fp = col2 - col1 ;                  
        pr = Rp [f] ;                       
        fn = Rp [f+1] - pr ;               

        if (keepH)
        {
            Stair = HStair + pr ;           
            fm = Hm [f] ;                  
            h = 0 ;                        
        }

        // ---------------------------------------------------------------------
        // 用R的压缩上三角形部分求解
        // ---------------------------------------------------------------------

        rm = 0 ;                            
        for (k = 0 ; k < fp ; k++)
        {
            j = col1 + k ;
            
            if (j+n1cols >= n) return ;     
            if (keepH)
            {
                t = Stair [k] ;             
                
                if (t == 0)
                {
                    live = FALSE ;          
                    t = rm ;               
                    h = rm ;
                }
                else
                {
                    live = (rm < fm) ;      
                    h = rm + 1 ;            
                }
                
            }
            else
            {

                live = (!Rdead [j]) ;                                  
            }

            if (live)
            {
                
                X1 = X + row1 ;
                
                for (kk = 0 ; kk < nrhs ; kk++)
                {
                    xi = X1 [rm] ;     
                    for (i = 0 ; i < rm ; i++)
                    {
                        xi -=  R [i] * X1 [i] ;
                    }
                    X1 [rm] = qr_divide (xi, R [rm], cc) ;
                    X1 += m ;
                }
                // rm递增后，它表示到目前为止在R块中的行数
                rm++ ;
            }

            // 在R块中前进到R的下一列
            R += rm + (keepH ? (t-h) : 0) ;
        }

        // ---------------------------------------------------------------------
        // 用R的矩形部分求解
        // ---------------------------------------------------------------------

        for ( ; k < fn ; k++)
        {
            j = Rj [pr + k] ;
            
            jj = j + n1cols ;
            if (jj >= n) break ;            

            Long ii = Rmap ? Rmap [jj] : jj ;
            
            if (ii < rank)
            {
                X2 = X + row1 ;                   
                X1 = X ;
                for (kk = 0 ; kk < nrhs ; kk++)
                {
                    xi = X1 [ii] ;                  
                    for (i = 0 ; i < rm ; i++)
                    {
                        xi -= R [i] * X2 [i] ;
                    }
                    X1 [ii] = xi ;
                    X1 += m ;
                    X2 += m ;
                }
            }

            // 到第 R 行
            R += rm ;
            if (keepH)
            {
                t = Stair [k] ;             
                
                h = MIN (h+1, fm) ;         
                
                R += (t-h) ;
            }
        }

        row1 += rm ;    
    }
    
}