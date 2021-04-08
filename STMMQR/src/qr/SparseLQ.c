/******************************************************************************
 * VERSION: 1.1
 * DATE:    2020年9月27日
 * FILE:    SparseLQ.c
 * BRIEF:   LQ 分解函数
 * FUNCTION: SparseLQ函数，用来计算矩阵的LQ分解，并保留Household变换
 *****************************************************************************/
/*******************************
 *         INCLUDE
 ******************************/ 
#include "SparseQR.h"
#include <float.h>
#include <string.h>

/*******************************
 *         FUNCTION
 ******************************/ 
void PRINT_CSR( Sparse_long nrows, Sparse_long ncols, Sparse_long nnz, 
    Sparse_long *rowIndex, Sparse_long *columns, double *values)
{
    int i;
    if ( nrows <= 0  || ncols <= 0 || nnz <= 0 ) 
    {
        printf(" PRINT Error: Matrix value is wrong! \n");
        return ;
    }

    printf ("\n rowindex = ");
    for ( i = 0; i < nrows + 1; i++)
    {
        printf ("%d ", rowIndex[i]);
    }

    printf ("\n columns = ");
    for ( i = 0; i < nnz; i++)
    {
        printf ("%d ", columns[i]);
    }
    
    printf ("\n values = ");
    for ( i = 0; i < nnz; i++)
    {
        printf ("%lf ", values[i]);
    }
    printf ("\n");
}

int CSR_transpose( Sparse_long* nrows, Sparse_long* ncols, Sparse_long nnz,
    Sparse_long** RowIndex, Sparse_long** columns, double** values)
{
    if ( *nrows <= 0  || *ncols <= 0 || nnz <= 0 ) 
    {
        printf(" Error: Matrix values wrong! \n");
        return 0;
    }
    int i, j;
    // 交换行列
    Sparse_long temp;
    temp = *nrows;
    *nrows = *ncols;
    *ncols = temp;

    Sparse_long *rowIndex_out, *columns_out; // 更新后的行索引、列和值
    double *values_out;

    rowIndex_out = (Sparse_long*)malloc( (*nrows + 2) * sizeof(Sparse_long));
    columns_out = (Sparse_long*)malloc( nnz * sizeof(Sparse_long));
    values_out = (double *)malloc( nnz * sizeof(double));

    memset(rowIndex_out, 0, (*nrows + 2) * sizeof(Sparse_long));
    memset(columns_out, 0, nnz * sizeof(Sparse_long));
    // count per column
    for ( i = 0; i < nnz; i++)
    {
        ++rowIndex_out[ (*columns)[i] + 2];
    }
    
    // from count per column generate new rowPtr (but shifted)
    for ( i = 2; i < (*nrows + 2); ++i) 
    {
        // create incremental sum
        rowIndex_out[i] += rowIndex_out[i - 1];
    }

    // perform the main part
    for ( i = 0; i < *ncols; ++i) {
        for ( j = (*RowIndex)[i]; j < (*RowIndex)[i + 1]; ++j) {
            // calculate index to transposed matrix at which we should place current element, and at the same time build final rowPtr
            const int new_index = rowIndex_out[ (*columns)[j] + 1]++;
            values_out[new_index] = *values[j];
            columns_out[new_index] = i;
        }
    }
    Sparse_long *ptr_row_t = *RowIndex;
    *RowIndex = rowIndex_out;
    free(ptr_row_t);

    Sparse_long *ptr_col_t = *columns;
    *columns = columns_out;
    free(ptr_col_t);

    double *ptr_val_t = *values;
    *values = values_out;
    free(ptr_val_t);

}

void qr_rcount
(
    // inputs, not modified
    qr_symbolic *QRsym,
    qr_numeric *QRnum,

    Long n1rows,        
    Long econ,          
    Long n2,            
    int getT,          

    // input/output
    Long *Ra,           
    Long *Rb,           

    Long *H2p,         

    Long *p_nh          
)
{
    double **Rblock, *R, *Tau, *HTau ;
    Long *Rp, *Rj, *Super, *HStair, *Stair, *Hm ;
    char *Rdead ;
    Long nf, j, f, col1, fp, pr, fn, rm, k, i, t, fm, h, getRa, getRb, nh,
        row1, keepH, getH, hnz ;

    // -------------------------------------------------------------------------
    // 获取QRsym和QRnum对象的内容
    // -------------------------------------------------------------------------

    keepH = QRnum->keepH ; 

    getRa = (Ra != NULL) ;
    getRb = (Rb != NULL) ;
    getH  = (H2p != NULL && p_nh != NULL) && keepH ;
    if (!(getRa || getRb || getH))
    {
        return ;
    }

    nf = QRsym->nf ;
    // n = QRsym->n ;
    Rblock = QRnum->Rblock ;
    Rp = QRsym->Rp ;
    Rj = QRsym->Rj ;
    Super = QRsym->Super ;
    Rdead = QRnum->Rdead ;

    HStair = QRnum->HStair ;
    HTau = QRnum->HTau ;
    Hm = QRnum->Hm ;
    Stair = NULL ;
    Tau = NULL ;
    fm = 0 ;
    h = 0 ;
    t = 0 ;
    nh = 0 ;
    hnz = 0 ;

    // -------------------------------------------------------------------------
    // 检查每个波前阵的压缩块
    // -------------------------------------------------------------------------

    row1 = n1rows ;
    for (f = 0 ; f < nf ; f ++)
    {
        R = Rblock [f] ;
        col1 = Super [f] ;                  
        fp = Super [f+1] - col1 ;          
        pr = Rp [f] ;                       
        fn = Rp [f+1] - pr ;                

        if (keepH)
        {
            Stair = HStair + pr ;           
            Tau = HTau + pr ;               
            fm = Hm [f] ;                   
            h = 0 ;                         
        }

        rm = 0 ;                            
        for (k = 0 ; k < fn ; k++)
        {
            // -----------------------------------------------------------------
            // 得到列和它的阶梯staircase
            // -----------------------------------------------------------------

            if (k < fp)
            {
                
                j = col1 + k ;
                
                if (keepH)
                {
                    t = Stair [k] ;             
                    if (t == 0)
                    {
                        t = rm ;               
                    }
                    else if (rm < fm)
                    {
                        rm++ ;                 
                    }
                    h = rm ;                    
                }
                else
                {
                    if (!Rdead [j])
                    {
                        rm++ ;                 
                    }
                }
            }
            else
            {
                //波前阵 F 的非关键列
                j = Rj [pr + k] ;
                
                if (keepH)
                {
                    t = Stair [k] ;             
                    h = MIN (h+1, fm) ;         
                }
            }

            // -----------------------------------------------------------------
            // 对这个R块计数nnz (R(0:econ-1,j))
            // -----------------------------------------------------------------

            for (i = 0 ; i < rm ; i++)
            {
                double rij = *(R++) ;
                if (rij != (double) 0)
                {
                    if (j < n2)
                    {
                        if (getRa && row1 + i < econ)
                        {
                            Ra [j]++ ;
                        }
                    }
                    else
                    {
                        if (getRb && row1 + i < econ)
                        {
                            if (getT)
                            {
                                Rb [row1+i]++ ;
                            }
                            else
                            {
                                Rb [j-n2]++ ;
                            }
                        }
                    }
                }
            }

            // -----------------------------------------------------------------
            //  计算 nnz (H (:,pr+k))
            // -----------------------------------------------------------------

            if (keepH && t >= h)
            {
                // Household 变换非空
                if (getH && Tau [k] != (double) 0)
                {
                    H2p [nh++] = hnz++ ;    
                    for (i = h ; i < t ; i++)
                    {
                        double hij = *(R++) ;
                        if (hij != (double) 0)
                        {
                            hnz++ ;         
                        }
                    }
                }
                else
                {
                    R += (t-h) ;            
                }
            }
        }
        row1 += rm ;                       
    }

    // -------------------------------------------------------------------------
    // 最后确定H的列指针
    // -------------------------------------------------------------------------

    if (getH)
    {
        H2p [nh] = hnz ;
        *p_nh = nh ;
    }
}

void qr_rconvert
(
    // inputs, not modified
    qr_symbolic *QRsym,
    qr_numeric *QRnum,

    Long n1rows,        
    Long econ,         
    Long n2,           
    int getT,          

    // input/output
    Long *Rap,         

    // output, not defined on input
    Long *Rai,          
    double *Rax,        

    // input/output
    Long *Rbp,          

    // output, not defined on input
    Long *Rbi,          
    double *Rbx,         

    // input
    Long *H2p,          

    // output, not defined on input
    Long *H2i,          
    double *H2x,        

    double *H2Tau        
)
{
    double rij, hij ;
    double **Rblock, *R, *Tau, *HTau ;
    Long *Rp, *Rj, *Super, *HStair, *Hii, *Stair, *Hip, *Hm, *Hi ;
    char *Rdead ;
    Long nf, j, f, col1, fp, pr, fn, rm, k, i, p, getRa, getRb, row1, fm,
        h, getH, keepH, ph, t, nh ;

    // -------------------------------------------------------------------------
    // 获取QRsym和QRnum 结构体的内容
    // -------------------------------------------------------------------------

    keepH = QRnum->keepH ;
    getRa = (Rap != NULL && Rai != NULL && Rax != NULL) ;
    getRb = (Rbp != NULL && Rbi != NULL && Rbx != NULL) ;
    getH  = (H2p != NULL && H2i != NULL && H2x != NULL && H2Tau != NULL)
            && keepH ;
    if (!(getRa || getRb || getH))
    {
        
        return ;
    }

    nf = QRsym->nf ;
    Rblock = QRnum->Rblock ;
    Rp = QRsym->Rp ;
    Rj = QRsym->Rj ;
    Super = QRsym->Super ;
    Rdead = QRnum->Rdead ;

    HStair = QRnum->HStair ;
    HTau = QRnum->HTau ;
    Hm = QRnum->Hm ;
    Hii = QRnum->Hii ;
    Hip = QRsym->Hip ;
    Stair = NULL ;
    Hi = NULL ;
    Tau = NULL ;
    fm = 0 ;
    h = 0 ;
    t = 0 ;
    nh = 0 ;

    // -------------------------------------------------------------------------
    // 为每个波前阵F转换压缩块
    // -------------------------------------------------------------------------

    row1 = n1rows ;
    ph = 0 ;                                
    for (f = 0 ; f < nf ; f ++)
    {
        R = Rblock [f] ;
        col1 = Super [f] ;                 
        fp = Super [f+1] - col1 ;           
        pr = Rp [f] ;                       
        fn = Rp [f+1] - pr ;               

        if (keepH)
        {
            Stair = HStair + pr ;           
            Tau = HTau + pr ;               
            Hi = &Hii [Hip [f]] ;           
            fm = Hm [f] ;                   

            h = 0 ;                         
        }

        // ---------------------------------------------------------------------
        // 提取R或R+H块的每一列
        // ---------------------------------------------------------------------

        rm = 0 ;                            
        for (k = 0 ; k < fn ; k++)
        {
            // -----------------------------------------------------------------
            // 获取列和它的staircase
            // -----------------------------------------------------------------

            if (k < fp)
            {
                j = col1 + k ;

                if (keepH)
                {
                    t = Stair [k] ;             

                    if (t == 0)
                    {
                        t = rm ;                
                    }
                    else if (rm < fm)
                    {
                        rm++ ;                  
                    }
                    h = rm ;                    
                }
                else
                {
                    
                    if (!Rdead [j])
                    {
                        rm++ ;                  
                    }

                }
            }
            else
            {

                j = Rj [pr + k] ;
                
                
                if (keepH)
                {
                    t = Stair [k] ;             
                    h = MIN (h+1, fm) ;      
                }
            }

            // -----------------------------------------------------------------
            // 提取R的列
            // -----------------------------------------------------------------

            for (i = 0 ; i < rm ; i++)
            {
                rij = *(R++) ;

                if (rij != (double) 0)
                {
                    if (j < n2)
                    {
                        if (getRa && row1 + i < econ)
                        {
                            p = Rap [j]++ ;
                            Rai [p] = row1 + i ;
                            Rax [p] = rij ;
                        }
                    }
                    else
                    {
                        if (getRb && row1 + i < econ)
                        {
                            
                                p = Rbp [j-n2]++ ;
                                Rbi [p] = row1 + i ;
                                Rbx [p] = rij ;
                        }
                    }
                }
            }

            // -----------------------------------------------------------------
            // 从 H 中提取列
            // -----------------------------------------------------------------

            if (keepH && t >= h)
            {

                if (getH && Tau [k] != (double) 0)
                {
                    H2Tau [nh++] = Tau [k] ;
                    H2i [ph] = Hi [h-1] + n1rows ;  
                    H2x [ph] = 1 ;
                    ph++ ;
                    for (i = h ; i < t ; i++)
                    {
                        hij = *(R++) ;
                        if (hij != (double) 0)
                        {
                            H2i [ph] = Hi [i] + n1rows ;
                            H2x [ph] = hij ;
                            ph++ ;
                        }
                    }
                }
                else
                {
                    R += (t-h) ;            
                }
            }
        }
        row1 += rm ;                        
    }
}

// 返回 R 的rank， 失败则为EMPTY
Long qr_trapezoidal 
(
    // inputs, not modified

    Long n,         
    Long *Rp,       
    Long *Ri,       
    double *Rx,      

    Long bncols,    

    Long *Qfill,    

    int skip_if_trapezoidal,        

    // outputs, not allocated on input
    Long **p_Tp,    
    Long **p_Ti,    
    double **p_Tx,   

    Long **p_Qtrap, 
    // workspace and parameters
    sparse_common *cc
)
{
    double *Tx ;
    Long *Tp, *Ti, *Qtrap ;
    Long rnz, i, rank, k, p, pend, len, t1nz, t2nz, k1, k2, p1, p2, found_dead,
        is_trapezoidal ;

    // -------------------------------------------------------------------------
    // 求R、nnz(T1)、nnz(T2)的秩
    // -------------------------------------------------------------------------

    rank = 0 ;              
    t1nz = 0 ;              
    t2nz = 0 ;              
    found_dead = FALSE ;    
    is_trapezoidal = TRUE ; 

    *p_Tp = NULL ;
    *p_Ti = NULL ;
    *p_Tx = NULL ;
    *p_Qtrap = NULL ;

    for (k = 0 ; k < n ; k++)
    {
        
        p = Rp [k] ;
        pend = Rp [k+1] ;
        len = pend - p ;
        i = (len > 0) ? Ri [pend - 1] : EMPTY ;

        if (i > rank)
        {
            
            return (EMPTY) ;
        }
        else if (i == rank)
        {
            rank++ ;
            t1nz += len ;
            if (found_dead)
            {

                is_trapezoidal = FALSE ;
            }
        }
        else
        {
            found_dead = TRUE ;
            t2nz += len ;
        }
    }

    // -------------------------------------------------------------------------
    // 若已梯形，可快速返回
    // -------------------------------------------------------------------------

    if (is_trapezoidal)
    {
        if (skip_if_trapezoidal)
        {
            return (rank) ;
        }
    }

    // -------------------------------------------------------------------------
    // 分配结果(T和Qtrap)
    // -------------------------------------------------------------------------

    rnz = Rp [n] ;

    Tp    = (Long  *) SparseCore_malloc (n+1,      sizeof (Long),  cc) ;
    Ti    = (Long  *) SparseCore_malloc (rnz,      sizeof (Long),  cc) ;
    Tx    = (double *) SparseCore_malloc (rnz,      sizeof (double), cc) ;
    Qtrap = (Long  *) SparseCore_malloc (n+bncols, sizeof (Long),  cc) ;

    if (cc->status < SPARSE_OK)
    {
        //  内存溢出
        SparseCore_free (n+1,      sizeof (Long),  Tp,    cc) ;
        SparseCore_free (rnz,      sizeof (Long),  Ti,    cc) ;
        SparseCore_free (rnz,      sizeof (double), Tx,    cc) ;
        SparseCore_free (n+bncols, sizeof (Long),  Qtrap, cc) ;
        return (EMPTY) ;
    }


    // -------------------------------------------------------------------------
    // 找到列指针Tp并排列Qtrap
    // -------------------------------------------------------------------------

    k1 = 0 ;                
    k2 = rank ;             
    p1 = 0 ;                
    p2 = t1nz ;             
    rank = 0 ;              

    for (k = 0 ; k < n ; k++)
    {
        p = Rp [k] ;
        pend = Rp [k+1] ;
        len = pend - p ;
        i = (len > 0) ? Ri [pend - 1] : EMPTY ;

        if (i == rank)
        {
            rank++ ;
            Tp [k1] = p1 ;
            Qtrap [k1] = Qfill ? Qfill [k] : k ;
            k1++ ;
            for ( ; p < pend ; p++)
            {
                Ti [p1] = Ri [p] ;
                Tx [p1] = Rx [p] ;
                p1++ ;
            }
        }
        else
        {
            Tp [k2] = p2 ;
            Qtrap [k2] = Qfill ? Qfill [k] : k ;
            k2++ ;
            for ( ; p < pend ; p++)
            {
                Ti [p2] = Ri [p] ;
                Tx [p2] = Rx [p] ;
                p2++ ;
            }
        }
    }

    for ( ; k < n+bncols ; k++)
    {
        Qtrap [k] = Qfill ? Qfill [k] : k ;
    }

    // -------------------------------------------------------------------------
    // 最后确定列指针并返回结果
    // -------------------------------------------------------------------------

    
    Tp [n] = rnz ;
    *p_Tp = Tp ;
    *p_Ti = Ti ;
    *p_Tx = Tx ;
    *p_Qtrap = Qtrap ;
    return (rank) ;
}


SparseQR_factorization *SparseLQ(
    int ordering, double tol,
    sparse_csc *A, sparse_common *cc 
)
{
    Long i, j, k, p, m, n, n2, econ, xtype, rank, getR, getE, nh, rnz, pr,
        getT, getH, n1cols, rjsize, n1rows, bncols, k2;
    Long *R1p, *R1j, *Rp, *Ri, *Rap, *H2p, *Hp, *Hi, *E, *HP1inv, *Q1fill;
    double *R1x, *Rx, *Hx;

    qr_symbolic *LQsym ;
    qr_numeric *LQnum ;
    SparseQR_factorization *LQ ;
    sparse_csc *AT, *R, *H;  // A转置后做QR分解，得到的R也需要转置获得矩阵L
    dense_array *HTau;
    nh = 0;
    rnz = 0;
    bncols = 0;
    xtype = 1 ;
    Rap = NULL ;

    R = NULL;
    H = NULL;
    H2p = NULL;
    Rp = NULL ;
    Ri = NULL ;
    Rx = NULL ;

    H = NULL ;
    Hp = NULL ;
    Hi = NULL ;
    Hx = NULL ;
    H2p = NULL ;

    AT = SparseCore_transpose (A, 1, cc) ; 
    m = AT->nrow;
    n = AT->ncol;
    // printf("AT:  %d-by-%d \n", m, n);
    econ = AT->ncol;
    
    LQ = SparseQR (QR_ORDERING_DEFAULT, tol, AT, cc, "LQresult.txt") ;
    // 本质上还是QR分解的sym和num 
    
    return (LQ) ;
}