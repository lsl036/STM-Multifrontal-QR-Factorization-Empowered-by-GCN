#ifndef SparseQR_H
#define SparseQR_H


/*******************************
 *         INCLUDE
 ******************************/ 
#include <string.h>
#include "SparseQR_struct.h"
#include "SparseQR_internal.h"
#include"tpsm.h"

/*******************************
 *         DEFINE
 ******************************/ 
size_t FCHUNK;
size_t SMALL;
size_t MINCHUNK;
size_t MINCHUNK_RATIO;

// =============================================
// === SparseQR  functions =====================
// =============================================

SparseQR_factorization  *SparseQR(
    int ordering,          
    double tol,            
    sparse_csc *A,      
    sparse_common *cc,
    char *result_name
);

int SparseQR_free(
    SparseQR_factorization **QR, 
    sparse_common *cc
);

void qr_freefac(
    SparseQR_factorization **QR_handle,
    sparse_common *cc
);

void qr_freenum(
    qr_numeric **QRnum_handle,
    sparse_common *cc
);

void qr_freesym(
    qr_symbolic **QRsym_handle,
    sparse_common *cc
);

int qr_1colamd ( 
    int ordering,           
    double tol,            
    sparse_csc *A,     
    Long **p_Q1fill,        
    Long **p_R1p,          
    Long **p_P1inv,        
    sparse_csc **p_Y,   
    Long *p_n1cols,         
    Long *p_n1rows,         
    sparse_common *cc,
    char *result_name
);

double qr_tol(
    sparse_csc *A,
    sparse_common *cc
);

void qr_shift(
    Long n,
    Long *X                   
);

Long qr_cumsum(            
    Long n,
    Long *X                
);

double qr_maxcolnorm
(
    sparse_csc *A,
    sparse_common *cc
);

int qr_rmap
(
    SparseQR_factorization *QR,
    sparse_common *cc
);

/*************************************************/
// ********* End of SparseQR Function ********** //
/*************************************************/


/*************************************************/
// ********* qr_analyze Function ************ //
/*************************************************/
qr_symbolic *qr_analyze(
    sparse_csc *A,
    int ordering,           
    Long *Quser,           
    int do_rank_detection,  
    sparse_common *cc,
    char *result_name
);

void qr_stranspose1( 
    sparse_csc *A,  
    Long *Qfill,       
    Long *Sp,           
    Long *Sj,          
    Long *PLinv,       
    Long *Sleft,       
    Long *W           
) ;
/*************************************************/
// ********* End of qr_analyze  ************* //
/*************************************************/

/*************************************************/
// ********* qr_factorize Function ********** //
/*************************************************/
qr_numeric *qr_factorize
(
    sparse_csc **Ahandle,
    Long freeA,                   
    double tol,                    
    Long ntol,                     
    qr_symbolic *QRsym,
    sparse_common *cc
) ;

int chunk_getSettings
(
    size_t const _FCHUNK_size,
    size_t const _SMALL_size,
    size_t const _MINCHUNK_size,
    size_t const _MINCHUNK_RATIO
) ;

void qr_stranspose2( 
    sparse_csc *A, 
    Long *Qfill,        
    Long *Sp,          
    Long *PLinv,       
    double *Sx,         
    Long *W           
);

void qr_kernel(
    Long task,
    qr_blob *Blob
);

void qr_hpinv(
    qr_symbolic *QRsym,
    qr_numeric *QRnum,
    Long *W             
);

Long qr_fsize  (   
    Long f,
    Long *Super,          
    Long *Rp,              
    Long *Rj,             
    Long *Sleft,           
    Long *Child,          
    Long *Childp,          
    Long *Cm,             
    Long *Fmap,            
    Long *Stair            
);

void qr_assemble(
    Long f,               
    Long fm,               
    int keepH,              
    Long *Super,
    Long *Rp,
    Long *Rj,
    Long *Sp,
    Long *Sj,
    Long *Sleft,
    Long *Child,
    Long *Childp,
    double *Sx,
    Long *Fmap,
    Long *Cm,
    double **Cblock,
    Long *Hr,
    Long *Stair,
    Long *Hii,             
    Long *Hip,
    double *F,
    Long *Cmap
);

Long qr_csize  (  
    Long c,                
    Long *Rp,               
    Long *Cm,              
    Long *Super            
);

Long qr_front
(
    Long m,           
    Long n,
    Long npiv,         
    double tol,         
    Long ntol,         
    Long fchunk,       
    double *F,          
    Long *Stair,       
    char *Rdead,       
    double *Tau,         
    double *W,           
    double *wscale,
    double *wssq,

    sparse_common *cc
);

Long qr_fcsize   ( 
    Long m,                
    Long n,               
    Long npiv,             
    Long g                 
);

Long qr_cpack  (  
    Long m,                 
    Long n,                
    Long npiv,              
    Long g,                 
    double *F,               
    double *C                
);

Long qr_rhpack  (  
    int keepH,             
    Long m,                 
    Long n,                
    Long npiv,              
    Long *Stair,           
    double *F,              
    double *R,              
    Long *p_rm            
);

void qr_larftb
(
    int method,     
    Long m,         
    Long n,
    Long k,        
    Long ldc,       
    Long ldv,       
    double *V,       
    double *Tau,     
    double *C,       
    double *W,       
    sparse_common *cc
);

void qr_multithreads(
    Long ntasks,
    qr_blob *Blob
);

/*************************************************/
// ********* End of  qr_factorize  ********** //
/*************************************************/
// LQ 相关
int CSR_transpose( Sparse_long* nrows, Sparse_long* ncols, Sparse_long nnz,
    Sparse_long** RowIndex, Sparse_long** columns, double** values);

SparseQR_factorization *SparseLQ(
    int ordering, double tol,
    sparse_csc *A, sparse_common *cc 
);

void qr_rcount
(
    qr_symbolic *QRsym,
    qr_numeric *QRnum,

    Long n1rows,        
    Long econ,          
    Long n2,           
    int getT,           
    Long *Ra,           
    Long *Rb,           
    Long *H2p,          
    Long *p_nh          
);

void qr_rconvert
(
    qr_symbolic *QRsym,
    qr_numeric *QRnum,

    Long n1rows,        
    Long econ,         
    Long n2,            
    int getT,         
    Long *Rap,          
    Long *Rai,          
    double *Rax,         
    Long *Rbp,        
    Long *Rbi,         
    double *Rbx,         
    Long *H2p,          
    Long *H2i,         
    double *H2x,        

    double *H2Tau       
);

Long qr_trapezoidal 
(
    Long n,        
    Long *Rp,       
    Long *Ri,      
    double *Rx,     

    Long bncols,   

    Long *Qfill,   
    int skip_if_trapezoidal,        
    Long **p_Tp,   
    Long **p_Ti,   
    double **p_Tx,   

    Long **p_Qtrap,  
    sparse_common *cc
);

/**************************************
 *           求解器验证需要
**************************************/
Long qr_private_get_H_vectors
(
    Long f,               
    SparseQR_factorization *QR,
    double *H_Tau,          
    Long *H_start,         
    Long *H_end,          
    sparse_common *cc
);

Long qr_private_load_H_vectors
(
    Long h1,            
    Long h2,
    Long *H_start,      
    Long *H_end,       
    double *R,         
    double *V,          
    sparse_common *cc
);

void qr_panel
(
    int method,        
    Long m,
    Long n,
    Long v,            
    Long h,             
    Long *Vi,         
    double *V,           
    double *Tau,        
    Long ldx,
    double *X,          
    double *C,           
    double *W,          

    sparse_common *cc
);

void qr_private_Happly
(
    int method,          
    SparseQR_factorization *QR,
    Long hchunk,           
    Long m,
    Long n,
    double *X,               
    double *H_Tau,          
    Long *H_start,          
    Long *H_end,            
    double *V,              
    double *C,               
    double *W,               
    sparse_common *cc
);

// 返回规模为m*n 的 Y ，或者在失败的情况下返回NULL。
dense_array *QR_qmult
(
    int method,            
    SparseQR_factorization *QR,
    dense_array *Xdense, 
    sparse_common *cc
);

dense_array *QR_solve    
(
    int system,                 
    SparseQR_factorization  *QR,    
    dense_array *B,          
    sparse_common *cc
);//

// 使用 SparseQR 得到的 QR 对象，求解 X = E*(R\B) 或 X=R\B  
void qr_rsolve
(
    SparseQR_factorization *QR,
    int use_Q1fill,       
    Long nrhs,           
    Long ldb,             
    double *B,               
    double *X,            
    double **Rcolp,         
    Long *Rlive,           
    double *W,             
    sparse_common *cc
);

void qr_private_rtsolve
(
    SparseQR_factorization  *QR,
    int use_Q1fill,
    Long nrhs,              
    Long ldb,               
    double *B,              
    double *X,              
    sparse_common *cc
);
#endif 