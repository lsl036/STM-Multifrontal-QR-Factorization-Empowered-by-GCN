
#include "SparseQR.h"
#include <sys/time.h>
#include <float.h>
#include <stdio.h>
#include "tpsm.h"
#include <math.h>
#define Long Sparse_long

double check_error (sparse_csc *A, SparseQR_factorization *QR, sparse_common *cc)
{
    int i, j;
    double one [2] = {1,0}, minusone [2] = {-1,0}, zero [2] = {0,0} ;
    dense_array *X, *Y, *B, *X_sol; // 原始X，AX=B， 求解得到X_sol

    Long n = A->ncol;
    // n = ones (n,1) 一个稠密的向量: 给稠密矩阵分配空间并置为0 
    X = SparseCore_zeros (n, 1, A->xtype, cc) ;
    B = SparseCore_zeros (n, 1, A->xtype, cc) ;
    //X_sol = SparseCore_zeros (n, 1, A->xtype, cc) ;

    // 初始化X = [0, 1, 2,..., n]
    double *x ;
    x = (double *) (X->x) ;
    for (i = 0 ; i < n ; i++)
    {
        x [i] = i ;
    }
    // 计算B = A*X
    SparseCore_sdmult (A, 0, one, zero, X, B, cc) ; 
    // 计算Y = Q'*B  (--> Rx = Y)
    Y = QR_qmult(QR_QTX, QR, B, cc);
    // 求解X = R\(E*Y)    (即 inv(R)*(E*Y) )
    X_sol = QR_solve (QR_RETX_EQUALS_B, QR, Y, cc) ;

    // 计算两个X的误差
    double diff_norm = 0.0;
    double xx = 0.0;
    double *x_check = (double *) (X_sol->x);
    for (j = 0; j < X->nrow; j++)
    {
        xx =  x_check [j] - j;
        diff_norm += xx * xx;
    }
    diff_norm = sqrt(diff_norm)/X->nrow;

    SparseCore_free_dense (&Y, cc) ;
    SparseCore_free_dense (&X, cc) ;
    SparseCore_free_dense (&B, cc) ;
    SparseCore_free_dense (&X_sol, cc) ;
    return diff_norm;
}

/**
 * @brief   根据矩阵文件所在根目录以及矩阵名称得到完整路径
 * 
 * @param file_path 
 * @param matrix_name 
 * @param matrix 
 * @return int 
 */

// =============================================================================
int main (int argc, char* argv[])
{
    Long rnk;
    double res;
    int cycleNum = 1;
    double timeStart, timeEnd;
    struct timeval tv;
    sparse_common Common, *cc ;
    sparse_csc *A ;
    int mtype ,j;
    Long m, n ;
    
    // start HNUCHOL
    cc = &Common ;
    SparseCore_start (cc) ;

    /* 矩阵路径 */
    char *fmatrix = argv[1];
    FILE *fp;
    fp = fopen(fmatrix,"r");
    if(fp == NULL) {
        printf("%s file is not exist!\n", fmatrix);
        return 0;
    }
    
    A = (sparse_csc *) SparseCore_read_matrix (fp, 1, &mtype, cc) ;

    /* 结果输出到文件中 */
    char *result_file = (char *)"./ARMQR.txt";
    FILE* fresult;
    fresult = fopen(result_file,"a+");
    fprintf(fresult, "%s\t", fmatrix);

    if (mtype != SPARSE_CSC)
    {
        printf ("input matrix must be sparse\n") ;
        exit (1) ;
    }

    // [m n] = size (A) ;
    m = A->nrow ;
    n = A->ncol ;
    
    printf ("Matrix %6ld-by-%-6ld nnz: %6ld\n", m, n, SparseCore_nnz (A, cc)) ;

    // 提前计算 tol
    double tol = QR_DEFAULT_TOL;
    double max2Norm = 0.0;
    max2Norm = qr_maxcolnorm (A, cc);
    if (max2Norm == 0)
    {
        max2Norm = 1;
    }
    tol = 20 * ((double) A->nrow + (double) A->ncol) * DBL_EPSILON * max2Norm;

    int CORE = 128;
    cc->SPQR_grain = (double) (CORE * 2);
    
    cc->status = SPARSE_OK ;
    // 开线程池
    TPSM_t *tpool = TPSM_init(128, 0 );  // 0 为均分线程， 1为集中线程

    if (A->xtype == SPARSE_REAL)
    {
        SparseQR_factorization *QR ;

        Long i ;
        gettimeofday(&tv, NULL);
        timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
        
        for (i = 0; i < cycleNum; ++i) {
            // factorize ORDERING_DEFAULT = 7
            QR = SparseQR (tpool, QR_ORDERING_DEFAULT, tol, A, cc) ;
        }
        gettimeofday(&tv, NULL);
        timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
        printf("HnuSparseQR run time: %f\n\n", (timeEnd-timeStart) / cycleNum);
        fprintf(fresult, "%lf\t", (timeEnd-timeStart) / cycleNum);
        // 等待全部完成,销毁线程池
        TPSM_destroy(tpool, 2);
        /***************************************
        * QR求解  Ax = b   ( QRx = b ) 误差检查
        ****************************************/
        res = check_error(A, QR, cc);
        printf ("res = %8.1e\n", res);
        fprintf (fresult, "%8.1e\n", res) ;
        // free QR
        SparseQR_free (&QR, cc) ;
    }
    // -------------------------------------------------------------------------
    // free everything that remains
    // -------------------------------------------------------------------------
    SparseCore_free_sparse (&A, cc) ;
    SparseCore_finish (cc) ;
    fclose(fp);
    fclose(fresult);
    return (0) ;
}
