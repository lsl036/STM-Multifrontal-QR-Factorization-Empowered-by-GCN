
#include "SparseQR.h"
#include <sys/time.h>
#include <float.h>
#include <stdio.h>
#include "tpsm.h"
#include <math.h>
#define Long Sparse_long

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

    /* 读源矩阵数据 的 路径 */
    char *fmatrix = argv[1];
    int graph_id = atoi(argv[2]);
    // 选择 QR内部使用的重排序方法 
    int method_selected = atoi(argv[3]);
    int ordering;
    switch (method_selected)
    {
    case 0:
        ordering = QR_ORDERING_AMD;
        //c.method [0].ordering = SPARSE_AMD;
        break;

    case 1:
        ordering = QR_ORDERING_COLAMD;  
        // c.method [0].ordering = SPARSE_COLAMD;
        break;

    case 2:
        ordering = QR_ORDERING_ONLYMETIS;
        // c.method [0].ordering = SPARSE_METIS;
        break;

    case 3:
        ordering = QR_ORDERING_NESDIS;
        // c.method [0].ordering = SPARSE_NESDIS;
        break;

    default:
        
        return -1;
    }

    FILE *fp;
    fp = fopen(fmatrix,"r");
    if(fp == NULL) {
        printf("%s file is not exist!\n", fmatrix);
        return 0;
    }
    /******************************
     * 对输入的名称字符串进行处理 *
     ******************************/
    int len = strlen(fmatrix), step = 0;
    // // 拼接目标文件的地址
    // char result_name[100]="./data_clean/";
    // 提取了矩阵名字
    char str2[50];
    for (int i = 23; i < 23+(len - 28)/2; i++)
    {
        str2[step++] = fmatrix[i];
    }
    
    /* Node 和Edge 文件 写入矩阵名字 */
    char *result_file = (char *)"./data_clean/QR_Node.txt";
    char *result_file_E = (char *)"./data_clean/QR_Edge.txt";
    FILE* fresult_node, *fresult_edge;
    fresult_node = fopen(result_file,"a+");
    fresult_edge = fopen(result_file_E,"a+");

    A = (sparse_csc *) SparseCore_read_matrix (fp, 1, &mtype, cc, fresult_node, fresult_edge, graph_id) ;

    if (mtype != SPARSE_CSC)
    {
        printf ("input matrix must be sparse\n") ;
        exit (1) ;
    }

    fclose(fresult_node);
    fclose(fresult_edge);
    // 但是没有写
    /**     矩阵对应图的边、结点信息在读矩阵的时候就已经写完了      **/

#ifndef write_graph
    // 重排序结果写在 result 中
    char *result = (char *)"./reorder_result/reorder.txt";
    FILE* fresult = fopen(result,"a+");
    // 先记录图标号
    fprintf(fresult,"%d\n", graph_id);
    fclose(fresult);

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

    TPSM_init(128, 2000, 3000, TPSM_NODE_AFFINITY );
    // 获取分块参数
    chunk_getSettings(32, 5000, 4, 4);
    Relaxfactor_setting (n, SparseCore_nnz (A, cc), RELAX_FOR_QR, cc);

    if (A->xtype == SPARSE_REAL)
    {
        SparseQR_factorization *QR ;
        
        
        QR = SparseQR (ordering, tol, A, cc, result) ;
        //QR = SparseQR (QR_ORDERING_BEST, tol, A, cc, str2) ;
        //QR = SparseQR (QR_ORDERING_GIVEN, tol, A, cc, str2) ;
        qr_symbolic *QRsym ;
        QRsym = QR->QRsym;
        Sparse_long *Qperm = QRsym->Qfill;
        
        fresult = fopen(result,"a+");
        for (int i = 0; i < n ; i++)
        {
            
            fprintf(fresult, "%d ",Qperm[i] );
        }
        fprintf(fresult,"\n");
        fclose(fresult);

        // 等待全部完成,销毁线程池
        TPSM_destroy(TPSM_SHUTDOWN_GENTLY);
        
        // free QR
        SparseQR_free (&QR, cc) ;
    }
    // -------------------------------------------------------------------------
    // free everything that remains
    // -------------------------------------------------------------------------
#endif
    SparseCore_free_sparse (&A, cc) ;
    SparseCore_finish (cc) ;
    fclose(fp);
    
    return (0) ;
}




