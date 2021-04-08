#include "Sparse.h"
#include "SparseLU.h"
#include <sys/time.h>
#include <math.h>
#define ABS(x) ((x) >= 0 ? (x) : -(x))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif
#define Long Sparse_long


double check_error
(
    sparse_csc *A,
    Long Ap [ ],
    Long Ai [ ],
    double Ax [ ],
    void *Numeric,
    sparse_common *cc
)
{
    double *null = (double *) NULL ;
    int i, j, status = 0;
    double one [2] = {1,0}, minusone [2] = {-1,0}, zero [2] = {0,0} ;
    dense_array *X, *B, *X_sol; // 原始X，AX=B， 求解得到X_sol
    int n = A->ncol;
    //  初始化
    X = SparseCore_zeros (n, 1, A->xtype, cc) ;
    B = SparseCore_zeros (n, 1, A->xtype, cc) ;
    X_sol = SparseCore_zeros (n, 1, A->xtype, cc) ;
    // 初始化X = [0, 1, 2,..., n]
    double *x ;
    x = (double *) (X->x) ;
    for (i = 0 ; i < n ; i++)
    {
        x [i] = i ;
    }
    // 计算B = A*X
    SparseCore_sdmult (A, 0, one, zero, X, B, cc) ; 

    //求解
    double diff_norm = 0.0;
    double xx = 0.0;
    double *x_check = (double *) (X_sol->x);
    double *b = (double *) (B->x);
    status = sparselu_solve (SparseLU_A, Ap, Ai, Ax, x_check, b, Numeric, null, null) ;
    if (status < 0)
    {
        printf ("sparselu_solve failed, status = %d\n", status) ;
        return -1;
    }
    // 计算两个X的误差
    for (j = 0; j < X->nrow; j++)
    {
        xx =  x_check [j] - j;
        diff_norm += xx * xx;
    }
    diff_norm = sqrt(diff_norm)/X->nrow;

    SparseCore_free_dense (&X, cc) ;
    SparseCore_free_dense (&B, cc) ;
    SparseCore_free_dense (&X_sol, cc) ;
    return diff_norm;
}

int main (int argc, char* argv[])
{
    int cycleNum = 1;
    double timeStart, timeEnd, analyze_time, factor_time, rnorm;
    struct timeval tv;
    double one [2] = {1,0}, m1 [2] = {-1,0} ;
    int mtype ;
    sparse_csc *A ;
    Long n_row, n_col, *Ap, *Ai, status, i ;
    double *Ax;         /* size nz, numerical values */
    double *null = (double *) NULL ;
    void *Symbolic, *Numeric ;
    sparse_common ch ;
    SparseCore_start (&ch) ;

    //char* matrix_name = argv[1];
    
    /* 矩阵路径 */

    char *fmatrix = argv[1];
    
    printf("%s\n", fmatrix);
    FILE *fp;
    fp = fopen(fmatrix,"r");
    if(fp == NULL) {
        printf("%s file is not exist!\n", fmatrix);
        return 0;
    }
    // 读矩阵，强制转换为General形式
    A = (sparse_csc *) SparseCore_read_matrix (fp, 1, &mtype, &ch) ;
    /* 结果输出到文件中 */
    char *result_file = (char *)"./ARMLU.txt";
    FILE* fresult;
    fresult = fopen(result_file,"a+");
    fprintf(fresult, "%s\t", fmatrix);

    n_row = A->nrow;
    n_col = A->ncol;
    
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;
    if (!Ap || !Ai || !Ax) return (0) ;
    printf("A: %d x %d\n", n_row, n_col);
    
    gettimeofday(&tv, NULL);
    timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
    (void) sparselu_symbolic (n_row, n_col, Ap, Ai, Ax, &Symbolic, null, null) ;      // symbolic ordering and analysis
    gettimeofday(&tv, NULL);
    timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
    analyze_time = timeEnd-timeStart;

    printf("LU analyze   time: %f\n", analyze_time);
    fprintf(fresult, "%f\t", analyze_time);

    gettimeofday(&tv, NULL);
    timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
    (void) sparselu_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;    // factorization
    gettimeofday(&tv, NULL);
    timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
    factor_time = timeEnd-timeStart;

    printf("LU factor    time: %f\n", factor_time);
    printf("LU total run time: %f\n\n", analyze_time + factor_time);
    fprintf(fresult, "%lf\t", factor_time);
    fprintf(fresult, "%lf\t", analyze_time + factor_time);

    // solve
    if (n_row != n_col)
    {
        printf("ERROR : not a square matrix ! \n");
        fclose(fp);
        fclose(fresult);
        return 0;
    }
    double res = 0.;
    res = check_error(A, Ap, Ai, Ax, Numeric, &ch);
    printf ("res = %8.1e \n ", res) ;
    fprintf (fresult, "%8.1e\n ", res) ;

    sparselu_free_symbolic (&Symbolic) ;
    sparselu_free_numeric (&Numeric) ;
    SparseCore_free_sparse (&A, &ch) ;
    SparseCore_finish (&ch) ;
    fclose(fp);
    fclose(fresult);

    return (0) ;
}

