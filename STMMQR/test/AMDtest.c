#include"Sparse.h"

int main (int argc, char* argv[])
{
    double one [2] = {1,0}, m1 [2] = {-1,0} ;
    sparse_csc *A;
    sparse_factor *L;
    sparse_common c;
    char *fmatrix = argv[1];
    int graph_id = atoi(argv[2]);

    printf("%s\n",fmatrix);
    SparseCore_start (&c) ;
    FILE *fp;
    fp = fopen(fmatrix,"r");
    int mtype;
    /******************************
     * 对输入的名称字符串进行处理 *
     ******************************/
    int len = strlen(fmatrix), step = 0;
    // 拼接目标文件的地址
    char result_name[100]="./reorder_result/reorder.txt";
    //char str2[50];
    // for (int i = 23; i < 23+(len - 28)/2; i++)
    // {
    //     str2[step++] = fmatrix[i];
    // }
    // str2[step++]='.';str2[step++]='t';str2[step++]='x';str2[step++]='t';
    // strcat(result_name, str2);
    // puts(result_name);
    char *result_file = (char *)"./data_clean/QR_Node.txt";
    char *result_file_E = (char *)"./data_clean/QR_Edge.txt";
    FILE* fresult_node, *fresult_edge;
    fresult_node = fopen(result_file,"a+");
    fresult_edge = fopen(result_file_E,"a+");

    /* 结果输出到文件中 */
    FILE* fresult;
    fresult = fopen(result_name,"a+");
    fprintf(fresult,"%d\n", graph_id);
    fclose(fresult);
    /******************************
     *     读  矩  阵  信  息     *
     ******************************/
    #ifdef SYM
    A = SparseCore_read_sparse(fp, &c) ;
    fprintf (fresult, "0\n");
    fclose(fresult);
    #else
    A = SparseCore_read_matrix(fp, 1, &mtype, &c, fresult_node, fresult_edge, graph_id);
    //fprintf (fresult, "1\n");
    //fclose(fresult);
    fclose(fresult_node);
    fclose(fresult_edge);
    #endif
    //SparseCore_print_sparse (A, "A", &c) ;		    /* print the matrix */

    int i;
    #ifdef SYM
    // 这里必须是对称矩阵
    if (A == NULL || A->stype == 0)		    /* A must be symmetric */
    {
        printf("Matrix must be symmetric !\n");
	    SparseCore_free_sparse (&A, &c) ;  
	    SparseCore_finish (&c) ;
	return (0) ;
    }
    #endif
    // [m n] = size (A) ;
    Int m = A->nrow ;
    Int n = A->ncol ;
    printf ("Matrix %6ld-by-%-6ld nnz: %6ld\n", m, n, SparseCore_nnz (A, &c)) ;

    c.final_asis = 1;                   
    c.supernodal = SPARSE_AUTO; 

    L = SparseChol_analyze (A, &c, result_name) ;           /* analyze */

    // for (int i = 0; i < 3; i++)
    // {
    //     printf("method %d lnz = %lf \n", i, c.method[i].lnz);
    // }
    
    // fresult = fopen(result_name,"a+");
    // fprintf(fresult,"%d \n\n", c.selected);
    // fclose(fresult);
    // Int *Lperm;
    // Lperm = L->Perm ;
    //printf("Perm array: \n");
    //fprintf(fresult, "Perm array: \n", fmatrix);
    // for (int i = 0; i < L->n ; i++)
    // {
    //     //printf( "%d ", Lperm[i] );
    //     fprintf(fresult, "%d ",Lperm[i] );
    // }
    // //printf("\n\n");
    // fprintf(fresult,"\n\n");
    /* free matrices */
    SparseCore_free_factor (&L, &c) ;		    
    SparseCore_free_sparse (&A, &c) ;
    fclose(fp);
    
    return 0;
}