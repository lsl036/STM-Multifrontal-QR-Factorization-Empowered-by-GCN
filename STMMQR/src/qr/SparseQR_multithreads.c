/******************************************************************************
 * VERSION: 1.1
 * DATE:    2020年9月27日
 * FILE:    SparseQR_multithreads.c
 * BRIEF:   并行函数
 * FUNCTION: multithreads函数，用于将ntasks加入线程池中并行执行，实现任务级并行
 *****************************************************************************/

#include"SparseQR.h"
#include "SparseQR_internal.h"
#include"tpsm.h"
#include<stdio.h>
#include<stdlib.h>
#ifdef PRINT_TIME
#include <sys/time.h>
#endif
#define NUMANODE_DEFAULT 0
/*********************************************
            需要并行的 任务
*********************************************/
void *task_execute(
    void* arg 
)
{   
    int i = 0;
    
    tasktodo *TASKROOT = (tasktodo *)arg;
    
    Long *TaskChildp = TASKROOT->Blob->QRsym->TaskChildp ;
    Long *TaskChild  = TASKROOT->Blob->QRsym->TaskChild ;
    Long pfirst = TaskChildp [TASKROOT->id] ;
    Long plast  = TaskChildp [TASKROOT->id+1] ;
    Long nchildren = plast - pfirst ;   //计算 子任务总数
    
    #ifdef NUMA_ALLOC
        Long s = TASKROOT->Blob->QRsym->TaskStack [TASKROOT->id]; // 使用 s 记一下task对应的栈号
        size_t numa_node_rank = TASKROOT->Blob->Work[s].numa_node_rank;
        //printf("stack = %d, numa_node_rank = %d, stack mod 4 = %d\n", s, numa_node_rank, s%4);
    #endif
	//printf("TPSM_MYRANK = %d", TPSM_get_myrank(tpool) );
    int checkcode;
    if (nchildren > 0)
    {
        
        // 创建子任务空间，每个子任务一个
        tasktodo *TASK_Child = malloc(nchildren*sizeof(*TASK_Child)) ;
        if (TASK_Child == NULL)
        {
            printf(" TASK %d's Child ERROR! \n", TASKROOT->id);
        }
        
        // 添加子任务
        for (Long i = 0; i < nchildren; i++)
        {
            Long child = TaskChild [pfirst + i] ;
            
            // 创建新的属性 传递给线程
            TASK_Child[i].id = child ;
            TASK_Child[i].Blob = TASKROOT->Blob ;
            // TASK_Child[i].tpool =  TASKROOT->tpool;
            TASK_Child[i].lock_id = TASKROOT->lock_id + 1; // 下一层了，锁的id + 1
            
            // 添加子任务，并行执行
            #ifdef NUMA_ALLOC
                checkcode = TPSM_addTaskOnNode( task_execute, &(TASK_Child[i]), TASKROOT->lock_id + 1, numa_node_rank) ;
            #else
                checkcode = TPSM_addTask( task_execute, &(TASK_Child[i]), TASKROOT->lock_id + 1) ;
            #endif
            if(checkcode!=0) printf("添加任务失败了，返回了非0值\n");
        }

        TPSM_barrier_tag( TASKROOT-> lock_id + 1) ;
    }
    
    //  孩子节点都做完了or 没有子节点，自己开始做任务
    #ifdef PRINT_TIME
    // double timeStart, timeEnd;
    // struct timeval tv;
    // gettimeofday(&tv, NULL);
    // timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
    #endif
    qr_kernel (TASKROOT->id, TASKROOT->Blob);
    #ifdef PRINT_TIME
    // gettimeofday(&tv, NULL);
    // timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
    // printf ("TASK %d 's kernel Function time: %lf\n\n", TASKROOT->id, timeEnd - timeStart);
    #endif

    return (NULL);
}

// 并行执行kernel运算函数
void qr_multithreads
(
    Long ntasks,
    qr_blob *Blob
)
{
    int checkcode;
    //printf("HERE WORKED POOL\n");

    tasktodo *TASK_ROOT = malloc(sizeof(*TASK_ROOT)) ;

    // 给根节点的结构体赋值
    TASK_ROOT->id = ntasks-1;
    TASK_ROOT->Blob = Blob;
    // TASK_ROOT->tpool = tpool;
    TASK_ROOT->lock_id = 1;
    // 创建一个根任务
    
    #ifdef NUMA_ALLOC
        Long s = Blob->QRsym->TaskStack [TASK_ROOT->id]; // 使用 s 记一下task对应的栈号
        size_t numa_node_rank = Blob->Work[s].numa_node_rank;
        //printf("stack = %d, numa_node_rank = %d, stack mod 4 = %d\n", s, numa_node_rank, s%4);
        checkcode = TPSM_addTaskOnNode(task_execute, TASK_ROOT, TASK_ROOT->lock_id, numa_node_rank);
    #else
        checkcode = TPSM_addTask(task_execute, TASK_ROOT, TASK_ROOT->lock_id);
	#endif
    if(checkcode!=0) printf("ROOT添加任务失败了，返回了非0值\n");
    TPSM_barrier_tag(TASK_ROOT->lock_id) ;
    
}
