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
    
    TPSM_t *tpool = TASKROOT->tpool;
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
            TASK_Child[i].tpool =  TASKROOT->tpool;
            TASK_Child[i].lock_id = TASKROOT->lock_id + 1; // 下一层了，锁的id + 1
            
            // 添加子任务，并行执行
            checkcode = TPSM_addTask( tpool, task_execute, &(TASK_Child[i]), TASKROOT->lock_id + 1, 0) ;
			if(checkcode!=0) printf("添加任务失败了，返回了非0值\n");
        
        }

        TPSM_barrier( tpool, TASKROOT-> lock_id + 1) ;
    }
    
    //  孩子节点都做完了or 没有子节点，自己开始做任务

    qr_kernel (TASKROOT->id, TASKROOT->Blob);

    return (NULL);
}

// 并行执行kernel运算函数
void qr_multithreads
(
    Long ntasks,
    qr_blob *Blob,
    TPSM_t *tpool
)
{
    int checkcode;
    //printf("HERE WORKED POOL\n");

    tasktodo *TASK_ROOT = malloc(sizeof(*TASK_ROOT)) ;

    // 给根节点的结构体赋值
    TASK_ROOT->id = ntasks-1;
    TASK_ROOT->Blob = Blob;
    TASK_ROOT->tpool = tpool;
    TASK_ROOT->lock_id = 1;
    
    // 创建一个根任务

    checkcode = TPSM_addTask(tpool, task_execute, TASK_ROOT, TASK_ROOT->lock_id, 0);
	if(checkcode!=0) printf("ROOT添加任务失败了，返回了非0值\n");
    TPSM_barrier( tpool, TASK_ROOT->lock_id) ;
    
}
