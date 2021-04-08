/******************************************************************************
 * VERSION: 1.1
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年9月24日
 * FILE:    tpsm.h
 * BRIEF:   线程池
 * CHLOG:	v1.1:增加了子线程初始化后向主线程的通知,使得在TPSM_init函数结束后，所有子线程必定处于可工作状态
 *****************************************************************************/

#ifndef TPSM_H
#define TPSM_H

/**************************************************************
* INCLUDE
**************************************************************/
#include "tpsm_base.h"

#include <pthread.h>


/**************************************************************
* DEFINE
**************************************************************/

#define TASK_BUFFER_SIZE (10000)

#define SYNCHRONIZATION_MAX_SIZE (300)
/**************************************************************
* STRUCTURE
**************************************************************/
typedef enum
{   
    TPSM_KEEP_RUNNING = 0,
    TPSM_SHUTDOWN_IMMEDIATELY = 1,
    TPSM_SHUTDOWN_GENTLY = 2

}TPSM_SHUTDOWN_t;

typedef enum
{   
	TPSM_EVEN_DISTRIBUTION = 0,
	TPSM_CENTRALIZED_DISTRIBUTION = 1

}TPSM_CREATE_STRATEGY_t;

typedef struct
{	
	void * (*function_name)(void*);
	void * parameters;
	size_t synchronization_tag;
	pthread_cond_t  cond_signal; //48 Byte
	pthread_mutex_t lock; //48 Byte
	int state; //4 Byte
}TPSM_TCB_t;


typedef struct
{
    void * (*function_name)(void*);
    void * parameters;
    size_t synchronization_tag;
	size_t numa_rank;
}TPSM_TASK_NODE_t;


typedef struct
{
	TPSM_TASK_NODE_t tasks[TASK_BUFFER_SIZE];
	int length;
    int head;
    int tail;
    int protected[TASK_BUFFER_SIZE];
}TPSM_TASK_BUFFER_t;


//线程池总结构体
typedef struct
{
	//系统参数
	size_t  sys_cores; //系统总核心数
	size_t  numa_nodes; //系统numa节点数
	size_t  *  numa_cores; //numa节点的核心数
	cpu_set_t * cpu_sets; //CPU集
	
	//线程池参数
	size_t pool_max_threads; //线程池最大常驻线程数
	size_t pool_max_workers; //线程池最大常驻工作者线程数
	size_t numa_nodes_utilization;//线程节点利用数
	
	int shutdown;
	pthread_mutex_t task_lock; //48 Byte
	TPSM_TASK_BUFFER_t * task_buffer; //总任务buffer
	
	pthread_attr_t * thrd_attr_onnode;//对应每个numa节点线程属性，用来绑定线程
	pthread_attr_t quicker_attr; //不绑定核，全系统竞争
	
	size_t * thrd_nr;		//线程对应绑定的numa节点号
	pthread_t * thrd_ids;	//线程对应的pthread id
	size_t * thrd_lwp;  	//线程对应的LWP RANK
	size_t * numa_node_workers; //每个numa节点的工作者线程数
	size_t ** numa_node_worker_ranks;//开在manger上


	TPSM_TCB_t ** tcbs;

	pthread_mutex_t pool_lock; //48 Byte
	pthread_cond_t pool_cond_signal;
	int thread_available_nbr; //线程达到可使用的数目


	//同步任务相关
    pthread_mutex_t * synchronization_locks;   //同步锁
    pthread_cond_t * synchronization_cond_signals;  //同步信号;
    int * synchronization_tasks;

}TPSM_t;


typedef struct
{	
	TPSM_t * tpool;
	void * (*function_name)(void*);
	void * parameters;
	size_t synchronization_tag;
}TPSM_QUICKER_TCB_t;


/**************************************************************
* PUBLIC FUNTIONS
**************************************************************/
TPSM_t * TPSM_init(size_t const threads_nbr, int const strategy);
size_t TPSM_get_myrank(TPSM_t * tpool);
size_t TPSM_get_mylwp(TPSM_t * tpool);
int TPSM_addTask(TPSM_t * const tpool, void* (*function_name)(void*), void * const parameters, int const synchronization_tag, size_t const numarank);
int TPSM_barrier(TPSM_t * const tpool, int const synchronization_tag);
int TPSM_destroy(TPSM_t * const tpool, int const shutdown_option);

#endif