/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_buffer.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

/******************************************************************************
 * INCLUDE
 ******************************************************************************/
#include "tpsm_barrier.h"

/******************************************************************************
 * EXTERN GLOBAL
 ******************************************************************************/
extern int 							tpsm_available;
extern int 							available_original_workers;
extern int 					* 		available_workers_on_node;

extern int 							pool_basic_size;
extern int  						synchronization_size;
extern int 							pool_total_size;

extern pthread_mutex_t 				pool_lock; 						//线程池总锁
extern pthread_cond_t 				pool_cond;						//线程池条件信号

extern pthread_t 			* 		members_pid;
extern int 					* 		members_nrk;

extern TPSM_TASK_BUFFER_t 	** 		task_original_buffers;

extern pthread_mutex_t 		* 		synchronization_locks;   		//同步锁
extern pthread_cond_t 		* 		synchronization_cond_signals;  	//同步信号;
extern int 					* 		synchronization_tasks; 			//同步任务数组

/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/
/**
* @brief TPSM并行框架：同步指定同步阻塞表标识的线程
* @param synchronization_tag 同步阻塞表标识
* @return 0
*/
int TPSM_barrier_tag(int const synchronization_tag)
{	
	
	#ifdef PREVENTION
	if(!tpsm_available) return 1;
	if(synchronization_tag<0 || synchronization_tag>=synchronization_size) return 1;
	#endif

	pthread_mutex_lock(&synchronization_locks[synchronization_tag]);
	while(synchronization_tasks[synchronization_tag] != 0)
	{
		pthread_cond_wait(&synchronization_cond_signals[synchronization_tag], &synchronization_locks[synchronization_tag]);
	}
	pthread_mutex_unlock(&synchronization_locks[synchronization_tag]);

	#ifdef DEBUG
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t结束\n",worker_rank,__FUNCTION__);
    #endif

	return 0;
}

