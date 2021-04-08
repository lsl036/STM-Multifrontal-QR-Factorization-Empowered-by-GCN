/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_synchronization.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本	
 ********************************************************************************************************************************/

/******************************************************************************
 * INCLUDE
 ******************************************************************************/
#include "tpsm_synchronization.h"

/******************************************************************************
 * EXTERN
 ******************************************************************************/
extern int 						synchronization_size;				//同步数组大小
extern pthread_mutex_t 	* 		synchronization_locks;   			//同步锁
extern pthread_cond_t 	* 		synchronization_cond_signals;  		//同步信号;
extern int 				* 		synchronization_tasks; 				//同步任务数组

/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/
/**
* @brief TPSM并行框架：初始化同步阻塞表
* @return 成功:0;失败:1
*/
int tpsm_initializeSynchronization(void)
{	
	if(synchronization_locks!=NULL) return 1;
	if(synchronization_cond_signals!=NULL) return 1;
	if(synchronization_tasks!=NULL) return 1;

	synchronization_tasks = TPSM_Malloc_Align(synchronization_size*sizeof(*synchronization_tasks));
	TPSM_assert(synchronization_tasks,1);
	
	synchronization_locks = TPSM_Malloc_Align(synchronization_size*sizeof(*synchronization_locks));
	TPSM_assert(synchronization_locks,1);
	
	synchronization_cond_signals = TPSM_Malloc_Align(synchronization_size*sizeof(*synchronization_cond_signals));
	TPSM_assert(synchronization_cond_signals,1);
	
	for(int i=0; i<synchronization_size; ++i)
	{
		pthread_mutex_init(&synchronization_locks[i],NULL);
		pthread_cond_init(&synchronization_cond_signals[i],NULL);
	}
	
	return 0;
}