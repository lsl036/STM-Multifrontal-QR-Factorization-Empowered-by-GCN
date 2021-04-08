/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_auxiliary.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

/******************************************************************************
 * INCLUDE
 ******************************************************************************/
#include "tpsm_auxiliary.h"

/*****************************************************************
 * 外部全局变量
 *****************************************************************/

//线程pthread_id
extern pthread_t 	* 	members_pid;

extern int  			pool_total_size;

extern int 				tpsm_available;

extern int 			* 	synchronization_tasks; //同步任务数组
/*****************************************************************
 * PUBLIC FUNCTION
 *****************************************************************/

/**
* @brief TPSM并行框架：获取线程序号
* @return 成功:0; 失败:1
*/
int TPSM_get_myrank(void)
{	

    pthread_t const my_id = pthread_self();
	//遍历线程id数组
	for(int i=0; i< pool_total_size; ++i)
	{
		if( pthread_equal(members_pid[i], my_id) ) return i;
	}
    return -1;
}

/**
* @brief TPSM并行框架：获取指定同步阻塞数表序号中的同步阻塞数表的值
* @param idx 同步阻塞数表序号
* @return 成功:0; 失败:1
*/
int TPSM_get_synValueOnTag(int const idx)
{	
    return synchronization_tasks[idx];
}