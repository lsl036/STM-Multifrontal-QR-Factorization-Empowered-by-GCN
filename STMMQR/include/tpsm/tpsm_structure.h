/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_structure.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

#ifndef TPSM_STRUCTURE_H
#define TPSM_STRUCTURE_H

#include "pthread.h"
/******************************************************************************
 * STRUCTURE
 ******************************************************************************/

//任务队列中的节点结构体定义
typedef struct
{
    void * (*function_name)(void*);
    void * parameters;
    int synchronization_tag;
}TPSM_TASK_NODE_t;

//任务队列结构体定义
typedef struct
{	
	pthread_mutex_t lock;
	pthread_cond_t cond_signal;
	TPSM_TASK_NODE_t * tasks;
	int size;
	int volatile length;
    int volatile head;
    int volatile tail;
    int * protected;
	int numa_rank;
}TPSM_TASK_BUFFER_t;

//TCB(Thread Control Block)结构体定义
typedef struct
{	
	void * (*function_name)(void*);
	void * parameters;
	int synchronization_tag;
	pthread_cond_t  cond_signal; //48 Byte
	pthread_mutex_t lock; //48 Byte
	int volatile state; //4 Byte
}TPSM_TCB_t;

#endif