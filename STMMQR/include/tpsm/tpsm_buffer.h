/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_buffer.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

#ifndef TPSM_BUFFER_H
#define TPSM_BUFFER_H


/******************************************************************************
 * INCLUDES
 ******************************************************************************/
#include "tpsm_base.h"
#include "tpsm_sysinfo.h"
#include "tpsm_structure.h"
#include "tpsm_default.h"
#include "pthread.h"

/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/

/**
* @brief TPSM并行框架：初始化任务缓冲队列(NUMA节点亲和，每个NUMA节点一个缓冲队列)
* @return 成功:0;失败:1
*/
int tpsm_initializeBuffer(void);


/**
* @brief TPSM并行框架：往并行框架中添加单个任务
* @param function_name 任务函数名称
* @param parameters 输入任务函数的参数
* @param synchronization_tag 同步阻塞表标识
* @return 成功:0; 失败:1
*/
int TPSM_addTask
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag
);


/**
* @brief TPSM并行框架：向NUMA节点中的添加单个任务
* @param function_name 任务函数名称
* @param parameters 输入任务函数的参数
* @param synchronization_tag 同步阻塞表标识
* @param numa_node_rank NUMA节点序号
* @return 成功:0; 失败:1
*/
int TPSM_addTaskOnNode
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag, 
	int const numa_node_rank
);


/**
* @brief TPSM并行框架：向指定线程添加单个任务
* @param function_name 任务函数名称
* @param parameters 输入任务函数的参数
* @param synchronization_tag 同步阻塞表标识
* @param worker_rank 工作者线程序号
* @return 成功:0; 失败:1
*/
int TPSM_addTaskOnWorker
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag, 
	int const worker_rank
);

#endif