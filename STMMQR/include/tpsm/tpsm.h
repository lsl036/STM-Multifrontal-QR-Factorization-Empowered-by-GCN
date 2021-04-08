/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

#ifndef TPSM_H
#define TPSM_H

/******************************************************************************
 * INCLUDE
 ******************************************************************************/

#include "tpsm_main.h"

/******************************************************************************
 * PUBLIC FUNCTION API
 ******************************************************************************/

/**
* @brief TPSM线程池初始化
* @param _pool_basic_size 线程池常驻线程数
* @param _buffer_size 线程池原始任务缓存队列尺寸
* @param _synchronization_size 同步阻塞数组尺寸
* @param _affinity_mode 线程亲和模式
* @return 0
*/
int TPSM_init
(	
	int const _pool_basic_size,
	int const _buffer_size,
	int const _synchronization_size,
	int const _affinity_mode
);


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


/**
* @brief TPSM并行框架：同步指定同步阻塞表标识的线程
* @param synchronization_tag 同步阻塞表标识
* @return 0
*/
int TPSM_barrier_tag(int const synchronization_tag);


/**
* @brief TPSM并行框架：销毁该框架
* @param shutdown_option 关闭选项
* @return 成功返回0,失败返回1
*/
int TPSM_destroy(int const shutdown_option);

#endif