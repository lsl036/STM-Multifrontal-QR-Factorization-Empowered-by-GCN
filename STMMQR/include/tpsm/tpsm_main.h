/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_main.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

#ifndef TPSM_MAIN_H
#define TPSM_MAIN_H

/******************************************************************************
 * INCLUDE
 ******************************************************************************/
#include "tpsm_base.h"
#include "tpsm_sysinfo.h"
#include "tpsm_default.h"
#include "tpsm_structure.h"

#include "tpsm_distribution.h"
#include "tpsm_buffer.h"
#include "tpsm_synchronization.h"
#include "tpsm_threads.h"
#include "tpsm_tcb.h"
#include "tpsm_barrier.h"
#include "tpsm_auxiliary.h"

#include "pthread.h"
/******************************************************************************
 * 默认基本设置
 ******************************************************************************/

//线程池默认最大线程数
#define TPSM_MAX_SIZE 								TPSM_SYSCORES

//线程池默认常驻线程数
#define TPSM_BASIC_DEFAULT_SIZE						(128)

//线程池NUMA默认的任务Buffer大小
#define TPSM_BUFFER_DEFAULT_SIZE					(1000)

//线程池同步数组大小
#define TPSM_SYNCHRONIZATION_DEFAULT_SIZE 			(1000)

//线程池默认线程分布模式
#define TPSM_DISTRIBUTE_DEFAULT_MODE				TPSM_EVEN_DISTRIBUTION			
 
//线程池默认线程亲和模式
#define TPSM_AFFINITY_DEFAULT_MODE					TPSM_NODE_AFFINITY

//线程池默认的任务缓冲队列模式
#define TPSM_BUFFER_DEFAULT_MODE					TPSM_NODE_BUFFER

//备用线程数默认值
#define TPSM_BACKUP_DEFAULT_SIZE					(TPSM_SYSCORES/TPSM_NUMANODES)
/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/ 

/**
* @brief TPSM并行框架：初始化
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
* @brief TPSM并行框架：销毁该框架
* @param shutdown_option 关闭选项
* @return 成功返回0,失败返回1
*/
int TPSM_destroy(int const shutdown_option);

#endif