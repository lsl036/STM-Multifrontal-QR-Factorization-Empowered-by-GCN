/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_default.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

#ifndef TPSM_DEFAULT_H
#define TPSM_DEFAULT_H

/******************************************************************************
 * 基本定义
 ******************************************************************************/

//线程池线程亲和模式定义
#define TPSM_NODE_AFFINITY					(0)
#define TPSM_CORE_AFFINITY					(1)

//线程池关闭选项定义
#define TPSM_KEEP_RUNNING 					(0)
// #define TPSM_SHUTDOWN_QUICKLY 				(1)
#define TPSM_SHUTDOWN_GENTLY 				(2)

//线程状态定义
#define TPSM_WORKER_LEISURE					(0) //空闲（挂起）
#define TPSM_WORKER_RESERVED				(1)	//被预定
#define TPSM_WORKER_EXCUTABLE  				(2)	//可执行TCB中的任务
#define TPSM_WORKER_INSPECT					(3)	//已执行完TCB中的任务，检查buffer中的任务

#endif