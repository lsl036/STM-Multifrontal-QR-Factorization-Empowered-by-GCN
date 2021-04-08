/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_threads.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

#ifndef TPSM_THREADS_H
#define TPSM_THREADS_H

/******************************************************************************
 * INCLUDES
 ******************************************************************************/
#include "tpsm_base.h"
#include "tpsm_sysinfo.h"
#include "tpsm_structure.h"
#include "tpsm_default.h"
#include "tpsm_auxiliary.h"
#include <unistd.h>
#include <sched.h>
#include "pthread.h"
/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/

/**
* @brief TPSM并行框架：初始化线程
* @return 成功:0
*/
int tpsm_initializeThreads(void);

#endif