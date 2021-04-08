/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_barrier.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

#ifndef TPSM_BARRIER_H
#define TPSM_BARRIER_H

/******************************************************************************
 * INCLUDES
 ******************************************************************************/
#include "tpsm_base.h"
#include "tpsm_sysinfo.h"
#include "tpsm_structure.h"
#include "tpsm_default.h"
#include "tpsm_threads.h"
#include "tpsm_auxiliary.h"
#include "pthread.h"

/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/

/**
* @brief TPSM并行框架：同步指定同步阻塞表标识的线程
* @param synchronization_tag 同步阻塞表标识
* @return 0
*/
int TPSM_barrier_tag(int const synchronization_tag);


#endif