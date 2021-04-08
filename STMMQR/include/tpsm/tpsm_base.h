/******************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_base.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 *****************************************************************************/

#ifndef TPSM_BASE_H
#define TPSM_BASE_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

/******************************************************************************
 * INCLUDES:
 *****************************************************************************/
#include <stdlib.h>	//malloc
#include <stddef.h> //size_t头文件
#include <stdio.h>
#include <errno.h>
#include <string.h> //memset
#include <stdarg.h>//主要目的为让函数能够接收可变参数
#include <sched.h>

//numa亲和相关
//#include <unistd.h>
#include <numa.h> //numa_alloc_xxx
#include <sys/syscall.h>
#include<unistd.h>


/******************************************************************************
 * PUBLIC FUNCTIONS: API
 *****************************************************************************/

/**
* @brief 接受可变参数并输出信息后退出
* @param form 语句输出标准格式
*/
void TPSM_xerbla(char *form, ...);

/**
* @brief 读环境变量值
* @param env 环境变量名称.
* @return 环境变量值(int类型)
*/
int TPSM_Read_Env(char *env);


/**
* @brief 申请 'bytes' 个内存空间,并以32-bit对齐.
* @param bytes 申请空间的字节数.
* @return 返回空间指针.
*/
void * TPSM_Malloc_Align(size_t const bytes);


/**
* @brief 释放由TPSM_Malloc_Align所开辟的空间
* @param ptr 指向要free的空间的指针
*/
void TPSM_Malloc_Free(void * ptr);

/**
* @brief 获得linux的LWP号
* @return 返回LWP号
*/
size_t TPSM_Get_Lwp(void);

/**
* @brief 获得调用该函数的进程所运行的cpu核序号
* @return cpu核序号
*/
size_t TPSM_Get_CpuRankSelf(void);

/**
* @brief 检测系统是否可以使用numa
* @return 0表示不可用，1表示可用
*/
int TPSM_Numa_Available(void);


/**
* @brief 获取系统可用核数
* @return size_t类型 系统可用核数
*/
size_t TPSM_Sys_GetCores(void);

/**
* @brief 获取系统numa节点数
* @return size_t类型 系统numa节点数
*/
size_t TPSM_Numa_GetNodes(void);

/**
* @brief 获取系统指定numa节点序号的核数
* @param numa_node_rank numa节点序号
* @return size_t类型 系统numa核数
*/
int TPSM_Numa_GetNodeCores(size_t const numa_node_rank);

/**
* @brief 获取系统NUMA结点的最大序号
* @return 得到最大numa结点序号
*/
size_t TPSM_Numa_GetMaxNodeRank(void);

/**
* @brief 获取CPU_rank 对应的 numa结点rank
* @return numa结点序号
*/
size_t TPSM_Numa_GetNodeRankOfCpu(size_t const cpu_rank);

/**
* @brief 获取调用该函数的进程所在numa节点序号
* @return numa结点序号
*/
size_t TPSM_Numa_GetNodeRankSelf(void);

/**
* @brief 获取numa节点开始的cpu号
* @param numa_node_rank numa节点序号
* @return cpu序号
*/
size_t TPSM_Numa_GetStartCpuRankOfNode(size_t const numa_node_rank);

/**
* @brief 在指定的numa节点上申请空间并初始化置0
* @param bytes 申请空间的字节数
* @param numa_node_rank numa节点序号
* @return 返回空间指针
*/
void *  TPSM_Numa_MallocOnNode(size_t const bytes, size_t const numa_node_rank);

/**
* @brief 在当前的numa节点上申请空间并初始化置0
* @param bytes 申请空间的字节数
* @return 返回空间指针
*/
void *  TPSM_Numa_MallocLocal(size_t const bytes);

/**
* @brief 重新分配空间
* @param old_addr 指向旧的空间指针
* @param old_size 旧空间大小
* @param new_size 新空间大小
*/
void* TPSM_Numa_Realloc(void *old_addr, size_t old_size, size_t new_size);


/**
* @brief 在全局开启NUMA NODE相互交织的内存空间
* @param bytes 申请空间的字节数
* @return 返回空间指针
*/
void *  TPSM_Numa_MallocInterleaved(size_t const bytes);


/**
* @brief 在释放由TPSM_NUMA_Mallocxxxx申请的空间
* @param ptr 指向要free的空间的指针
*  @param bytes 指向要free的空间的大小
*/
void TPSM_Numa_Free(void * ptr, size_t const bytes);


/**
* @brief 获取指定numa节点上与其他numa节点的距离排序
* @param numa_node_rank numa节点序号
* @return 一个具有对numa_node_rank亲和的数组
*/
int * TPSM_Numa_SearchNodeSequence(int const numa_node_rank);


/******************************************************************************
 * DEFINE FUNCTION:
 *****************************************************************************/
#define TPSM_TRUE (1)
#define TPSM_FALSE (0)

#define TPSM_SUCCESS (0)
#define TPSM_FAIL (1)

#define TPSMMakeStr(_x) # _x

#define TPSM_assert(_x, _right) \
{\
   switch((_right)) \
   { \
      case TPSM_SUCCESS: \
         if((_x)) \
         { \
             TPSM_xerbla("assertion %s FAILED, line %d of file %s\n", \
                  TPSMMakeStr(_x), __LINE__, __FILE__); \
         } \
         break; \
      case 1: \
         if(!(_x)) \
         { \
            TPSM_xerbla("assertion %s FAILED, line %d of file %s\n", \
                  TPSMMakeStr(_x), __LINE__, __FILE__); \
         } \
         break;\
      default: \
         TPSM_xerbla("assertion arg %s set WRONG, line %d of file %s\n", \
                  TPSMMakeStr(_right), __LINE__, __FILE__); \
         break;\
   } \
}



#define TPSM_Atomic_addone(_x) __sync_fetch_and_add(&(_x),1)
#define TPSM_Atomic_subone(_x) __sync_fetch_and_sub(&(_x),1)
#define TPSM_Atomic_set(_x,_y) __sync_lock_test_and_set(&(_x),(_y))

#endif