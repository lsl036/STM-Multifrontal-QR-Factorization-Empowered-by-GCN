/******************************************************************************
 * VERSION: 2.0
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年9月24日
 * FILE:    tpsm_base.h
 * BRIEF:   基本函数：包括申请空间功能,判断功能
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


//numa亲和相关
//#include <unistd.h>
#include <numa.h> //numa_alloc_xxx
#include <sys/syscall.h>
#include<unistd.h>



/******************************************************************************
 * PUBLIC FUNCTIONS: API
 *****************************************************************************/

/**
* @brief 申请 'bytes' 个内存空间,并以32-bit对齐.
* @param bytes 申请空间的字节数.
* @return 返回空间指针.
*/
void * TPSM_malloc_align(size_t const bytes);

/**
* @brief 释放由HNU_malloc所开辟的空间
* @param ptr 指向要free的空间的指针
*/
void TPSM_align_free(void * ptr);

/**
* @brief 接受可变参数并输出信息后退出
* @param form 语句输出标准格式
*/
void TPSM_xerbla(char *form, ...);

/**
* @brief 在指定的numa节点上申请空间
* @param bytes 申请空间的字节数
* @param numa_node_rank numa节点序号
* @return 返回空间指针
*/
void *  TPSM_numalloc_onnode(size_t const bytes, size_t const numa_node_rank);

/**
* @brief 在当前的numa节点上申请空间
* @param bytes 申请空间的字节数
* @return 返回空间指针
*/
void *  TPSM_numalloc_local(size_t const bytes);

/**
* @brief 在释放由TPSM_numalloc_xxxx申请的空间
* @param ptr 指向要free的空间的指针
*/
void TPSM_numa_free(void * ptr, size_t const bytes);

/**
* @brief 获取系统可用核数
* @return size_t类型 系统可用核数
*/
size_t TPSM_get_syscores(void);

/**
* @brief 获取系统numa节点数
* @return size_t类型 系统numa节点数
*/
size_t TPSM_get_numanodes(void);

/**
* @brief 获取系统指定numa节点序号的核数
* @param numa_node_rank numa节点序号
* @return size_t类型 系统numa核数
*/
size_t TPSM_get_numacores(size_t numa_node_rank);

/******************************************************************************
 * DEFINE FUNCTION:
 *****************************************************************************/
#define TPSMMakeStr(_x) # _x

#define TPSM_assert(_x, _right) \
{\
   switch((_right)) \
   { \
      case 0: \
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


#endif