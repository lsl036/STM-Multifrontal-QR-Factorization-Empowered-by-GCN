/******************************************************************************
 * VERSION: 2.0
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年9月24日
 * FILE:    tpsm_base.c
 * BRIEF:   基本函数：包括申请空间功能,判断功能
 *****************************************************************************/


/******************************************************************************
 * INCLUDES:
 *****************************************************************************/
#include "tpsm_base.h"

/******************************************************************************
 * FUNCTIONS:
 *****************************************************************************/

/**
* @brief 申请 'bytes' 个内存空间,并以32-bit对齐.
* @param bytes 申请空间的字节数.
* @return 返回空间指针.
*/
void * TPSM_malloc_align(size_t const bytes) 
{
    void * ptr;
    //以32bit对齐的方式开bytes个字节的内存空间
    int const success = posix_memalign(&ptr, 32, bytes);
    if(success != 0) {//success等于0就是成功，否则出错
        switch(success) {
        case ENOMEM:
            fprintf(stderr, "base: posix_memalign() returned ENOMEM. "
                            "Insufficient memory.\n");
            break;
        case EINVAL:
            fprintf(stderr, "base: posix_memalign() returned EINVAL. "
                            "Alignment must be power of two.\n");
            break;
        default:
            fprintf(stderr, "base: posix_memalign() returned '%d'.\n", success);
            break;
        }
        ptr = NULL;
    }else{
        memset(ptr,0,bytes);
    }
    return ptr;
}


/**
* @brief 释放由TPSM_malloc_align所开辟的空间
* @param ptr 指向要free的空间的指针
*/
void TPSM_align_free(void * ptr){
    free(ptr);
}


/**
* @brief 接受可变参数并输出信息后退出
* @param form 语句输出标准格式
*/
void TPSM_xerbla(char *form, ...)
{
   va_list argptr; //声明可接受变参变量

   va_start(argptr, form); //将form后的参数全部赋给argptr（不包括form）

   vfprintf(stderr, form, argptr);

   va_end(argptr); //注销这个参数

   exit(-1);
}




/**
* @brief 在指定的numa节点上申请空间
* @param bytes 申请空间的字节数
* @param numa_node_rank numa节点序号
* @return 返回空间指针
*/
void *  TPSM_numalloc_onnode(size_t const bytes, size_t const numa_node_rank)
{	
	if(numa_node_rank>numa_max_node()) return NULL;
	void * ptr = numa_alloc_onnode(bytes,numa_node_rank);
	if(ptr!=NULL) memset(ptr,0,bytes);
	return ptr;
}



/**
* @brief 在当前的numa节点上申请空间
* @param bytes 申请空间的字节数
* @return 返回空间指针
*/
void *  TPSM_numalloc_local(size_t const bytes)
{	
	void * ptr = numa_alloc_local(bytes);
	if(ptr!=NULL) memset(ptr,0,bytes);
	return ptr;
}

/**
* @brief 在释放由TPSM_numalloc_xxxx申请的空间
* @param ptr 指向要free的空间的指针
*  @param bytes 指向要free的空间的大小
*/
void TPSM_numa_free(void * ptr, size_t const bytes)
{
	numa_free(ptr,bytes);
}


/**
* @brief 获取系统可用核数
* @return size_t类型 系统可用核数
*/
size_t TPSM_get_syscores(void)
{	
	//sysconf(_SC_NPROCESSORS_ONLN)的返回值真正的代表了系统当前可用的核数
	return sysconf(_SC_NPROCESSORS_ONLN);
}



/**
* @brief 获取系统numa节点数
* @return size_t类型 系统numa节点数
*/
size_t TPSM_get_numanodes(void)
{	
	return 1+numa_max_node();
}


/**
* @brief 获取系统指定numa节点序号的核数
* @param numa_node_rank numa节点序号
* @return size_t类型 系统numa核数
*/
size_t TPSM_get_numacores(size_t numa_node_rank)
{	
	if(numa_node_rank>numa_max_node()) return -1;
	
	size_t const syscores = TPSM_get_syscores();
	size_t numa_cores = 0;
	
	for(int i = 0; i< syscores; ++i)
	{
		if(numa_node_of_cpu(i)==numa_node_rank)
		{
			numa_cores++;
		}
	}
	
	return numa_cores;
}


