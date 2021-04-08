/******************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_base.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 *****************************************************************************/


/******************************************************************************
 * INCLUDES:
 *****************************************************************************/
#include "tpsm_base.h"

/******************************************************************************
 * FUNCTIONS:
 *****************************************************************************/

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
* @brief 读环境变量值
* @param env 环境变量名称.
* @return 环境变量值(int类型)
*/
int TPSM_Read_Env(char *env)
{
  char *p;
  //getenv函数是gcc自带的，再stdlib.h中
  if (( p = getenv(env)))
  	return (atoi(p));
  else
	return(0);
}

/**
* @brief 申请 'bytes' 个内存空间,并以32-bit对齐.
* @param bytes 申请空间的字节数.
* @return 返回空间指针.
*/
void * TPSM_Malloc_Align(size_t const bytes) 
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
* @brief 释放由TPSM_Malloc_Align所开辟的空间
* @param ptr 指向要free的空间的指针
*/
void TPSM_Malloc_Free(void * ptr){
    free(ptr);
}


/**
* @brief 获得linux的LWP号
* @return 返回LWP号
*/
size_t TPSM_Get_Lwp(void)
{
    return syscall(SYS_gettid);
}



/**
* @brief 获得调用该函数的进程所运行的cpu核序号
* @return cpu核序号
*/
size_t TPSM_Get_CpuRankSelf(void)
{
    return sched_getcpu();
}


/**
* @brief 检测系统是否可以使用numa
* @return 0表示不可用，1表示可用
*/
int TPSM_Numa_Available(void)
{
	if(numa_available()<0) return TPSM_FALSE;
	else return TPSM_TRUE;
}


/**
* @brief 获取系统可用核数
* @return size_t类型 系统可用核数
*/
size_t TPSM_Sys_GetCores(void)
{	
	//sysconf(_SC_NPROCESSORS_ONLN)的返回值真正的代表了系统当前可用的核数
	return sysconf(_SC_NPROCESSORS_ONLN);
}



/**
* @brief 获取系统numa节点数
* @return size_t类型 系统numa节点数
*/
size_t TPSM_Numa_GetNodes(void)
{	
	return 1+numa_max_node();
}


/**
* @brief 获取系统指定numa节点序号的核数
* @param numa_node_rank numa节点序号
* @return size_t类型 系统numa核数
*/
int TPSM_Numa_GetNodeCores(size_t const numa_node_rank)
{	
	if(numa_node_rank>numa_max_node()) return -1;
	
	size_t const syscores = TPSM_Sys_GetCores();
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


/**
* @brief 获取系统NUMA结点的最大序号
* @return 得到最大numa结点序号
*/
size_t TPSM_Numa_GetMaxNodeRank(void)
{
	return numa_max_node();
}


/**
* @brief 获取CPU_rank 对应的 numa结点rank
* @return numa结点序号
*/
size_t TPSM_Numa_GetNodeRankOfCpu(size_t const cpu_rank)
{
	return numa_node_of_cpu(cpu_rank);
}


/**
* @brief 获取调用该函数的进程所在numa节点序号
* @return numa结点序号
*/
size_t TPSM_Numa_GetNodeRankSelf(void)
{
	return TPSM_Numa_GetNodeRankOfCpu(TPSM_Get_CpuRankSelf());
}

/**
* @brief 获取numa节点开始的cpu号
* @param numa_node_rank numa节点序号
* @return cpu序号
*/
size_t TPSM_Numa_GetStartCpuRankOfNode(size_t const numa_node_rank)
{	
	//检查输入合法性
	if(numa_node_rank>TPSM_Numa_GetMaxNodeRank())
		return -1;
	
	for(int i=0;i<TPSM_Sys_GetCores();++i)
	{
		if(numa_node_rank==numa_node_of_cpu(i))
			return i;
	}

	return -1;
}

/**
* @brief 在指定的numa节点上申请空间并初始化置0
* @param bytes 申请空间的字节数
* @param numa_node_rank numa节点序号
* @return 返回空间指针
*/
void *  TPSM_Numa_MallocOnNode(size_t const bytes, size_t const numa_node_rank)
{	
	if(numa_node_rank>numa_max_node()) return NULL;
	void * ptr = numa_alloc_onnode(bytes,numa_node_rank);
	if(ptr!=NULL) memset(ptr,0,bytes);
	return ptr;
}



/**
* @brief 在当前的numa节点上申请空间并初始化置0
* @param bytes 申请空间的字节数
* @return 返回空间指针
*/
void *  TPSM_Numa_MallocLocal(size_t const bytes)
{	
	void * ptr = numa_alloc_local(bytes);
	if(ptr!=NULL) memset(ptr,0,bytes);
	return ptr;
}



/**
* @brief 重新分配空间
* @param old_addr 指向旧的空间指针
* @param old_size 旧空间大小
* @param new_size 新空间大小
*/
void * TPSM_Numa_Realloc(void *old_addr, size_t old_size, size_t new_size)
{
	numa_realloc(old_addr,old_size,new_size);
}

/**
* @brief 在全局开启NUMA NODE相互交织的内存空间
* @param bytes 申请空间的字节数
* @return 返回空间指针
*/
void *  TPSM_Numa_MallocInterleaved(size_t const bytes)
{
	return numa_alloc_interleaved(bytes);
}

/**
* @brief 在释放由TPSM_NUMA_Mallocxxxx申请的空间
* @param ptr 指向要free的空间的指针
*  @param bytes 指向要free的空间的大小
*/
void TPSM_Numa_Free(void * ptr, size_t const bytes)
{
	numa_free(ptr,bytes);
}

/**
* @brief 获取指定numa节点上与其他numa节点的距离排序
* @param numa_node_rank numa节点序号
* @return 一个具有对numa_node_rank亲和的数组
*/
int * TPSM_Numa_SearchNodeSequence(int const numa_node_rank)
{	
	int const nodes_nbr = TPSM_Numa_GetNodes();
	
	int * sequence = TPSM_Numa_MallocOnNode( nodes_nbr*sizeof(*sequence), numa_node_rank);
	memset(sequence,0,nodes_nbr*sizeof(*sequence));

	//桶号表示：目标节点序号； 值表示：目标节点与numa_node_rank间的距离
	int * distances = TPSM_Malloc_Align(nodes_nbr*sizeof(*distances));
	memset(distances,0,nodes_nbr*sizeof(*distances));

	for(int i = 0;i<nodes_nbr;++i)
	{
		distances[i] = numa_distance(numa_node_rank,i);
	}

	int * protected = TPSM_Malloc_Align(nodes_nbr*sizeof(*protected));
	memset(protected,0,nodes_nbr*sizeof(*protected));
	
	//从小到大排序
	int max_distance;
	int max_distance_idx;
	for(int i=nodes_nbr-1; i>=0; --i)
	{	
		max_distance = -1;
		max_distance_idx = -1;
		for(int j=0;j<nodes_nbr;++j)
		{
			if(!protected[j])
			{
				if(distances[j] > max_distance )
				{
					max_distance = distances[j];
					max_distance_idx = j;
				}
			}
		}
		sequence[i] = max_distance_idx;
		protected[max_distance_idx] = 1;
	}

	TPSM_Malloc_Free(distances);
	TPSM_Malloc_Free(protected);

	// for(int i=0; i<nodes_nbr; ++i)
	// {
	// 	printf("numa_node_rank = %d :: sequence[%d] = %d \n",numa_node_rank,i,sequence[i]);
	// }
	return  sequence;
}