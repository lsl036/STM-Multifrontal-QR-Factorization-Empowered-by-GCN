/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_distribution.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

/******************************************************************************
 * INCLUDE
 ******************************************************************************/
#include "tpsm_distribution.h"

/******************************************************************************
 * GLOBAL
 ******************************************************************************/
extern int 		* 			available_workers_on_node;		//节点内可用工作者线程数

extern int  				pool_basic_size;				//常驻线程数

extern int 					utilized_node_nbr;
extern int 		* 			workers_nbr_on_node;
extern int 		** 			workers_rank_on_node; 			//每个numa节点下，线程的rank

/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/

/**
* @brief TPSM并行框架：初始化线程分布(平均分配线程到所有NUMA节点)
* @return 成功:0;失败:1
*/
int initializeThreadsDistribution(void)
{	
	//输入参数判断
	if(workers_nbr_on_node != NULL) return 1;
	if(workers_rank_on_node != NULL) return 1;

	//开空间
	workers_nbr_on_node = TPSM_Malloc_Align(TPSM_NUMANODES*sizeof(*workers_nbr_on_node));
	TPSM_assert(workers_nbr_on_node,1);

	//开空间
	available_workers_on_node = TPSM_Malloc_Align(TPSM_NUMANODES*sizeof(*available_workers_on_node));
	TPSM_assert(available_workers_on_node,1);

	int j = 0;
	int count = 0;
	for(int i=0;i<pool_basic_size;++i)
	{
		workers_nbr_on_node[j]++;
		available_workers_on_node[j]++;
		j++;
		j = (j==TPSM_NUMANODES) ? 0 : j;
	}
	utilized_node_nbr = TPSM_NUMANODES;

	workers_rank_on_node = TPSM_Malloc_Align(utilized_node_nbr*sizeof(*workers_rank_on_node));
	TPSM_assert(workers_rank_on_node,1);
	
	//0号线程是主线程，工作者线程号从1开始
	j=1;
	for(int i=0; i<utilized_node_nbr; ++i)
	{
		workers_rank_on_node[i] = TPSM_Numa_MallocOnNode(workers_nbr_on_node[i]*sizeof(*workers_rank_on_node[i]), i);
		TPSM_assert(workers_rank_on_node[i], 1);
		count = 0;
		while(count<workers_nbr_on_node[i])
		{
			workers_rank_on_node[i][count] = j;
			j++;
			count++;
		}
	}
	
	return 0;
}

