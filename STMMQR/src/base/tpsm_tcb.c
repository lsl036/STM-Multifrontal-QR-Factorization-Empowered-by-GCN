/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_tcb.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/


/******************************************************************************
 * INCLUDE
 ******************************************************************************/

#include "tpsm_tcb.h"

/******************************************************************************
 * EXTERN GLOBAL
 ******************************************************************************/

extern TPSM_TCB_t 	** 			tcbs;					//线程TCB(Thread Control Block)
extern int 						pool_total_size;		//线程的总线程数（就是多包含一个主线程）
/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/
/**
* @brief TPSM并行框架：初始化TCB(Thread Control Block)
* @return 成功:0
*/
int tpsm_initializeTcbs(void)
{	
	if(tcbs != NULL) return 0;
	
	tcbs  = TPSM_Malloc_Align(pool_total_size*sizeof(*tcbs));
	TPSM_assert(tcbs,1);

	return 0;	
}