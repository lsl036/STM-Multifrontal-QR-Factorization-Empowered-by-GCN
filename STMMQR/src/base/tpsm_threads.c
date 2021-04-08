/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_threads.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

/******************************************************************************
 * INCLUDE
 ******************************************************************************/
#include "tpsm_threads.h"

/******************************************************************************
 * DEFINE
 ******************************************************************************/
#define TPSM_THREAD_STACK_SIZE (16UL*1024UL*1024UL)

/******************************************************************************
 * EXTERN
 ******************************************************************************/
extern int 							available_original_workers;		//最初可用工作者线程数
extern int 					* 		available_workers_on_node; 		//节点内可用工作者线程数

extern int 							pool_basic_size;				//常驻线程数
extern int 							pool_total_size;				//线程的总线程数（就是多包含一个主线程）

extern int 							buffer_size;					//任务缓冲队列大小
extern int 							affinity_mode;					//亲和模式

extern int 							pool_shutdown; 					//线程池关闭参数: 0不关，1立马关，2缓慢关

extern pthread_mutex_t 				pool_lock;						//线程池总锁
extern pthread_cond_t 				pool_cond_signal;				//线程池条件信号


extern int 					* 		workers_nbr_on_node;			//根据线程分配模式（集中或均分）与线程总数获得每个节点的线程数
extern int 					** 		workers_rank_on_node; 			//每个numa节点下，工作者线程的rank



extern pthread_t 			* 		members_pid;					//线程pthread_id

#ifdef DEBUG
extern int 					* 		members_lid;					//线程lwp_id
#endif

extern int 					* 		members_original_crk;			//线程对应的cpu_rank--仅在core亲和时一一对应
extern int 					* 		members_nrk;					//线程的numa_node_rank

extern TPSM_TASK_BUFFER_t 	** 		task_original_buffers;			//原始任务缓冲队列指针数组


extern TPSM_TCB_t 			** 		tcbs;							//线程TCB(Thread Control Block)

extern pthread_mutex_t 		* 		synchronization_locks;   		//同步锁
extern pthread_cond_t 		* 		synchronization_cond_signals;  	//同步信号;
extern int 					* 		synchronization_tasks; 			//同步任务数组

/******************************************************************************
 * DENOTE
 ******************************************************************************/

int initializeOriginalCpuRankAndNodeRank(void);

cpu_set_t *  private_initializeCpuSet(void);

pthread_attr_t *  private_initializePthreadAttr(cpu_set_t * cpusets);

int  private_initializeWorkers(void);

int  private_initializeBackups(void);

void * TPSM_worker(void * arg);

/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/

/**
* @brief TPSM并行框架：初始化线程
* @return 成功:0
*/
int  tpsm_initializeThreads(void)
{	
	int checkcode;

	checkcode = initializeOriginalCpuRankAndNodeRank();
	TPSM_assert(checkcode,0);

	checkcode = private_initializeWorkers();
	TPSM_assert(checkcode,0);

	return  0;
}



/***********************************************************
 * PRIVATE FUNCTION 
 ***********************************************************/
/**
* @brief TPSM并行框架：初始化CPU集（根据亲和模式）
* @return 成功:0; 失败:1
*/
cpu_set_t *  private_initializeCpuSet(void)
{	

	//输出参数(亲和模式，affinity_mode)的合法性判断
	if(affinity_mode!=0 && affinity_mode!=1) return NULL;

	cpu_set_t * cpusets = TPSM_Malloc_Align(TPSM_SYSCORES*sizeof(*cpusets));
	TPSM_assert(cpusets,1);
	
	//置零
	for(int i=0;i<TPSM_SYSCORES;++i)
	{
		CPU_ZERO(&cpusets[i]);
	}
		

	int start_idx;
	int current_node;
	int i;
	int j;
	//affinity_mode = 0: numa亲和；affinity_mode = 1: core亲和；
	switch(affinity_mode)
	{
		case TPSM_NODE_AFFINITY :
			start_idx = 0;
			current_node = TPSM_Numa_GetNodeRankOfCpu(0);
			for(i=0; i<TPSM_SYSCORES; ++i)
			{	
				if(TPSM_Numa_GetNodeRankOfCpu(i)!= current_node)
				{
					current_node = TPSM_Numa_GetNodeRankOfCpu(i);
					start_idx = i;
				}	
				for(j=start_idx; j < TPSM_SYSCORES; ++j)
				{	
					if(TPSM_Numa_GetNodeRankOfCpu(j)!=TPSM_Numa_GetNodeRankOfCpu(i)) 
						break;
					else
						CPU_SET(j, &cpusets[i]);
				}
			}
			break;
		case TPSM_CORE_AFFINITY :
			//赋值
			for(i=0; i<TPSM_SYSCORES; ++i)
			{	
				CPU_SET(i,&cpusets[i]);
			} 
			break;
		default:
			TPSM_Malloc_Free(cpusets);
			return NULL;
	}
	return cpusets;
}


/**
* @brief TPSM并行框架：根据线程数与分布模式设置线程对应的cpuset的idx
* @return 成功:0; 失败:1
*/
int initializeOriginalCpuRankAndNodeRank(void)
{
	
	if(workers_nbr_on_node == NULL) return 1;
	if(members_original_crk != NULL) return 1;
	if(members_nrk != NULL) return 1;
	//按线程池总线程数创建线程的cpu_rank记录表(仅用于记录最初的线程cpu核心的对应关系)
	members_original_crk = TPSM_Malloc_Align(pool_total_size*sizeof(*members_original_crk));
	TPSM_assert(members_original_crk, 1);
	
	members_nrk = TPSM_Malloc_Align(pool_total_size*sizeof(*members_nrk));
	TPSM_assert(members_nrk,1);
	
	int numa_node_rank = 0;//numa结点序号 0-4 (kunpeng 920下)
	int count = 0;//每个结点的线程设置数
	int cpu_core_rank = 0; //cpu核心的序号 0-127 (kunpeng 920下)
	//从1号工作者线程开始(0号为主线程)
	for(int i=1; i<pool_total_size; ++i)
	{	
		members_original_crk[i] = cpu_core_rank;
		members_nrk[i] = numa_node_rank;
		count++;

		//如果一个节点内的的线程数大于的了节点的cpu个数，那么就又回到起始位置来分配cpu
		if(count>=TPSM_SYSCORES/TPSM_NUMANODES)
		{
			cpu_core_rank = TPSM_Numa_GetStartCpuRankOfNode(numa_node_rank)+count%(TPSM_SYSCORES/TPSM_NUMANODES);
		}
		else
		{
			cpu_core_rank++;
		}
		
		if(count==workers_nbr_on_node[numa_node_rank])
		{
			count = 0;
			numa_node_rank++;
			cpu_core_rank = TPSM_Numa_GetStartCpuRankOfNode(numa_node_rank);
		}

	}
	
	return 0;
}

/**
* @brief TPSM并行框架：根据cpusets与members_original_crk设置线程属性attr
* @return 成功:0; 失败:1
*/
pthread_attr_t * private_initializePthreadAttr( cpu_set_t * cpusets )
{
	if(members_original_crk == NULL) return NULL;
	if(cpusets == NULL) return NULL;

	int checkcode;
	
	pthread_attr_t * attrs = TPSM_Malloc_Align( pool_total_size * sizeof(*attrs));
	TPSM_assert(attrs,1);

	int idx;
	//从1号工作者线程开始设置
	for(int i=1; i<pool_total_size; ++i)
	{	
		idx = members_original_crk[i];
		checkcode = pthread_attr_init(&attrs[i]);//初始化
        TPSM_assert(checkcode,0);
		checkcode = pthread_attr_setdetachstate(&attrs[i], PTHREAD_CREATE_JOINABLE); //可join
        TPSM_assert(checkcode,0);
		#ifdef SETSTACKSIZE
		checkcode = pthread_attr_setstacksize(&attrs[i],TPSM_THREAD_STACK_SIZE); //设置线程栈的大小
        TPSM_assert(checkcode,0);
		#endif
		checkcode = pthread_attr_setaffinity_np(&attrs[i], sizeof(cpusets[idx]), &cpusets[idx]); //绑定CPU集
		TPSM_assert(checkcode,0);
	}

	return  attrs;
}


/**
* @brief TPSM并行框架：初始化常规线程
* @return 成功:0; 失败:1
*/
int  private_initializeWorkers(void)
{	
	if(members_pid!=NULL) return 1;

	int checkcode;
	
	cpu_set_t * cpusets = private_initializeCpuSet();
	TPSM_assert(cpusets,1);

	pthread_attr_t * attrs =  private_initializePthreadAttr(cpusets);
	TPSM_assert(attrs,1);

	#ifdef DEBUG
	members_lid = TPSM_Malloc_Align(pool_total_size*sizeof(*members_lid));
	TPSM_assert(members_lid,1);
	#endif

	members_pid = TPSM_Malloc_Align(pool_total_size*sizeof(*members_pid));
	TPSM_assert(members_pid,1);

	for(int i=1;i<pool_total_size;++i)
	{
		checkcode = pthread_create(&members_pid[i], &attrs[i], TPSM_worker, NULL);
		TPSM_assert(checkcode,0);
	}

	//等待常规线程初始化完成
	pthread_mutex_lock(&pool_lock);
	while( available_original_workers < pool_basic_size )
	{
		pthread_cond_wait(&pool_cond_signal, &pool_lock);
	}
	pthread_mutex_unlock(&pool_lock);

	TPSM_Malloc_Free(cpusets);
	TPSM_Malloc_Free(attrs);

	return 0;
}


/**
* @brief TPSM并行框架：工作者线程
* @param arg 线程输入参数
* @return 成功:0; 失败:1
*/
void * TPSM_worker(void * arg)
{	
	int checkcode;

	//***获得自身rank号***********************************
	int const my_rank = TPSM_get_myrank();
	
	//***获得自身的所在的numa节点的序号(nrk, my node rank)********
	int const my_nrk = members_nrk[my_rank];

	#ifdef DEBUG
	//***获得自身的lwp号(lid,lwp id)并写入全局数组**********************
	int const my_lid = TPSM_Get_Lwp();
	members_lid[my_rank] = my_lid;
	#endif
	
	//**创建自己的本地tcb并初始化**************************

	TPSM_TCB_t * const my_tcb_p = TPSM_Numa_MallocLocal(sizeof(*my_tcb_p));
	TPSM_assert(my_tcb_p,1);
	
	//初始化本地tcb中的内容
	pthread_mutex_init(&my_tcb_p->lock,NULL); //本地TCB的锁初始化
	pthread_cond_init(&my_tcb_p->cond_signal,NULL); //本地TCB的条件信号量初始化

	my_tcb_p->state = TPSM_WORKER_LEISURE; //线程状态初始化
	tcbs[my_rank] =  my_tcb_p; //将创建好的tcb指针存入全局tcbs数组中

	//**创建一个指向任务缓冲队列的指针**************************
	TPSM_TASK_BUFFER_t * const  my_taskqueue = task_original_buffers[my_nrk];

	//可用线程数增加，并向发送信号（为了同步到TPSM_init）
	pthread_mutex_lock(&pool_lock);
	available_original_workers++;
	if( available_original_workers == pool_basic_size )
		pthread_cond_signal(&pool_cond_signal);
	pthread_mutex_unlock(&pool_lock);

	#ifdef DEBUG
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t工作者线程初始化完成。运行在cpu=%d\n",my_rank,__FUNCTION__,TPSM_Get_CpuRankSelf());
	#endif

	//一些可以保存到本地的常变量
	int register state;  //本地状态
	int register queidx; //任务队列索引
	int register syntag; //同步数组索引
	while(1) //FSM有限状态机
	{	
		// printf("%d pool_shutdown =%d \n",my_rank, pool_shutdown);
		pthread_mutex_lock(&my_tcb_p->lock);
		my_tcb_p->state = state;
		while( (state = my_tcb_p->state)==TPSM_WORKER_LEISURE && pool_shutdown == TPSM_KEEP_RUNNING )
		{   
			pthread_cond_wait(&my_tcb_p->cond_signal, &my_tcb_p->lock);
		}
		pthread_mutex_unlock(&my_tcb_p->lock);

		// if(pool_shutdown == TPSM_SHUTDOWN_IMMEDIATELY)
        // {	
        //     break;
        // }

        if(pool_shutdown == TPSM_SHUTDOWN_GENTLY && my_taskqueue->length==0 )
        {	
            break;
        }

		switch(state)
		{
			case TPSM_WORKER_EXCUTABLE:
				#ifdef DEBUG
				printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t将执行TCB中的任务\n",my_rank,__FUNCTION__);
				#endif

				syntag = my_tcb_p->synchronization_tag;

				(my_tcb_p->function_name)(my_tcb_p->parameters);
								
				#ifdef DEBUG
				printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t已执行TCB中的任务\n",my_rank,__FUNCTION__);
				#endif

				//更新同步数组
				pthread_mutex_lock(&synchronization_locks[syntag]);
				synchronization_tasks[syntag]--;
				if(synchronization_tasks[syntag]==0)
				{	
					pthread_cond_broadcast(&synchronization_cond_signals[syntag]);
				}
				pthread_mutex_unlock(&synchronization_locks[syntag]);
				
				//修改状态
				state = TPSM_WORKER_INSPECT;

			case TPSM_WORKER_INSPECT:
				queidx = -1;
				pthread_mutex_lock(&my_taskqueue->lock);
				if(my_taskqueue->length > 0)
				{	
					//取信息
					queidx = my_taskqueue->head;
					syntag = my_taskqueue->tasks[queidx].synchronization_tag;
					
					//更新队列
					my_taskqueue->length--;
					my_taskqueue->head++;
					my_taskqueue->head = (my_taskqueue->head == buffer_size) ? 0 : my_taskqueue->head;
				}
				pthread_mutex_unlock(&my_taskqueue->lock);
			
				//如果取到了任务 则 执行任务
				if(queidx >= 0)
				{	

					(my_taskqueue->tasks[queidx].function_name)(my_taskqueue->tasks[queidx].parameters);
				
					//任务执行完后，更新proctected
					//TPSM_Atomic_set(my_taskqueue->protected[queidx],0);
					pthread_mutex_lock(&my_taskqueue->lock);
					my_taskqueue->protected[queidx]=0;
					pthread_mutex_unlock(&my_taskqueue->lock);
					
					//更新同步数组
					pthread_mutex_lock(&synchronization_locks[syntag]);
					synchronization_tasks[syntag]--;
					if(synchronization_tasks[syntag]==0)
					{	
						pthread_cond_broadcast(&synchronization_cond_signals[syntag]);
					}
					pthread_mutex_unlock(&synchronization_locks[syntag]);		
				}
				else
				{
					state = TPSM_WORKER_LEISURE;
				}
				
				break;

			default:
				state = TPSM_WORKER_LEISURE;
		}

	}

	//free自己的TCB
	pthread_mutex_destroy(&my_tcb_p->lock); //锁销毁
	pthread_cond_destroy(&my_tcb_p->cond_signal); //条件信号量销毁
	
	TPSM_Numa_Free(my_tcb_p, sizeof(*my_tcb_p));


	#ifdef DEBUG
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t退出线程\n", my_rank,__FUNCTION__);
	#endif
	
	//退出线程
	pthread_exit(NULL);

	return ((void*)0);
}
