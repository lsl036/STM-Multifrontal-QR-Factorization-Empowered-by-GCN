/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_buffer.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/

/******************************************************************************
 * INCLUDE
 ******************************************************************************/
#include "tpsm_buffer.h"

/******************************************************************************
 * EXTERN GLOBAL
 ******************************************************************************/

extern int 						tpsm_available;			//通过本地变量tpsm_available来检查线程池初始化是否完成
extern int  					pool_basic_size; 		//常驻线程数
extern int  					buffer_size;

extern int 					* 	workers_nbr_on_node;	//根据线程分配模式（集中或均分）与线程总数获得每个节点的线程数
extern int 					** 	workers_rank_on_node; 	//每个numa节点下，工作者线程的rank

extern TPSM_TASK_BUFFER_t 	** 	task_original_buffers;

extern TPSM_TCB_t 			** 	tcbs;

extern pthread_mutex_t 		* 	synchronization_locks;
extern int 					* 	synchronization_tasks;

extern int 					** 	search_node_sequences;

/******************************************************************************
 * DENOTE
 ******************************************************************************/
int private_findIdleWorkerOnNode(int const numa_node_rank);

int private_addTaskOnWorker
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag, 
	int const worker_rank
);

int private_addTaskOnOriginalBuffer
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag, 
	int const numa_node_rank
);

/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/

/**
* @brief TPSM并行框架：初始化任务缓冲队列(NUMA节点亲和，每个NUMA节点一个缓冲队列)
* @return 成功:0;失败:1
*/
int tpsm_initializeBuffer(void)
{	
	if(task_original_buffers!=NULL) return 1;
	if(workers_nbr_on_node==NULL) return 1;

	//开4个指针，每个都指向一个的task_buffer
	task_original_buffers = TPSM_Malloc_Align(TPSM_NUMANODES*sizeof(*task_original_buffers));
	TPSM_assert(task_original_buffers,1);

	//每个numa节点都开一个task_buffer,对每个结点的task_buffer进行初始化
	for(int i=0; i<TPSM_NUMANODES; ++i)
	{
		if(workers_nbr_on_node[i]>0)
		{
			task_original_buffers[i] = TPSM_Numa_MallocOnNode(sizeof(*task_original_buffers[i]),i);
			TPSM_assert(task_original_buffers[i],1);

			pthread_mutex_init(&task_original_buffers[i]->lock,NULL);
			pthread_cond_init(&task_original_buffers[i]->cond_signal,NULL);
			
			task_original_buffers[i]->size = buffer_size;
			task_original_buffers[i]->length = 0;
			task_original_buffers[i]->head = 0;
			task_original_buffers[i]->tail = task_original_buffers[i]->head;
			task_original_buffers[i]->numa_rank = i;
			
			task_original_buffers[i]->tasks = TPSM_Numa_MallocOnNode(buffer_size*sizeof(*task_original_buffers[i]->tasks),i);
			TPSM_assert(task_original_buffers[i]->tasks,1);
			
			task_original_buffers[i]->protected = TPSM_Numa_MallocOnNode(buffer_size*sizeof(*task_original_buffers[i]->protected),i);
			TPSM_assert(task_original_buffers[i]->protected,1);
		}	
	}

	return 0;
}



/**
* @brief TPSM并行框架：向指定NUMA节点中添加单个任务
* @param function_name 任务函数名称
* @param parameters 输入任务函数的参数
* @param synchronization_tag 同步阻塞表标识
* @param numa_node_rank NUMA节点序号
* @return 成功:0; 失败:1
*/
int TPSM_addTaskOnNode
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag, 
	int const numa_node_rank
)
{   
	#ifdef DEBUG
	int const worker_rank = TPSM_get_myrank();
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t开始\n",worker_rank, __FUNCTION__);
	#endif

	#ifdef PREVENTION
	if(numa_node_rank >= TPSM_NUMANODES) return 1;
	
	//判断该numa节点下是否有常驻线程
	if(workers_nbr_on_node[numa_node_rank]<1) return 1;

	if(!tpsm_available) return 1;
	#endif

	int checkcode;
	//*********一些临时声明**************************
	TPSM_TASK_BUFFER_t * taskqueue  = task_original_buffers[numa_node_rank];
	int idx; 
	
	//找对应numa_node_rank中的空闲的线程并预定线程
	int const active_rank = private_findIdleWorkerOnNode(numa_node_rank);
	//加入任务队列后，需要再唤醒一下所有的线程，以免有线程刚好沉睡
	int broadcast_rank;
	if(active_rank >= 0) //如果找到了,直接放到tcb中
	{	
		checkcode = private_addTaskOnWorker(function_name, parameters, synchronization_tag, active_rank);
		pthread_cond_signal(&tcbs[active_rank]->cond_signal);
	}
	else //如果没找到空闲的线程，那么就直接加入任务队列末尾
	{	
		checkcode = private_addTaskOnOriginalBuffer(function_name, parameters, synchronization_tag,numa_node_rank);
	}

	#ifdef DEBUG
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t结束\n",worker_rank, __FUNCTION__);
	#endif

	return checkcode;
}


/**
* @brief TPSM并行框架：往并行框架中添加单个任务
* @param function_name 任务函数名称
* @param parameters 输入任务函数的参数
* @param synchronization_tag 同步阻塞表标识
* @return 成功:0; 失败:1
*/
int TPSM_addTask
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag
)
{	
	#ifdef DEBUG
	int worker_rank = TPSM_get_myrank();
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t开始\n",worker_rank,__FUNCTION__);
	#endif
	
	#ifdef PREVENTION
	if(!tpsm_available) return 1;
	#endif

	int checkcode;
	
	int const my_nrk = TPSM_Numa_GetNodeRankSelf();//当前运行在哪个节点上
	int idx;
	int active_rank;

	//按节点查找次序(TPSM_SEARCH_NODE_SEQUENCE)来依次查找节点内是否有空闲线程，保证节点间距离由小到大
	for(int i=0; i<TPSM_NUMANODES; ++i)
	{	
		idx = search_node_sequences[my_nrk][i];
		if(workers_nbr_on_node[idx] > 0)
		{
			active_rank = private_findIdleWorkerOnNode(idx);//找看看有没有空闲线程,找到就会预定该线程
			if(active_rank >= 0)//如果找到了空闲线程
			{
				checkcode = private_addTaskOnWorker(function_name,parameters,synchronization_tag,active_rank);
				pthread_cond_signal(&tcbs[active_rank]->cond_signal);
				return checkcode;
			}
		}
	}

	#ifdef DEBUG
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t没有找到空闲线程，将添加到node=%d的任务队列中\n",worker_rank,__FUNCTION__,my_nrk);
	#endif

	//能出while循环证明没找到
	//那么就还是加入到my_nrk的numa节点对应的任务队列中
	checkcode =  TPSM_addTaskOnNode(function_name, parameters, synchronization_tag, my_nrk);

	#ifdef DEBUG
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t结束\n",worker_rank,__FUNCTION__,my_nrk);
	#endif

	return checkcode;
}

/**
* @brief TPSM并行框架：向指定线程添加单个任务
* @param function_name 任务函数名称
* @param parameters 输入任务函数的参数
* @param synchronization_tag 同步阻塞表标识
* @param worker_rank 工作者线程序号
* @return 成功:0; 失败:1
*/
int TPSM_addTaskOnWorker
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag, 
	int const worker_rank
)
{
	#ifdef PREVENTION
	if(worker_rank < 1 || worker_rank > pool_basic_size) return 1;
	#endif

	//拿预定的线程锁
	pthread_mutex_lock(&tcbs[worker_rank]->lock);
	
	//0表示tcb空闲，1表示线程被预定，但tcb中没有任务
	if(tcbs[worker_rank]->state!=TPSM_WORKER_LEISURE && tcbs[worker_rank]->state!=TPSM_WORKER_RESERVED)
	{	
		pthread_mutex_unlock(&tcbs[worker_rank]->lock);
		return 1;
	}
	
	//将任务信息直接写到其tcb中
	tcbs[worker_rank]->function_name = function_name;
	tcbs[worker_rank]->parameters = parameters; 
	tcbs[worker_rank]->synchronization_tag = synchronization_tag; 

	tcbs[worker_rank]->state = TPSM_WORKER_EXCUTABLE;//TPSM_WORKER_EXCUTABLE = 2 表示tcb中有未执行的任务

	//更新同步队列
	//TPSM_Atomic_addone(synchronization_tasks[synchronization_tag]);
	pthread_mutex_lock(&synchronization_locks[synchronization_tag]);
	synchronization_tasks[synchronization_tag]++;
	pthread_mutex_unlock(&synchronization_locks[synchronization_tag]);
	
	pthread_cond_signal(&tcbs[worker_rank]->cond_signal);
	
	pthread_mutex_unlock(&tcbs[worker_rank]->lock);

	return 0;
}

/******************************************************************************
 * PRIVATE FUNCTION
 ******************************************************************************/

//找到指定NUMA结点上的空闲线程并预定
int private_findIdleWorkerOnNode(int const numa_node_rank)
{	
	int active_rank;
	int const workers_nbr = workers_nbr_on_node[numa_node_rank];
	int * ptr = workers_rank_on_node[numa_node_rank];
	for(int i=0; i<workers_nbr; ++i)
	{
		active_rank = ptr[i];
		//printf("%d在private_findStandbyWorkerOnNode[%d]中找active_rank=%d---i=%d\n",TPSM_get_myrank(),numa_node_rank,active_rank,i);
		pthread_mutex_lock(&tcbs[active_rank]->lock);
		if(tcbs[active_rank]->state == TPSM_WORKER_LEISURE)
		{	
			tcbs[active_rank]->state = TPSM_WORKER_RESERVED; //状态TPSM_WORKER_RESERVED = 1表示已经预定
			pthread_mutex_unlock(&tcbs[active_rank]->lock);
			return active_rank;
		}
		pthread_mutex_unlock(&tcbs[active_rank]->lock);
	}
	return -1;
}

//将任务添加到worker的TCB中
int private_addTaskOnWorker
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag, 
	int const worker_rank
)
{	
	#ifdef PREVENTION
	if(worker_rank < 1 || worker_rank > pool_basic_size) return 1;
	#endif

	//拿预定的线程锁
	pthread_mutex_lock(&tcbs[worker_rank]->lock);
	
	//0表示tcb空闲，1表示线程被预定，但tcb中没有任务
	if(tcbs[worker_rank]->state!=TPSM_WORKER_LEISURE && tcbs[worker_rank]->state!=TPSM_WORKER_RESERVED)
	{	
		pthread_mutex_unlock(&tcbs[worker_rank]->lock);
		return 1;
	}
	
	//将任务信息直接写到其tcb中
	tcbs[worker_rank]->function_name = function_name;
	tcbs[worker_rank]->parameters = parameters; 
	tcbs[worker_rank]->synchronization_tag = synchronization_tag; 

	tcbs[worker_rank]->state = TPSM_WORKER_EXCUTABLE;//TPSM_WORKER_EXCUTABLE = 2 表示tcb中有未执行的任务

	//更新同步队列
	//TPSM_Atomic_addone(synchronization_tasks[synchronization_tag]);
	pthread_mutex_lock(&synchronization_locks[synchronization_tag]);
	synchronization_tasks[synchronization_tag]++;
	pthread_mutex_unlock(&synchronization_locks[synchronization_tag]);
	
	pthread_mutex_unlock(&tcbs[worker_rank]->lock);

	return 0;
}

//将任务添加到原始任务缓冲队列中
int private_addTaskOnOriginalBuffer
(	
	void* (*function_name)(void*), 
	void * const parameters, 
	int const synchronization_tag, 
	int const numa_node_rank
)
{	
	TPSM_TASK_BUFFER_t * const taskqueue  = task_original_buffers[numa_node_rank];
	//拿任务锁
	pthread_mutex_lock(&taskqueue->lock);
	
	int const idx = taskqueue->tail; //更新任务buffer之前保存索引
	
	//如果发现该处已被保护，则无法添加新任务
	if(taskqueue->protected[idx]) 
	{
		pthread_mutex_unlock(&taskqueue->lock);
		return 1;
	}

	//将任务信息添加到任务buffer中
	taskqueue->tasks[idx].function_name = function_name;
	taskqueue->tasks[idx].parameters = parameters; 
	taskqueue->tasks[idx].synchronization_tag = synchronization_tag; 
	
	//添加保护
	taskqueue->protected[idx] = 1;

	//更新任务buffer队列
	taskqueue->length++;
	taskqueue->tail++;
	taskqueue->tail = (taskqueue->tail == buffer_size) ? 0 : taskqueue->tail;
	
	//更新同步队列
	//TPSM_Atomic_addone(synchronization_tasks[synchronization_tag]);
	pthread_mutex_lock(&synchronization_locks[synchronization_tag]);
	synchronization_tasks[synchronization_tag]++;
	pthread_mutex_unlock(&synchronization_locks[synchronization_tag]);

	pthread_mutex_unlock(&taskqueue->lock);
	
	//唤醒该node的所有线程	
	// int broadcast_rank;
	// for(int i=0; i<workers_nbr_on_node[numa_node_rank]; ++i)
	// {	
	// 	broadcast_rank = workers_rank_on_node[numa_node_rank][i];
	// 	pthread_cond_signal(&tcbs[broadcast_rank]->cond_signal);
	// }

	return 0;
}


