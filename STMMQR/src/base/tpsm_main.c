/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_main.c
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
 ********************************************************************************************************************************/


/******************************************************************************
 * INCLUDE
 ******************************************************************************/

#include "tpsm_main.h"

/******************************************************************************
 * GLOBAL 
 ******************************************************************************/

//通过本地变量tpsm_available来检查线程池初始化是否完成
int 						tpsm_available = 0;
int 				** 		search_node_sequences;

//可用工作者线程数
int 						available_original_workers = 0;	//最初可用工作者线程数
int 				* 		available_workers_on_node; 		//节点内可用工作者线程数

//线程池基本设置
int  						pool_basic_size; 				//常驻线程数
int  						buffer_size;					//任务缓冲队列大小
int  						synchronization_size;			//同步数组大小
int  						affinity_mode; 					//亲和模式

//线程的总线程数（就是多包含一个主线程）
int 						pool_total_size;				//pool_total_size = pool_basic_size + 1

//线程池关闭参数: 0不关，1立马关，2缓慢关
int 						pool_shutdown;

//线程池全局锁
pthread_mutex_t 			pool_lock; 						//线程池总锁
pthread_cond_t 				pool_cond_signal;				//线程池条件信号

//线程分布相关
int 						utilized_node_nbr; 				//根据线程分配模式（均分）与线程总数获得被利用的numa节点数
int 				* 		workers_nbr_on_node;			//根据线程分配模式（均分）与线程总数获得每个节点的线程数
int 				**	 	workers_rank_on_node; 			//每个numa节点下，工作者线程的rank


//原始任务缓冲队列指针数组
TPSM_TASK_BUFFER_t 	** 		task_original_buffers;

//线程TCB(Thread Control Block)
TPSM_TCB_t 			** 		tcbs;

//线程pthread_id
pthread_t 			* 		members_pid;

#ifdef DEBUG
//线程lwp_id
int 				* 		members_lid;
#endif

//线程对应的cpu_rank--仅在core亲和时一一对应
int 				* 		members_original_crk;

//线程的numa_node_rank
int 				* 		members_nrk;

//同步任务相关
pthread_mutex_t 	* 		synchronization_locks  __attribute__((aligned(TPSM_CACHELINE)));   //同步锁(chcheline对齐)
pthread_cond_t 		* 		synchronization_cond_signals  __attribute__((aligned(TPSM_CACHELINE)));  //同步信号(chcheline对齐)
int 				* 		synchronization_tasks __attribute__((aligned(TPSM_CACHELINE))); //同步任务数组(chcheline对齐)


/******************************************************************************
 * PRIVATE DENOTE
 ******************************************************************************/

int tpsm_getSettings
(	
	int const _pool_basic_size,
	int const _buffer_size,
	int const _synchronization_size,
	int const _affinity_mode
);

int tpsm_checkSettings(void);

/******************************************************************************
 * PUBLIC FUNCTION
 ******************************************************************************/

/**
* @brief TPSM并行框架：初始化
* @param _pool_basic_size 线程池常驻线程数
* @param _buffer_size 线程池原始任务缓存队列尺寸
* @param _synchronization_size 同步阻塞数组尺寸
* @param _affinity_mode 线程亲和模式
* @return 0
*/
int TPSM_init
(	
	int const _pool_basic_size,
	int const _buffer_size,
	int const _synchronization_size,
	int const _affinity_mode
)
{
	//*************如果线程池已经初始化了，那么直接快速返回************************
	if( tpsm_available ) return 0;

	//*************checkcode变量：用来接收返回值*********************************
	int checkcode;

	//*************获得线程池基本参数********************************************
	checkcode = tpsm_getSettings(_pool_basic_size, _buffer_size,\
					_synchronization_size, _affinity_mode);
	TPSM_assert(checkcode,0);

	//*************检查线程池基本参数是否合法*************************************
	checkcode = tpsm_checkSettings();
	TPSM_assert(checkcode,0);

	//*************创建numa节点搜索队列*************************************
	search_node_sequences = TPSM_Malloc_Align(TPSM_NUMANODES*sizeof(*search_node_sequences));
	TPSM_assert(search_node_sequences,1);
	for(int i=0;i<TPSM_NUMANODES;++i)
	{
		search_node_sequences[i] = TPSM_Numa_SearchNodeSequence(i);
	}

	//*************设置线程池基本的固定属性***************************************
	pool_total_size = pool_basic_size+1;
	pool_shutdown = TPSM_KEEP_RUNNING;

	//***********线程池全局锁与全局条件信号初始化**********************************
	pthread_mutex_init(&pool_lock,NULL);
	pthread_cond_init(&pool_cond_signal,NULL);

	//****************设置线程分布***********************************************
	checkcode = initializeThreadsDistribution();
	TPSM_assert(checkcode,0);
	
	//****************任务buffer初始化*******************************************
	checkcode = tpsm_initializeBuffer();
	TPSM_assert(checkcode,0);
	
	//****************同步数组初始化*********************************************
	checkcode = tpsm_initializeSynchronization();
	TPSM_assert(checkcode,0);
	
	//****************TCB设置****************************************************
	checkcode = tpsm_initializeTcbs();
	TPSM_assert(checkcode,0);

	//****************线程初始化*************************************************
	checkcode = tpsm_initializeThreads();
	TPSM_assert(checkcode,0);

	//****************主线程加入*************************************************
	members_pid[0] = pthread_self();

	#ifdef DEBUG
	members_lid[0] = TPSM_Get_Lwp();
	#endif
	
	members_original_crk[0] = TPSM_Get_CpuRankSelf();//得到当前主线程运行在哪个CPU核上
	members_nrk[0] = 0; //表示未定。主线程没有亲和，由操作系统调度

	//**************TPSM线程池完成初始化*****************************************
	tpsm_available = 1; 
	
	#ifdef DEBUG
	printf("TPSM DEBUG:\tTPSM初始化完成\n");
	#endif
	
	return 0;
}

/**
* @brief TPSM并行框架：销毁该框架
* @param shutdown_option 关闭选项
* @return 成功返回0,失败返回1
*/
int TPSM_destroy(int const shutdown_option)
{	
	//如果线程池没有初始化,也就不存在销毁线程池
	if(tpsm_available !=1 ) return 1;

	int checkcode;

    //对参数进行判断
    switch(shutdown_option)
    {
        case TPSM_KEEP_RUNNING:
            return 1;
        // case TPSM_SHUTDOWN_IMMEDIATELY:
        //     pool_shutdown = TPSM_SHUTDOWN_IMMEDIATELY;
        //     break;
        case TPSM_SHUTDOWN_GENTLY:
            pool_shutdown = TPSM_SHUTDOWN_GENTLY;
            break;
        default:
            pool_shutdown = TPSM_SHUTDOWN_GENTLY;
            break;
    }
	
	//唤醒所有线程
    for(int i=1; i<pool_total_size; ++i)
    {	
        pthread_cond_signal(&tcbs[i]->cond_signal);
    }


	//回收工作者线程
    for(int i=1; i<pool_total_size; ++i)
    {
        pthread_join(members_pid[i], NULL);
    }

	
	TPSM_Malloc_Free(members_pid);
	TPSM_Malloc_Free(members_original_crk);
	
	#ifdef DEBUG
	TPSM_Malloc_Free(members_lid);
	#endif

	TPSM_Malloc_Free(members_nrk);

	//线程都join后，都已释放自己的tcb，那么这时就可以再释放tcbs了
	TPSM_Malloc_Free(tcbs);

	for(int i=0;i<TPSM_NUMANODES; ++i)
	{
		TPSM_Numa_Free(workers_rank_on_node[i],workers_nbr_on_node[i]*sizeof(*workers_rank_on_node[i]));
	}
	TPSM_Malloc_Free(workers_rank_on_node);
	TPSM_Malloc_Free(workers_nbr_on_node);

	//释放任务队列资源
	for(int i=0;i<TPSM_NUMANODES; ++i)
	{	
		if(task_original_buffers[i]!=NULL)
		{	
			pthread_cond_destroy(&task_original_buffers[i]->cond_signal);
			pthread_mutex_destroy(&task_original_buffers[i]->lock);
			TPSM_Numa_Free(task_original_buffers[i]->tasks, buffer_size*sizeof(*(task_original_buffers[i]->tasks)));
			TPSM_Numa_Free(task_original_buffers[i]->protected, buffer_size*sizeof(*(task_original_buffers[i]->protected)));		
			TPSM_Numa_Free(task_original_buffers[i], sizeof(*task_original_buffers[i]));	
		}
		
	}

	TPSM_Malloc_Free(task_original_buffers);

	//释放同步相关资源
	for(int i=0; i<synchronization_size; ++i)
    {
        pthread_cond_destroy(&synchronization_cond_signals[i]);
        pthread_mutex_destroy(&synchronization_locks[i]);
    }
	TPSM_Malloc_Free(synchronization_tasks);
	TPSM_Malloc_Free(synchronization_locks);
	TPSM_Malloc_Free(synchronization_cond_signals);

	//销毁线程池锁
	pthread_mutex_destroy(&pool_lock);
	pthread_cond_destroy(&pool_cond_signal);
	

	//释放搜索队列空间
	for(int i=0;i<TPSM_NUMANODES;++i)
	{
		TPSM_Numa_Free(search_node_sequences[i],TPSM_NUMANODES*sizeof(*search_node_sequences[i])); 
	}
	TPSM_Malloc_Free(search_node_sequences);


	tpsm_available = 0;
	
	#ifdef DEBUG
	printf("TPSM DEBUG:\t当前线程(%d)\t函数:%s\t成功destory了线程池\n",TPSM_get_myrank(),__FUNCTION__);
	#endif

	return 0;
}



/******************************************************************************
 * PRIVATE FUNCTION
 ******************************************************************************/


/**
* @brief TPSM并行框架:环境变量获取
* @param _pool_basic_size 线程池常驻线程数
* @param _buffer_size 线程池原始任务缓存队列尺寸
* @param _synchronization_size 同步阻塞数组尺寸
* @param _affinity_mode 线程亲和模式
* @return 成功为0 失败为1
*/
int tpsm_getSettings
(	
	int const _pool_basic_size,
	int const _buffer_size,
	int const _synchronization_size,
	int const _affinity_mode
)
{	
	/*获得基本设定:
	*	环境变量输入关时(默认是关闭：env_on = 0)
	*		基本设置从函数输入设置
	*	环境变量输入开时(默认是关闭：env_on = 0)
	*		基本设置与函数输入无关，有环境变量输入则设置环境变量，否则设置默认宏定义
	*/

	//已环境变量为准： env_on = 0 表示不以环境变量为准
	int env_on = 0;

	//如果在环境变量中设置了env_on
	if( TPSM_Read_Env("TPSM_ENV_ON") ) 
		env_on = TPSM_Read_Env("TPSM_ENV_ON");
	
	switch(env_on)
	{	
		//如果以输入参数为准：则以输入参数为准
		case 0:
			pool_basic_size = _pool_basic_size;
			buffer_size = _buffer_size;
			synchronization_size = _synchronization_size;
			affinity_mode = _affinity_mode;
			break;

		//如果以环境变量为准：优先以环境变量为准，环境变量未设置，则以默认define为准
		case 1:
			//pool_basic_size设置
			if( TPSM_Read_Env("TPSM_BASIC_SIZE") ) 
				pool_basic_size = TPSM_Read_Env("TPSM_BASIC_SIZE");
			else
				pool_basic_size = TPSM_BASIC_DEFAULT_SIZE;
			
			//buffer_size设置
			if( TPSM_Read_Env("TPSM_BUFFER_SIZE") ) 
				buffer_size = TPSM_Read_Env("TPSM_BUFFER_SIZE");
			else
				buffer_size = TPSM_BUFFER_DEFAULT_SIZE;
			
			//synchronization_size设置
			if( TPSM_Read_Env("TPSM_SYNCHRONIZATION_SIZE") ) 
				synchronization_size = TPSM_Read_Env("TPSM_SYNCHRONIZATION_SIZE");
			else
				synchronization_size = TPSM_SYNCHRONIZATION_DEFAULT_SIZE;
			
			//affinity_mode设置
			if( TPSM_Read_Env("TPSM_AFFINITY_MODE") ) 
				affinity_mode = TPSM_Read_Env("TPSM_AFFINITY_MODE");
			else
				affinity_mode = TPSM_AFFINITY_DEFAULT_MODE;

			break;

		default:
			return 1;
	}

	return 0;
}


/**
* @brief TPSM并行框架:基本设置参数是否合法
* @return 成功为0 失败为正整数
*/
int tpsm_checkSettings(void)
{		
	
	//基本线程数必须要大于1
	if(pool_basic_size < 1)
	{
		fprintf(stderr, "WRONG : 线程池的常驻线程数设置错误，其值必须要大于等于1!\n");
		return 1;
	}

	//基本线程数必须能被numa结点数整除
	if(pool_basic_size%TPSM_NUMANODES != 0)
	{
		fprintf(stderr, "WRONG : 线程池的常驻线程数设置错误，其值必须是NUMA节点数的倍数!\n");
		return 2;
	}

	//任务缓冲区好歹也得来个10以上
	if(buffer_size < 10)
	{
		fprintf(stderr, "WRONG : 线程池的任务缓冲队列大小设置错误，其值好歹也来个10以上吧!\n");
		return 3;
	}
		
	//任务缓冲区好歹也得来个10以上
	if(synchronization_size < 10)
	{
		fprintf(stderr, "WRONG : 线程池的同步数组大小设置错误，其值好歹也来个10以上吧!\n");
		return 4;
	}
		
	if(affinity_mode != 0 && affinity_mode != 1)
	{
		fprintf(stderr, "WRONG : 线程池线程亲和模式设置错误，0(Numa Node亲和) or 1(Core亲和)!\n");
		return 5;
	}
		
	return 0;
}



