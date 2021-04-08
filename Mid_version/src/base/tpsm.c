/******************************************************************************
 * VERSION: 1.1
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年9月24日
 * FILE:    tpsm.h
 * BRIEF:   线程池
 * CHLOG:	v1.1:增加了子线程初始化后向主线程的通知,使得在TPSM_init函数结束后，所有子线程必定处于可工作状态
 *****************************************************************************/

/**************************************************************
* INCLUDE
**************************************************************/
#include "tpsm.h"

/**************************************************************
* DENOTES
**************************************************************/
int TPSM_initialize_numacores(TPSM_t * tpool);
int TPSM_initialize_cpusets(TPSM_t * tpool);
int TPSM_initialize_taskbuffer(TPSM_t * tpool);
int TPSM_initialize_pthreadattr(TPSM_t * tpool);

size_t TPSM_numa_nodes_utilization(size_t const threads_nbr, size_t const numa_nodes, size_t  * const numa_cores, int const strategy);
int TPSM_initialize_distribution(TPSM_t * tpool, int const strategy);
int TPSM_initialize_threads(TPSM_t * tpool, int const strategy);

int TPSM_synchronization_init(TPSM_t * tpool);

void TPSM_manager(void * arg);
void TPSM_worker(void * arg);
void TPSM_quicker(void * arg);
int TPSM_create_quicker
(
	TPSM_t * tpool,
	void * (*function_name)(void*),
    void * parameters,
	size_t synchronization_tag
);
int TPSM_free(TPSM_t * const tpool);
/**************************************************************
* PUBLIC FUNCTIONS
**************************************************************/
TPSM_t * TPSM_init(size_t const threads_nbr, int const strategy)
{
	int checkcode;

	//*******创建线程池结构体********
	TPSM_t * tpool = TPSM_malloc_align(sizeof(*tpool));
	TPSM_assert(tpool,1);
	
	//*******初始化关闭参数********
	tpool->shutdown = TPSM_KEEP_RUNNING;

	//*******获取CPU系统参数*********
	tpool->sys_cores = TPSM_get_syscores();
	tpool->numa_nodes = TPSM_get_numanodes();

	//*******初始化numa节点核数的数据存入TPSM中*********
	checkcode = TPSM_initialize_numacores(tpool);
	TPSM_assert(checkcode,0);

	//*******初始化CPU集*********
	checkcode = TPSM_initialize_cpusets(tpool);
	TPSM_assert(checkcode,0);
	
	//******获取线程池最大线程数*********
	tpool->pool_max_workers = threads_nbr;
	tpool->pool_max_threads = tpool->pool_max_workers+1;

	//******初始化线程池锁*********
	checkcode = pthread_mutex_init(&tpool->pool_lock,NULL);
	TPSM_assert(checkcode,0);

	//******初始化线程条件信号量*********
	checkcode = pthread_cond_init(&tpool->pool_cond_signal,NULL);
	TPSM_assert(checkcode,0);

	tpool->thread_available_nbr = 0;

	//*****创建总任务buffer*********
	checkcode = TPSM_initialize_taskbuffer(tpool);
	TPSM_assert(checkcode,0);

	//创建同步
	checkcode = TPSM_synchronization_init(tpool);
	TPSM_assert(checkcode,0);

	//*******初始化线程池中的pthread线程属性，用来绑定线程*******
	checkcode = TPSM_initialize_pthreadattr(tpool);
	TPSM_assert(checkcode,0);
	
	//*******线程池线程初始化*******
	checkcode = TPSM_initialize_threads(tpool, strategy);
	TPSM_assert(checkcode,0);


	//********等待所有线程达到可用状态*******
	pthread_mutex_lock(&tpool->pool_lock);
	while(tpool->thread_available_nbr!=tpool->pool_max_threads)
	{
		pthread_cond_wait(&tpool->pool_cond_signal, &tpool->pool_lock);
	}
	pthread_mutex_unlock(&tpool->pool_lock);

	return tpool;
}

size_t TPSM_get_myrank(TPSM_t * tpool)
{
    pthread_t const my_id = pthread_self();
	//遍历线程id数组
	for(size_t i=0; i<tpool->pool_max_threads; ++i)
	{
		if( pthread_equal(tpool->thrd_ids[i], my_id) ) return i;
	}
    return -1;
}

size_t TPSM_get_mylwp(TPSM_t * tpool)
{
    return syscall(SYS_gettid);
}



int TPSM_addTask(TPSM_t * tpool, void* (*function_name)(void*), void * const parameters, int const synchronization_tag, size_t const numa_rank)
{   
	if(tpool==NULL) return -1;

    if(tpool->shutdown != TPSM_KEEP_RUNNING )
    {	

        return 1;
    }

	TPSM_TASK_BUFFER_t *  taskqueue  = tpool->task_buffer;

    // 要修改tpool，必须进入临界区，则必须先取得互斥锁
    pthread_mutex_lock(&(tpool->task_lock));
    
    //检测任务队列是否已满
    if(taskqueue->length == TASK_BUFFER_SIZE){
        pthread_mutex_unlock(&(tpool->task_lock));
        return 1;
    }
	
    //在tail处加入任务
    int idx = taskqueue->tail;
    //在添加任务前查看tail处的任务是否是用过的
    if(taskqueue->protected[idx]==1)//如果任务不可添加
    {   
        pthread_mutex_unlock(&(tpool->task_lock));
        return 1;
    }

    (taskqueue->tasks[idx]).function_name = function_name;
    (taskqueue->tasks[idx]).parameters = parameters;
    (taskqueue->tasks[idx]).synchronization_tag = synchronization_tag;
	(taskqueue->tasks[idx]).numa_rank = numa_rank;
    taskqueue->protected[idx]=1;//该处不可添加任务

    pthread_mutex_lock(&tpool->synchronization_locks[synchronization_tag]);
    tpool->synchronization_tasks[synchronization_tag]++;
    pthread_mutex_unlock(&tpool->synchronization_locks[synchronization_tag]);
    
    //更新任务队列
    taskqueue->length++;
    taskqueue->tail++;
    taskqueue->tail = (taskqueue->tail==TASK_BUFFER_SIZE) ? 0 : taskqueue->tail;
    //放开锁之后退出临界区
    pthread_mutex_unlock(&tpool->task_lock);
    //唤醒一个管理者线程（已经做好准备了的线程）
    pthread_cond_signal(&((tpool->tcbs[0])->cond_signal));
    return 0;
}


int TPSM_barrier(TPSM_t * const tpool, int const synchronization_tag)
{	
	if(tpool==NULL) return -1;
    pthread_mutex_lock(&tpool->synchronization_locks[synchronization_tag]);
    while( (tpool->synchronization_tasks[synchronization_tag]!=0) )
    {
        pthread_cond_wait(&tpool->synchronization_cond_signals[synchronization_tag], &tpool->synchronization_locks[synchronization_tag]);
    }
    pthread_mutex_unlock(&tpool->synchronization_locks[synchronization_tag]);
    return 0;
}



int TPSM_destroy(TPSM_t * const tpool, int const shutdown_option)
{	
	int checkcode;
	int shutdown;
    //对参数进行判断
    switch(shutdown_option)
    {
        case TPSM_KEEP_RUNNING:
            return 1;
        case TPSM_SHUTDOWN_IMMEDIATELY:
            shutdown = shutdown_option;
            break;
        case TPSM_SHUTDOWN_GENTLY:
            shutdown = shutdown_option;
            break;
        default:
            shutdown = TPSM_SHUTDOWN_GENTLY;
            break;
    }
	tpool->shutdown = shutdown;
	
	//唤醒所有线程
    for(int i=0;i<tpool->pool_max_threads;++i)
    {
        pthread_cond_signal(&(tpool->tcbs[i]->cond_signal));
    }
	//优先回收工作者线程
    for(int i=1;i<tpool->pool_max_threads;++i)
    {
        pthread_join(tpool->thrd_ids[i], NULL);
    }

	//最后回收管理者线程
    pthread_join(tpool->thrd_ids[0], NULL);
	//释放同步相关资源
	for(int i=0;i<SYNCHRONIZATION_MAX_SIZE;++i)
    {
        pthread_cond_destroy(&tpool->synchronization_cond_signals[i]);
        pthread_mutex_destroy(&tpool->synchronization_locks[i]);
    }
	TPSM_align_free(tpool->synchronization_tasks);
	//释放task_buffer资源
	pthread_mutex_destroy(&tpool->task_lock);
	TPSM_align_free(tpool->task_buffer);
	
	//释放线程相关资源
	TPSM_align_free(tpool->thrd_attr_onnode);

	for(int i=0;i<tpool->numa_nodes;++i)
	{
		TPSM_numa_free(tpool->numa_node_worker_ranks[i],(tpool->numa_node_workers[i])*sizeof(*tpool->numa_node_worker_ranks[i]));
	}
	TPSM_numa_free(tpool->numa_node_worker_ranks,tpool->numa_nodes_utilization*sizeof(*tpool->numa_node_worker_ranks));	
	TPSM_align_free(tpool->thrd_nr);
	TPSM_align_free(tpool->thrd_lwp);
	TPSM_numa_free(tpool->numa_node_workers, tpool->numa_nodes_utilization*sizeof(*tpool->numa_node_workers));
	TPSM_numa_free(tpool->tcbs,tpool->pool_max_threads*sizeof(*tpool->tcbs));
	
	pthread_mutex_destroy(&tpool->pool_lock);
	
	free(tpool);
	return 0;
}

/**************************************************************
* PRIVATE FUNCTIONS
**************************************************************/

//初始化numa节点核数的数据存入TPSM中
int TPSM_initialize_numacores(TPSM_t * tpool)
{	
	if(tpool==NULL) return -1;
	tpool->numa_cores = TPSM_malloc_align(tpool->numa_nodes*sizeof(*tpool->numa_cores));
	TPSM_assert(tpool->numa_cores,1);
	for(int i=0;i<tpool->numa_nodes;++i)
	{
		tpool->numa_cores[i] = TPSM_get_numacores(i);
	}
	return 0;
}

//初始化CPU集
int TPSM_initialize_cpusets(TPSM_t * tpool)
{
	if(tpool==NULL) return -1;
	if(tpool->sys_cores<=0) return -1;
	if(tpool->numa_nodes<=0) return -1;

	tpool->cpu_sets = TPSM_malloc_align(tpool->numa_nodes*sizeof(*tpool->cpu_sets));
	for(int i=0;i<tpool->numa_nodes;++i)
	{
		CPU_ZERO(&tpool->cpu_sets[i]);
	}

	size_t numa_rank;
	for(int i=0;i<tpool->sys_cores;++i)
	{
		numa_rank = numa_node_of_cpu(i);
		CPU_SET(i,&tpool->cpu_sets[numa_rank]);
	}

	return 0;
}

//初始化总任务buffer
int TPSM_initialize_taskbuffer(TPSM_t * tpool)
{
	if(tpool==NULL) return -1;
	if(tpool->sys_cores<=0) return -1;
	if(tpool->numa_nodes<=0) return -1;
	if(tpool->numa_cores==NULL) return -1;
	if(tpool->cpu_sets==NULL) return -1;

	tpool->task_buffer = TPSM_malloc_align(sizeof(*tpool->task_buffer));
	TPSM_assert(tpool->task_buffer,1);

	pthread_mutex_init(&tpool->task_lock,NULL); //初始化任务buffer 互斥锁

	return 0;
}

//初始化同步任务相关
int TPSM_synchronization_init(TPSM_t * tpool)
{	
	if(tpool==NULL) return -1;
    tpool->synchronization_locks = TPSM_malloc_align(SYNCHRONIZATION_MAX_SIZE*sizeof(*tpool->synchronization_locks));
    tpool->synchronization_cond_signals = TPSM_malloc_align(SYNCHRONIZATION_MAX_SIZE*sizeof(*tpool->synchronization_cond_signals));
    tpool->synchronization_tasks = TPSM_malloc_align(SYNCHRONIZATION_MAX_SIZE*sizeof(*tpool->synchronization_tasks));
	for(int i=0; i<SYNCHRONIZATION_MAX_SIZE; ++i)
	{
		pthread_mutex_init(&tpool->synchronization_locks[i],NULL);
    	pthread_cond_init(&tpool->synchronization_cond_signals[i],NULL);
	}
	return 0;
}
//初始化线程pthread_attr_t属性
int TPSM_initialize_pthreadattr(TPSM_t * tpool)
{
	if(tpool==NULL) return -1;
	if(tpool->sys_cores<=0) return -1;
	if(tpool->numa_nodes<=0) return -1;
	if(tpool->numa_cores==NULL) return -1;
	if(tpool->cpu_sets==NULL) return -1;

	int checkcode;

	size_t stacksize = 16UL*1024UL*1024UL; //设置线程栈的大小
	
	tpool->thrd_attr_onnode = TPSM_malloc_align(tpool->numa_nodes*sizeof(*tpool->thrd_attr_onnode));
	
	for(int i=0; i<tpool->numa_nodes; ++i)
	{	
		checkcode = pthread_attr_init(&tpool->thrd_attr_onnode[i]);//初始化
        TPSM_assert(checkcode,0);
        checkcode = pthread_attr_setdetachstate(&tpool->thrd_attr_onnode[i],PTHREAD_CREATE_JOINABLE); //可join
        TPSM_assert(checkcode,0);
		checkcode = pthread_attr_setaffinity_np(&tpool->thrd_attr_onnode[i], sizeof(*tpool->cpu_sets), &tpool->cpu_sets[i]); //绑定CPU集
		TPSM_assert(checkcode,0);
		checkcode = pthread_attr_setstacksize(&tpool->thrd_attr_onnode[i],stacksize); //设置线程栈的大小
        TPSM_assert(checkcode,0);
	}

	//设置quicker_attr：不绑定核，全系统竞争
	checkcode = pthread_attr_init(&tpool->quicker_attr);//初始化
	TPSM_assert(checkcode,0);
	checkcode = pthread_attr_setdetachstate(&tpool->quicker_attr,PTHREAD_CREATE_JOINABLE); //可join
	TPSM_assert(checkcode,0);
	checkcode = pthread_attr_setstacksize(&tpool->quicker_attr,stacksize); //设置线程栈的大小
	TPSM_assert(checkcode,0);
	checkcode = pthread_attr_setscope(&tpool->quicker_attr,PTHREAD_SCOPE_SYSTEM);
	TPSM_assert(checkcode,0);

	return 0;
}

//根据策略获取利用到的节点数
size_t TPSM_numa_nodes_utilization(size_t const threads_nbr, size_t const numa_nodes, size_t  * const numa_cores, int const strategy)
{	
	#ifdef DEBUG
	printf("FUNCTION: %s\n",__FUNCTION__);
	#endif
	size_t nodes_utilization = 1;
	size_t temp = threads_nbr;
	int i = 0;
	switch(strategy)
	{
		case TPSM_EVEN_DISTRIBUTION://均分线程
			if(threads_nbr%numa_nodes!=0)
			{
				return -1;
			}
			else
			{	
				nodes_utilization = numa_nodes;
				break;
			}
		case TPSM_CENTRALIZED_DISTRIBUTION://集中分线程
			
			for(i=0; i<numa_nodes; ++i)
			{
				if(temp>numa_cores[i])
				{
					temp-=numa_cores[i];
					nodes_utilization++;
				}
				else
				{
					break;
				}
			}

			if(threads_nbr%nodes_utilization!=0)
			{
				return -1;
			}
			else
			{
				break;
			}

		default:
			return -1;
	}
	
	return nodes_utilization;
}


//根据策略初始化线程分配
int TPSM_initialize_distribution(TPSM_t * tpool, int const strategy)
{	
	if(tpool==NULL) return -1;
	if(tpool->thrd_nr==NULL) return -1;
	if(tpool->pool_max_workers%2!=0) return -1; //简单的判断，被分配的线程数必须是偶数
	
	//计算出要利用的节点数:这里算出来的值必定能整除pool_max_workers
	size_t const numa_nodes_utilization = TPSM_numa_nodes_utilization(tpool->pool_max_workers, tpool->numa_nodes, tpool->numa_cores, strategy);
	tpool->numa_nodes_utilization = numa_nodes_utilization;

	//管理线程的numa节点设置为1
	tpool->thrd_nr[0] = 1; 

	tpool->numa_node_workers = TPSM_numalloc_onnode(numa_nodes_utilization*sizeof(*tpool->numa_node_workers), tpool->thrd_nr[0]);
	TPSM_assert(tpool->numa_node_workers,1);
	
	//计算出每个numa节点的工作者线程个数
	size_t numa_rank = 0;
	for(int i=0;i<tpool->pool_max_workers;++i)
	{	
		tpool->numa_node_workers[numa_rank]++;
		numa_rank++;
		numa_rank = (numa_rank==numa_nodes_utilization) ? 0 : numa_rank; 
	}

	tpool->numa_node_worker_ranks = TPSM_numalloc_onnode(numa_nodes_utilization*sizeof(*tpool->numa_node_worker_ranks), tpool->thrd_nr[0]);
	TPSM_assert(tpool->numa_node_worker_ranks,1);

	for(int i=0;i<numa_nodes_utilization;++i)
	{	
		tpool->numa_node_worker_ranks[i] = TPSM_numalloc_onnode((tpool->numa_node_workers[i])*sizeof(*tpool->numa_node_worker_ranks[i]), i);
		TPSM_assert(tpool->numa_node_worker_ranks[i],1);
	}

	numa_rank = 0;
	size_t count = 0;
	for(int i=1;i<tpool->pool_max_threads;++i)
	{
		tpool->thrd_nr[i] = numa_rank;
		tpool->numa_node_worker_ranks[numa_rank][count] = i;
		count++;
		if(count==tpool->numa_node_workers[numa_rank])
		{
			numa_rank++;
			count = 0;
		} 
	}

	#ifdef DEBUG
	printf("pool_max_threads = %d\npool_max_workers=%d\nnuma_nodes_utilization=%d\n",tpool->pool_max_threads,tpool->pool_max_workers,numa_nodes_utilization);
	for(int i=0;i<numa_nodes_utilization;++i)
	{	
		printf("numa_node (%d) : ",i);
		for(int j=0;j<tpool->numa_node_workers[i];++j)
		{
			printf("%d ",tpool->numa_node_worker_ranks[i][j]);
		}
		printf("\n");
	}

	printf("thrd_nr = ");
	for(int i=0;i<tpool->pool_max_threads;++i)
	{
		printf("%d ",tpool->thrd_nr[i]);
	}
	printf("\n");
	#endif

	return 0;
}


//初始化线程(创建管理者线程与工作者线程)
int TPSM_initialize_threads(TPSM_t * tpool, int const strategy)
{	
	if(tpool==NULL) return -1;
	if(tpool->thrd_attr_onnode==NULL) return -1;
	
	int checkcode;
	
	tpool->thrd_nr = TPSM_malloc_align(tpool->pool_max_threads*sizeof(*tpool->thrd_nr));
	TPSM_assert(tpool->thrd_nr,1);

	tpool->thrd_ids = TPSM_malloc_align(tpool->pool_max_threads*sizeof(*tpool->thrd_ids));
	TPSM_assert(tpool->thrd_ids,1);
	
	tpool->thrd_lwp = TPSM_malloc_align(tpool->pool_max_threads*sizeof(*tpool->thrd_lwp));
	TPSM_assert(tpool->thrd_lwp,1);

	//根据策略初始化线程分配
	checkcode = TPSM_initialize_distribution(tpool, strategy);
	TPSM_assert(checkcode,0);

	//创建tcb指针数组，开空间到管理者线程所在的numa节点下
	tpool->tcbs = TPSM_numalloc_onnode(tpool->pool_max_threads*sizeof(*tpool->tcbs),tpool->thrd_nr[0]);
	TPSM_assert(tpool->tcbs,1);

	//创建管理者线程
	checkcode = pthread_create(&tpool->thrd_ids[0],&(tpool->thrd_attr_onnode[tpool->thrd_nr[0]]),&TPSM_manager,tpool);
	TPSM_assert(checkcode,0);

	//创建工作者线程
	for(int i=1;i<tpool->pool_max_threads;++i)
	{
		checkcode = pthread_create(&tpool->thrd_ids[i],&(tpool->thrd_attr_onnode[tpool->thrd_nr[i]]),&TPSM_worker,tpool);
		TPSM_assert(checkcode,0);
	}

	return 0;
}


int TPSM_find_workers_rank(TPSM_t * tpool,size_t const task_numa_rank)
{
	int active_rank;
	int numa_idx = task_numa_rank;
	int count = 0;
	//遍历该numa节点的所有线程
	while(count < tpool->numa_nodes)
	{	
		
		for(int i=0; i<tpool->numa_node_workers[numa_idx]; ++i)
		{	
			active_rank = tpool->numa_node_worker_ranks[numa_idx][i];
			if( (tpool->tcbs[active_rank])->state==0 )
			{
				break;
			}
			else
			{
				active_rank = -1;
			}
		}
		if(active_rank <0 )
		{
			count++;
			numa_idx++;
			numa_idx = (numa_idx == tpool->numa_nodes) ? 0 : numa_idx;
		}
		else
		{
			break;
		}

	}
	
	return active_rank;
}

void TPSM_manager(void * arg)
{	
	int checkcode;
	TPSM_t * tpool = (TPSM_t*)arg;
	size_t const my_rank = TPSM_get_myrank(tpool);
	size_t const my_numa_node_rank = tpool->thrd_nr[my_rank];
	size_t const my_lwp = TPSM_get_mylwp(tpool);
	tpool->thrd_lwp[my_rank] = my_lwp;

	//创建自己的本地tcb并初始化
	TPSM_TCB_t * my_tcb_p = TPSM_numalloc_local(sizeof(*my_tcb_p));
	//TPSM_TCB_t * my_tcb = TPSM_numalloc_onnode(sizeof(*my_tcb),my_numa_node_rank);
	TPSM_assert(my_tcb_p,1);
	pthread_mutex_init(&my_tcb_p->lock,NULL);
	pthread_cond_init(&my_tcb_p->cond_signal,NULL);
	my_tcb_p->state = 1;

	//将创建好的本地tcb指针存储tpool中
	tpool->tcbs[my_rank] = my_tcb_p;

	TPSM_TASK_BUFFER_t * const taskqueue = tpool->task_buffer;
	
	//可用线程数增加，并向发送信号（为了同步到TPSM_init）
	pthread_mutex_lock(&tpool->pool_lock);
	tpool->thread_available_nbr++;
	pthread_cond_signal(&tpool->pool_cond_signal);
	pthread_mutex_unlock(&tpool->pool_lock);


	#ifdef DEBUG
	printf("my_rank = %d\n",my_rank);
	#endif

	int active_rank = -1;
	size_t task_numa_rank;
	while(1)
    {   
        pthread_mutex_lock(&my_tcb_p->lock);
		// my_tcb_p->state = 0;
        while( taskqueue->length==0 && tpool->shutdown == TPSM_KEEP_RUNNING)
        {   
            pthread_cond_wait(&my_tcb_p->cond_signal,&my_tcb_p->lock);
        }
		// my_tcb_p->state = 1; //仅管理者线程有此步，工作者线程由管理者线程设置
        pthread_mutex_unlock(&my_tcb_p->lock);
		
		#ifdef DEBUG
		printf("管理者线程醒了\n");
		#endif
		
		if(tpool->shutdown==TPSM_SHUTDOWN_IMMEDIATELY)
        {
            break;
        }

        if(tpool->shutdown==TPSM_SHUTDOWN_GENTLY && taskqueue->length==0)
        {
            break;
        }

		//能走到此步，证明任务buffer里必定有一个任务
		
		//获得任务的numa_rank
		task_numa_rank = taskqueue->tasks[taskqueue->head].numa_rank;
		//根据任务的numa_rank给找到numa_rank中的线程

		active_rank = TPSM_find_workers_rank(tpool,task_numa_rank);
		

		#ifdef DEBUG
		printf("要激活%d号线程\n",active_rank);
		#endif

		pthread_mutex_lock(&tpool->task_lock);
		if(active_rank<0)//没找到可激活的线程
		{	
			checkcode = TPSM_create_quicker(tpool,(taskqueue->tasks[taskqueue->head]).function_name, (taskqueue->tasks[taskqueue->head]).parameters, (taskqueue->tasks[taskqueue->head]).synchronization_tag);
			TPSM_assert(checkcode,0);
		}
		else //找到了可激活的线程
		{
			//把任务放到这个rank的TCB中然后更新任务队列
			(tpool->tcbs[active_rank])->function_name = (taskqueue->tasks[taskqueue->head]).function_name;
			(tpool->tcbs[active_rank])->parameters = (taskqueue->tasks[taskqueue->head]).parameters;
			(tpool->tcbs[active_rank])->synchronization_tag = (taskqueue->tasks[taskqueue->head]).synchronization_tag;
			pthread_mutex_lock(&((tpool->tcbs[active_rank])->lock));
			(tpool->tcbs[active_rank])->state = 1;
			pthread_mutex_unlock(&((tpool->tcbs[active_rank])->lock));
			//激活这个线程
			pthread_cond_signal(&((tpool->tcbs[active_rank])->cond_signal));
		}
		
		//更新任务队列
        taskqueue->protected[taskqueue->head] = 0;
        taskqueue->head++;
        taskqueue->head = (taskqueue->head == TASK_BUFFER_SIZE) ? 0 : taskqueue->head;
        taskqueue->length--;
		pthread_mutex_unlock(&tpool->task_lock);

	}
	
	//将自己tcb中的锁和条件信号销毁
	TPSM_free(tpool);
	pthread_exit(NULL);
}



void TPSM_worker(void * arg)
{
	TPSM_t * tpool = (TPSM_t*)arg;
	size_t const my_rank = TPSM_get_myrank(tpool);
	size_t const my_numa_node_rank = tpool->thrd_nr[my_rank];
	size_t const my_lwp = TPSM_get_mylwp(tpool);
	tpool->thrd_lwp[my_rank] = my_lwp;

	//创建自己的本地tcb并初始化
	TPSM_TCB_t * my_tcb_p = TPSM_numalloc_local(sizeof(*my_tcb_p));
	//size_t const my_numa_node_rank = tpool->thrd_nr[my_rank];
	//TPSM_TCB_t * my_tcb = TPSM_numalloc_onnode(sizeof(*my_tcb),my_numa_node_rank);
	TPSM_assert(my_tcb_p,1);
	pthread_mutex_init(&my_tcb_p->lock,NULL);
	pthread_cond_init(&my_tcb_p->cond_signal,NULL);
	my_tcb_p->state = 0;

	//将创建好的本地tcb指针存储tpool中
	tpool->tcbs[my_rank] = my_tcb_p;

	TPSM_TASK_BUFFER_t * const taskqueue = tpool->task_buffer;

	//可用线程数增加，并向发送信号（为了同步到TPSM_init）
	pthread_mutex_lock(&tpool->pool_lock);
	tpool->thread_available_nbr++;
	pthread_cond_signal(&tpool->pool_cond_signal);
	pthread_mutex_unlock(&tpool->pool_lock);
	
	#ifdef DEBUG
	printf("my_rank = %d\n",my_rank);
	#endif

	while(1)
    {   

        pthread_mutex_lock(&my_tcb_p->lock);
		my_tcb_p->state = 0;
        while( my_tcb_p->state == 0 && tpool->shutdown== TPSM_KEEP_RUNNING)
        {   
            pthread_cond_wait(&my_tcb_p->cond_signal,&my_tcb_p->lock);
        }
        pthread_mutex_unlock(&my_tcb_p->lock);

		if(tpool->shutdown==TPSM_SHUTDOWN_IMMEDIATELY)
        {
            break;
        }

        if(tpool->shutdown==TPSM_SHUTDOWN_GENTLY && taskqueue->length==0)
        {
            break;
        }
		
		#ifdef DEBUG
		printf("工作线程%d醒了\n",my_rank);
		#endif

		(my_tcb_p->function_name)(my_tcb_p->parameters);

		pthread_mutex_lock(&tpool->synchronization_locks[my_tcb_p->synchronization_tag]);
        tpool->synchronization_tasks[my_tcb_p->synchronization_tag]--;
        if(tpool->synchronization_tasks[my_tcb_p->synchronization_tag]==0)
        {
            pthread_cond_broadcast(&tpool->synchronization_cond_signals[my_tcb_p->synchronization_tag]);
        }
        pthread_mutex_unlock(&tpool->synchronization_locks[my_tcb_p->synchronization_tag]);
	}

	//将自己tcb中的锁和条件信号销毁
	TPSM_free(tpool);
	//free自己的TCB
	pthread_exit(NULL);
}



void TPSM_quicker(void * arg)
{
	TPSM_QUICKER_TCB_t * my_tcb_p = (TPSM_QUICKER_TCB_t*)arg;

	TPSM_t * tpool = my_tcb_p->tpool;

	(my_tcb_p->function_name)(my_tcb_p->parameters);
	
	pthread_mutex_lock(&tpool->synchronization_locks[my_tcb_p->synchronization_tag]);
	tpool->synchronization_tasks[my_tcb_p->synchronization_tag]--;
	if(tpool->synchronization_tasks[my_tcb_p->synchronization_tag]==0)
	{
		pthread_cond_broadcast(&tpool->synchronization_cond_signals[my_tcb_p->synchronization_tag]);
	}
	pthread_mutex_unlock(&tpool->synchronization_locks[my_tcb_p->synchronization_tag]);

	//free自己的TCB
	TPSM_numa_free(my_tcb_p,sizeof(*my_tcb_p));
	pthread_exit(NULL);
}


//创建快速线程
int TPSM_create_quicker
(
	TPSM_t * tpool,
	void * (*function_name)(void*),
    void * parameters,
	size_t synchronization_tag
)
{
	if(tpool==NULL) return -1;
	TPSM_QUICKER_TCB_t * tcb = TPSM_malloc_align(sizeof(*tcb));
	TPSM_assert(tcb,1);
	tcb->tpool = tpool;
	tcb->function_name = function_name;
	tcb->parameters = parameters;
	tcb->synchronization_tag = synchronization_tag;
	pthread_t tid;
	pthread_create(&tid,&tpool->quicker_attr,&TPSM_quicker,tcb);
	return 0;
}


int TPSM_free(TPSM_t * const tpool)
{	
	size_t const my_rank = TPSM_get_myrank(tpool);
	size_t const my_numa_node_rank = tpool->thrd_nr[my_rank];
	TPSM_TCB_t * my_tcb_p =  tpool->tcbs[my_rank];
	
	//free TCB
	pthread_cond_destroy(&my_tcb_p->cond_signal);
    pthread_mutex_destroy(&my_tcb_p->lock);
	TPSM_numa_free(my_tcb_p,sizeof(*my_tcb_p));

	return 0;
}
