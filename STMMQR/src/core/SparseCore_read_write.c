/**
 * @file SparseCore_read_write.c
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-20
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef NCHECK

#include "Sparse_internal.h"
#include "SparseCore.h"
#include <string.h>
#include <ctype.h>

/* mtx格式规定最大行长度为1024 */
#define MAXLINE 1030

/**
 * @brief 读取文件的一行，如果成功返回TRUE，如果EOF返回FALSE。
 * 
 * @param f 
 * @param buf 
 * @return int 
 */
static int get_line (FILE *f, char *buf)
{
    buf [0] = '\0' ;
    buf [1] = '\0' ;
    buf [MAXLINE] = '\0' ;
    return (fgets (buf, MAXLINE, f) != NULL) ;
}

/**
 * @brief 用 +/-Inf 替换巨大的值，因为scanf和printf不能正确处理Inf。
 * 
 * @param x 
 * @return double 
 */
static double fix_inf (double x)
{
    if ((x >= HUGE_DOUBLE) || (x <= -HUGE_DOUBLE))
    {
	/* 假设2*x导致溢出，将其视为+/- Inf */
	x = 2*x ;
    }
    return (x) ;
}

/**
 * @brief 如果s是空行或注释，为真，否则为假
 * 
 * @param s 
 * @return int 
 */
static int is_blank_line
(
    char *s
)
{
    int c, k ;
    if (s [0] == '%')
    {
	/* 一个注释行 */
	return (TRUE) ;
    }
    for (k = 0 ; k <= MAXLINE ; k++)
    {
	c = s [k] ;
	if (c == '\0')
	{
	    /* 行结束 */
	    break ;
	}
	if (!isspace (c))
	{
	    /* non-space character */
	    return (FALSE) ;
	}
    }
    return (TRUE) ;
}


#define STYPE_UNKNOWN 999
#define STYPE_SYMMETRIC_UPPER 1
#define STYPE_UNSYMMETRIC 0     //“一般”存储(非对称矩阵)，两部分都存在
#define STYPE_SYMMETRIC_LOWER -1
#define STYPE_SKEW_SYMMETRIC -2

/**
 * @brief 	读抬头 ，这包含零个或多个注释行(空白，或在第一列以“%”开头)，
 *  		后跟一个最多包含四个数值的数据行。
 * 
 */
static int read_header	
(
    /* ---- input ---- */
    FILE *f,		/* 需要读取的文件 */
    /* ---- output --- */
    char *buf,		/* 读取行，长度为MAXLINE+1 */
    int *mtype,		/* SPARSE_TRIPLET三元组 或者 SPARSE_DENSE稠密 */
    size_t *nrow,	/* 矩阵的行 */
    size_t *ncol,	/* 矩阵的列 */
    size_t *nnz,	/* 三元组中的成员数量（稠密矩阵为0）*/
    int *stype		/*  */
)
{
    char *p ;
    int first = TRUE, got_mm_header = FALSE, c, c2, nitems ;
    double l1, l2, l3, l4 ;

    *mtype = SPARSE_TRIPLET ;
    *nrow = 0 ;
    *ncol = 0 ;
    *nnz = 0 ;
    *stype = STYPE_UNKNOWN ;

    for ( ; ; )
    {

	/* 获取下一行 */
	if (!get_line (f, buf))
	{
	    /* 文件提前读完 */
	    return (FALSE) ;
	}

	if (first && (strncmp (buf, "%%MatrixMarket", 14) == 0))  //读数据头
	{

	    /* -------------------------------------------------------------- */
	    /* 读矩阵的头 */
	    /* -------------------------------------------------------------- */

	    got_mm_header = TRUE ;
	    p = buf ;

	    /* -------------------------------------------------------------- */
	    /* 读头一行的matrix */
	    /* -------------------------------------------------------------- */

	    while (*p && !isspace (*p)) p++ ;
	    while (*p &&  isspace (*p)) p++ ;
	    c = tolower (*p) ;
	    if (c != 'm')
	    {
		/* bad format */
		return (FALSE) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* 区分第一行的coordinate和 array*/
	    /* -------------------------------------------------------------- */

	    while (*p && !isspace (*p)) p++ ;
	    while (*p &&  isspace (*p)) p++ ;
	    c = tolower (*p) ;
	    if (c == 'c')
	    {
		*mtype = SPARSE_TRIPLET  ; //三元组形式
	    }
	    else if (c == 'a')
	    {
		*mtype = SPARSE_DENSE  ; // 稠密形式
	    }
	    else
	    {
		/* 格式错误 */
		return (FALSE) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* 得到数据形式 (real, pattern, integer) */
	    /* -------------------------------------------------------------- */

	    while (*p && !isspace (*p)) p++ ;
	    while (*p &&  isspace (*p)) p++ ;
	    c = tolower (*p) ;
	    if (!(c == 'r' || c == 'p' || c == 'c' || c == 'i')) // real /pattern /complex /integer
	    {
		/* 错误的格式 */
		return (FALSE) ;
	    }
		// if (c == 'r')
		// {
		// 	printf("THIS IS A REAL MATRIX !\n");
		// }
		// else if( c == 'p')
		// {
		// 	printf("THIS IS A PATTERN MATRIX !\n");
		// }
		// else if( c == 'c')
		// {
		// 	printf("ERROR about complex MATRIX !\n");
		// }
		// else
		// {
		// 	printf("THIS IS A integer MATRIX !\n");
		// }
		
	    /* -------------------------------------------------------------- */
	    /* 获取矩阵类型(general, hermitian, symmetric, skew)  */
	    /* -------------------------------------------------------------- */

	    while (*p && !isspace (*p)) p++ ;
	    while (*p &&  isspace (*p)) p++ ;
	    c = tolower (*p) ;
	    c2 = tolower (*(p+1)) ;
		// 添加的写文件
		// char *result_file = (char *)"./ARMQR.txt";
		// FILE* fresult;
		// fresult = fopen(result_file,"a+");
		

	    if (c == 'g')
	    {
			/* 一般非对称矩阵（全部） */
			*stype = STYPE_UNSYMMETRIC ;
			//printf("UNSYMMETRIC WHEN READ\n");
			// fprintf(fresult, "UNSYMMETRIC WHEN READ\t" );
	    }
	    else if (c == 's' && c2 == 'y')
	    {
		    /* 实对称矩阵（下三角部分） */
		    *stype = STYPE_SYMMETRIC_LOWER ;
			//printf("SYMMETRIC_LOWER WHEN READ\n");
			// fprintf(fresult, "SYMMETRIC WHEN READ\t" );
	    }
	    else if (c == 'h')
	    {
		/* 埃米尔特矩阵（下三角部分） */
			*stype = STYPE_SYMMETRIC_LOWER ;
			//printf("SYMMETRIC_LOWER WHEN READ\n");
			// fprintf(fresult, "SYMMETRIC WHEN READ\t" );
	    }
	    else if (c == 's' && c2 == 'k')
	    {
			/* 反对称矩阵 A'+ A = 0（下三角部分） */
			*stype = STYPE_SKEW_SYMMETRIC ;
			//printf("SKEW_SYMMETRIC WHEN READ\n");
	    }
	    else
	    {
		/* 错误的格式 */
		return (FALSE) ;
	    }
		// fclose(fresult);
	}
	else if (is_blank_line (buf)) // 空行跳过
	{
	    continue ;
	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /*  读取第一行数据并返回 */
	    /* -------------------------------------------------------------- */

	    /* 格式: nrow ncol nnz stype */
	    l1 = EMPTY ;
	    l2 = EMPTY ;
	    l3 = 0 ;
	    l4 = 0 ;
	    nitems = sscanf (buf, "%lg %lg %lg %lg\n", &l1, &l2, &l3, &l4) ;
	    if (nitems < 2 || nitems > 4 || l1 > Int_max || l2 > Int_max)
	    {
		/* 检查错误 */
		return (FALSE) ;
	    }
	    *nrow = l1 ;
	    *ncol = l2 ;
	    if (nitems == 2)
	    {
		/* 稠密矩阵 */
		if (!got_mm_header)
		{
		    *mtype = SPARSE_DENSE ;
		    *stype = STYPE_UNSYMMETRIC ;
		}
	    }
	    if (nitems == 3 || nitems == 4)
	    {
		/* 稀疏矩阵 */
		*nnz = l3 ;
		if (!got_mm_header)
		{
		    *mtype = SPARSE_TRIPLET ;
		}
	    }

	    if (nitems == 4)
	    {
		/* 这里指定的样式只能是1、0或-1,设定对称矩阵存的位置 */
		if (l4 < 0)
		{
		    *stype = STYPE_SYMMETRIC_LOWER ;
		}
		else if (l4 > 0)
		{
		    *stype = STYPE_SYMMETRIC_UPPER ;
		}
		else
		{
		    *stype = STYPE_UNSYMMETRIC ;
		}
	    }
	    if (*nrow != *ncol)
	    {
		/* 方形矩阵一定是非对称的 */
		*stype = STYPE_UNSYMMETRIC ;
	    }
	    return (TRUE) ;
	}

	first = FALSE ;
    }
}


/**
 * @brief 头已经读完了，读三元组
 * 
 */
static sparse_triplet *read_triplet
(
    /* ---- input ---- */
    FILE *f,		    /* 需要读取的文件，必须已经打开 */
    size_t nrow,	    /* 行数 */
    size_t ncol,	    /* 列数 */
    size_t nnz,		    /* 需要读取的三元组的数目 */
    int stype,		    /* unknown或者从抬头读取到的数据类型 */
    int prefer_unsym,	/* 如果为TRUE，T->stype类型的zero */
    /* ---- workspace */
    char *buf,		    /* 大小为MAXLINE+1 */
    /* --------------- */
    sparse_common *Common
)
{
    double x, z ;
    double *Tx ;   //并不区分格式，统一使用double保存数据
    Int *Ti, *Tj, *Rdeg, *Cdeg ;
    sparse_triplet *T ;
    double l1, l2 ;
    Int nitems, xtype, unknown, k, nshould, is_lower, is_upper, one_based, i, j,
	imax, jmax, skew_symmetric, p ;
    size_t s, nnz2, extra ;
    int ok = TRUE ;

    /* 判断矩阵是否为空 */
    if (nrow == 0 || ncol == 0 || nnz == 0)
    {
	return (SparseCore_allocate_triplet (nrow, ncol, 0, 0, SPARSE_REAL,
		    Common)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 区分特殊格式 */
    /* ---------------------------------------------------------------------- */

    unknown = (stype == STYPE_UNKNOWN) ;
    skew_symmetric = (stype == STYPE_SKEW_SYMMETRIC) ;

    extra = 0 ;
    if (stype < STYPE_SYMMETRIC_LOWER
	|| (prefer_unsym && stype != STYPE_UNSYMMETRIC)) // 对称矩阵走的是这里
    {
	/* 999: unknown可能被转换为非对称 */
	/*  1:  如果prefer_unsym为真，对称上三角转换成非对称 */
	/* -1:  如果prefer_unsym为真，对称下三角转换成非对称 */
	/* -2:  实斜对称转换成非对称 */
	//printf("change to UNSYMMETRIC\n");
	stype = STYPE_UNSYMMETRIC ;
	extra = nnz ;
    }
    nnz2 = SparseCore_add_size_t (nnz, extra, &ok) ;

    /* ---------------------------------------------------------------------- */
    /* 分配内存空间 */
    /* ---------------------------------------------------------------------- */

    /* s = nrow + ncol */
    s = SparseCore_add_size_t (nrow, ncol, &ok) ;
    if (!ok || nrow > Int_max || ncol > Int_max || nnz > Int_max)
    {
	ERROR (SPARSE_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }

    SparseCore_allocate_work (0, s, 0, Common) ;
    Rdeg = Common->Iwork ;	/* 大小为nrow */
    Cdeg = Rdeg + nrow ;	/* 大小为ncol */

    /* ---------------------------------------------------------------------- */
    /* 读取三元组 */
    /* ---------------------------------------------------------------------- */

    is_lower = TRUE ;
    is_upper = TRUE ;
    one_based = TRUE ;
    imax = 0 ;
    jmax = 0 ;

    Tx = NULL ;
    Ti = NULL ;
    Tj = NULL ;
    xtype = 999 ;
    nshould = 0 ;

    for (k = 0 ; k < (Int) nnz ; k++)
    {

	/* ------------------------------------------------------------------ */
	/*        获取下一个三元组，跳过空行和注释行       */
	/* ------------------------------------------------------------------ */

	l1 = EMPTY ;
	l2 = EMPTY ;
	x = 0 ;
	z = 0 ;

	for ( ; ; )
	{
	    if (!get_line (f, buf)) // 没有行了返回
	    {
		/* 过早结束文件——没有读入足够的三胞胎 */
		ERROR (SPARSE_INVALID, "premature EOF") ;
		return (NULL) ;
	    }
	    if (is_blank_line (buf)) // 空行跳过
	    {
		/* 空行或者注释 */
		continue ;
	    }
	    nitems = sscanf (buf, "%lg %lg %lg %lg\n", &l1, &l2, &x, &z) ;
	    x = fix_inf (x) ;
	    z = fix_inf (z) ;
	    break ;
	}

	nitems = (nitems == EOF) ? 0 : nitems ;
	i = l1 ;
	j = l2 ;

	/* ------------------------------------------------------------------ */
	/*       对于第一个三元组: 确定类型并分配三重值矩阵       */
	/* ------------------------------------------------------------------ */

	if (k == 0)
	{
	    if (nitems < 2 || nitems > 4) // 报错
	    {
		/* 错误 */
		ERROR (SPARSE_INVALID, "invalid format") ;
		return (NULL) ;
	    }
	    else if (nitems == 2)
	    {
		/* 稍后将转换为实矩阵 */
		xtype = SPARSE_PATTERN ;/* 图片，没有数值 */
	    }
	    else if (nitems == 3)
	    {
		xtype = SPARSE_REAL ;   /* 实矩阵 */
	    }

		//   其余的行应该具有相同的项数
	    nshould = nitems ;

	    /* 给三元组矩阵分配空间 */
	    T = SparseCore_allocate_triplet (nrow, ncol, nnz2, stype,
		    (xtype == SPARSE_PATTERN ? SPARSE_REAL : xtype), Common) ;
	    if (Common->status < SPARSE_OK)
	    {
		/* out of memory */
		return (NULL) ;
	    }
	    Ti = T->i ;
	    Tj = T->j ;
	    Tx = T->x ;
	    T->nnz = nnz ;
	}

	/* ------------------------------------------------------------------ */
	/* 保存三元矩阵中的元素 */
	/* ------------------------------------------------------------------ */

	if (nitems != nshould || i < 0 || j < 0)
	{
	    /* 错误的格式、过早的文件结束或负索引 */
	    SparseCore_free_triplet (&T, Common) ;
	    ERROR (SPARSE_INVALID, "invalid matrix file") ;
	    return (NULL) ;
	}

	Ti [k] = i ;
	Tj [k] = j ;

	if (i < j)
	{
	    /* 这一项在上三角部分 */
	    is_lower = FALSE ;
	}
	if (i > j)
	{
	    /* 这一项在下三角部分 */
	    is_upper = FALSE ;
	}

	if (xtype == SPARSE_REAL) // 三元组形式
	{
	    Tx [k] = x ;
	}

	if (i == 0 || j == 0)
	{
	    one_based = FALSE ;
	}

	imax = MAX (i, imax) ;
	jmax = MAX (j, jmax) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 如果原来的矩阵基于1开始则 转换为基于零开始*/
    /* ---------------------------------------------------------------------- */

    if (one_based)   // 转换
    {
	/* 输入矩阵为从1开始的，转换为从0开始 */
	for (k = 0 ; k < (Int) nnz ; k++)
	{
	    Ti [k]-- ;
	    Tj [k]-- ;
	}
    }

    if (one_based ?
	(imax >  (Int) nrow || jmax >  (Int) ncol) :
	(imax >= (Int) nrow || jmax >= (Int) ncol)) // 检查是否超出size
    {
	/* 索引超出范围 */
	SparseCore_free_triplet (&T, Common) ;
	ERROR (SPARSE_INVALID, "indices out of range") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 确定数据类型 */
    /* ---------------------------------------------------------------------- */

    if (unknown)
    {
	if (is_lower && is_upper)
	{
	    /* 对角矩阵，上半部分对称 */
	    stype = STYPE_SYMMETRIC_UPPER ;
	}
	else if (is_lower && !is_upper)
	{
	    /* 对称矩阵，下半部分对称 */
	    stype = STYPE_SYMMETRIC_LOWER ;
	}
	else if (!is_lower && is_upper)
	{
	    /* 对称矩阵，上半部分对称 */
	    stype = STYPE_SYMMETRIC_UPPER ;
	}
	else
	{
	    /* 非对称 */
	    stype = STYPE_UNSYMMETRIC ;
	    extra = 0 ;
	}
    }

    /* ---------------------------------------------------------------------- */
	/*          将对称矩阵、斜对称矩阵或厄米矩阵的剩余部分相加             */
    /* ---------------------------------------------------------------------- */

    /* 注意，除非prefer_unsym为真，否则此步骤不会用于真实对称矩阵 */
    if (extra > 0)
    {
	p = nnz ;
	for (k = 0 ; k < (Int) nnz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i != j)
	    {
		Ti [p] = j ;
		Tj [p] = i ;
		if (xtype == SPARSE_REAL)
		{
		    if (skew_symmetric)
		    {
			Tx [p] = -Tx [k] ;
		    }
		    else
		    {
			Tx [p] =  Tx [k] ; // 再拷贝一份到对称位置
		    }
		}
		p++ ;
	    }
	}
	T->nnz = p ;
	nnz = p ;
    }

    T->stype = stype ;

    /* ---------------------------------------------------------------------- */
    /* 为只使用模式的矩阵创建值 */
    /* ---------------------------------------------------------------------- */

    if (xtype == SPARSE_PATTERN)
    {
	if (stype == STYPE_UNSYMMETRIC || Common->prefer_binary)
	{
	    /* 非对称情况，或者二进制情况 */
		//printf("PATTERN and UNSYMMETRIC\n");
		for (k = 0 ; k < (Int) nnz ; k++)
		{
		Tx [k] = 1 ;
		}
	}
	else  // 没有转换成general的 skew-sym pattern矩阵 
	{	//printf("PATTERN and SYMMETRIC\n");
		// 计算行和列的度数(不包括对角线)
	    for (i = 0 ; i < (Int) nrow ; i++) // 行初始化
	    {
		Rdeg [i] = 0 ;
	    }
	    for (j = 0 ; j < (Int) ncol ; j++) // 列初始化
	    {
		Cdeg [j] = 0 ;
	    }
	    for (k = 0 ; k < (Int) nnz ; k++) 
	    {
		i = Ti [k] ;
		j = Tj [k] ;
		if ((stype < 0 && i > j) || (stype > 0 && i < j))// stype < 0 是对称存在下方，>0是存在上方
		{
		    /* a(i,j)，a(j,i)都在矩阵中 */
		    Rdeg [i]++ ;
		    Cdeg [j]++ ;
		    Rdeg [j]++ ;
		    Cdeg [i]++ ;
		}
	    }
	    /* 赋值 */
	    for (k = 0 ; k < (Int) nnz ; k++)
	    {
		i = Ti [k] ;
		j = Tj [k] ;
		Tx [k] = (i == j) ? (1 + MAX (Rdeg [i], Cdeg [j])) : (-1) ;
		// 对角线上的元素为1 + MAX (Rdeg [i], Cdeg [j])
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 返回新的三元组矩阵 */
    /* ---------------------------------------------------------------------- */

    return (T) ;
}


/**
 * @brief 头已经被读入，包括第一行(nrow ncol)。读取稠密矩阵。
 * 
 */
static dense_array *read_dense
(
    /* ---- input ---- */
    FILE *f,		    /* 读取的文件，必须已经打开 */
    size_t nrow,	    /* 行数 */
    size_t ncol,	    /* 列数 */
    int stype,		    /* 抬头读到的数据类型 */
    /* ---- workspace */
    char *buf,		    /* 大小为MAXLINE+1 */
    /* --------------- */
    sparse_common *Common
)
{
    double x, z ;
    double *Xx = NULL ;
    dense_array *X ;
    Int nitems, xtype = -1, nshould = 0, i, j, k, kup, first ;

    /* ---------------------------------------------------------------------- */
    /* 判断矩阵是否为空 */
    /* ---------------------------------------------------------------------- */

    if (nrow == 0 || ncol == 0)
    {
	/* 返回一个空的稠密矩阵 */
	return (SparseCore_zeros (nrow, ncol, SPARSE_REAL, Common)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 读元素*/
    /* ---------------------------------------------------------------------- */

    first = TRUE ;

    for (j = 0 ; j < (Int) ncol ; j++)
    {

	/* ------------------------------------------------------------------ */
	// 获取文件中列j的第一个条目的行索引
	/* ------------------------------------------------------------------ */

	if (stype == STYPE_UNSYMMETRIC)
	{
	    i = 0 ;
	}
	else if (stype == STYPE_SKEW_SYMMETRIC)
	{
	    i = j+1 ;
	}
	else /* 读对称的下三角部分 */
	{
	    i = j ;
	}

	/* ------------------------------------------------------------------ */
	/* 得到j列 */
	/* ------------------------------------------------------------------ */

	for ( ; i < (Int) nrow ; i++)
	{

	    /* -------------------------------------------------------------- */
	    /* 获取下一个条目，跳过空白行和注释行 */
	    /* -------------------------------------------------------------- */

	    x = 0 ;
	    z = 0 ;
	    for ( ; ; )
	    {

		if (!get_line (f, buf))
		{
		    /* 文件过早结束——没有读入足够的条目 */
		    ERROR (SPARSE_INVALID, "premature EOF") ;
		    return (NULL) ;
		}

		if (is_blank_line (buf))
		{
		    /* 空行或者注释 */
		    continue ;
		}
		nitems = sscanf (buf, "%lg %lg\n", &x, &z) ;
		x = fix_inf (x) ;
		z = fix_inf (z) ;
		break ;
	    }

	    nitems = (nitems == EOF) ? 0 : nitems ;

	    /* -------------------------------------------------------------- */
	    /* 第一个条目:确定类型和给稠密矩阵分配空间 */
	    /* -------------------------------------------------------------- */

	    if (first)
	    {
		first = FALSE ;

		if (nitems < 1 || nitems > 2)
		{
		    /* 检查错误 */
		    ERROR (SPARSE_INVALID, "invalid format") ;
		    return (NULL) ;
		}
		else if (nitems == 1)
		{
		    /* 实矩阵 */
		    xtype = SPARSE_REAL ;
		}

		/* 其余的行应该有相同数量的条目 */
		nshould = nitems ;

		/* 给结果分配存储空间 */
		X = SparseCore_zeros (nrow, ncol, xtype, Common) ;
		if (Common->status < SPARSE_OK)
		{
		    /* 超出内存 */
		    return (NULL) ;
		}
		Xx = X->x ;
	    }

	    /* -------------------------------------------------------------- */
	    /* 存储在稠密矩阵中 */
	    /* -------------------------------------------------------------- */

	    if (nitems != nshould)
	    {
		/* 错误的格式或过早结束文件 */
		SparseCore_free_dense (&X, Common) ;
		ERROR (SPARSE_INVALID, "invalid matrix file") ;
		return (NULL) ;
	    }

	    k = i + j*nrow ;
	    kup = j + i*nrow ;

	    if (xtype == SPARSE_REAL)
	    {
		/* 实矩阵 */
		Xx [k] = x ;
		if (k != kup)
		{
		    if (stype == STYPE_SYMMETRIC_LOWER)
		    {
			/* 实对称矩阵 */
			Xx [kup] = x ;
		    }
		    else if (stype == STYPE_SKEW_SYMMETRIC)
		    {
			/* 实斜对称矩阵 */
			Xx [kup] = -x ;
		    }
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 返回新的稠密矩阵 */
    /* ---------------------------------------------------------------------- */

    return (X) ;
}


/**
 * @brief 从文件中读取一个三元组矩阵。
 * 
 */
sparse_triplet *SparseCore_read_triplet
(
    /* ---- input ---- */
    FILE *f,		/* 读取的文件，必须已经打开 */
    /* --------------- */
    sparse_common *Common
)
{
    char buf [MAXLINE+1] ;
    size_t nrow, ncol, nnz ;
    int stype, mtype ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (f, NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 读取标题和第一个数据行 */
    /* ---------------------------------------------------------------------- */

    if (!read_header (f, buf, &mtype, &nrow, &ncol, &nnz, &stype) ||
	mtype != SPARSE_TRIPLET)
    {
	/* 非法输入，这个函数只能读取一个三元组矩阵 */
	ERROR (SPARSE_INVALID, "invalid format") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 读取三元组矩阵 */
    /* ---------------------------------------------------------------------- */

    return (read_triplet (f, nrow, ncol, nnz, stype, FALSE, buf, Common)) ;
}


/**
 * @brief 从文件中读取稀疏矩阵。有关文件格式的讨论，请参阅sparse_read_triplet。
 * 
 */
sparse_csc *SparseCore_read_sparse
(
    /* ---- input ---- */
    FILE *f,		/* 读取的文件，必须已经打开 */
    /* --------------- */
    sparse_common *Common
)
{
    sparse_csc *A, *A2 ;
    sparse_triplet *T ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (f, NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 转成CSC格式*/
    /* ---------------------------------------------------------------------- */

    T = SparseCore_read_triplet (f, Common) ;
    A = SparseCore_triplet_to_sparse (T, 0, Common) ;
    SparseCore_free_triplet (&T, Common) ;

    if (Common->prefer_upper && A != NULL && A->stype == -1)
    {
	/* A=A' */
	A2 = SparseCore_transpose (A, 2, Common) ;
	SparseCore_free_sparse (&A, Common) ;
	A = A2 ;
    }
    return (A) ;
}


/**
 * @brief 从文件中读取稠密矩阵
 * 
 */
dense_array *SparseCore_read_dense
(
    /* ---- input ---- */
    FILE *f,		/* 读取的文件，必须已经打开 */
    /* --------------- */
    sparse_common *Common
)
{
    char buf [MAXLINE+1] ;
    size_t nrow, ncol, nnz ;
    int stype, mtype ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (f, NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 读取标题和第一个数据行 */
    /* ---------------------------------------------------------------------- */

    if (!read_header (f, buf, &mtype, &nrow, &ncol, &nnz, &stype) ||
	mtype != SPARSE_DENSE)
    {
	/* 非法输入，这个函数只能读取一个稠密矩阵 */
	ERROR (SPARSE_INVALID, "invalid format") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 读取稠密矩阵 */
    /* ---------------------------------------------------------------------- */

    return (read_dense (f, nrow, ncol, stype, buf, Common)) ;
}


/**
 * @brief 	从文件中读取三连矩阵、稀疏矩阵或稠密矩阵。
 * 			返回一个指向sparse_triplet、sparse_csc
 * 			或dense_array对象的空指针。对象的类型作为mtype
 * 			参数传递回调用者。
 * 
 */
void *SparseCore_read_matrix
(
    /* ---- input ---- */
    FILE *f,		/* 要读取的文件，必须已经打开 */
    int prefer,		/* 如果为0，则稀疏矩阵总是以sparse_triplet形式返回。
				     * 它可以有任何样式(对称-下、不对称或对称-上)。
					 * 
					 * 如果1，则返回一个非对称sparse_csc形式(a ->stype == 0)的
					 * 稀疏矩阵，其中包含上三角和下三角部分。
					 * 这是MATLAB mread mexFunction所做的，因为MATLAB没有样式。
					 * 
					 * 如果是2，则返回一个样式为0或1的稀疏矩阵(不对称，或存储上半部分
					 * 的对称)。这个参数对稠密矩阵没有影响。
			 		 */
    /* ---- output---- */
    int *mtype,		/* SPARSE_TRIPLET, SPARSE_CSC, SPARSE_DENSE */
    /* --------------- */
    sparse_common *Common,
	FILE *fresult_node,
	FILE *fresult_edge,
	int graph_id
)
{
    void *G = NULL ;
    sparse_csc *A, *A2 ;
    sparse_triplet *T ;
    char buf [MAXLINE+1] ;
    size_t nrow, ncol, nnz ;
    int stype ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (f, NULL) ;
    RETURN_IF_NULL (mtype, NULL) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 读取标题以确定mtype */
    /* ---------------------------------------------------------------------- */

    if (!read_header (f, buf, mtype, &nrow, &ncol, &nnz, &stype))
    {
	/* 非法矩阵 */
	ERROR (SPARSE_INVALID, "invalid format") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 读取矩阵 */
    /* ---------------------------------------------------------------------- */

    if (*mtype == SPARSE_TRIPLET)
    {
	/* 如果 prefer == 1，读取三重矩阵转换为非对称格式 */
	T = read_triplet (f, nrow, ncol, nnz, stype, prefer == 1, buf, Common) ;
	#ifdef write_graph
	// 写入矩阵信息 到result_Node 文件中
	/* 行号  */
	//fprintf(fresult_node, "%d ", nrow );
	// 获取三元组中的 i,j, 以及真实的非零元数量real_nnz;
	double *Tx;
	Int *Ti, *Tj;
	size_t real_nnz = T->nnz;
	//printf("nnz = %ld, real_nnz = %d\n",nnz, real_nnz);
	Ti = T->i; Tj = T->j;
	Tx = T->x;
	// 初始化出入度
	Int *Rdeg, *Cdeg, max_Rdeg = -1, min_Rdeg = 999, max_Cdeg = -1, min_Cdeg = 999;
	Rdeg = (Int *)malloc(nrow * sizeof(Int));
	Cdeg = (Int *)malloc(ncol * sizeof(Int));
	for (size_t k = 0; k < nrow; k++)
	{
		Rdeg[k] = 0;
	}
	for (size_t k = 0; k < ncol; k++)
	{
		Cdeg[k] = 0;
	}
	// 写入边信息，并且计算各个结点的出入度
	for (size_t k = 0; k < real_nnz; k++)
	{
		fprintf(fresult_edge, "%d %ld %ld %lg\n",graph_id ,Ti[k], Tj[k], Tx[k] );
		if (Ti[k] == Tj[k])
			continue;
		Rdeg[ Ti[k] ]++;
		Cdeg[ Tj[k] ]++;
	}
	// 计算度信息
	double avg_Rdegree, MaxRdeg_percent, MaxCdeg_percent;
	Int Rdeg_sum = 0, MaxRdeg_sum = 0, MaxCdeg_sum = 0;
	for (size_t k = 0; k < ncol; k++)
	{	max_Rdeg = MAX(max_Rdeg, Rdeg[k]);
		min_Rdeg = MIN(min_Rdeg, Rdeg[k]);
		max_Cdeg = MAX(max_Cdeg, Cdeg[k]);
		min_Cdeg = MIN(min_Cdeg, Cdeg[k]);
		Rdeg_sum += Rdeg[k]; 
		// fprintf(fresult_node, "%d %ld %ld\n", graph_id, k, Rdeg[k]);
		// fprintf(fresult_node, "%d %ld %ld %ld \n", graph_id, k, Rdeg[k], Cdeg[k] );
	}
	// 度的频率信息
	for (size_t k = 0; k < ncol; k++)
	{
		if ( Rdeg[k]== max_Rdeg )     // 出度为0的节点比例
		{
			MaxRdeg_sum++;
		}
		if ( Cdeg[k]== max_Cdeg )
		{
			MaxCdeg_sum++;
		}
	}
	

	// 2.28 - 新增的消去属性计算 - 
	Int* Eli_add = (Int *)malloc(nrow * sizeof(Int)); // 最终的写入结果，保存消去这个节点后会新增多少边
	Int* Nei_Rstore = (Int *)malloc(max_Rdeg * sizeof(Int)); // 用节点间最大的度来开辟空间，保存 节点k 的邻居
	Int* Nei_Cstore = (Int *)malloc(max_Cdeg * sizeof(Int)); // 用节点间最大的度来开辟空间，保存 节点k 的邻居
	size_t h_c, h_r;
	for (size_t k = 0; k < ncol; k++)  // 获取各个节点的消去属性
	{
		Eli_add[k] = Cdeg[k] * Rdeg[k]; // 有向图: Cdeg * Rdeg - 
		// Eli_add[k] = Cdeg[k] * ( Cdeg[k] - 1 ) / 2; // 无向图: n*(n-1)/2 - 邻居的相互连接数量.
		h_c = 0, h_r = 0;
		for (size_t kk = 0; kk < real_nnz; kk++) // kk 表示专门找一个节点的遍历符
		{
			if ( Tj[kk] == k && Ti[kk] != k ) // 只保留(t_R,k), (k,k)不属于邻居节点, 这里存的是指向k的节点
			{
				Nei_Cstore[h_c++] = Ti[kk];   // 保存指向k的节点, 规模 Cdeg[k]
			}
			
			if ( Ti[kk] == k && Tj[kk] != k) // 只保留(k,t_c), (k,k)不属于邻居节点, 这里存的是指向k的节点
			{
				Nei_Rstore[h_r++] = Tj[kk];   // 保存k指向的节点, 规模 Rdeg[k]
			}
			// 无向图的情况
			// if (Tj[kk] == k && Ti[kk] != k) // 只保留(:,k), (k,k)不属于邻居节点
			// {
			// 	Nei_store[h] = Ti[kk];   //将 节点k的邻居保存到Nei_store.
			// 	++h;
			// }
		}

		for (size_t kk = 0; kk < Cdeg[k]; kk++)
		{
			for (size_t search_c = 0; search_c < real_nnz; search_c++)
			{
				if ( Ti[search_c] == Nei_Cstore[kk] )  // 找到(t_R, :), 且排除(t_R, t_R)
				{
					for (size_t kkp = 0; kkp < Rdeg[k]; kkp++)
					{
						if (Tj[search_c] == Nei_Rstore[kkp])   // 如果(t_R, t_c) 存在， 减一条新增边
							Eli_add[k]--;
					}
				}
				
			}
			

			// for (size_t search_s = 0; search_s < nnz; search_s++)   // 对称、一条边只需要减一次，所以只遍历一半图
			// {
			// 	if ( Tj[search_s] == Nei_store[kk] && Ti[search_s] != Nei_store[kk])  
			// 	{
			// 		for (size_t kkp = 0; kkp < Cdeg[k]; kkp++)
			// 		{
			// 			if (Ti[search_s] == Nei_store[kkp])   // 如果邻居的邻居包含在 自己的邻居里了, 消去属性减少
			// 				Eli_add[k]--;
			// 		}
					
			// 	}
			// }
		}
	}
	// 写入节点属性
	for (size_t k = 0; k < ncol; k++)
	{
		// 无向图
		// fprintf(fresult_node, "%d %ld %ld %ld\n", graph_id, k, Rdeg[k], Eli_add[k] );
		// 有向图
		fprintf(fresult_node, "%d %ld %ld %ld %ld\n", graph_id, k, Rdeg[k], Cdeg[k], Eli_add[k] );
	}
	free(Eli_add);
	free(Nei_Rstore);
	free(Nei_Cstore);
	/*------------------------------*/
	/*        记录额外信息          */
	/*------------------------------*/
	double density = (double) real_nnz / (nrow * ncol);  // 矩阵稠密度m_dense
	avg_Rdegree = (double) Rdeg_sum / nrow;                 // 图的平均度avg_d 
	MaxRdeg_percent = (double) MaxRdeg_sum / nrow;       // 出度为0的节点所占的比例 zero_per
	MaxCdeg_percent = (double) MaxCdeg_sum / nrow; 
	FILE *fresult_extinfo;
	char *extinfo_file = (char *)"./Results/QR_extinfo.txt";
	// char *extinfo_file = (char *)"./Results/Chol_extinfo.txt";
	fresult_extinfo = fopen(extinfo_file ,"a+");
	// 节点数、边数、稠密度、 平均出入度、度为max的节点所占比例。
	fprintf(fresult_extinfo, "%d %ld %ld %lg %lg %lg %lg ", 
		graph_id, nrow, real_nnz, density, avg_Rdegree, MaxRdeg_percent, MaxCdeg_percent);
	// QR 最大出度 最小出度； 最大入度， 最小入度！ 有向图
	fprintf(fresult_extinfo, "%ld %ld %ld %ld\n", max_Rdeg, min_Rdeg, max_Cdeg, min_Cdeg);
	// 读cholesky 对称正定矩阵，只需要一个度就行( 对称的QR矩阵也是只需要一个度就行.)
	// fprintf(fresult_extinfo, "%ld %ld\n", max_Rdeg, min_Rdeg);
	fclose(fresult_extinfo);
	
	free(Rdeg);
	free(Cdeg);
	#endif
	if (prefer == 0)
	{
	    /* 返回原来的三元组形式的矩阵 */
	    G = T ;
	}
	else
	{
	    /* 返回CSC格式的矩阵 */
	    A = SparseCore_triplet_to_sparse (T, 0, Common) ;
	    SparseCore_free_triplet (&T, Common) ;
	    if (A != NULL && prefer == 2 && A->stype == -1)
	    { // 对称的统一转换为上三角 
		A2 = SparseCore_transpose (A, 2, Common) ;
		SparseCore_free_sparse (&A, Common) ;
		A = A2 ;
	    }
	    *mtype = SPARSE_CSC ;
	    G = A ;
	}
    }
    else if (*mtype == SPARSE_DENSE)
    {
	/* 返回一个稠密矩阵 */
	G = read_dense (f, nrow, ncol, stype, buf, Common) ;
    }
    return (G) ;
}



/* 以mtx格式向文件写入一个矩阵 */

#define MMLEN 1024

/**
 * @brief 	读取注释文件(如果存在)，并将其复制到mtx文件中。
 * 			在每行前面加上一个“%”。如果成功返回TRUE，否则返回FALSE。
 * 
 * @param f 
 * @param comments 
 * @return int 
 */
static int include_comments (FILE *f, const char *comments)
{
    FILE *cf = NULL ;
    char buffer [MAXLINE] ;
    int ok = TRUE ;
    if (comments != NULL && comments [0] != '\0')
    {
	cf = fopen (comments, "r") ;
	if (cf == NULL)
	{
	    return (FALSE) ;
	}
	while (ok && fgets (buffer, MAXLINE, cf) != NULL)
	{
	    /* 确保行不要太长 */
	    buffer [MMLEN-1] = '\0' ;
	    buffer [MMLEN-2] = '\n' ;
	    ok = ok && (fprintf (f, "%%%s", buffer) > 0) ;
	}
	fclose (cf) ;
    }
    return (ok) ;
}


/**
 * @brief 得到矩阵中的第p个值。
 * 
 */
static void get_value
(
    double *Ax,	    /* SPARSE_COMPLEX类型的实数值或者实数的值 */
    double *Az,	    /* SPARSE_ZOMPLEX的虚数值 */
    Int p,	    	/* 得到第p个值 */
    Int xtype,	    /* A->xtype: pattern, real */
    double *x,	    /* 实部 */
    double *z	    /* 虚部 */
)
{
    switch (xtype)
    {
	case SPARSE_PATTERN:
	    *x = 1 ;
	    *z = 0 ;
	    break ;

	case SPARSE_REAL:
	    *x = Ax [p] ;
	    *z = 0 ;
	    break ;
    }
}


/**
 * @brief 将数值打印到文件中，使用最短的格式确保精确地写入值。如果成功返回TRUE，否则返回FALSE。
 * 
 */
static int print_value
(
    FILE *f,	    /* 要写入的文件 */
    double x,	    /* 写入的值 */
    Int is_integer  /* 如果写入值为整数则为true */
)
{
    double y ;
    char s [MAXLINE], *p ;
    Int i, dest = 0, src = 0 ;
    int width, ok ;

    if (is_integer)
    {
	i = (Int) x ;
	ok = (fprintf (f, ID, i) > 0) ;
	return (ok) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 处理无穷值和nan值 */
    /* ---------------------------------------------------------------------- */

    /* 将-inf改为-HUGE_DOUBLE，将+inf和nan改为+HUGE_DOUBLE */
    if (SPARSE_IS_NAN (x) || x >= HUGE_DOUBLE)
    {
	x = HUGE_DOUBLE ;
    }
    else if (x <= -HUGE_DOUBLE)
    {
	x = -HUGE_DOUBLE ;
    }

    /* ---------------------------------------------------------------------- */
    /* 找到最小的可接受的精度 */
    /* ---------------------------------------------------------------------- */

    for (width = 6 ; width < 20 ; width++)
    {
	sprintf (s, "%.*g", width, x) ;
	sscanf (s, "%lg", &y) ;
	if (x == y) break ;
    }

    /* ---------------------------------------------------------------------- */
    /* 缩短字符串 */
    /* ---------------------------------------------------------------------- */

    /* 把“e+0”换成“e”，把“e+”换成“e”，把“e-”换成“e-” */
    for (i = 0 ; i < MAXLINE && s [i] != '\0' ; i++)
    {
	if (s [i] == 'e')
	{
	    if (s [i+1] == '+')
	    {
		dest = i+1 ;
		if (s [i+2] == '0')
		{
		    /* 删除字符s[i+1]和s[i+2] */
		    src = i+3 ;
		}
		else
		{
		    /* 删除字符s[i+1] */
		    src = i+2 ;
		}
	    }
	    else if (s [i+1] == '-')
	    {
		dest = i+2 ;
		if (s [i+2] == '0')
		{
		    /* 删除字符s[i+2] */
		    src = i+3 ;
		}
		else
		{
		    /* 无变化 */
		    break ;
		}
	    }
	    while (s [src] != '\0')
	    {
		s [dest++] = s [src++] ;
	    }
	    s [dest] = '\0' ;
	    break ;
	}
    }

    /* 如果没有必要，请删除前导“0” */
    p = s ;
    s [MAXLINE-1] = '\0' ;
    i = strlen (s) ;
    if (i > 2 && s [0] == '0' && s [1] == '.')
    {
	/* "0.x"变为".x" */
	p = s + 1 ;
    }
    else if (i > 3 && s [0] == '-' && s [1] == '0' && s [2] == '.')
    {
	/* "-0.x"变为"-.x" */
	s [1] = '-' ;
	p = s + 1 ;
    }

#if 0
    /* 重复检查 */
    i = sscanf (p, "%lg", &z) ;
    if (i != 1 || y != z)
    {
	/* 在上面的“e+0”编辑中出现了一些错误。*/
	sprintf (s, "%.*g", width, x) ;
	p = s ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* 将值写入到文件中 */
    /* ---------------------------------------------------------------------- */

    ok = (fprintf (f, "%s", p) > 0) ;
    return (ok) ;
}


/**
 * @brief 打印一个三元组，将其转换为基于一的。如果成功返回TRUE，否则返回FALSE
 * 
 */
static int print_triplet
(
    FILE *f,		/* 写入的文件 */
    Int is_binary,	/* 文件为时"pattern"为TRUE */
    Int is_integer,	/* 未见为时"integer"为TRUE */
    Int i,			/* 行索引 0开始 */
    Int j,			/* 列索引 0开始 */
    double x,		/* 实部 */
    double z		/* 虚部 */
)
{
    int ok ; 
    ok = (fprintf (f, ID " " ID, 1+i, 1+j) > 0) ;
    if (!is_binary)
    {
	fprintf (f, " ") ;
	ok = ok && print_value (f, x, is_integer) ;
    }
    ok = ok && (fprintf (f, "\n") > 0) ;
    return (ok) ;
}


/**
 * @brief 计算将从矩阵A打印到文件的三元组的数目
 * 
 */
static Int ntriplets
(
    sparse_csc *A,	    /* 将被写入的矩阵 */
    Int is_sym		    	/* 如果文件是对称的，为真(仅下三角)*/
)
{
    Int *Ap, *Ai, *Anz, packed, i, j, p, pend, ncol, stype, nz = 0 ;
    if (A == NULL)
    {
	/* 矩阵Z为空 */
	return (0) ;
    }
    stype = A->stype ;
    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;
    ncol = A->ncol ;
    for (j = 0 ; j < ncol ; j++)
    {
	p = Ap [j] ;
	pend = (packed) ? Ap [j+1] : p + Anz [j] ;
	for ( ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if ((stype < 0 && i >= j) || (stype == 0 && (i >= j || !is_sym)))
	    {
		/* HNUCHOL矩阵下三角对称(文件也是如此);
		 * 或者HNUCHOL矩阵不对称，要么A(i,j)在下面，要么文件不对称。 
		 */
		nz++ ;
	    }
	    else if (stype > 0 && i <= j)
	    {
		/* HNUCHOL矩阵是上三角对称，但是文件是下三角对称。需要转置这个元素。 */
		nz++ ;
	    }
	}
    }
    return (nz) ;
}


/**
 * @brief 用mtx格式将一个稀疏矩阵写入文件
 * 
 */
int SparseCore_write_sparse
(
    /* ---- input ---- */
    FILE *f,		    	/* 要写入的文件，必须已经打开 */
    sparse_csc *A,	    /* 要写入的矩阵 */
    sparse_csc *Z,	    /* 可选的矩阵与明确的零的模式 */
    const char *comments,   /* 要包含的注释的可选文件名 */
    /* --------------- */
    sparse_common *Common
)
{
    double x = 0, z = 0 ;
    double *Ax, *Az ;
    Int *Ap, *Ai, *Anz, *Zp, *Zi, *Znz ;
    Int nrow, ncol, symmetry, i, j, q, iz, p, nz, is_binary, stype,
	is_integer, asym, is_sym, xtype, apacked, zpacked, pend, qend, zsym ;
    int ok ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (f, EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, SPARSE_PATTERN, SPARSE_REAL, EMPTY) ;
    if (Z != NULL && (Z->nrow == 0 || Z->ncol == 0))
    {
	/* Z是非空的，但是是空的，所以把它当作一个空矩阵 */
	Z = NULL ;
    }
    if (Z != NULL)
    {
	RETURN_IF_XTYPE_INVALID (Z, SPARSE_PATTERN, SPARSE_REAL, EMPTY) ;
	if (Z->nrow != A->nrow || Z->ncol != A->ncol || Z->stype != A->stype)
	{
	    ERROR (SPARSE_INVALID, "dimension or type of A and Z mismatch") ;
	    return (EMPTY) ;
	}
    }
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 得到矩阵A */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    xtype = A->xtype ;
    apacked = A->packed ;

    if (xtype == SPARSE_PATTERN)
    {
	/* HNUCHOL模式矩阵作为“pattern”在文件中写入 */
	is_binary = TRUE ;
	is_integer = FALSE ;
    }
    else if (xtype == SPARSE_REAL)
    {
	/* 确定一个实矩阵实际上是二进制还是整数 */
	is_binary = TRUE ;
	is_integer = TRUE ;
	for (j = 0 ; (is_binary || is_integer) && j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (apacked) ? Ap [j+1] : p + Anz [j] ;
	    for ( ; (is_binary || is_integer) && p < pend ; p++)
	    {
		x = Ax [p] ;
		if (x != 1)
		{
		    is_binary = FALSE ;
		}
		/* 转换为Int，然后返回到double */
		i = (Int) x ;
		z = (double) i ;
		if (z != x)
		{
		    is_integer = FALSE ;
		}
	    }
	}
    }
    else
    {
	is_binary = FALSE ;
	is_integer = FALSE ;
    }

    /* ---------------------------------------------------------------------- */
    /* 得到Z矩阵(只考虑pattern) */
    /* ---------------------------------------------------------------------- */

    Zp = NULL ;
    Zi = NULL ;
    Znz = NULL ;
    zpacked = TRUE ;
    if (Z != NULL)
    {
	Zp = Z->p ;
	Zi = Z->i ;
	Znz = Z->nz ;
	zpacked = Z->packed ;
    }

    /* ---------------------------------------------------------------------- */
    /* 确定A和Z的对称性 */
    /* ---------------------------------------------------------------------- */

    stype = A->stype ;
    if (A->nrow != A->ncol)
    {
	asym = SPARSE_MM_RECTANGULAR ;
    }
    else if (stype != 0)
    {
	/* HNUCHOL的A和Z矩阵具有对称(和匹配)的样式。注意，对角线是不勾选的。 */
	asym = SPARSE_MM_SYMMETRIC ;
    }
    else if (!A->sorted)
    {
	/* A是不对称存储，但未排序 */
	asym = SPARSE_MM_UNSYMMETRIC ;
    }
    else
    {
	/* HNUCHOL的stype为零(以不对称形式存储) */
	asym = EMPTY ;
	zsym = EMPTY ;

	if (asym == EMPTY || zsym <= SPARSE_MM_UNSYMMETRIC)
	{
	    /* 没有计算，内存不足，或者Z不对称 */
	    asym = SPARSE_MM_UNSYMMETRIC ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* 写mtx文件抬头说明 */
    /* ---------------------------------------------------------------------- */

    ok = fprintf (f, "%%%%MatrixMarket matrix coordinate") > 0 ;

    if (is_binary)
    {
	ok = ok && (fprintf (f, " pattern") > 0) ;
    }
    else if (is_integer)
    {
	ok = ok && (fprintf (f, " integer") > 0) ;
    }
    else
    {
	ok = ok && (fprintf (f, " real") > 0) ;
    }

    is_sym = FALSE ;

    switch (asym)
    {
	case SPARSE_MM_RECTANGULAR:
	case SPARSE_MM_UNSYMMETRIC:
	    /* A是矩形的或不对称的 */
	    ok = ok && (fprintf (f, " general\n") > 0) ;
	    is_sym = FALSE ;
	    symmetry = SPARSE_MM_UNSYMMETRIC ;
	    break ;

	case SPARSE_MM_SYMMETRIC:
	case SPARSE_MM_SYMMETRIC_POSDIAG:
	    /* A是对称的 */
	    ok = ok && (fprintf (f, " symmetric\n") > 0) ;
	    is_sym = TRUE ;
	    symmetry = SPARSE_MM_SYMMETRIC ;
	    break ;

	case SPARSE_MM_HERMITIAN:
	case SPARSE_MM_HERMITIAN_POSDIAG:
	    /* A是Hermitian矩阵 */
	    ok = ok && (fprintf (f, " Hermitian\n") > 0) ;
	    is_sym = TRUE ;
	    symmetry = SPARSE_MM_HERMITIAN ;
	    break ;

	case SPARSE_MM_SKEW_SYMMETRIC:
	    /* A是斜对称的 */
	    ok = ok && (fprintf (f, " skew-symmetric\n") > 0) ;
	    is_sym = TRUE ;
	    symmetry = SPARSE_MM_SKEW_SYMMETRIC ;
	    break ;
    }

    /* ---------------------------------------------------------------------- */
    /* 如果有的话，包括评论 */
    /* ---------------------------------------------------------------------- */

    ok = ok && include_comments (f, comments) ;

    /* ---------------------------------------------------------------------- */
    /* 写入稀疏矩阵(A和Z) */
    /* ---------------------------------------------------------------------- */

    nz = ntriplets (A, is_sym) + ntriplets (Z, is_sym) ;

    /* 使用三元组的nrow、ncol和#写入第一行数据 */
    ok = ok && (fprintf (f, ID " " ID " " ID "\n", nrow, ncol, nz) > 0) ;

    for (j = 0 ; ok && j < ncol ; j++)
    {
	/* 合并A和Z的列 */
	p = Ap [j] ;
	pend = (apacked) ? Ap [j+1] : p + Anz [j] ;
	q = (Z == NULL) ? 0 : Zp [j] ;
	qend = (Z == NULL) ? 0 : ((zpacked) ? Zp [j+1] : q + Znz [j]) ;
	while (ok)
	{
	    /* 从A和Z中获取下一行索引 */
	    i  = (p < pend) ? Ai [p] : (nrow+1) ;
	    iz = (q < qend) ? Zi [q] : (nrow+2) ;
	    if (i <= iz)
	    {
		/* 得到A(i,j),或者在A和Z都完成的时候退出 */
		if (i == nrow+1) break ;
		get_value (Ax, Az, p, xtype, &x, &z) ;
		p++ ;
	    }
	    else
	    {
		/* 得到Z(i,j) */
		i = iz ;
		x = 0 ;
		z = 0 ;
		q++ ;
	    }
	    if ((stype < 0 && i >= j) || (stype == 0 && (i >= j || !is_sym)))
	    {
		/* HNUCHOL矩阵对称下三角(文件也是如此);
		 * 或者HNUCHOL矩阵不对称，要么A(i,j)在下三角，要么文件不对称。 
		 */
		ok = ok && print_triplet (f, is_binary, is_integer,
		    i,j, x,z) ;
	    }
	    else if (stype > 0 && i <= j)
	    {
		/* HNUCHOL矩阵是对称上三角的，但是文件是对称下三角的。
		 * 如果矩阵是实的，需要转置这个元素。
		 */
		if (z != 0)
		{
		    z = -z ;
		}
		ok = ok && print_triplet (f, is_binary, is_integer,
		    j,i, x,z) ;
	    }
	}
    }

    if (!ok)
    {
	ERROR (SPARSE_INVALID, "error reading/writing file") ;
	return (EMPTY) ;
    }

    return (asym) ;
}


/**
 * @brief 以mtx格式将稠密矩阵写入文件。
 * 
 */
int SparseCore_write_dense
(
    /* ---- input ---- */
    FILE *f,		    	/* 要写入的文件，必须已经打开 */
    dense_array *X,	    /* 要被写入的矩阵 */
    const char *comments,   /* 要包含的注释的可选文件名 */
    /* --------------- */
    sparse_common *Common
)
{
    double x = 0, z = 0 ;
    double *Xx, *Xz ;
    Int nrow, ncol, i, j, xtype, p ;
    int ok ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入 */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (f, EMPTY) ;
    RETURN_IF_NULL (X, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (X, SPARSE_REAL, SPARSE_REAL, EMPTY) ;
    Common->status = SPARSE_OK ;

    /* ---------------------------------------------------------------------- */
    /* 得到矩阵X */
    /* ---------------------------------------------------------------------- */

    Xx = X->x ;
    Xz = X->z ;
    nrow = X->nrow ;
    ncol = X->ncol ;
    xtype = X->xtype ;

    /* ---------------------------------------------------------------------- */
    /* 写mtx抬头说明 */
    /* ---------------------------------------------------------------------- */

    ok = (fprintf (f, "%%%%MatrixMarket matrix array") > 0) ;
	ok = ok && (fprintf (f, " real general\n") > 0) ;

    /* ---------------------------------------------------------------------- */
    /* 如果有的话，包括注释 */
    /* ---------------------------------------------------------------------- */

    ok = ok && include_comments (f, comments) ;

    /* ---------------------------------------------------------------------- */
    /* 写入稠密矩阵*/
    /* ---------------------------------------------------------------------- */

    /* 用nrow和ncol写入第一行数据 */
    ok = ok && (fprintf (f, ID " " ID "\n", nrow, ncol) > 0) ;

    Xx = X->x ;
    Xz = X->z ;
    for (j = 0 ; ok && j < ncol ; j++)
    {
	for (i = 0 ; ok && i < nrow ; i++)
	{
	    p = i + j*nrow ;
	    get_value (Xx, Xz, p, xtype, &x, &z) ;
	    ok = ok && print_value (f, x, FALSE) ;
	    ok = ok && (fprintf (f, "\n") > 0) ;
	}
    }

    if (!ok)
    {
	ERROR (SPARSE_INVALID, "error reading/writing file") ;
	return (EMPTY) ;
    }

    return ((nrow == ncol) ? SPARSE_MM_UNSYMMETRIC : SPARSE_MM_RECTANGULAR) ;
}
#endif
