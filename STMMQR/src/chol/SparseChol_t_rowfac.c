/**
 * @file t_SparseChol_rowfac.c
 * @author your name (you@domain.com)
 * @brief 	SparseChol_rowfac模板历程
 * @version 0.1
 * @date 2020-09-22
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Sparse_template.h"

#ifdef MASK
static int TEMPLATE (SparseChol_rowfac_mask)
#else
static int TEMPLATE (SparseChol_rowfac)
#endif
(
    /* ---- input ---- */
    sparse_csc *A,	/* 分解的稀疏矩阵 */
    sparse_csc *F,	/* 仅用于A*A'. F=A'或者A(:,f)' */
    double beta [2],	/* 分解beta*I+A或者beta*I+AA' (仅beta [0]) */
    size_t kstart,		/* 分解的第一行 */
    size_t kend,		/* 分解的最后一行为kend-1 */
#ifdef MASK
    /* 这些输入仅用于SparseChol_rowfac_mask */
    Int *mask,		/* 大小为A->nrow. 如果mask[i] >= maskmark，W(i)置为0 */
    Int maskmark,
    Int *RLinkUp,	/* 大小为A->nrow. 链接要计算的行列表 */
#endif
    /* ---- in/out --- */
    sparse_factor *L,
    /* --------------- */
    sparse_common *Common
)
{
    double yx [2], lx [2], fx [2], dk [1], di [1], fl = 0 ;
    double *Ax, *Az, *Lx, *Lz, *Wx, *Wz, *Fx, *Fz ;
    Int *Ap, *Anz, *Ai, *Lp, *Lnz, *Li, *Lnext, *Flag, *Stack, *Fp, *Fi, *Fnz,
	*Iwork ;
    Int i, p, k, t, pf, pfend, top, s, mark, pend, n, lnz, is_ll, multadds,
	use_dbound, packed, stype, Fpacked, sorted, nzmax, len, parent ;
#ifndef REAL
    Int dk_imaginary ;
#endif

    /* ---------------------------------------------------------------------- */
    /* 得到输入 */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    stype = A->stype ;

    if (stype > 0)
    {
	/* 对称上三角情况: 不需要F.它可能为空 */
	Fp = NULL ;
	Fi = NULL ;
	Fx = NULL ;
	Fz = NULL ;
	Fnz = NULL ;
	Fpacked = TRUE ;
    }
    else
    {
	/* 非对称情况: 需要F */
	Fp = F->p ;
	Fi = F->i ;
	Fx = F->x ;
	Fz = F->z ;
	Fnz = F->nz ;
	Fpacked = F->packed ;
    }

    Ap = A->p ;		/* A的列指针 大小为A->ncol+1 */
    Ai = A->i ;		/* A的行索引，大小为nz = Ap [A->ncol] */
    Ax = A->x ;		/* A的数值，大小为nz */
    Az = A->z ;
    Anz = A->nz ;
    packed = A->packed ;
    sorted = A->sorted ;

    use_dbound = IS_GT_ZERO (Common->dbound) ;

    /* 得到当前因子L(和LDL'的D);如果需要的话，再分配空间 */
    is_ll = L->is_ll ;
    if (L->xtype == SPARSE_PATTERN)
    {
		/* ------------------------------------------------------------------ */
		/* L只是符号，指定并初始化L (以及LDL'的D) */
		/* ------------------------------------------------------------------ */

		/* 不需要额外的工作空间 */
		SparseCore_change_factor (A->xtype, is_ll, FALSE, FALSE, TRUE, L, Common);
		if (Common->status < SPARSE_OK)
		{
			/* 内存溢出 */
			return (FALSE) ;
		}
    }
    else if (kstart == 0 && kend == (size_t) n)
    {
	/* ------------------------------------------------------------------ */
	/* 重置L->nz和L->minor重新进行因数分解 */
	/* ------------------------------------------------------------------ */

		L->minor = n ;
		Lnz = L->nz ;
		for (k = 0 ; k < n ; k++)
		{
			Lnz [k] = 1 ;
		}
    }

    /* 输入，可以在输出修改: */
    Lp = L->p ;		/* 大小为n+1 */

    /* 输出，仅对增量情况下的输入定义的内容: */
    Lnz = L->nz ;		/* 大小为n */
    Lnext = L->next ;	/* 大小为n+2 */
    Li = L->i ;			/* 大小为L->nzmax, 能够修改 */
    Lx = L->x ;			/* 大小为L->nzmax or 2*L->nzmax, 能够修改 */
    Lz = L->z ;			/* 大小为L->nzmax for zomplex case, 能够修改 */
    nzmax = L->nzmax ;

    /* ---------------------------------------------------------------------- */
    /* 获取工作空间e */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Stack = Iwork ;			/* 大小为n (i/i/l) */
    Flag = Common->Flag ;	/* 大小为n, 必须保持Flag [i] < mark */
    Wx = Common->Xwork ;	/* 大小为n， 必须保持Xwork [i] == 0 */
    Wz = Wx + n ;			/* 仅在复数情况下大小为n */
    mark = Common->mark ;

    /* ---------------------------------------------------------------------- */
    /* 按行计算LDL'或LL'因数分解 */
    /* ---------------------------------------------------------------------- */

#ifdef MASK
#define NEXT(k) k = RLinkUp [k]
#else
#define NEXT(k) k++
#endif

    for (k = kstart ; k < ((Int) kend) ; NEXT(k))
    {

		/* ------------------------------------------------------------------ */
		/* 计算L的第k行pattern，分散第k个输入列 */
		/* ------------------------------------------------------------------ */

		/* L的第k列目前为空 */

		top = n ;			/* 栈为空 */
		Flag [k] = mark ;	/* 堆栈中不包括对角条目 */

		/* 使用Li [Lp [i]+1]表示etree */
	#define PARENT(i) (Lnz [i] > 1) ? (Li [Lp [i] + 1]) : EMPTY

		if (stype > 0)
		{
			/* 分发triu(beta*I+AA')的第k列, 得到pattern L(k,:) */
			p = Ap [k] ;
			pend = (packed) ? (Ap [k+1]) : (p + Anz [k]) ;
			/* W [i] = Ax [i] ; scatter column of A */
		#define SCATTER ASSIGN(Wx,Wz,i, Ax,Az,p)
			SUBTREE ;
		#undef SCATTER
		}
		else
		{
			/* 分发triu (beta*I+AA')的第k列, 得到pattern L(k,:) */
			pf = Fp [k] ;
			pfend = (Fpacked) ? (Fp [k+1]) : (pf + Fnz [k]) ;
			for ( ; pf < pfend ; pf++)
			{
				/* 获取F (t,k)的非零元 */
				t = Fi [pf] ;
				/* fk = Fx [pf] */
				ASSIGN (fx, fz, 0, Fx, Fz, pf) ;
				p = Ap [t] ;
				pend = (packed) ? (Ap [t+1]) : (p + Anz [t]) ;
				multadds = 0 ;
				/* W [i] += Ax [p] * fx ; 分发A*A'的列 */
			#define SCATTER MULTADD (Wx,Wz,i, Ax,Az,p, fx,fz,0) ; multadds++  ;
				SUBTREE ;
			#undef SCATTER
				fl += 2 * ((double) multadds) ;
			}
		}

	#undef PARENT

	/* ------------------------------------------------------------------ */
	/* 如果定义mask，则将W中的对应项设置为零 */
	/* ------------------------------------------------------------------ */

	#ifdef MASK
			/* 删除Wx的无用项 */
			if (mask != NULL)
			{

	#if 0
				for (p = n; p > top;)
				{
					i = Stack [--p] ;
					if ( mask [i] >= 0 )
			{
				CLEAR (Wx,Wz,i) ;	/* 将W(i)置为0 */
			}
				}
	#endif

				for (s = top ; s < n ; s++)
				{
					i = Stack [s] ;
					if (mask [i] >= maskmark)
			{
				CLEAR (Wx,Wz,i) ;	/* 将W(i)置为0 */
			}
				}

			}
	#endif

	/* L的第k行非零pattern现在在Stack[top..n-1]中.
	 * Flag [Stack [top..n-1]]==mark,但不再需要 */

	/* mark = SparseCore_clear_flag (Common) ; */
	SPARSE_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	/* ------------------------------------------------------------------ */
	/* 计算L的第k行，并以列的形式存储 */
	/* ------------------------------------------------------------------ */

	/* 求解L (0:k-1, 0:k-1) * y (0:k-1) = b (0:k-1) 
	 * 其中b (0:k) = A (0:k,k)或者A(0:k,:) * F(:,k) 在W和Stack中.
	 *
	 * LDL':
	 * L (k, 0:k-1) = y (0:k-1) ./ D (0:k-1)
	 * D (k) = b (k) - L (k, 0:k-1) * y (0:k-1)
	 *
	 * LL':
	 * L (k, 0:k-1) = y (0:k-1)
	 * L (k,k) = sqrt (b (k) - L (k, 0:k-1) * L (0:k-1, k))
	 */

	/* dk = W [k] + beta */
	ADD_REAL (dk,0, Wx,k, beta,0) ;

#ifndef REAL
	/* 在非对称情况下，W[k]的虚部必须是实的，因为假设F是a的复共轭转置。
	 * 在对称情况下，W[k]是a的对角，如果W[k]的虚部非零，则无法计算Cholesky分解;
	 * A不是正定的 */
	dk_imaginary = (stype > 0) ? (IMAG_IS_NONZERO (Wx,Wz,k)) : FALSE ;
#endif

	/* W [k] = 0.0 ; */
	CLEAR (Wx,Wz,k) ;

	for (s = top ; s < n ; s++)
	{
	    /* 对每个非零项L(k,i)取i */
	    i = Stack [s] ;

	    /* y = W [i] ; */
	    ASSIGN (yx,yz,0, Wx,Wz,i) ;

	    /* W [i] = 0.0 ; */
	    CLEAR (Wx,Wz,i) ;

	    lnz = Lnz [i] ;
	    p = Lp [i] ;
	    pend = p + lnz ;

	    /* di = Lx [p] ; 实数对角项L或D(i,i) */
	    ASSIGN_REAL (di,0, Lx,p) ;

	    if (i >= (Int) L->minor || IS_ZERO (di [0]))
	    {
			/* LL'：L(i,i)为0.  
			* LDL'：D(i,i)为0.  
			* 跳过L的第i列,并且设置L(k,i) = 0. */
			CLEAR (lx,lz,0) ;
			p = pend ;
	    }
	    else if (is_ll)
	    {
			fl += 2 * ((double) (pend - p - 1)) + 3 ;

			/* 使用L (i:(k-1),i)向前求解 */
			/* 除以L(i,i)它必须是实数且非零 */
			/* y /= di [0] */
			DIV_REAL (yx,yz,0, yx,yz,0, di,0) ;
			for (p++ ; p < pend ; p++)
			{
				/* W [Li [p]] -= Lx [p] * y ; */
				MULTSUB (Wx,Wz,Li[p], Lx,Lz,p, yx,yz,0) ;
			}
			/* 不要缩放L;计算L(k,k)的点积 */
			/* L(k,i) = conj(y) ; */
			ASSIGN_CONJ (lx,lz,0, yx,yz,0) ;
			/* d -= conj(y) * y ; */
			LLDOT (dk,0, yx,yz,0) ;
	    }
	    else
	    {

			fl += 2 * ((double) (pend - p - 1)) + 3 ;

			/* 使用D (i,i)和L ((i+1):(k-1),i)向前求解 */
			for (p++ ; p < pend ; p++)
			{
				/* W [Li [p]] -= Lx [p] * y ; */
				MULTSUB (Wx,Wz,Li[p], Lx,Lz,p, yx,yz,0) ;
			}
			/* 为LDL'分解缩放L (k,0:k-1)计算D (k,k)*/

			/* L(k,i) = y/d */
			lx [0] = yx [0] / di [0] ;
			/* d -= L(k,i) * y */
			dk [0] -= lx [0] * yx [0] ;

	    }

	    /* 确定L的第i列是否可以保存新的L(k,i)项 */
	    if (p >= Lp [Lnext [i]])
	    {
			/* 第i列需要grow */
			if (!SparseCore_reallocate_column (i, lnz + 1, L, Common))
			{
				/* 内存溢出，L现在为简单符号 */
				for (i = 0 ; i < n ; i++)
				{
					/* W [i] = 0 ; */
					CLEAR (Wx,Wz,i) ;
				}
				return (FALSE) ;
			}
			Li = L->i ;		/* L->i, L->x, L->z可能移动了*/
			Lx = L->x ;
			Lz = L->z ;
			p = Lp [i] + lnz ;	/* L->p发生变化 */
	    }

	    /* 将L (k,i)存储在列存储的L矩阵中 */
	    Li [p] = k ;
	    /* Lx [p] = L(k,i) ; */
	    ASSIGN (Lx,Lz,p, lx,lz,0) ;
	    Lnz [i]++ ;
	}

	/* ------------------------------------------------------------------ */
	/* 如果给出了dbound，确保abs (d) >= dbound并且将其存入L中 */
	/* ------------------------------------------------------------------ */

	p = Lp [k] ;
	Li [p] = k ;

	if (k >= (Int) L->minor)
	{
	    /* 矩阵已经不是正定的了 */
	    dk [0] = 0 ;
	}
	else if (use_dbound)
	{
	    /* 修改对角线使LL'或LDL'存在 */
	    dk [0] = SparseCore_dbound (is_ll ? fabs (dk [0]) : dk [0], Common) ;
	}
	else if ((is_ll ? (IS_LE_ZERO (dk [0])) : (IS_ZERO (dk [0])))
#ifndef REAL
		|| dk_imaginary
#endif
		)
	{
	    /* 这个矩阵不是正定的 */
	    dk [0] = 0 ;
	    L->minor = k ;
	    ERROR (SPARSE_NOT_POSDEF, "not positive definite") ;
	}

	if (is_ll)
	{
	    /* 视为以此计算 */
	    dk [0] = sqrt (dk [0]) ;
	}

	/* Lx [p] = D(k,k) = d 实部 */
	ASSIGN_REAL (Lx,p, dk,0) ;
	CLEAR_IMAG (Lx,Lz,p) ;
    }

#undef NEXT

    if (is_ll) fl += MAX ((Int) kend - (Int) kstart, 0) ;  
    Common->rowfacfl = fl ;
    return (TRUE) ;
}
#undef PATTERN
#undef REAL
