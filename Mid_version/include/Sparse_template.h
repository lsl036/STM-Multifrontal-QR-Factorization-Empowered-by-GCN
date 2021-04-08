/* ========================================================================== */
/* === Include/Sparse_template.h =========================================== */
/* ========================================================================== */

#undef TEMPLATE
#undef TEMPLATE2
#undef XTYPE
#undef XTYPE2
#undef XTYPE_OK
#undef ENTRY_IS_NONZERO
#undef ENTRY_IS_ZERO
#undef ENTRY_IS_ONE
#undef IMAG_IS_NONZERO

#undef ASSEMBLE
#undef ASSIGN
#undef ASSIGN_CONJ
#undef ASSIGN2
#undef ASSIGN2_CONJ
#undef ASSIGN_REAL
#undef MULT
#undef MULTADD
#undef ADD
#undef ADD_REAL
#undef MULTSUB
#undef MULTADDCONJ
#undef MULTSUBCONJ
#undef LLDOT
#undef CLEAR
#undef DIV
#undef DIV_REAL
#undef MULT_REAL
#undef CLEAR_IMAG
#undef LDLDOT
#undef PREFIX

#undef ENTRY_SIZE

#undef XPRINT0
#undef XPRINT1
#undef XPRINT2
#undef XPRINT3

/* -------------------------------------------------------------------------- */
/* pattern */
/* -------------------------------------------------------------------------- */


#ifdef PATTERN

#define PREFIX				    p_
#define TEMPLATE(name)			    P_TEMPLATE(name)
#define TEMPLATE2(name)			    P_TEMPLATE(name)
#define XTYPE				    SPARSE_PATTERN
#define XTYPE2				    SPARSE_REAL
#define XTYPE_OK(type)			    (TRUE)
#define ENTRY_IS_NONZERO(ax,az,q)	    (TRUE)
#define ENTRY_IS_ZERO(ax,az,q)		    (FALSE)
#define ENTRY_IS_ONE(ax,az,q)		    (TRUE)
#define IMAG_IS_NONZERO(ax,az,q)	    (FALSE)
#define ENTRY_SIZE			    0

#define ASSEMBLE(x,z,p,ax,az,q)
#define ASSIGN(x,z,p,ax,az,q)
#define ASSIGN_CONJ(x,z,p,ax,az,q)
#define ASSIGN2(x,z,p,ax,az,q)		    P_ASSIGN2(x,z,p,ax,az,q)
#define ASSIGN2_CONJ(x,z,p,ax,az,q)	    P_ASSIGN2(x,z,p,ax,az,q)
#define ASSIGN_REAL(x,p,ax,q)
#define MULT(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD_REAL(x,p, ax,q, bx,r)
#define MULTSUB(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define LLDOT(x,p,ax,az,q)
#define CLEAR(x,z,p)
#define CLEAR_IMAG(x,z,p)
#define DIV(x,z,p,ax,az,q)
#define DIV_REAL(x,z,p, ax,az,q, bx,r)
#define MULT_REAL(x,z,p, ax,az,q, bx,r)
#define LDLDOT(x,p, ax,az,q, bx,r)

#define XPRINT0(x,z,p)			    P_PRINT(0,x,z,p)
#define XPRINT1(x,z,p)			    P_PRINT(1,x,z,p)
#define XPRINT2(x,z,p)			    P_PRINT(2,x,z,p)
#define XPRINT3(x,z,p)			    P_PRINT(3,x,z,p)

/* -------------------------------------------------------------------------- */
/* real */
/* -------------------------------------------------------------------------- */

#elif defined (REAL)

#define PREFIX				    r_
#define TEMPLATE(name)			    R_TEMPLATE(name)
#define TEMPLATE2(name)			    R_TEMPLATE(name)
#define XTYPE				    SPARSE_REAL
#define XTYPE2				    SPARSE_REAL
#define XTYPE_OK(type)			    R_XTYPE_OK(type)
#define ENTRY_IS_NONZERO(ax,az,q)	    R_IS_NONZERO(ax,az,q)
#define ENTRY_IS_ZERO(ax,az,q)		    R_IS_ZERO(ax,az,q)
#define ENTRY_IS_ONE(ax,az,q)		    R_IS_ONE(ax,az,q)
#define IMAG_IS_NONZERO(ax,az,q)	    (FALSE)
#define ENTRY_SIZE			    1

#define ASSEMBLE(x,z,p,ax,az,q)		    R_ASSEMBLE(x,z,p,ax,az,q) 
#define ASSIGN(x,z,p,ax,az,q)		    R_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN_CONJ(x,z,p,ax,az,q)	    R_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN2(x,z,p,ax,az,q)		    R_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN2_CONJ(x,z,p,ax,az,q)	    R_ASSIGN(x,z,p,ax,az,q)
#define ASSIGN_REAL(x,p,ax,q)		    R_ASSIGN_REAL(x,p,ax,q)
#define MULT(x,z,p,ax,az,q,bx,bz,pb)	    R_MULT(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADD(x,z,p,ax,az,q,bx,bz,pb)     R_MULTADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD(x,z,p,ax,az,q,bx,bz,pb)	    R_ADD(x,z,p,ax,az,q,bx,bz,pb)
#define ADD_REAL(x,p, ax,q, bx,r)	    R_ADD_REAL(x,p, ax,q, bx,r)
#define MULTSUB(x,z,p,ax,az,q,bx,bz,pb)     R_MULTSUB(x,z,p,ax,az,q,bx,bz,pb)
#define MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb) \
    R_MULTADDCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb) \
    R_MULTSUBCONJ(x,z,p,ax,az,q,bx,bz,pb)
#define LLDOT(x,p,ax,az,q)		    R_LLDOT(x,p,ax,az,q)
#define CLEAR(x,z,p)			    R_CLEAR(x,z,p) 
#define CLEAR_IMAG(x,z,p)		    R_CLEAR_IMAG(x,z,p) 
#define DIV(x,z,p,ax,az,q)		    R_DIV(x,z,p,ax,az,q)
#define DIV_REAL(x,z,p, ax,az,q, bx,r)	    R_DIV_REAL(x,z,p, ax,az,q, bx,r)
#define MULT_REAL(x,z,p, ax,az,q, bx,r)	    R_MULT_REAL(x,z,p, ax,az,q, bx,r)
#define LDLDOT(x,p, ax,az,q, bx,r)	    R_LDLDOT(x,p, ax,az,q, bx,r)

#define XPRINT0(x,z,p)			    R_PRINT(0,x,z,p)
#define XPRINT1(x,z,p)			    R_PRINT(1,x,z,p)
#define XPRINT2(x,z,p)			    R_PRINT(2,x,z,p)
#define XPRINT3(x,z,p)			    R_PRINT(3,x,z,p)

#endif
