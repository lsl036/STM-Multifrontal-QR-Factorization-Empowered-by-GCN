/* ========================================================================== */
/* === SparseBase_config.h ==================================================== */
/* ========================================================================== */

#ifndef SparseBase_CONFIG_H
#define SparseBase_CONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <limits.h>
#include <stdlib.h>

/* ========================================================================== */
/* === Sparse_long ===================================================== */
/* ========================================================================== */

#ifndef Sparse_long

#ifdef _WIN64

#define Sparse_long __int64
#define Sparse_long_max _I64_MAX
#define Sparse_long_idd "I64d"

#else

#define Sparse_long long
#define Sparse_long_max LONG_MAX
#define Sparse_long_idd "ld"

#endif
#define Sparse_long_id "%" Sparse_long_idd
#endif

struct SparseBase_config_struct
{
    void *(*malloc_func) (size_t) ;            
    void *(*calloc_func) (size_t, size_t) ;     
    void *(*realloc_func) (void *, size_t) ;    
    void (*free_func) (void *) ;               
    int (*printf_func) (const char *, ...) ;   
    double (*hypot_func) (double, double) ;    
    int (*divcomplex_func) (double, double, double, double, double *, double *);
} ;

extern struct SparseBase_config_struct SparseBase_config ;

void SparseBase_start ( void ) ;   

void SparseBase_finish ( void ) ;  

void *SparseBase_malloc    
(
    size_t nitems,          
    size_t size_of_item     
) ;

void *SparseBase_calloc    
(
    size_t nitems,          
    size_t size_of_item     
) ;

void *SparseBase_realloc   
(
    size_t nitems_new,      
    size_t nitems_old,      
    size_t size_of_item,    
    void *p,                
    int *ok                 
) ;

void *SparseBase_free      
(
    void *p                 
) ;

void SparseBase_tic    
(
    double tic [2]      
) ;

double SparseBase_toc  
(
    double tic [2]      
) ;

double SparseBase_time  
(
    void
) ;


double SparseBase_hypot (double x, double y) ;


int SparseBase_divcomplex
(
    double ar, double ai,	
    double br, double bi,	
    double *cr, double *ci	
) ;


#ifndef NTIMER
#ifdef _POSIX_C_SOURCE
#if    _POSIX_C_SOURCE >= 199309L
#define SPARSE_TIMER_ENABLED
#endif
#endif
#endif


#define SPARSE_PRINTF(params) \
{ \
    if (SparseBase_config.printf_func != NULL) \
    { \
        (void) (SparseBase_config.printf_func) params ; \
    } \
}

#ifdef __cplusplus
}
#endif
#endif
