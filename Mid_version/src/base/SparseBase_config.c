/* ========================================================================== */
/* === SparseBase_config =================================================== */
/* ========================================================================== */

/* SparseBase configuration : memory manager and printf functions. */

#include <math.h>
#include <stdlib.h>

#ifndef NPRINT
#include <stdio.h>
#endif

#ifndef NULL
#define NULL ((void *) 0)
#endif

#include "SparseBase_config.h"

/* -------------------------------------------------------------------------- */
/* SparseBase_config : a global extern struct */
/* -------------------------------------------------------------------------- */

struct SparseBase_config_struct SparseBase_config =
{

    /* memory management functions */
    #ifndef NMALLOC
            /* standard ANSI C: */
            malloc, calloc, realloc, free,
    #else
        /* no memory manager defined; you must define one at run-time: */
        NULL, NULL, NULL, NULL,
    #endif

    /* printf function */
    #ifndef NPRINT
            /* standard ANSI C: */
            printf,
    #else
        /* printf is disabled */
        NULL,
    #endif

    SparseBase_hypot,
    SparseBase_divcomplex

} ;

/* -------------------------------------------------------------------------- */
/* SparseBase_start */
/* -------------------------------------------------------------------------- */

/* All applications that use SparseBase should call SparseBase_start prior
   to using any SparseBase function.  Only a single thread should call this
   function, in a multithreaded application.  Currently, this function is
   optional, since all this function currently does is to set the four memory
   function pointers to NULL (which tells SparseBase to use the default
   functions).  In a multi- threaded application, only a single thread should
   call this function.
 */

void SparseBase_start ( void )
{

    /* memory management functions */
    #ifndef NMALLOC
            /* standard ANSI C: */
            SparseBase_config.malloc_func  = malloc ;
            SparseBase_config.calloc_func  = calloc ;
            SparseBase_config.realloc_func = realloc ;
            SparseBase_config.free_func    = free ;
    #else
        /* no memory manager defined; you must define one after calling
           SparseBase_start */
        SparseBase_config.malloc_func  = NULL ;
        SparseBase_config.calloc_func  = NULL ;
        SparseBase_config.realloc_func = NULL ;
        SparseBase_config.free_func    = NULL ;
    #endif

    /* printf function */
    #ifndef NPRINT
            /* standard ANSI C: */
            SparseBase_config.printf_func = printf ;
    #else
        /* printf is disabled */
        SparseBase_config.printf_func = NULL ;
    #endif

    /* math functions */
    SparseBase_config.hypot_func = SparseBase_hypot ;
    SparseBase_config.divcomplex_func = SparseBase_divcomplex ;
}

/* -------------------------------------------------------------------------- */
/* SparseBase_finish */
/* -------------------------------------------------------------------------- */

/* This currently does nothing, but in the future, applications should call
   SparseBase_start before calling any SparseBase function, and then
   SparseBase_finish after calling the last SparseBase function, just before
   exiting.  In a multithreaded application, only a single thread should call
   this function.
 */

void SparseBase_finish ( void )
{
    /* do nothing */ ;
}

/* -------------------------------------------------------------------------- */
/* SparseBase_malloc: malloc wrapper */
/* -------------------------------------------------------------------------- */

void *SparseBase_malloc    /* pointer to allocated block of memory */
(
    size_t nitems,          /* number of items to malloc */
    size_t size_of_item     /* sizeof each item */
)
{
    void *p ;
    size_t size ;
    if (nitems < 1) nitems = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems * size_of_item  ;

    if (size != ((double) nitems) * size_of_item)
    {
        /* size_t overflow */
        p = NULL ;
    }
    else
    {
        p = (void *) (SparseBase_config.malloc_func) (size) ;
    }
    return (p) ;
}


/* -------------------------------------------------------------------------- */
/* SparseBase_calloc: calloc wrapper */
/* -------------------------------------------------------------------------- */

void *SparseBase_calloc    /* pointer to allocated block of memory */
(
    size_t nitems,          /* number of items to calloc */
    size_t size_of_item     /* sizeof each item */
)
{
    void *p ;
    size_t size ;
    if (nitems < 1) nitems = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems * size_of_item  ;

    if (size != ((double) nitems) * size_of_item)
    {
        /* size_t overflow */
        p = NULL ;
    }
    else
    {
        p = (void *) (SparseBase_config.calloc_func) (nitems, size_of_item) ;
    }
    return (p) ;
}

/* -------------------------------------------------------------------------- */
/* SparseBase_realloc: realloc wrapper */
/* -------------------------------------------------------------------------- */

/* If p is non-NULL on input, it points to a previously allocated object of
   size nitems_old * size_of_item.  The object is reallocated to be of size
   nitems_new * size_of_item.  If p is NULL on input, then a new object of that
   size is allocated.  On success, a pointer to the new object is returned,
   and ok is returned as 1.  If the allocation fails, ok is set to 0 and a
   pointer to the old (unmodified) object is returned.
 */

void *SparseBase_realloc   /* pointer to reallocated block of memory, or
                               to original block if the realloc failed. */
(
    size_t nitems_new,      /* new number of items in the object */
    size_t nitems_old,      /* old number of items in the object */
    size_t size_of_item,    /* sizeof each item */
    void *p,                /* old object to reallocate */
    int *ok                 /* 1 if successful, 0 otherwise */
)
{
    size_t size ;
    if (nitems_old < 1) nitems_old = 1 ;
    if (nitems_new < 1) nitems_new = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems_new * size_of_item  ;

    if (size != ((double) nitems_new) * size_of_item)
    {
        /* size_t overflow */
        (*ok) = 0 ;
    }
    else if (p == NULL)
    {
        /* a fresh object is being allocated */
        p = SparseBase_malloc (nitems_new, size_of_item) ;
        (*ok) = (p != NULL) ;
    }
    else if (nitems_old == nitems_new)
    {
        /* the object does not change; do nothing */
        (*ok) = 1 ;
    }
    else
    {
        /* change the size of the object from nitems_old to nitems_new */
        void *pnew ;
        pnew = (void *) (SparseBase_config.realloc_func) (p, size) ;
        if (pnew == NULL)
        {
            if (nitems_new < nitems_old)
            {
                /* the attempt to reduce the size of the block failed, but
                   the old block is unchanged.  So pretend to succeed. */
                (*ok) = 1 ;
            }
            else
            {
                /* out of memory */
                (*ok) = 0 ;
            }
        }
        else
        {
            /* success */
            p = pnew ;
            (*ok) = 1 ;
        }
    }
    return (p) ;
}

/* -------------------------------------------------------------------------- */
/* SparseBase_free: free wrapper */
/* -------------------------------------------------------------------------- */

void *SparseBase_free      /* always returns NULL */
(
    void *p                 /* block to free */
)
{
    if (p)
    {
        (SparseBase_config.free_func) (p) ;
    }
    return (NULL) ;
}


/* -------------------------------------------------------------------------- */
/* SparseBase_tic: return current wall clock time */
/* -------------------------------------------------------------------------- */

/* Returns the number of seconds (tic [0]) and nanoseconds (tic [1]) since some
 * unspecified but fixed time in the past.  If no timer is installed, zero is
 * returned.  A scalar double precision value for 'tic' could be used, but this
 * might cause loss of precision because clock_getttime returns the time from
 * some distant time in the past.  Thus, an array of size 2 is used.
 *
 * The timer is enabled by default.  To disable the timer, compile with
 * -DNTIMER.  If enabled on a POSIX C 1993 system, the timer requires linking
 * with the -lrt library.
 *
 * example:
 *
 *      double tic [2], r, s, t ;
 *      SparseBase_tic (tic) ;     // start the timer
 *      // do some work A
 *      t = SparseBase_toc (tic) ; // t is time for work A, in seconds
 *      // do some work B
 *      s = SparseBase_toc (tic) ; // s is time for work A and B, in seconds
 *      SparseBase_tic (tic) ;     // restart the timer
 *      // do some work C
 *      r = SparseBase_toc (tic) ; // s is time for work C, in seconds
 *
 * A double array of size 2 is used so that this routine can be more easily
 * ported to non-POSIX systems.  The caller does not rely on the POSIX
 * <time.h> include file.
 */

#ifdef SPARSE_TIMER_ENABLED

#include <time.h>

void SparseBase_tic
(
    double tic [2]      /* output, contents undefined on input */
)
{
    /* POSIX C 1993 timer, requires -librt */
    struct timespec t ;
    clock_gettime (CLOCK_MONOTONIC, &t) ;
    tic [0] = (double) (t.tv_sec) ;
    tic [1] = (double) (t.tv_nsec) ;
}

#else

void SparseBase_tic
(
    double tic [2]      /* output, contents undefined on input */
)
{
    /* no timer installed */
    tic [0] = 0 ;
    tic [1] = 0 ;
}

#endif


/* -------------------------------------------------------------------------- */
/* SparseBase_toc: return time since last tic */
/* -------------------------------------------------------------------------- */

/* Assuming SparseBase_tic is accurate to the nanosecond, this function is
 * accurate down to the nanosecond for 2^53 nanoseconds since the last call to
 * SparseBase_tic, which is sufficient for SparseBase (about 104 days).  If
 * additional accuracy is required, the caller can use two calls to
 * SparseBase_tic and do the calculations differently.
 */

double SparseBase_toc  /* returns time in seconds since last tic */
(
    double tic [2]  /* input, not modified from last call to SparseBase_tic */
)
{
    double toc [2] ;
    SparseBase_tic (toc) ;
    return ((toc [0] - tic [0]) + 1e-9 * (toc [1] - tic [1])) ;
}


/* -------------------------------------------------------------------------- */
/* SparseBase_time: return current wallclock time in seconds */
/* -------------------------------------------------------------------------- */

/* This function might not be accurate down to the nanosecond. */

double SparseBase_time  /* returns current wall clock time in seconds */
(
    void
)
{
    double toc [2] ;
    SparseBase_tic (toc) ;
    return (toc [0] + 1e-9 * toc [1]) ;
}

/* -------------------------------------------------------------------------- */
/* SparseBase_hypot */
/* -------------------------------------------------------------------------- */

/* There is an equivalent routine called hypot in <math.h>, which conforms
 * to ANSI C99.  However, SparseBase does not assume that ANSI C99 is
 * available.  You can use the ANSI C99 hypot routine with:
 *
 *      #include <math.h>
 *i     SparseBase_config.hypot_func = hypot ;
 *
 * Default value of the SparseBase_config.hypot_func pointer is
 * SparseBase_hypot, defined below.
 *
 * s = hypot (x,y) computes s = sqrt (x*x + y*y) but does so more accurately.
 * The NaN cases for the double relops x >= y and x+y == x are safely ignored.
 * 
 * Source: Algorithm 312, "Absolute value and square root of a complex number,"
 * P. Friedland, Comm. ACM, vol 10, no 10, October 1967, page 665.
 */

double SparseBase_hypot (double x, double y)
{
    double s, r ;
    x = fabs (x) ;
    y = fabs (y) ;
    if (x >= y)
    {
        if (x + y == x)
        {
            s = x ;
        }
        else
        {
            r = y / x ;
            s = x * sqrt (1.0 + r*r) ;
        }
    }
    else
    {
        if (y + x == y)
        {
            s = y ;
        }
        else
        {
            r = x / y ;
            s = y * sqrt (1.0 + r*r) ;
        }
    } 
    return (s) ;
}

/* -------------------------------------------------------------------------- */
/* SparseBase_divcomplex */
/* -------------------------------------------------------------------------- */

/* c = a/b where c, a, and b are complex.  The real and imaginary parts are
 * passed as separate arguments to this routine.  The NaN case is ignored
 * for the double relop br >= bi.  Returns 1 if the denominator is zero,
 * 0 otherwise.
 *
 * This uses ACM Algo 116, by R. L. Smith, 1962, which tries to avoid
 * underflow and overflow.
 *
 * c can be the same variable as a or b.
 *
 * Default value of the SparseBase_config.divcomplex_func pointer is
 * SparseBase_divcomplex.
 */

int SparseBase_divcomplex
(
    double ar, double ai,       /* real and imaginary parts of a */
    double br, double bi,       /* real and imaginary parts of b */
    double *cr, double *ci      /* real and imaginary parts of c */
)
{
    double tr, ti, r, den ;
    if (fabs (br) >= fabs (bi))
    {
        r = bi / br ;
        den = br + r * bi ;
        tr = (ar + ai * r) / den ;
        ti = (ai - ar * r) / den ;
    }
    else
    {
        r = br / bi ;
        den = r * br + bi ;
        tr = (ar * r + ai) / den ;
        ti = (ai * r - ar) / den ;
    }
    *cr = tr ;
    *ci = ti ;
    return (den == 0.) ;
}
