#include<assert.h>

#include "debug_tools.h"
#include "sputil.h"

/* If this has a positive value k, then it is guaranteed that the k-th call to
 * malloc results in NULL.  If the value is non-positive, no such pretended OOM
 * condition will be triggered. */
static long GLOBAL_failcounter = 0;

void FGNSRdbg_setmallocfail(long next_fail) {
    assert(next_fail >= 0);
#pragma omp critical (MALLOCDEBUG)
    {
    GLOBAL_failcounter = next_fail;
    }
}

void * FGNSRdbg_malloc(size_t s) {
    void *retval;
#pragma omp critical (MALLOCDEBUG)
    {
    if ( (GLOBAL_failcounter > 0) && (--GLOBAL_failcounter == 0) )
        retval = NULL;
    else
        retval = malloc(s);
    }

    return retval;

}

void * FGNSRdbg_calloc(size_t c, size_t s) {
    void *retval;
#pragma omp critical (MALLOCDEBUG)
    {
    if ( (GLOBAL_failcounter > 0) && (--GLOBAL_failcounter == 0) )
        retval = NULL;
    else
        retval = calloc(c,s);
    }

    return retval;

}

void FGNSRdbg_free(void * p) {
    free(p);
}

