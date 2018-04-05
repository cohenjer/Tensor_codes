#ifndef DEBUG_TOOLS_H_INCDEF
#define DEBUG_TOOLS_H_INCDEF

#include<stdlib.h>
#include "fgnsrdef.h"

void FGNSRdbg_setmallocfail(long next_fail);
void * FGNSRdbg_malloc(size_t size);
void * FGNSRdbg_calloc(size_t count, size_t size);
void FGNSRdbg_free(void * p);


#endif /* DEBUG_TOOLS_H_INCDEF */

