#include<stdlib.h>

#include "CuTest.h"

#include "sputil.h"
#include "fgnsrdef.h"
#include "testing_tools.h"

/* These tests are written with the intention that they be run with valgrind */

void TestSpAlloc1by1(CuTest *tc) {
    int numvec=1, nnz=1;
    FGNSRSPNZ *begin=NULL;
    FGNSRSPDIM *index=NULL;
    double *value=NULL;
    int status = 0;

    status = FGNSRalloc_spmatrix(numvec, nnz, &begin, &index, &value);
    CuAssertIntEquals(tc, 0, status);
    CuAssertPtrNotNull(tc, begin);
    CuAssertPtrNotNull(tc, value);
    CuAssertPtrNotNull(tc, index);
            
    /* No elements sitting in the matrix, hence nz count zero */
    CuAssertIntEquals(tc, 0, (int)begin[0]);
    CuAssertIntEquals(tc, 0, (int)begin[1]);

    /* Touch all of the memory for valgrind, nothing to test  */
    begin[0] = 0;
    begin[1] = 1;
    value[0] = 235.22677e12;
    index[0] = 21;

    FGNSRfree((void **) &begin);
    FGNSRfree((void **) &index);
    FGNSRfree((void **) &value);
}

#ifdef FGNSRDEBUG
static int fun_TestSpAllocOOM(void *data) {
    int numvec = 4;
    int nnz = 17;
    FGNSRSPNZ *begin=NULL;
    FGNSRSPDIM *index=NULL;
    double *value=NULL;
    int status = 0;

    status = FGNSRalloc_spmatrix(numvec, nnz, &begin, &index, &value);
    if (status == FGNSRSTAT_OK) {
        /* Allocation succeeded, free memory */
        FGNSRfree((void **) &begin);
        FGNSRfree((void **) &index);
        FGNSRfree((void **) &value);
    }

    /* else:  Some allocation failed, but the allocating function should free
     * succeeded mallocs by itself, so we do not call free ourselves.  If it
     * doesn't, valgrind would detect it. */

    return status;
}

void TestSpAlloOOM(CuTest *tc) {
    FGNSRoom_loop(tc, 10, fun_TestSpAllocOOM, NULL);
}
#endif /* FGNSRDEBUG */

CuSuite* GetSpUtilSuite() {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, TestSpAlloc1by1);
#ifdef FGNSRDEBUG
    SUITE_ADD_TEST(suite, TestSpAlloOOM);
#endif /* FGNSRDEBUG */
    return suite;
}

