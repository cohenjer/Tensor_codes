#include<stdlib.h>
#include<assert.h>

#include "CuTest.h"
#include "debug_tools.h"
#include "util.h"

void TestDbgMallocFail(CuTest *tc) {
    int i;
    size_t msize = 8;
    void *p = NULL;

    FGNSRdbg_setmallocfail(0l);
    for (i=0; i<5; i++) {
        p = FGNSRdbg_malloc(msize);
        CuAssertPtrNotNull(tc, p);
        FGNSRdbg_free(p);
    }

    FGNSRdbg_setmallocfail(1l);
    p = FGNSRdbg_malloc(msize);
    CuAssertPtrEquals(tc, NULL, p);

    p = FGNSRdbg_malloc(msize);
    CuAssertPtrNotNull(tc, p);
    FGNSRdbg_free(p);

    p = FGNSRdbg_malloc(msize);
    CuAssertPtrNotNull(tc, p);
    FGNSRdbg_free(p);

    FGNSRdbg_setmallocfail(6l);
    for (i=0; i<5; i++) {
        p = FGNSRdbg_malloc(msize);
        CuAssertPtrNotNull(tc, p);
        FGNSRdbg_free(p);
    }
    FGNSRdbg_setmallocfail(1l);
    p = FGNSRdbg_malloc(msize);
    CuAssertPtrEquals(tc, NULL, p);
}

CuSuite* GetDebugToolsSuite() {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, TestDbgMallocFail);
    return suite;
}

