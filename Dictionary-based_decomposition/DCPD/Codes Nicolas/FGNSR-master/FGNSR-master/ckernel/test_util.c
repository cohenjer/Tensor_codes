#include<stdlib.h>
#include<assert.h>
#include<unistd.h>

#include "util.h"

#include "CuTest.h"

/* These tests are written with the intention that they be run with valgrind */

void TestFreeNull(CuTest *tc) {
    double *dblptr = NULL;
    dblptr = (double *) malloc(sizeof(double));
    assert(dblptr != NULL);

    FGNSRfree((void **) &dblptr);
    CuAssertPtrEquals(tc, NULL, dblptr);
}

void TestFGNSRmalloc(CuTest *tc) {
    int *intptr = NULL;
    intptr = (int *) FGNSRmalloc(3 * sizeof(int));
    CuAssertPtrNotNull(tc, intptr);

    intptr[0] = 9;
    intptr[1] = 8;
    intptr[2] = 7;

    FGNSRfree((void **) &intptr);
    CuAssertPtrEquals(tc, NULL, intptr);
}

void TestFGNSRcalloc(CuTest *tc) {
    int *intptr = NULL;
    intptr = (int *) FGNSRcalloc(3, sizeof(int));
    CuAssertPtrNotNull(tc, intptr);

    CuAssertIntEquals(tc, 0, intptr[0]);
    CuAssertIntEquals(tc, 0, intptr[1]);
    CuAssertIntEquals(tc, 0, intptr[2]);
    intptr[0] = 9;
    intptr[1] = 8;
    intptr[2] = 7;

    FGNSRfree((void **) &intptr);
    CuAssertPtrEquals(tc, NULL, intptr);
}

void TestGetTimestamp(CuTest *tc) {
    DATUM stamp;

    FGNSRtimestamp(&stamp);
    /* Nothing to test here, just make sure this function runs clean with
     * valgrind */
    CuAssertIntEquals(tc, 0, 0);
}

void TestWtime(CuTest *tc) {
    DATUM start, end;
    double wc_time = 0.0;

    FGNSRtimestamp(&start);
    sleep(1.0);
    FGNSRtimestamp(&end);

    /* Elapsed wallclock time should be about 1s */
    wc_time = FGNSRwctime(start, end);
    CuAssertDblEquals(tc, 1.0, wc_time, 0.1);

}

void TestCPUtime(CuTest *tc) {
    DATUM start, end;
    double cpu_time = 0.0;

    FGNSRtimestamp(&start);
    sleep(1.0);
    FGNSRtimestamp(&end);

    /* Elapsed CPU time should be about 0s */
    cpu_time = FGNSRcputime(start, end);
    CuAssertDblEquals(tc, 0.0, cpu_time, 0.1);
}

CuSuite* GetUtilSuite() {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, TestFreeNull);
    SUITE_ADD_TEST(suite, TestFGNSRmalloc);
    SUITE_ADD_TEST(suite, TestFGNSRcalloc);
    SUITE_ADD_TEST(suite, TestGetTimestamp);
    SUITE_ADD_TEST(suite, TestWtime);
    SUITE_ADD_TEST(suite, TestCPUtime);
    return suite;
}

