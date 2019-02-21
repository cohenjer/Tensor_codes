#include<stdlib.h>
#include<assert.h>

#include "testing_tools.h"
#include "util.h"

/* Are equal arrays compared equal? */
void TestIntArrCompare(CuTest *tc) {
    int a[3] = {1,2,3};
    int b[3] = {1,2,3};
    int diff_pos;

    diff_pos = FGNSRassert_intarr_eql(tc, 3, &a[0], &b[0]);
    CuAssertIntEquals(tc, 3, diff_pos);
}

void TestDenseToSp1by1(CuTest *tc) {
    int numvec=1, lenvec=1;
    double dmatrix[1] = {4.29348577};
    FGNSRSPNZ *begin=NULL;
    FGNSRSPDIM *index=NULL;
    double *value=NULL;
    int status=0;

    status = FGNSRdtosp(numvec, lenvec, dmatrix, &begin, &index, &value);
    CuAssertIntEquals(tc, 0, status);

    CuAssertPtrNotNull(tc, begin);
    CuAssertPtrNotNull(tc, index);
    CuAssertPtrNotNull(tc, value);
    CuAssertIntEquals(tc, 0, begin[0]);
    CuAssertIntEquals(tc, 1, begin[1]);

    FGNSRassert_dblarr_abseql(tc, 1, &dmatrix[0], value, 1e-8);

    FGNSRfree((void **) &begin);
    FGNSRfree((void **) &index);
    FGNSRfree((void **) &value);
}

void TestSpToDense1by1(CuTest *tc) {
    int numvec=1, lenvec=1;
    double value[1] = {3.14159267};
    FGNSRSPNZ begin[2] = {0,1};
    FGNSRSPDIM index[1] = {0};
    double *dmatrix=NULL;
    int status=0;

    status = FGNSRsptod(&begin[0], &index[0], &value[0], numvec, lenvec,
            &dmatrix);

    CuAssertIntEquals(tc, 0, status);
    CuAssertPtrNotNull(tc, dmatrix);
    FGNSRassert_dblarr_abseql(tc, 1, value, dmatrix, 1e-8);

    FGNSRfree((void **) &dmatrix);
}

void TestDenseToSpEmpty(CuTest *tc) {
    int numvec=2, lenvec=2;
    double dmatrix[4] = {0.,0.,0.,0.};
    FGNSRSPNZ *begin=NULL;
    FGNSRSPDIM *index=NULL;
    double *value=NULL;
    FGNSRSPNZ ex_begin[3] = {0,0,0};
    int status=0;

    status = FGNSRdtosp(numvec, lenvec, dmatrix, &begin, &index, &value);
    CuAssertIntEquals(tc, 0, status);

    CuAssertPtrNotNull(tc, begin);
    CuAssertPtrNotNull(tc, index);
    CuAssertPtrNotNull(tc, value);
    FGNSRassert_nzarr_eql(tc, 3, &ex_begin[0], begin);

    FGNSRfree((void **) &begin);
    FGNSRfree((void **) &index);
    FGNSRfree((void **) &value);
}

void TestSpToDenseEmpty(CuTest *tc) {
    int numvec=2, lenvec=2;
    double value[1] = {3.14159267}; /* Should not be touched */
    FGNSRSPNZ begin[3] = {0,0,0};
    FGNSRSPDIM index[1] = {1}; /* Should not be touched */
    double *dmatrix=NULL;

    double ex_dmatrix[4] = {0.,0.,0.,0.};
    int status=0;

    status = FGNSRsptod(&begin[0], &index[0], &value[0], numvec, lenvec,
            &dmatrix);

    CuAssertIntEquals(tc, 0, status);
    CuAssertPtrNotNull(tc, dmatrix);
    FGNSRassert_dblarr_abseql(tc, 4, ex_dmatrix, dmatrix, 0.0);

    FGNSRfree((void **) &dmatrix);
}

void TestDenseToSp(CuTest *tc) {
    int numvec=3, lenvec=4;
    double dmatrix[12] = {
        .5, 1.,   0., 1.,
        0., 8.,   2., 0.,
        1., 0.25, 0., 0.,
    };
    FGNSRSPNZ *begin=NULL;
    FGNSRSPNZ ex_begin[4] = {0, 3, 5, 7};
    FGNSRSPDIM *index=NULL;
    FGNSRSPDIM ex_index[7] = {0, 1, 3, 1, 2, 0, 1};
    double *value=NULL;
    double ex_value[7] = {.5, 1., 1., 8., 2., 1., 0.25};
    int status=0;

    status = FGNSRdtosp(numvec, lenvec, dmatrix, &begin, &index, &value);
    CuAssertIntEquals(tc, 0, status);

    FGNSRassert_nzarr_eql(tc, 4, &ex_begin[0], begin);
    FGNSRassert_dimarr_eql(tc, 7, &ex_index[0], index);
    FGNSRassert_dblarr_abseql(tc, 7, &ex_value[0], value, 0.0);

    FGNSRfree((void**) &begin);
    FGNSRfree((void**) &index);
    FGNSRfree((void**) &value);
}

void TestSpToDense(CuTest *tc) {
    int numvec=3, lenvec=4;
    double *dmatrix = NULL;
    double ex_dmatrix[12] = {
        .5, 1.,   0., 1.,
        0., 8.,   2., 0.,
        1., 0.25, 0., 0.,
    };
    FGNSRSPNZ begin[4] = {0, 3, 5, 7};
    FGNSRSPDIM index[7] = {0, 1, 3, 1, 2, 0, 1};
    double value[7] = {.5, 1., 1., 8., 2., 1., 0.25};
    int status=0;

    status = FGNSRsptod(&begin[0], &index[0], &value[0], numvec, lenvec,
            &dmatrix);
    CuAssertIntEquals(tc, 0, status);

    FGNSRassert_dblarr_abseql(tc, 7, &ex_dmatrix[0], dmatrix, 0.0);

    FGNSRfree( (void **) &dmatrix);
}

void TestCopyTwice(CuTest *tc) {
    int numvec=3, lenvec=4;
    double dmatrix[12] = {
        0, 1, 2, 3,
        0, 0, 0, 0,
        4, 0, 6, 0,
    };

    FGNSRSPNZ *begin=NULL;
    FGNSRSPDIM *index=NULL;
    double *value=NULL;
    double *dmatrix_copy=NULL;
    int status = 0;

    status = FGNSRdtosp(numvec, lenvec, dmatrix, &begin, &index, &value);
    CuAssertIntEquals(tc, 0, status);
    status = FGNSRsptod(&begin[0], &index[0], &value[0], numvec, lenvec,
            &dmatrix_copy);
    CuAssertIntEquals(tc, 0, status);

    FGNSRassert_dblarr_abseql(tc, numvec*lenvec, dmatrix, dmatrix_copy, 0.0);
    FGNSRfree((void **) &dmatrix_copy);
    FGNSRfree((void **) &begin);
    FGNSRfree((void **) &index);
    FGNSRfree((void **) &value);
}


CuSuite* GetTestingToolsSuite() {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, TestIntArrCompare);
    SUITE_ADD_TEST(suite, TestSpToDense1by1);
    SUITE_ADD_TEST(suite, TestDenseToSp1by1);
    SUITE_ADD_TEST(suite, TestSpToDenseEmpty);
    SUITE_ADD_TEST(suite, TestDenseToSpEmpty);
    SUITE_ADD_TEST(suite, TestDenseToSp);
    SUITE_ADD_TEST(suite, TestSpToDense);
    SUITE_ADD_TEST(suite, TestCopyTwice);
    return suite;
}

