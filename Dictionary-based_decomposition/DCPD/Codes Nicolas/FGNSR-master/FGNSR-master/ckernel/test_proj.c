#include<stdlib.h>
#include<assert.h>
#include<string.h>

#include "proj.h"
#include "util.h"
#include "fgnsrdef.h"

#include "testing_tools.h"

#ifdef FGNSRDEBUG
#include "debug_tools.h"
#endif

/* Data driven test routine:  Input data along with expected result in dense
 * format.  Data is converted automatically to sparse data and rerun. */
static void datatest_dense_sparse(CuTest *tc, int numvec, int lenvec,
        double *weights, double *dmatrix, double ub, int *diag_idx,
        int verbose, double *ex_dmatrix, double tol)
{
    FGNSRSPNZ  *begin=NULL;
    FGNSRSPDIM *index=NULL;
    FGNSRSPDIM *sp_diag_idx=NULL;
    double   *value=NULL;
    double *sp_dmatrix=NULL;

    PROJDATAptr pdata=NULL;
    SPPROJDATAptr spdata=NULL;
    int status;

    /* Obtain copy of the input data as sparse matrix */
    status = FGNSRdtosp(numvec, lenvec, dmatrix, &begin, &index,
            &value);
    assert(status==0);

    status = FGNSRpromote_spdim(numvec, diag_idx, &sp_diag_idx);
    assert(status==0);

    /* Test dense projection */
    status = FGNSRproj_init(lenvec, weights, verbose, &pdata);
    CuAssertIntEquals(tc, 0, status);

    status = FGNSRproj_project_dmatrix(pdata, numvec, dmatrix, ub, diag_idx);
    assert(status == 0);

    FGNSRassert_dblarr_abseql(tc, numvec*lenvec, ex_dmatrix, dmatrix, tol);

    FGNSRproj_free(&pdata);

    /* Test sparse projection */
    status = FGNSRproj_spinit(lenvec, weights, verbose, &spdata);
    assert(status == 0);

    status = FGNSRproj_project_spmatrix(spdata, numvec, begin, index, value, ub,
            sp_diag_idx);
    assert(status == 0);

    status = FGNSRsptod(begin, index, value, numvec, lenvec, &sp_dmatrix);
    assert(status==0);

    FGNSRassert_dblarr_abseql(tc, numvec*lenvec, ex_dmatrix, sp_dmatrix, tol);

    FGNSRproj_spfree(&spdata);


    /* Memory cleanup */
    FGNSRfree((void **) &begin);
    FGNSRfree((void **) &index);
    FGNSRfree((void **) &value);
    FGNSRfree((void **) &sp_diag_idx);
    FGNSRfree((void **) &sp_dmatrix);

}

static void datatest_dense(CuTest *tc, int numvec, int lenvec,
        double *weights, double *dmatrix, double ub, int *diag_idx,
        int verbose, double *ex_dmatrix, double tol)
{
    PROJDATAptr pdata=NULL;
    int status;

    /* Test dense projection */
    status = FGNSRproj_init(lenvec, weights, verbose, &pdata);
    assert(status == 0);

    status = FGNSRproj_project_dmatrix(pdata, numvec, dmatrix, ub, diag_idx);
    CuAssertIntEquals(tc, 0, status);

    FGNSRassert_dblarr_abseql(tc, numvec*lenvec, ex_dmatrix, dmatrix, tol);

    FGNSRproj_free(&pdata);
}



/* Are pointers set from NULL and back ?*/
void TestInitFree(CuTest *tc) {
    int status = 0;
    int lenvec=3, verbose=0;
    double weights[3] = {1., .5, .2};
    PROJDATAptr pdata=NULL;

    status = FGNSRproj_init(lenvec, weights, verbose, &pdata);

    /* Given the small size, OOM conditions are virtually impossible */
    CuAssertIntEquals(tc, FGNSRSTAT_OK, status);
    CuAssertPtrNotNull(tc, pdata);

    FGNSRproj_free(&pdata);
    CuAssertPtrEquals(tc, NULL, pdata);
}

/* What happens for 1x1 input? */
void TestOneByOne1(CuTest *tc) {
    int lenvec = 1;
    int numvec = 1;
    double weights[1] = {1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[1] = {0};
    double matrix[1] = {1.};
    double ex_matrix[1] = {1.0};
    double tol = 0.0;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* What happens for 1x1 input with upper bound? */
void TestOneByOne2(CuTest *tc) {
    int lenvec = 1;
    int numvec = 1;
    double weights[1] = {1.};
    int verbose = 0;
    double ub = .5;
    int diag_idx[1] = {0};
    double matrix[1] = {1.};
    double ex_matrix[1] = {0.5};
    double tol = 0.0;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* What happens for 1x1 input with negative entry? */
void TestOneByOne3(CuTest *tc) {
    int lenvec = 1;
    int numvec = 1;
    double weights[1] = {1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[1] = {0};
    double matrix[1] = {-3.14159267};
    double ex_matrix[1] = {0.0};
    double tol = 0.0;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);

}


/* Simple functional test */
void TestTypicalProjection(CuTest *tc) {
    int lenvec = 3;
    int numvec = 3;
    double weights[3] = {.5, 2., 3.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[3] = {1,0,2};
    double matrix[9] = {
        1, 4, 7, 
        2, 5, 8, 
        3, 6, 9,
    };
    double ex_matrix[9] = {
        .25,   1.,    1.5,
        1.,    4.,    6.,
        1./6., 2./3., 1.
    };
    double tol = 1e-8;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* Project twice */
void TestProjectTwiceDense(CuTest *tc) {
    int status = 0;
    int lenvec = 3;
    int numvec = 3;
    double weights[3] = {.5, 2., 3.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[3] = {1,0,2};
    double matrix[9] = {
        1, 4, 7, 
        2, 5, 8, 
        3, 6, 9,
    };
    double ex_matrix[9] = {
        .25,   1.,    1.5,
        1.,    4.,    6.,
        1./6., 2./3., 1.
    };
    PROJDATAptr pdata = NULL;

    status = FGNSRproj_init(lenvec, weights, verbose, &pdata);
    assert(status == 0);

    /* This projection changes the data */
    status = FGNSRproj_project_dmatrix(pdata, numvec, matrix, ub, diag_idx);
    CuAssertIntEquals(tc, 0, status);

    /* This should not alter the data */
    status = FGNSRproj_project_dmatrix(pdata, numvec, matrix, ub, diag_idx);
    CuAssertIntEquals(tc, 0, status);

    FGNSRassert_dblarr_abseql(tc, lenvec*numvec, ex_matrix, matrix, 1e-8);

    FGNSRproj_free(&pdata);
}

void TestZeroDiagWeight(CuTest *tc) {
    int lenvec = 3;
    int numvec = 1;
    double weights[3] = {.0, 1., 1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[3] = {0, 1, 2};
    double matrix[3] = {2., 2., -1.};
    double ex_matrix[3] = {1., 2., 0.};
    double tol = 1e-8;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);

}

void TestJustClip(CuTest *tc) {
    int lenvec = 3;
    int numvec = 1;
    double weights[3] = {1., 1., 1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[1] = {0};
    double matrix[3] = {1, -1, 1};
    double ex_matrix[3] = {1, 0, 1};
    double tol = 1e-8;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* Simple 2x2 test for permutations */
void TestPermutation1(CuTest *tc) {
    int lenvec = 2;
    int numvec = 2;
    double weights[2] = {1., 1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[2] = {1,0};
    double matrix[4] = {
        1, 0,
        0, 1,
    };
    double ex_matrix[4] = {
        .5, .5,
        .5, .5,
    };
    double tol = 1e-8;

    /* Test only dense as diag values are zero */
    datatest_dense(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* Simple 2x2 test for permutations */
void TestPermutation2(CuTest *tc) {
    int lenvec = 2;
    int numvec = 2;
    double weights[2] = {1., 1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[2] = {1,0};
    double matrix[4] = {
        .75, .25,
        .25, .75,
    };
    double ex_matrix[4] = {
        .5, .5,
        .5, .5,
    };
    double tol = 1e-8;

    /* Test only dense as diag values are zero */
    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* 2x2 ID matrix */
void TestIdentity(CuTest *tc) {
    int lenvec = 2;
    int numvec = 2;
    double weights[2] = {1., 1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[2] = {0,1};
    double matrix[4] = {
        1, 0,
        0, 1,
    };
    double ex_matrix[4] = {
        1, 0,
        0, 1,
    };
    double tol = 1e-8;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* 1x1 corner case with negative diagonal entry of negative weight */
void TestNegValZeroWeight(CuTest *tc) {
    int lenvec = 1;
    int numvec = 1;
    double weights[1] = {0.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[1] = {0};
    double matrix[1] = {-1.};
    double ex_matrix[1] = {0.};
    double tol = 1e-8;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* Negative diagonal entry is a bit special, since it must not be clipped */
void TestNegDiag(CuTest *tc) {
    int lenvec = 2;
    int numvec = 1;
    double weights[2] = {1., 1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[1] = {0};
    double matrix[2] = {-.25, 1};
    double ex_matrix[2] = {3./8., 3./8.};
    double tol = 1e-8;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* Trivial zero assignments */
void TestZeroSetting(CuTest *tc) {
    int lenvec = 3;
    int numvec = 1;
    double weights[3] = {0., 1., 1.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[1] = {1};
    double matrix[3] = {2., .9, -1.};
    double ex_matrix[3] = {0.,.9,0.};
    double tol = 1e-8;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* Sample data where one really needs to sort violations */
void TestSorting(CuTest *tc) {
    int lenvec = 4;
    int numvec = 1;
    double weights[4] = {.75, .5, .5, .25};
    double ub = 100; /* This is not tight */
    int diag_idx[1] = {1};
    double matrix[4] = {4., 1., 16., 64.,};
    double ex_matrix[4] = {4., 26.4, 16., 13.2};
    int verbose = 0;
    double tol = 1e-8;
    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* Sample data where thresholding can be terminated early because upper bound
 * has been reached */
void TestCutUb(CuTest *tc) {
    int lenvec = 4;
    int numvec = 1;
    double weights[4] = {.75, .5, .5, .25};
    double ub = 16; /* This is tight */
    int diag_idx[1] = {1};
    double matrix[4] = {4., 1., 16., 64.,};
    double ex_matrix[4] = {4., 16., 16., 8.};
    int verbose = 0;
    double tol = 1e-8;
    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* Sample test with no particular intention; it just enables verbose output */
void TestVerboseSwitch(CuTest *tc) {
    int lenvec = 4;
    int numvec = 4;
    double weights[4] = {1.,1.,1.,1.};
    int verbose = 1;
    double ub = 100;
    int diag_idx[4] = {0,1,2,3};
    double matrix[16] = {
        5., 10.,  3., 6.,
        6.,  2., 10., 9.,
        9.,  2.,  6., 4.,
        2.,  8.,  4., 5.,
    };
    double ex_matrix[16] = {
        7.5, 7.5, 3.,  6.,
        6.,  7.,  7.,  7.,
        7.5, 2.,  7.5, 4.,
        2.,  6.5, 4.0, 6.5,
    };
    double tol = 1e-8;

    datatest_dense_sparse(tc, numvec, lenvec, &weights[0], &matrix[0], ub,
            &diag_idx[0], verbose, &ex_matrix[0], tol);
}

/* An error must be triggered if there is no diagonal element in a sparse
 * vector */
void TestDiagElementMissing(CuTest *tc) {
    int status = 0;
    int lenvec = 3;
    int numvec = 1;
    double weights[3] = {1., 1., 1.};
    int verbose = 0;
    double ub = 1.0;
    FGNSRSPDIM diag_idx[1] = {2};
    FGNSRSPNZ begin[2] = {0, 2};
    FGNSRSPDIM index[2] = {0,1};
    double value[2] = {.5, .5};
    SPPROJDATAptr pdata = NULL;

    status = FGNSRproj_spinit(lenvec, weights, verbose, &pdata);
    assert(status == 0);

    status = FGNSRproj_project_spmatrix(pdata, numvec, begin, index, value, ub,
            diag_idx);
    CuAssertIntEquals(tc, FGNSRSTAT_DATAERR, status);

    FGNSRproj_spfree(&pdata);
}


#ifdef FGNSRDEBUG
/* Enforce OOM conditions */
void TestOOMConditionsDense(CuTest *tc) {
    int status = 0;
    int lenvec = 3;
    int numvec = 3;
    double weights[3] = {.5, 2., 3.};
    int verbose = 0;
    double ub = 1.;
    int diag_idx[3] = {1,0,2};
    double matrix[9] = {
        1, 4, 7, 
        2, 5, 8, 
        3, 6, 9,
    };
    double ex_matrix[9] = {
        .25,   1.,    1.5,
        1.,    4.,    6.,
        1./6., 2./3., 1.
    };
    double ac_matrix[9];
    PROJDATAptr pdata = NULL;
    const int maxmallocs = 100;
    long failcount;
    int clean_pass = 0;

    for (failcount=1; failcount<maxmallocs; failcount++) {
        clean_pass = 0;
        memcpy(ac_matrix, matrix, 9 * sizeof(double));
        FGNSRdbg_setmallocfail(failcount);

        status = FGNSRproj_init(lenvec, weights, verbose, &pdata);
        if (status != FGNSRSTAT_OK) {
            CuAssertIntEquals(tc, FGNSRSTAT_OOM, status);
            continue;
        }

        /* No memory allocated in projection loop at time of writing */
        status = FGNSRproj_project_dmatrix(pdata,numvec, ac_matrix, ub, diag_idx);
        assert(status == 0);

        /* Projection has succeeded, no OOM conditions have happened */
        clean_pass = 1;
        FGNSRassert_dblarr_abseql(tc,lenvec*numvec, ex_matrix, ac_matrix, 1e-8);
        FGNSRproj_free(&pdata);
        break;
    }

    /* If this doesn't hold, one needs to increase maxmallocs.  It means that
     * not every malloc has been tested for OOM cleanliness */
    assert(clean_pass==1);

}

void TestOOMConditionsSparse(CuTest *tc) {
    int status = 0;
    int lenvec = 3;
    int numvec = 3;
    double weights[3] = {.5, 2., 3.};
    int verbose = 0;
    double ub = 1.;
    FGNSRSPDIM diag_idx[3] = {1,0,2};
    FGNSRSPNZ begin[4] = {0, 3, 6, 9};
    FGNSRSPDIM index[9] = {
        0, 1, 2,
        0, 1, 2,
        0, 1, 2,
    };
    double value[9] = {
        1, 4, 7, 
        2, 5, 8, 
        3, 6, 9,
    };
    double ex_value[9] = {
        .25,   1.,    1.5,
        1.,    4.,    6.,
        1./6., 2./3., 1.
    };
    double ac_value[9];
    SPPROJDATAptr pdata = NULL;
    const int maxmallocs = 100;
    long failcount;
    int clean_pass = 0;

    for (failcount=1; failcount<maxmallocs; failcount++) {
        clean_pass = 0;
        memcpy(ac_value, value, 9 * sizeof(double));
        FGNSRdbg_setmallocfail(failcount);

        status = FGNSRproj_spinit(lenvec, weights, verbose, &pdata);
        if (status != FGNSRSTAT_OK) {
            CuAssertIntEquals(tc, FGNSRSTAT_OOM, status);
            continue;
        }

        status = FGNSRproj_project_spmatrix(pdata, numvec, begin, index,
                ac_value, ub, diag_idx);
        /* No memory allocation in projection code at time of writing */
        assert(status == 0);

        /* Projection has succeeded, no OOM conditions have happened */
        clean_pass = 1;
        FGNSRassert_dblarr_abseql(tc,lenvec*numvec, ex_value, ac_value, 1e-8);
        FGNSRproj_spfree(&pdata);
        break;
    }

    /* If this doesn't hold, one needs to increase maxmallocs.  It means that
     * not every malloc has been tested for OOM cleanliness */
    assert(clean_pass==1);

}
#endif /* FGNSRDEBUG */

CuSuite* GetProjSuite() {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, TestInitFree);
    SUITE_ADD_TEST(suite, TestOneByOne1);
    SUITE_ADD_TEST(suite, TestOneByOne2);
    SUITE_ADD_TEST(suite, TestOneByOne3);
    SUITE_ADD_TEST(suite, TestTypicalProjection);
    SUITE_ADD_TEST(suite, TestProjectTwiceDense);
    SUITE_ADD_TEST(suite, TestZeroDiagWeight);
    SUITE_ADD_TEST(suite, TestJustClip);
    SUITE_ADD_TEST(suite, TestPermutation1);
    SUITE_ADD_TEST(suite, TestPermutation2);
    SUITE_ADD_TEST(suite, TestIdentity);
    SUITE_ADD_TEST(suite, TestNegValZeroWeight);
    SUITE_ADD_TEST(suite, TestNegDiag);
    SUITE_ADD_TEST(suite, TestZeroSetting);
    SUITE_ADD_TEST(suite, TestSorting);
    SUITE_ADD_TEST(suite, TestCutUb);
    SUITE_ADD_TEST(suite, TestVerboseSwitch);
    SUITE_ADD_TEST(suite, TestDiagElementMissing);
#ifdef FGNSRDEBUG
    SUITE_ADD_TEST(suite, TestOOMConditionsDense);
    SUITE_ADD_TEST(suite, TestOOMConditionsSparse);
#endif
    return suite;
}
