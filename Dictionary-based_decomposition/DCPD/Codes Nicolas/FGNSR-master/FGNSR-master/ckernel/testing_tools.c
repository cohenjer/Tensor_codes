#include<math.h>
#include<stdlib.h>

#include "testing_tools.h"
#include "sputil.h"

/* Check that two integer arrays are equal, return first position where they
 * differ  */
int FGNSRassert_intarr_eql(CuTest *tc, int len, int *a, int *b) {
    int k, diff=0;

    for (k=0; k<len; k++) {
        diff = a[k] - b[k];
        if (diff != 0)
            break;
    }

    CuAssertIntEquals(tc, 0, diff);
    return k;
}

void FGNSRassert_nzarr_eql(CuTest *tc, int len, FGNSRSPNZ *a, FGNSRSPNZ *b) {
    int k, diff=0;

    for (k=0; k<len; k++) {
        diff = a[k] - b[k];
        if (diff != 0)
            break;
    }

    CuAssertIntEquals(tc, 0, diff);
}

void FGNSRassert_dimarr_eql(CuTest *tc, int len, FGNSRSPDIM *a, FGNSRSPDIM *b) {
    int k, diff=0;

    for (k=0; k<len; k++) {
        diff = a[k] - b[k];
        if (diff != 0)
            break;
    }

    CuAssertIntEquals(tc, 0, diff);
}

/* Assert that max(abs(a-b)) <= atol */
void FGNSRassert_dblarr_abseql(CuTest *tc, int len, double *a, double *b,
        double atol)
{
    int k;
    double diff = 0.0, maxdiff = 0.0;

    for (k=0; k<len; k++) {
        diff = fabs(a[k] - b[k]);
        maxdiff = (diff > maxdiff) ? diff : maxdiff;
    }

    CuAssertDblEquals(tc, 0.0, maxdiff, atol);
}

int FGNSRdtosp(int numvec, int lenvec, double *dmatrix,
        FGNSRSPNZ **begin, FGNSRSPDIM **index, double **value)
{
    int nnz = 0;
    int k;
    int status = 0;
    FGNSRSPNZ *beg;
    FGNSRSPDIM *ind;
    double *val;
    FGNSRSPNZ which;
    FGNSRSPDIM cidx, ridx;
    double this_value;

    /* How much space do we need? */
    for (k=0; k<numvec*lenvec; k++) {
        if (dmatrix[k] != 0.0)
            nnz++;
    }

    /* Allocate */
    status = FGNSRalloc_spmatrix(numvec, nnz, begin, index, value);
    if (status)
        goto TERMINATE;

    /* Allocation succeeded, copy data */
    beg = *begin;
    ind = *index;
    val = *value;

    which = 0;
    for (cidx=0; cidx<numvec; cidx++) {
        beg[cidx] = which;
        for (ridx=0; ridx<lenvec; ridx++) {
            this_value = dmatrix[cidx*lenvec + ridx];
            if (this_value != 0.0) {
                ind[which] = ridx;
                val[which++] = this_value;
            }
        }
    }
    beg[numvec] = which;

TERMINATE:

    return status;
}

int FGNSRsptod(FGNSRSPNZ *begin, FGNSRSPDIM *index, double *value,
        int numvec, int lenvec, double **dmatrix)
{
    double *a = NULL;
    int status = 0;
    FGNSRSPDIM ridx, cidx;
    FGNSRSPNZ which;

    a = (double *) FGNSRcalloc(lenvec*numvec, sizeof(double));
    if (a==NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }
    *dmatrix = a;


    /* Copy into dense matrix, assume memory is set to zero */
    for (cidx=0; cidx<numvec; cidx++) {
        for (which=begin[cidx]; which<begin[cidx+1]; which++) {
            ridx = index[which];
            a[cidx*lenvec + ridx] = value[which];
        }
    }
    

TERMINATE:
    return status;
}

int FGNSRpromote_spdim(int n, int *a, FGNSRSPDIM **b) {
   FGNSRSPDIM *c=NULL;
   int k, status = 0;

   c = (FGNSRSPDIM *) malloc(n * sizeof(FGNSRSPDIM));
   if (c==NULL) {
       status = FGNSRSTAT_OOM;
       goto TERMINATE;
   }

   for (k=0; k<n; k++) {
       c[k] = (FGNSRSPDIM) a[k];
   }
   
   *b = c;

TERMINATE:
   return status;

}

#ifdef FGNSRDEBUG

/* Call a given function repeatedly until it returns NMSTAT_OK.  Initially, the
 * first malloc is guaranteed to result in a NULL pointer.  Each time the given
 * function is called, it is checked whether the returned value is FGNSRSTAT_OK
 * or FGNSRSTAT_OOM.  In the latter case, the next OOM condition is triggered at
 * the next call to malloc. v*/

void FGNSRoom_loop(CuTest *tc, long maxmallocs, int (*fun)(void *fun_data),
        void *fun_data) {
    int failcount;
    int status;
    int clean_pass = 0;

    for (failcount=1; failcount<maxmallocs; failcount++) {
        clean_pass = 0;
        FGNSRdbg_setmallocfail(failcount);
        
        status = fun(fun_data);
        if (status == 0) {
            clean_pass = 1;
        }
        else if (status == FGNSRSTAT_OOM) {
            clean_pass = 0;
        }
        else {
            /* Invalid status */
            CuAssertIntEquals(tc, 0, 1);
        }

        /* Projection has succeeded, no OOM conditions have happened */
        if ( clean_pass == 1) {
            break;
        }
    }

    /* If this doesn't hold, one needs to increase maxmallocs.  It means that
     * not every malloc has been tested for OOM cleanliness */
    CuAssertIntEquals(tc, 1, clean_pass);

    /* Reset fail counter */
    FGNSRdbg_setmallocfail(0);
}
#endif /* FGNSRDEBUG */
