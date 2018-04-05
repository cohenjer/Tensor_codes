#ifndef TESTING_TOOLS_H_INCDEF
#define TESTING_TOOLS_H_INCDEF

#include "fgnsrdef.h"
#include "CuTest.h"

#ifdef FGNSRDEBUG
#include "debug_tools.h"
#endif /* FGNSRDEBUG */

/* testing_tools.h
 *
 * Some extensions to the CuTest framework
 */

/* Check that two integer arrays are equal */
int FGNSRassert_intarr_eql(CuTest *tc, int len, int *a, int *b);
void FGNSRassert_dimarr_eql(CuTest *tc, int len, FGNSRSPDIM *a, FGNSRSPDIM *b);
void FGNSRassert_nzarr_eql(CuTest *tc, int len, FGNSRSPNZ *a, FGNSRSPNZ *b);

/* Assert that max(abs(a-b)) <= atol */
void FGNSRassert_dblarr_abseql(CuTest *tc, int len, double *a, double *b,
        double atol);

int FGNSRdtosp(int numvec, int lenvec, double *dmatrix,
        FGNSRSPNZ **begin, FGNSRSPDIM **index, double **value); 

int FGNSRsptod(FGNSRSPNZ *begin, FGNSRSPDIM *index, double *value,
        int numvec, int lenvec, double **dmatrix);

int FGNSRpromote_spdim(int n, int *a, FGNSRSPDIM **b);

#ifdef FGNSRDEBUG
/* Run code in a loop so as to trigger all possible OOM conditions */
void FGNSRoom_loop(CuTest *tc, long maxmallocs, int (*fun)(void *fun_data),
        void *fun_data);
#endif /*  FGNSRDEBUG */


#endif /* TESTING_TOOLS_H_INCDEF */
