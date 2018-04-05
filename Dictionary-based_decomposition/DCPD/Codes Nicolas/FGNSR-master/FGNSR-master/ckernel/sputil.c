#include<stdlib.h>
#include "sputil.h"

int FGNSRalloc_spmatrix(FGNSRSPDIM ncols, FGNSRSPNZ maxnnz, FGNSRSPNZ **begin,
        FGNSRSPDIM **index, double **value)
{
    FGNSRSPDIM *ind=NULL;
    FGNSRSPNZ *beg=NULL;
    double *val=NULL;
    int status = 0;

    beg = (FGNSRSPNZ *)  FGNSRcalloc((ncols+1), sizeof(FGNSRSPNZ));
    if (beg==NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    ind = (FGNSRSPDIM *) FGNSRmalloc(maxnnz   * sizeof(FGNSRSPDIM));
    if (ind==NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }
    val = (double *)   FGNSRmalloc(maxnnz   * sizeof(double));
    if (val==NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    *begin = beg;
    *index = ind;
    *value = val;

TERMINATE:
    if (status == FGNSRSTAT_OOM) {
        if (beg!=NULL) FGNSRfree((void **) &beg);
        if (ind!=NULL) FGNSRfree((void **) &ind);
        if (val!=NULL) FGNSRfree((void **) &val);
    }

    return status;
}
