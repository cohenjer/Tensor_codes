#ifndef SPUTIL_H_INCDEF
#define SPUTIL_H_INCDEF

#include "fgnsrdef.h"
#include "util.h"

int FGNSRalloc_spmatrix(FGNSRSPDIM ncols, FGNSRSPNZ maxnnz, FGNSRSPNZ **begin,
        FGNSRSPDIM **index, double **value);

#endif /* SPUTIL_H_INCDEF */

