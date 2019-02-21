#ifndef PROJ_H_INCDEF
#define PROJ_H_INCDEF

#include "fgnsrdef.h"

struct projdata;
typedef struct projdata *PROJDATAptr;

struct spprojdata;
typedef struct spprojdata *SPPROJDATAptr;

/* Obtain projdata pointer for dense data */
int FGNSRproj_init(int lenvec, double * weights, int verbose,
        PROJDATAptr *pdata_ptr);

/* Obtain projdata pointer for sparse data */
int FGNSRproj_spinit(FGNSRSPDIM lenvec, double * weights, int verbose,
        SPPROJDATAptr *pdata_ptr);

/* Release projdata pointer for dense data */
void FGNSRproj_free(PROJDATAptr *pdata);

/* Release projdata pointer for sparse data */
void FGNSRproj_spfree(SPPROJDATAptr *pdata);

/* The entries in 'matrix' will be overwritten
 * matrix is flattened matrix of numvec vectors each of length lenvec */
int FGNSRproj_project_dmatrix(PROJDATAptr pdata, int numvec, double *matrix,
        double ub, int *diag_idx);

int FGNSRproj_project_spmatrix(SPPROJDATAptr pdata, FGNSRSPDIM numvec,
        FGNSRSPNZ *begin, FGNSRSPDIM *index, double *value, double ub,
        FGNSRSPDIM *diag_idx);

#endif /* PROJ_H_INCDEF */
