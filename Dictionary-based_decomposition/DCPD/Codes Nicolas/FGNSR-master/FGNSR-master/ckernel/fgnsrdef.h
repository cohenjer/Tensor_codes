#ifndef FGNSRDEF_H_INCDEF
#define FGNSRDEF_H_INCDEF

/* fgnsrdef.h -- status codes etc */

/*************** Sparse matrix API setup *******************/

#if defined(FGNSRAPI_SPARSEMAT_MODEL_LI)
typedef long FGNSRSPNZ;
typedef int FGNSRSPDIM;
#elif defined(FGNSRAPI_SPARSEMAT_MODEL_LL)
typedef long FGNSRSPNZ;
typedef long FGNSRSPDIM;
#elif defined(FGNSRAPI_SPARSEMAT_MODEL_SS)
typedef size_t FGNSRSPNZ;
typedef size_t FGNSRSPDIM;
#else
#ifndef FGNSRAPI_SPARSEMAT_MODEL_II
#    define FGNSRAPI_SPARSEMAT_MODEL_II 1
#endif
typedef int FGNSRSPNZ;
typedef int FGNSRSPDIM;
#endif

/*************** Error codes *******************************/

#define FGNSRSTAT_OK        0
#define FGNSRSTAT_OOM       1
#define FGNSRSTAT_DATAERR   2
#define FGNSRSTAT_LASTID    2

/*************** Global constants **************************/

#define FGNSR_NUM_INF       1.0e101
#define FGNSR_NUM_INVINF    1.0e-101

#endif /* FGNSRDEF_H_INCDEF */
