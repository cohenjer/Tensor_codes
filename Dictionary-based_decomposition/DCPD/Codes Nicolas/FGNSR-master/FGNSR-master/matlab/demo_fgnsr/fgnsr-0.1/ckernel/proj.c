#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#include "proj.h"
#include "fgnsrdef.h"
#include "util.h"

#define XMAX(X, Y)  ((X) < (Y) ? Y : X) 
#define XMIN(X, Y)  ((X) < (Y) ? X : Y) 

/* Epsilon for relative feasibility of w_j * x_i <= w_i * x_j */
#define EP_FEAS_REL 1.e-8
#define EP_FEAS_ABS 1.e-8

/* TODO IDEA FOR OPTIMIZATION The values of the to-be-projected matrix do not
 * change at random.  We add a scaled gradient to a feasible point. Since the
 * gradient may have structure (sparse, only a few non-zero rows, sign of
 * diagonal entries known, etc), it may be of advantage to keep sorted arrays
 * throughout the iteration and process only the changes that destroy the
 * current ordering. */

/* TODO IDEA FOR OPTIMIZATION Depending on the length of a sparse column to
 * have a local copy of only the weights actually needed instead of all
 * weights. This would avoid all the indirect addressing over the row index
 * array in the CSC structure at the cost of some thread local memory and
 * copies.
 */

/************************************************************************/
/* Data structures */
/************************************************************************/

/* This is used to sort two different arrays w.r.t only one of the arrays */
struct dpair {
    int index;    /* Index k */
    double value; /* Value x_k / w_k */
};
typedef struct dpair DPAIR, *DPAIRptr;

struct sppair {
    FGNSRSPDIM pos; /* Position in compressed column corresponding to 'value' */
    double value;   /* Value x_k / w_k */
};
typedef struct sppair SPPAIR, *SPPAIRptr;

struct projdata {
    int lenvec;
    int verbose;
    int numthreads;
    double *weights;
    DPAIRptr *ratios; /* private data array for each thread */
};
typedef struct projdata PROJDATA;
/* pointer typedefs for this struct placed in proj.h */

struct spprojdata {
    int verbose;
    int numthreads;
    double *weights;
    SPPAIRptr *ratios; /* private data array for each thread */
    int *flags; /* thread private flag space, e.g. statuses */
};
typedef struct spprojdata SPPROJDATA;
/* pointer typedefs for this struct placed in proj.h */

/************************************************************************/
/* Prototypes for static functions */
/************************************************************************/

#ifdef _OPENMP
/* Setup openmp environment */
static int setup_openmp(int *numthreads);
#endif

/* Allocation function */
static int alloc_projdata(int lenvec, int numthreads, int verbose,
        PROJDATAptr *pdata_ptr);

static int alloc_projdata_sparse(FGNSRSPDIM lenvec, int numthreads, int verbose,
        SPPROJDATAptr *pdata_ptr);

/* Project a given dense vector.
 *
 * n        -- length of the vector
 * vec      -- the vector
 * diag_idx -- position of the diagonal element in the vector
 * ub       -- upper bound for the projected diagonal element 
 *
 * TODO Describe the polytope and the type of projection
 */
static void project_dvec(int n, double *vec, double *weights, int diag_idx,
        double ub, DPAIR* ratios);

/* veclen -- number of nonzeros in this column
 * rind   -- row indices of the vector
 * value  -- nz values of the vector
 * weights
 * diag_idx -- row index of the diagonal element in the vector
 * ub
 * ratios -- unpopulated work array
 */
static void project_spvec(FGNSRSPDIM veclen, FGNSRSPDIM *rind, double *vec,
        double *weights, FGNSRSPDIM diag_idx, double ub, SPPAIR *ratios,
        int *status);

/* Push vector to feasible orthant, copy infeasible entries to DPAIR array for
 * thresholding phase.  Returns the number of violations found and stored to
 * DPAIR. */
static int preclean_dvec(int n, double *vec, double *weights, int diag_idx,
        double ub, DPAIR* ratios);
static int preclean_spvec(FGNSRSPDIM n, FGNSRSPDIM *rind, double *vec,
        double *weights, FGNSRSPDIM diag_pos, double ub, SPPAIR* ratios);

/* Determine optimal clipped value for diagonal element by thresholding
 *
 * Return number of tight upper bounds
 */
static void iterthres_dvec(int numelm, DPAIR* ratios, double *vec,
        double *weights, int diag_idx, double ub);
static void iterthres_spvec(FGNSRSPDIM numelem, SPPAIR* ratios, FGNSRSPDIM *rind,
        double *vec, double *weights, FGNSRSPDIM diag_pos, double ub);

/* Clip whole vector at given value v for the diagonal element.  Violated
 * bounds will be corrected as well. */
static void clipatval_dvec(int n, double *vec, double *weights, int diag_idx,
        double v);
static void clipatval_spvec(FGNSRSPDIM n, FGNSRSPDIM *rind, double *vec,
        double *weights, FGNSRSPDIM diag_pos, double v);

/* Move vector to feasible orthant, that is:
 *
 * vec[diag_idx] <= ub
 * vec >= 0
 *
 * This function is dedicated to the case where weight[diag_idx] == 0 (or tiny).
 */
static void enforce_bounds(int n, double *vec, int diag_idx, double ub);

static void enforce_bounds_spvec(FGNSRSPDIM nnz_vec, double *vec,
        FGNSRSPDIM diag_pos, double ub);

/* Initialize heap order on 'ratios' by siftdown in reverse order */
static void init_heaporder(DPAIR *ratios, int num_elem);

/* NOTE: can be made inline */
/* Sift down element from node 'index'*/
static void siftdown(DPAIR *ratios, int index, int numelem);

static void init_heaporder_sppair (SPPAIR* store, FGNSRSPDIM numelem);

/* NOTE: can be made inline */
static void siftdown_sppair(SPPAIR* s, FGNSRSPDIM index, FGNSRSPDIM numelem);


/************************************************************************/
/* API functions */
/************************************************************************/

/* Obtain projdata pointer */
int FGNSRproj_init(int lenvec, double * weights, int verbose,
        PROJDATAptr *pdata_ptr)
{
    int numthreads = -1;
    int status = 0;

#ifdef _OPENMP
    status = setup_openmp(&numthreads);
    if (status != 0) {
        goto TERMINATE;
    }
#else
    numthreads = 1;
#endif

    status = alloc_projdata(lenvec, numthreads, verbose, pdata_ptr);
    if (status != 0) {
        goto TERMINATE;
    }

    (*pdata_ptr)->weights = weights;

TERMINATE:
    if (status != 0) {
        /* *pdata_ptr == NULL is possible but safe */
        FGNSRproj_free(pdata_ptr);
    }
    return status;

}

int FGNSRproj_spinit(FGNSRSPDIM lenvec, double * weights, int verbose, SPPROJDATAptr *pdata_ptr){
    int numthreads = -1;
    int status = 0;

#ifdef _OPENMP
    status = setup_openmp(&numthreads);
    if (status != 0) {
        goto TERMINATE;
    }
#else
    numthreads = 1;
#endif

    status = alloc_projdata_sparse(lenvec, numthreads, verbose, pdata_ptr);
    if (status != 0) {
        goto TERMINATE;
    }

    (*pdata_ptr)->weights = weights;

TERMINATE:
    if (status != 0) {
        /* *pdata_ptr == NULL is possible but safe */
        FGNSRproj_spfree(pdata_ptr);
    }
    return status;

}


void FGNSRproj_free(PROJDATAptr *pdata_ptr)
{
    PROJDATAptr pdata = *pdata_ptr;

    if (pdata != NULL) {
        pdata->weights = NULL; /* External memory, no free */
        if (pdata->ratios != NULL) {
            /* Only first thread's data points to malloc'ed memory */
            FGNSRfree((void **) pdata->ratios);
            FGNSRfree((void **) &(pdata->ratios));
        }
        FGNSRfree((void **) pdata_ptr);
    }
    else {
        /* Nothing to free, which is OK */
    }
}

void FGNSRproj_spfree(SPPROJDATAptr *pdata_ptr)
{
    SPPROJDATAptr pdata = *pdata_ptr;

    if (pdata != NULL) {
        pdata->weights = NULL; /* External memory, no free */
        if (pdata->ratios != NULL) {
            /* Only first thread's data points to malloc'ed memory */
            FGNSRfree((void **) pdata->ratios);
            FGNSRfree((void **) &(pdata->ratios));
            FGNSRfree((void **) &(pdata->flags));
        }
        FGNSRfree((void **) pdata_ptr);
    }
    else {
        /* Nothing to free, which is OK */
    }
}


int FGNSRproj_project_dmatrix(PROJDATAptr pdata, int numvec, double *matrix,
        double ub, int *diag_idx)
{
    int tid = 0;
    int i;
    DATUM t_start, t_end;
    int status = FGNSRSTAT_OK;

    int lenvec = pdata->lenvec;
    double *weights = pdata->weights;

    FGNSRtimestamp(&t_start);
#pragma omp parallel private(tid, i)
    {
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif

#pragma omp for nowait
    for (i=0; i<numvec; i++) {
        project_dvec(lenvec, &(matrix[i*lenvec]), weights, diag_idx[i], ub,
                (pdata->ratios)[tid]);
    }
    }
    FGNSRtimestamp(&t_end);

    if (pdata->verbose) {
        printf("Elapsed CPU time %.2f,  WC time %.2f (%d threads)\n",
                FGNSRcputime(t_start,t_end), FGNSRwctime(t_start, t_end),
                pdata->numthreads);
    }

    return status;
}

int FGNSRproj_project_spmatrix(SPPROJDATAptr pdata, FGNSRSPDIM numvec,
        FGNSRSPNZ *begin, FGNSRSPDIM *index, double *value, double ub,
        FGNSRSPDIM *diag_idx)
{
    int tid = 0;
    FGNSRSPDIM i;
    DATUM t_start, t_end;
    int *tidstatus = pdata->flags; /* Use flags as per-thread status */
    int status = FGNSRSTAT_OK;

    double *weights = pdata->weights;

    FGNSRtimestamp(&t_start);
#pragma omp parallel private(tid, i)
    {
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    tidstatus[tid] = 0;

#pragma omp for nowait
    for (i=0; i<numvec; i++) {
        /* Safe cast: no column can be longer than SPDIM, although begin[k] may
         * be larger */
        FGNSRSPDIM veclen = (FGNSRSPDIM) (begin[i+1] - begin[i]);
        project_spvec(veclen, &index[begin[i]], &value[begin[i]],
                weights, diag_idx[i], ub, (pdata->ratios)[tid],
                &tidstatus[tid]);
    }
    }
    FGNSRtimestamp(&t_end);

    for (i=0; i<pdata->numthreads; i++) {
        if (tidstatus[i] != FGNSRSTAT_OK) {
            /* Some thread has signaled error, we only signal tid-first error */
            status = tidstatus[i];
            break;
        }

    }

    if (pdata->verbose) {
        printf("Elapsed CPU time %.2f,  WC time %.2f (%d threads)\n",
                FGNSRcputime(t_start,t_end), FGNSRwctime(t_start, t_end),
                pdata->numthreads);
    }
    return status;
}

/************************************************************************/
/* Static functions */
/************************************************************************/

#ifdef _OPENMP
/* Setup openmp environment */
static int setup_openmp(int *numthreads) {
    int status = FGNSRSTAT_OK;
    int maxthreads;

    /* No dynamic thread assignments */
    omp_set_dynamic(0);

    /* No nested parallel regions */
    omp_set_nested(0);

    maxthreads = omp_get_max_threads();
    if (*numthreads <= 0 || *numthreads > maxthreads) {
        /* Automatic setting, use max threads */
        *numthreads = maxthreads;
    }
    omp_set_num_threads(*numthreads);

    return status;
}
#endif

/* Obtain projdata pointer */
static int alloc_projdata(int lenvec, int numthreads, int verbose,
        PROJDATAptr *pdata_ptr) {
    int status, k;
    DPAIRptr *ratios = NULL, tmpchunk = NULL;
    PROJDATAptr pdata = NULL;

    status = FGNSRSTAT_OK;
    
    pdata = (PROJDATAptr) FGNSRmalloc(sizeof(PROJDATA));
    if (pdata == NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    ratios = (DPAIRptr *) FGNSRmalloc (numthreads * sizeof(DPAIRptr));
    if (ratios == NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    tmpchunk = (DPAIR *) FGNSRmalloc (numthreads * (lenvec - 1) * sizeof(DPAIR));
    if (tmpchunk == NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    /* Split chunk apart into per-thread data */
    for (k=0; k<numthreads; k++) {
        ratios[k] = & ( tmpchunk[k * (lenvec - 1)] );
    }

    /* All mallocs have succeeded at this point, set up data structure */
    pdata->lenvec     = lenvec;
    pdata->numthreads = numthreads;
    pdata->verbose    = verbose;
    pdata->ratios     = ratios;
    *pdata_ptr        = pdata;

TERMINATE:
    if (status == FGNSRSTAT_OOM) {
        FGNSRfree((void **) &pdata);
        FGNSRfree((void **) &ratios);
        FGNSRfree((void **) &tmpchunk);
    }

    return status;
}

static int alloc_projdata_sparse(FGNSRSPDIM lenvec, int numthreads, int verbose,
        SPPROJDATAptr *pdata_ptr) {
    int status, k;
    SPPAIRptr *ratios = NULL, tmpchunk = NULL;
    int *flags = NULL;
    SPPROJDATAptr pdata = NULL;

    status = FGNSRSTAT_OK;
    
    pdata = (SPPROJDATAptr) FGNSRmalloc(sizeof(SPPROJDATA));
    if (pdata == NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    ratios = (SPPAIRptr *) FGNSRmalloc (numthreads * sizeof(SPPAIRptr));
    if (ratios == NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    flags = (int *) FGNSRmalloc (numthreads * sizeof(int));
    if (flags == NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    tmpchunk = (SPPAIR *) FGNSRmalloc (numthreads * (lenvec - 1) * sizeof(SPPAIR));
    if (tmpchunk == NULL) {
        status = FGNSRSTAT_OOM;
        goto TERMINATE;
    }

    /* Split chunk apart into per-thread data */
    for (k=0; k<numthreads; k++) {
        ratios[k] = & ( tmpchunk[k * (lenvec - 1)] );
    }

    /* All mallocs have succeeded at this point, set up data structure */
    pdata->numthreads = numthreads;
    pdata->verbose    = verbose;
    pdata->ratios     = ratios;
    pdata->flags      = flags;
    *pdata_ptr        = pdata;

TERMINATE:
    if (status == FGNSRSTAT_OOM) {
        FGNSRfree((void **) &pdata);
        FGNSRfree((void **) &ratios);
        FGNSRfree((void **) &flags);
        FGNSRfree((void **) &tmpchunk);
    }

    return status;
}

static void project_dvec(int n, double *vec, double *weights, int diag_idx,
        double ub, DPAIR* ratios)
{
    int num_vio = 0; /* Number of violations saved in 'ratios' */
    double diag_weight, diag_value;

    diag_weight = weights[diag_idx];
    if (diag_weight < FGNSR_NUM_INVINF) {
        
        /* Diagonal weight is essentially zero.  We only need to assert bound
         * feasibility of the vector */
        enforce_bounds(n, vec, diag_idx, ub);
        return;
    }

    diag_value = vec[diag_idx];
    if (diag_value >= (1. - EP_FEAS_REL) * ub) {
        /* In the optimal projection the diagonal element will be sitting at
         * ub. We clip the vector at ub and enforce feasibility constraints */
        clipatval_dvec(n, vec, weights, diag_idx, ub);
        return;
    }

    /* At this point we cleared out the trivial cases and assume diag_weight !=
     * 0.0.  Now we have to identify the violated feasibility constraints.
     * Since we move vec to memory anyway, we enforce bounds as we search for
     * violations.*/
        
    num_vio = preclean_dvec(n, vec, weights, diag_idx, ub, ratios);
    if (num_vio == 0) {
        /* This vector is feasible except for possibly the diagonal element.
         * Since we are asserting above that the upper bound is satisfied, the
         * only source of infeasibility is negativity.  By setting the diagonal
         * value to zero in this case, we do not introduce new infeasibilities,
         * as 0 is a lower bound for all components. */
        
        vec[diag_idx] = XMAX(vec[diag_idx], 0.0);
        return;
    }

    iterthres_dvec(num_vio, ratios, vec, weights, diag_idx, ub);
}

static void project_spvec(FGNSRSPDIM n, FGNSRSPDIM *rind, double *vec,
        double *weights, FGNSRSPDIM diag_idx, double ub, SPPAIR *ratios,
        int *status)
{
    FGNSRSPDIM num_vio = 0; /* Number of violations saved in 'ratios' */
    FGNSRSPDIM k, diag_pos;
    double diag_weight, diag_value;
    int diag_found = 0;

    /* Find diagonal entry in sparse vector */
    for (k=0; k<n; k++) {
        if (rind[k] == diag_idx) {
            diag_pos = k;
            diag_found = 1;
            break;
        }
    }
    if (!diag_found) {
        /* This column has no diagonal element, signal error and quit */
        *status = FGNSRSTAT_DATAERR;
        return;
    }

    diag_weight = weights[diag_idx];
    if (diag_weight < FGNSR_NUM_INVINF) {
        /* Diagonal weight is essentially zero.  We only need to assert bound
         * feasibility of the vector */
        enforce_bounds_spvec(n, vec, diag_pos, ub);
        return;
    }

    diag_value = vec[diag_pos];
    if (diag_value >= (1. - EP_FEAS_REL) * ub) {
        /* In the optimal projection the diagonal element will be sitting at
         * ub. We clip the vector at ub and enforce feasibility constraints */
        clipatval_spvec(n, rind, vec, weights, diag_pos, ub);
        return;
    }

    /* At this point we cleared out the trivial cases and assume diag_weight !=
     * 0.0.  Now we have to identify the violated feasibility constraints.
     * Since we move vec to memory anyway, we enforce bounds as we search for
     * violations.*/
        
    num_vio = preclean_spvec(n, rind, vec, weights, diag_pos, ub, ratios);
    if (num_vio == 0) {
        /* This vector is feasible except for possibly the diagonal element.
         * Since we are asserting above that the upper bound is satisfied, the
         * only source of infeasibility is negativity.  By setting the diagonal
         * value to zero in this case, we do not introduce new infeasibilities,
         * as 0 is a lower bound for all components. */
        
        vec[diag_pos] = XMAX(vec[diag_pos], 0.0);
        return;
    }

    iterthres_spvec(num_vio, ratios, rind, vec, weights, diag_pos, ub);
    return;
}

static int preclean_dvec(int n, double *vec, double *weights, int diag_idx,
        double ub, DPAIR* ratios)
{
    int numvio = 0;
    double factor;
    int k;

    /* We need to identify those indices k such that
     *
     *   weights[diag_idx] * vec[k] > weights[k] * vec[diag_idx]
     *
     * so we check vec[k] > vec[diag_idx] * weights[k] / weights[diag_idx]
     * since we can assert that weights[diag_dix] is sufficiently large.
     */

    assert(weights[diag_idx] > FGNSR_NUM_INVINF);
    assert(vec[diag_idx] < ub);

    factor = vec[diag_idx] / weights[diag_idx];

    for (k=0; k<n; k++) {
        if (k==diag_idx)
            continue;

        if (weights[k] == 0.0) {
            /* Corner case, no restrictions in this case TODO Treatment of this
             * case could be optimized since *all* values for this k-th row of
             * the matrix will need to go to zero. */

            vec[k] = 0.0;
        }
        else if (vec[k] <= 0.0) {
            /* This component will always be projected on 0.  No need to
             * consider it for thresholding, since the contribution to the
             * weighted sum will always be negative. */
            
            vec[k] = 0.0;
        }
        else if (vec[k] > weights[k] * factor + EP_FEAS_ABS) {
            /* k-th constraint for is violated. */
            assert(weights[k] > 0);
            ratios[numvio].value = vec[k]/weights[k];
            ratios[numvio].index = k;
            numvio++;
        }
        assert(k != diag_idx); /* Should have been treated separately */
    }

    return numvio;
}

static int preclean_spvec(FGNSRSPDIM n, FGNSRSPDIM *rind, double *vec,
        double *weights, FGNSRSPDIM diag_pos, double ub, SPPAIR* ratios)
{
    int numvio = 0;
    double factor;
    double this_weight;
    FGNSRSPDIM k;

    /* We need to identify those indices k such that
     *
     *   weights[diag_idx] * vec[k] > weights[k] * vec[diag_idx]
     *
     * so we check vec[k] > vec[diag_idx] * weights[k] / weights[diag_idx]
     * since we can assert that weights[diag_dix] is sufficiently large.
     */

    assert(weights[rind[diag_pos]] > FGNSR_NUM_INVINF);
    assert(vec[diag_pos] < ub);

    factor = vec[diag_pos] / weights[rind[diag_pos]];

    for (k=0; k<n; k++) {
        if (k==diag_pos)
            continue;
        this_weight = weights[rind[k]];

        if (this_weight == 0.0) {
            /* Corner case, no restrictions in this case TODO Treatment of this
             * case could be optimized since *all* values for this k-th row of
             * the matrix will need to go to zero. */

            vec[k] = 0.0;
        }
        else if (vec[k] <= 0.0) {
            /* This component will always be projected on 0.  No need to
             * consider it for thresholding, since the contribution to the
             * weighted sum will always be negative. */
            
            vec[k] = 0.0;
        }
        else if (vec[k] > this_weight * factor + EP_FEAS_ABS) {
            /* rind[k]-th constraint is violated. */
            assert(this_weight > 0);
            ratios[numvio].value = vec[k]/this_weight;
            ratios[numvio].pos = k;
            numvio++;
        }
        assert(k != diag_pos); /* Should have been treated separately */
    }

    return numvio;
}


/* iterthres_dvec
 *
 * Compute optimal clipping value for a dense vector.
 *
 * Assumptions:
 *
 * - Reference point is sitting in leading position of 'vec' and 'weights'.
 * - 'ratios' is descendingly ordered array of (n-1) ratios 
 *     vec[k] * weights[0] / weights[k],  1<= k <= n
 *
 * TODO proper treatment of inf / 0.0 / NaN cases
 */
static void iterthres_dvec(int numelem, DPAIR* ratios, double *vec,
        double *weights, int diag_idx, double ub)
{
    double nom, den, clipval, ratio, diag_weight, diag_value, diag_iweight;
    double diag_ratio;
    int k, index;
    int len_heap;
    DPAIR tmp;

    /* Initialize heap of ratios x_k/w_k decreasingly (for violated constraints
     * of index k).  Following this operation, ratios[0] holds minimal node. */

    len_heap = numelem;
    init_heaporder(ratios, len_heap);

    diag_weight  = weights[diag_idx];
    diag_iweight = 1.0 / diag_weight;
    diag_value  = vec[diag_idx];

    nom = diag_value  * diag_weight;
    den = diag_weight * diag_weight;
    clipval = diag_weight * nom/den; /* Initially == vec[diag_idx] */

    /* We know that there is some space to move towards the ub.  */
    assert(clipval < ub);

    /* First pass through the data: determine clipping value. If the loop
     * condition fails, we have pushed clipvalue far enough so that all
     * remaining constraints become feasible. This cannot happen in the first
     * iteration. */
    while ( (len_heap > 0) && (clipval < ub) &&
            (ratios[0].value > clipval * diag_iweight) )
    {
        /* Update clipval by consideration of one more violated constraint */
        index = ratios[0].index;
        nom += weights[index] * vec[index];
        den += weights[index] * weights[index];
        clipval = diag_weight * nom/den;

        /* Move current max to tail (aka delete) and siftdown last leave to
         * restore heap order */
        len_heap--;
        tmp = ratios[0];
        ratios[0] = ratios[len_heap];
        ratios[len_heap] = tmp;
        siftdown(ratios, 0, len_heap);
    }

    /* The determined clipping value may well be negative or exceed the bound */
    clipval = XMIN(clipval, ub);
    clipval = XMAX(clipval, 0.0);

    /* At this point we have determined an optimal clipping value */
    assert( (0.0 <= clipval) && (clipval <= ub) );

    vec[diag_idx] = clipval;
    diag_ratio = clipval * diag_iweight;

    /* Second pass through data in order to clip all entries of vec that exceed
     * the weighted optimal clipping value */

    /* First chunk: These constraints *may* be violated by our choice of
     * clipval*/
    for (k=0; k<len_heap; k++) {
        ratio = ratios[k].value;
        index = ratios[k].index;
        if (ratio > diag_ratio) {
            vec[index] = weights[index] * diag_ratio;
        }
    }

    /* Second chunk: These constraints *are* violated by our choice of
     * clipval*/
    for (k=len_heap; k<numelem; k++) {
        ratio = ratios[k].value;
        index = ratios[k].index;
        vec[index] = weights[index] * diag_ratio;
    }

    return;
}

static void iterthres_spvec(FGNSRSPDIM numelem, SPPAIR* ratios, FGNSRSPDIM *rind,
        double *vec, double *weights, FGNSRSPDIM diag_pos, double ub)
{
    double nom, den, clipval, ratio, diag_weight, diag_value, diag_iweight;
    double diag_ratio;
    FGNSRSPDIM k, index, pos;
    FGNSRSPDIM len_heap;
    SPPAIR tmp;

    /* Initialize heap of ratios x_k/w_k decreasingly (for violated constraints
     * of index k).  Following this operation, ratios[0] holds minimal node. */

    len_heap = numelem;
    init_heaporder_sppair(ratios, len_heap);

    diag_weight  = weights[rind[diag_pos]];
    diag_iweight = 1.0 / diag_weight;
    diag_value  = vec[diag_pos];

    nom = diag_value  * diag_weight;
    den = diag_weight * diag_weight;
    clipval = diag_weight * nom/den; /* Initially == vec[diag_idx] */

    /* We know that there is some space to move towards the ub.  */
    assert(clipval < ub);

    /* First pass through the data: determine clipping value. If the loop
     * condition fails, we have pushed clipvalue far enough so that all
     * remaining constraints become feasible. This cannot happen in the first
     * iteration. */
    while ( (len_heap > 0) && (clipval < ub) &&
            (ratios[0].value > clipval * diag_iweight) )
    {
        /* Update clipval by consideration of one more violated constraint */
        pos = ratios[0].pos;
        index = rind[pos];
        nom += weights[index] * vec[pos];
        den += weights[index] * weights[index];
        clipval = diag_weight * nom/den;

        /* Move current max to tail (aka delete) and siftdown last leave to
         * restore heap order */
        len_heap--;
        tmp = ratios[0];
        ratios[0] = ratios[len_heap];
        ratios[len_heap] = tmp;
        siftdown_sppair(ratios, 0, len_heap);
    }

    /* The determined clipping value may well be negative or exceed the bound */
    clipval = XMIN(clipval, ub);
    clipval = XMAX(clipval, 0.0);

    /* At this point we have determined an optimal clipping value */
    assert( (0.0 <= clipval) && (clipval <= ub) );

    vec[diag_pos] = clipval;
    diag_ratio = clipval * diag_iweight;

    /* Second pass through data in order to clip all entries of vec that exceed
     * the weighted optimal clipping value */

    /* First chunk: These constraints *may* be violated by our choice of
     * clipval*/
    for (k=0; k<len_heap; k++) {
        ratio = ratios[k].value;
        pos = ratios[k].pos;
        index = rind[pos];

        if (ratio > diag_ratio) {
            vec[pos] = weights[index] * diag_ratio;
        }
    }

    /* Second chunk: These constraints *are* violated by our choice of
     * clipval*/
    for (k=len_heap; k<numelem; k++) {
        ratio = ratios[k].value;
        pos = ratios[k].pos;
        index = rind[pos];
        vec[pos] = weights[index] * diag_ratio;
    }

    return;
}


static void clipatval_dvec(int n, double *vec, double *weights, int diag_idx,
        double v)
{
    int k;
    double value1, value2;

    assert(v >= 0.0);
    assert(weights[diag_idx] > FGNSR_NUM_INVINF);

    /* Set diagonal value */
    vec[diag_idx] = v;

    /* Precompute ratio to save divisions */
    v /= weights[diag_idx];

    for (k=0; k<n; k++) {
        /* For k==diag_idx value1==v*weights[diag_idx]==value2, so the loop
         * pass is harmless and has no effect but simplifies the logic
         */
        value1 = XMAX(0.0, vec[k]);
        value2 = v * weights[k];
        vec[k] = XMIN(value1, value2);
    }
}

static void clipatval_spvec(FGNSRSPDIM n, FGNSRSPDIM *rind, double *vec,
        double *weights, FGNSRSPDIM diag_pos, double v)
{
    FGNSRSPDIM k;
    double value1, value2;

    assert(v >= 0.0);
    assert(weights[rind[diag_pos]] > FGNSR_NUM_INVINF);

    /* Set diagonal value */
    vec[diag_pos] = v;

    /* Precompute ratio to save divisions */
    v /= weights[rind[diag_pos]];

    for (k=0; k<n; k++) {
        /* For k==diag_pos value1==v*weights[diag_idx]==value2, so the loop
         * pass is harmless and has no effect but simplifies the logic
         */
        value1 = XMAX(0.0, vec[k]);
        value2 = v * weights[rind[k]];
        vec[k] = XMIN(value1, value2);
    }
}

static void enforce_bounds(int n, double *vec, int diag_idx, double ub) {
    int k;

    vec[diag_idx] = XMIN(vec[diag_idx], ub);
    for (k=0; k<n; k++) {
        vec[k] = XMAX(vec[k], 0.0);
    }
}

static void enforce_bounds_spvec(FGNSRSPDIM nnz_vec, double *vec,
        FGNSRSPDIM diag_pos, double ub) {
    FGNSRSPDIM k;

    vec[diag_pos] = XMIN(vec[diag_pos], ub);
    for (k=0; k<nnz_vec; k++) {
        vec[k] = XMAX(vec[k], 0.0);
    }
}



/**********************************************************************/


static void
init_heaporder  (DPAIR* store, int numelem) {
    int i;

    for (i=numelem/2-1; i>=0; i--) {
        siftdown (store, i, numelem);
    }
}



static void siftdown (DPAIR* s, int index, int numelem) {
    int cur_node = index;
    int child_node;
    while ((child_node = 2*cur_node+1) < numelem) {
        /* Then cur_node has at least one child */
        int max_index=cur_node;
        DPAIR tmp;
        if (s[child_node].value > s[max_index].value) {
            max_index = child_node;
        }

        if (++child_node < numelem &&
            s[child_node].value > s[max_index].value)
        {
            /* Then cur_node has two children and we select the max child */
            max_index = child_node;
        }

        if (max_index == cur_node) {
            /* Element is no smaller than all children */
            break;
        }

        tmp = s[max_index];
        s[max_index] = s[cur_node];
        s[cur_node] = tmp;

        cur_node = max_index;
    }
}

static void
init_heaporder_sppair (SPPAIR* store, FGNSRSPDIM numelem) {
    FGNSRSPDIM i;

    for (i=numelem/2; i>0; i--) {
        siftdown_sppair (store, i-1, numelem);
    }
}



static void siftdown_sppair (SPPAIR* s, FGNSRSPDIM index, FGNSRSPDIM numelem) {
    FGNSRSPDIM cur_node = index;
    FGNSRSPDIM child_node;
    while ((child_node = 2*cur_node+1) < numelem) {
        /* Then cur_node has at least one child */
        int max_index=cur_node;
        if (s[child_node].value > s[max_index].value) {
            max_index = child_node;
        }

        if (++child_node < numelem &&
            s[child_node].value > s[max_index].value)
        {
            /* Then cur_node has two children and we select the max child */
            max_index = child_node;
        }

        if (max_index == cur_node) {
            /* Element is no smaller than all children */
            break;
        }

        {
            SPPAIR tmp;
            tmp = s[max_index];
            s[max_index] = s[cur_node];
            s[cur_node] = tmp;
        }

        cur_node = max_index;
    }
}
