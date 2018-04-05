/* parproj_priv.c  --  Matlab connector for parallel fast projection */

#include <assert.h>

#include "mex.h"
#include "fgnsr_ss.h"

/* matlab prototype
function X = parproj(Z, w, ub, diag_idx, verbose)
*/

static void do_project_dense(int nrows, int ncols, double *weights,
        double ub, int verbose, int *diag_idx, double *dmatrix);

static void do_project_sparse(int nrows, int ncols, double *weights,
        double ub, int verbose, FGNSRSPDIM *diag_idx, FGNSRSPNZ *begin,
        FGNSRSPDIM *index, double *value);

/* Gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    const mxArray *Z = NULL;        /* Input matrix to be projected */
    const mxArray *w = NULL;        /* Weights, usually some column norms */
    double ub;                      /* Upper bound input argument */
    const mxArray *diag_idx = NULL; /* Positions of the diagonal elements*/
    int verbose;                    /* Verbosity switch */

    mxArray *output = NULL;         /* The projected iterate */
    double *matrix;                 /* Pointer to the output matrix data */

    mwSize nrows, ncols;            /* Size of Z */
    mwSize tmp_m, tmp_n;
    int issparse;

    /************************************************************************/
    /* Setup and input checks ***********************************************/
    /************************************************************************/

    /* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("FGNSR:input", "Five inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("FGNSR:input", "One output required.");
    }

    /* First argument: Z matrix, iterate to be projected */
    Z = prhs[0];
    if( !mxIsDouble(Z) || mxIsComplex(Z) ) {
        mexErrMsgIdAndTxt("FGNSR:input",
                "Input matrix Z must be real double.");
    }

    /* NOTE:  Matlab stores 2d arrays flattened col major.  Thus if the
     * original Z matrix holds the to-be-projected points in its _rows_, the
     * wrapper parproj.m must pass the transpose of Z to this function. */
    nrows = mxGetM(Z);
    ncols = mxGetN(Z);
    issparse = mxIsSparse(Z);
    
    /* Second argument: weight vector (usually some column norms) */
    w = prhs[1];
    if( !mxIsDouble(w) || mxIsComplex(w) ){
        mexErrMsgIdAndTxt("FGNSR:input", "Weight vector must be real double");
    }

    /* Check for consistent length */
    tmp_m = mxGetM(w);
    tmp_n = mxGetN(w);
    if ( !(tmp_m==1 && tmp_n==nrows) && !(tmp_m==nrows && tmp_n==1) ) {
        mexErrMsgIdAndTxt("FGNSR:input",
                "Weight vector must be 1d of length %d", nrows);
    }

    /* Third argument */
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2])!=1 ) {
        mexErrMsgIdAndTxt("FGNSR:input","ub must be a real double scalar");
    }
    ub = mxGetScalar(prhs[2]);

    /* Fourth argument
     * NOTE: We assume that the indices are already base-0 shifted */
    diag_idx = prhs[3];
    if( !issparse && !mxIsInt32(diag_idx) ) {
        mexErrMsgIdAndTxt("FGNSR:input", "diag_idx vector must be int32");
    }

    if( issparse && !mxIsInt64(diag_idx) ) {
        mexErrMsgIdAndTxt("FGNSR:input", "diag_idx vector must be int64");
    }

    /* Check for consistent length */
    tmp_m = mxGetM(diag_idx);
    tmp_n = mxGetN(diag_idx);
    if ( !(tmp_m==1 && tmp_n==ncols) && !(tmp_m==ncols && tmp_n==1) ) {
        mexErrMsgIdAndTxt("FGNSR:input",
                "diag_idx vector must be 1d of length %d", nrows);
    }

    /* Fifth argument */
    if( !mxIsInt32(prhs[4]) || mxGetNumberOfElements(prhs[4])!=1 ) {
        mexErrMsgIdAndTxt("FGNSR:input","verbose must be a int scalar");
    }
    verbose = mxGetScalar(prhs[4]);

    /* create the output matrix */
    output = mxDuplicateArray(Z);
    /* If this call has failed, control is passed ot Matlab, hence no check */

    /************************************************************************/
    /* Dense projection case ************************************************/
    /************************************************************************/

    if (issparse) {
        mxAssert(sizeof(mwIndex) == sizeof(FGNSRSPNZ),  "Sparse API mismatach");
        mxAssert(sizeof(mwIndex) == sizeof(FGNSRSPDIM), "Sparse API mismatach");

        do_project_sparse(nrows, ncols, mxGetPr(w), ub, verbose,
                (FGNSRSPDIM *) mxGetData(diag_idx), mxGetJc(output), 
                mxGetIr(output), mxGetPr(output));

    }
    else {
        do_project_dense(nrows, ncols, mxGetPr(w), ub, verbose,
                (int *) mxGetData(diag_idx), mxGetPr(output));
    }
    
    plhs[0] = output;
}


static void do_project_dense(int nrows, int ncols, double *weights,
        double ub, int verbose, int *diag_idx, double *dmatrix)
{
    PROJDATAptr pdata = NULL;       /* Data structure for projection */
    int status = 0;

    /* CAUTION In this region:  make sure that mexErr* induced long jumps do
     * not cause memory leaks from the allocation of pdata */

    /* Allocate memory */
    status = FGNSRproj_init(nrows, weights, verbose, &pdata);
    if (status != 0) {
        mexErrMsgIdAndTxt("FGNSR:general", "Error %d\n", status);
        /* next statement is never reached */
        assert(0);
    }

    /* Do projection */
    status = FGNSRproj_project_dmatrix(pdata, ncols, dmatrix, ub, diag_idx);
    if (status != 0) {
        FGNSRproj_free(&pdata);
        mexErrMsgIdAndTxt("FGNSR:general", "Error %d\n", status);
        /* next statement is never reached */
        assert(0);
    }

    FGNSRproj_free(&pdata);
    /* *********** END CAUTION **************** */
}

static void do_project_sparse(int nrows, int ncols, double *weights,
        double ub, int verbose, FGNSRSPDIM *diag_idx, FGNSRSPNZ *begin,
        FGNSRSPDIM *index, double *value)
{
    SPPROJDATAptr pdata = NULL;       /* Data structure for projection */
    int status = 0;

    /* CAUTION In this region:  make sure that mexErr* induced long jumps do
     * not cause memory leaks from the allocation of pdata */

    /* Allocate memory */
    status = FGNSRproj_spinit(nrows, weights, verbose, &pdata);
    if (status != 0) {
        mexErrMsgIdAndTxt("FGNSR:general", "Error %d\n", status);
        /* next statement is never reached */
        assert(0);
    }

    /* Do projection */
    status = FGNSRproj_project_spmatrix(pdata, ncols, begin, index, value,
            ub, diag_idx);
    if (status != 0) {
        FGNSRproj_spfree(&pdata);
        mexErrMsgIdAndTxt("FGNSR:general", "Error %d\n", status);
        /* next statement is never reached */
        assert(0);
    }

    FGNSRproj_spfree(&pdata);
    /* *********** END CAUTION **************** */
}



