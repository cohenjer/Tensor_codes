#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "proj.h"

void fill_data(int numvec, int lenvec, double * matrix) {
    int k;
    const double shift = -.5;
    const double scale =  3.;

    for (k=0; k<numvec*lenvec; k++) {
        matrix[k] = (scale * drand48()) + shift;
    }
}

int main (int argc, char *argv[]) {
    int status = 0;
    int lenvec = 12000;
    int numvec = 12000;
    int data_size;

    PROJDATAptr pdata = NULL;
    int verbose = 1;
    double ub = 10.;

    double *weights = NULL;
    double *matrix = NULL;
    int *diag_idx = NULL;
    
    int k;
    unsigned int seed;

    if (argc < 2) {
        seed = 2341563;
    }
    else {
        seed = atoi(argv[1]);
    }

    srand48(seed);

    data_size = lenvec * numvec;

    weights = (double *) malloc(lenvec * sizeof(double));
    matrix = (double *) malloc(data_size * sizeof(double));
    diag_idx = (int *) malloc(numvec * sizeof(int));

    if (weights == NULL || matrix==NULL || diag_idx==NULL) {
        status = 1;
        goto TERMINATE;
    }

    for (k=0; k<lenvec; k++) {
        weights[k] = 1.0;
    }

    for (k=0; k<numvec; k++) {
        diag_idx[k] = k;
    }


    fill_data(numvec, lenvec, matrix);

    assert(lenvec >= 10);
    printf("The first values of the first column:\n");
    for (k=0; k<10; k++) {
        printf("%+.2f ", matrix[k]);
    }
    printf("\n");

    status = FGNSRproj_init(lenvec, weights, verbose, &pdata);
    if (status != 0)
        goto TERMINATE;

    status = FGNSRproj_project_dmatrix(pdata, numvec, matrix, ub, diag_idx);
    if (status != 0)
        goto TERMINATE;

    FGNSRproj_free(&pdata);

    for (k=0; k<10; k++) {
        printf("%+.2f ", matrix[k]);
    }
    printf("\n");

    free(weights);
    free(matrix);
    free(diag_idx);

TERMINATE:
    if (status != 0) {
        printf("Out of memory");
    }

    return status;
}
