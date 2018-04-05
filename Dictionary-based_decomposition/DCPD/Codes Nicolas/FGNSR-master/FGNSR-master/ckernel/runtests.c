#include <stdio.h>
#include <stdlib.h>

#include "CuTest.h"

#define NUMSUITS 5

CuSuite* GetUtilSuite();
CuSuite* GetSpUtilSuite();
CuSuite* GetProjSuite();
CuSuite* GetDebugToolsSuite();
CuSuite* GetTestingToolsSuite();

void RunAllTests(int mask[NUMSUITS]) {
    CuString *output = CuStringNew();
    int k;

    char *suitename[NUMSUITS] = {
        " 1 util.c --------------",
        " 2 sputil.c ------------",
        " 3 proj.c --------------",
        " 4 debug_tools.c -------",
        " 5 testing_tools.c -----",
    };

    CuSuite *suite[NUMSUITS];
    suite[0] = GetUtilSuite();
    suite[1] = GetSpUtilSuite();
    suite[2] = GetProjSuite();
    suite[3] = GetDebugToolsSuite();
    suite[4] = GetTestingToolsSuite();

    printf("The following %d suits are registered:\n", NUMSUITS);
    for (k=0; k<NUMSUITS; k++) {
        printf("%s\n", suitename[k]);
    }

    for (k=0; k<NUMSUITS; k++) {
        if (mask[k]) {
            /* Skip this suite */
            printf("%s\n    [SKIPPED]\n", suitename[k]);
            continue;
        }
        CuStringAppendFormat(output, "%s\n", suitename[k]);
        CuSuiteRun(suite[k]);
        CuSuiteSummary(suite[k], output);
        CuSuiteDetails(suite[k], output);
        CuSuiteDelete(suite[k]);
    }

    printf("%s\n", output->buffer);
    CuStringDelete(output);
}

int main(int argc, char *argv[]) {
    int mask[NUMSUITS] = {0,0,0,0};
    int k;
    int value;

    for (k=1; k<argc; k++) {
        value = atoi(argv[k]);
        if (value < 1 || value > NUMSUITS) {
            printf("No suite number %d, ignoring\n", value);
        }
        else {
            mask[--value] = 1;
        }
    }


    RunAllTests(mask);
    return 0;
}
