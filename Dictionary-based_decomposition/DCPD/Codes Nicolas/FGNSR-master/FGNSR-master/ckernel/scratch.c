#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

static void fun() {
    int tid;

    tid = omp_get_thread_num();

    printf("Hi from thread %d\n", tid);
}


int main (int argc, char *argv[]) {

#pragma omp parallel
    {
        fun();
    } /* All threads join master thread, barrier and disband */

    return 0;
}

