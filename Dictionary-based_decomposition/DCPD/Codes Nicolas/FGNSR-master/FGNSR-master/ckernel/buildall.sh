#!/bin/sh

make purge
MOD_OPT=dbg MOD_THREAD=seq make test
make clean
MOD_OPT=opt MOD_THREAD=seq make test
make clean
MOD_OPT=cov MOD_THREAD=seq make test
make clean
OMP_NUM_THREADS=2 MOD_OPT=dbg MOD_THREAD=par make test
make clean
OMP_NUM_THREADS=2 MOD_OPT=opt MOD_THREAD=par make test
make clean
OMP_NUM_THREADS=2 MOD_OPT=cov MOD_THREAD=par make test
make clean
