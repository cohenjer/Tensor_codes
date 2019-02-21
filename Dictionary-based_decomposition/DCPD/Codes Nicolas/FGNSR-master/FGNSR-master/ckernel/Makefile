ARCH := $(shell uname -s)
PROC := $(shell uname -m)
HOST := $(shell hostname -s)

ifeq ($(MOD_THREAD),)
# Sequential build or parallel build? 'seq' / 'par'
	# possible values seq par
	MOD_THREAD = seq
endif

ifeq ($(MOD_OPT),)
	# possible values: dbg, opt, prf, cov
	MOD_OPT = dbg
endif

ifeq ($(MOD_CC),)
	# possible values: gnu, intel
	MOD_CC = gnu
endif

ifeq ($(API_SPARSE),)
# Controls the type of API for sparse matrices.  Possible values are:
# II (int/int), IL (int/long), LL (long/long), SS (size_t/size_t)
    API_SPARSE = SS
endif

EXTSTR := $(ARCH).$(PROC).$(MOD_CC).$(MOD_THREAD).$(MOD_OPT)

# Default settings
CC = gcc
LD = gcc
CFLAGS = -g -O0 -Wall -pedantic -I. -fPIC
LDFLAGS = -shared
MEX = mex
MEXOPTS = -g -largeArrayDims
CUTESTPATH = ../devtools/cutest-1.5

SHAREDLIBEXT=so
MEXEXT=mexa64

# End of default settings

# Use platform specific setup
include make/make.$(EXTSTR)

# Use host specific setup if setup file exists
ifeq ($(wildcard make/make.host_$(HOST)),make/make.host_$(HOST)) 
	include make/make.host_$(HOST)
endif

################### END OF BUILD SETUP ##############################

OBJFILES = proj.o util.o debug_tools.o sputil.o
TESTBIN = runtests.$(EXTSTR).bin
TESTSRC = runtests.c testing_tools.c test_proj.c test_util.c \
test_debug_tools.c test_sputil.c test_testing_tools.c
MEXFILE = parproj_priv.$(MEXEXT)
SHAREDLIB = libfgnsr.$(EXTSTR).$(SHAREDLIBEXT)
STATICLIB = libfgnsr.$(EXTSTR).a

# TODO Preferred sparse API, shuold be some parameter somewhere
CFLAGS += -DFGNSRAPI_SPARSEMAT_MODEL_$(API_SPARSE)

all: $(STATICLIB) $(SHAREDLIB) $(TESTBIN)

matlab: $(MEXFILE)

test: $(TESTBIN)
	./$(TESTBIN)

links: $(SHAREDLIB)
	ln -s $(SHAREDLIB) libfgnsr.so

vtest: $(TESTBIN)
	valgrind --suppressions=vg_supp ./$(TESTBIN)

$(STATICLIB): $(OBJFILES)
	ar rsc $(STATICLIB) $(OBJFILES)
	ranlib $(STATICLIB)

$(SHAREDLIB): $(OBJFILES)
	$(LD) $(LDFLAGS) -o $(SHAREDLIB) $(OBJFILES) $(LDLIBS)

$(MEXFILE): parproj.m parproj_priv.c $(SHAREDLIB)
	$(MEX) $(MEXOPTS) parproj_priv.c $(SHAREDLIB) $(LDLIBS)

$(TESTBIN): $(TESTSRC) $(STATICLIB) testing_tools.h
	$(CC) $(CFLAGS) -I$(CUTESTPATH) -o $(TESTBIN) $(TESTSRC) \
$(CUTESTPATH)/CuTest.c $(STATICLIB) $(LDLIBS)

a.out: scratch.c
	$(CC) -fopenmp -Wall -pedantic -g scratch.c

perfproj: perfproj.c $(STATICLIB)
	$(CC) $(CFLAGS) -o perfproj perfproj.c $(STATICLIB)

clean:
	rm -f $(OBJFILES)
	rm -rf *.dSYM
	rm -rf a.out
	rm -f perfproj
	rm -f *.gcov
	rm -f *.gcda
	rm -f *.gcno
	rm -f gmon.out
	rm -f *.pyc

purge: clean
	rm -f runtests.*.bin
	rm -f libfgnsr.*.$(SHAREDLIBEXT)
	rm -f libfgnsr.*.a
	rm -f parproj_priv.o
	rm -f $(MEXFILE)
