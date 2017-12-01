SHELL=/bin/bash

# Matrix dimensions.
N=1000
M=1000

CC=mpicc

# Compiler optimization level.
OPT_LEVEL=-O3

##########################################
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
##########################################

OPENMP_FLAG=-fopenmp
# Intel wants a different flag.
CC_VERSION=$(shell $(CC) --version)
ifeq ($(findstring Intel,$(CC_VERSION)),Intel)
  OPENMP_FLAG=-qopenmp
endif

# Used for checking the results.
ERROR_THRESHOLD=1e-4

# The papi library location.
PAPI_LIB_DIR=/usr/local/apps/papi/5.4.1/lib
PAPI_INC_DIR=/usr/local/apps/papi/5.4.1/include

# Common program arguments.
COMMON_PROG_ARGS=-DN=$(N) -DM=$(M) -std=c99 $(OPT_LEVEL) -DERROR_THRESHOLD=$(ERROR_THRESHOLD) $(OPENMP_FLAG)

# Program arguments for cc.
cc_PROG_ARGS=$(COMMON_PROG_ARGS)

# Compilation command for cc, no PAPI.
cc_NOPAPI_CC=$(CC) -DNOPAPI $(cc_PROG_ARGS)

cc_PAPI_CC=$(CC) -I${PAPI_INC_DIR} -L$(PAPI_LIB_DIR) -lpapi $(cc_PROG_ARGS)

all: test_omp test_mpi

test_omp: *.c
	$(cc_PAPI_CC) *.c -o test_omp_$(N)_$(M) -DUSE_OMP

test_mpi: *.c
	$(cc_PAPI_CC) *.c -o test_mpi_$(N)_$(M) -DUSE_MPI

clean:
	rm -f *.out *.o *.optrpt test_omp_* test_mpi_*
