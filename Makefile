SHELL=/bin/bash

#To use this makefile type make -f Makefile.win
# Matrix dimensions.

CC=h5c++
# CC=CC

# Compiler optimization level.
OPT_LEVEL=-O3

##########################################
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
##########################################

#OPENMP_FLAG=-fopenmp
# Intel wants a different flag.
#CC_VERSION=$(shell $(CC) --version)
#ifeq ($(findstring Intel,$(CC_VERSION)),Intel)
#  OPENMP_FLAG=-qopenmp
#endif

# Used for checking the results.
#ERROR_THRESHOLD=1e-4

# The papi library location.
#We have linked the PAPI library by linker 

# Common program arguments.
COMMON_PROG_ARGS= $(OPT_LEVEL) 

# Program arguments for cc.
cc_PROG_ARGS=$(COMMON_PROG_ARGS)

# Compilation command for cc, no PAPI.
cc_NOPAPI_CC=$(CC) -DNOPAPI $(cc_PROG_ARGS)

cc_PAPI_CC=$(CC) -I $(GSL_INCLUDE_PATH) -L $(GSL_LIBRARY_PATH) $(cc_PROG_ARGS)
LIBS=-lgsl -lpapi -lgslcblas -lm 
OBJ = main.o hdf_output.o matrix.o math.o

all: LU.exe

LU.exe: $(OBJ)
	$(cc_PAPI_CC) -o LU.exe $^ $(LIBS) 

%.o: %.cpp
	$(cc_PAPI_CC) -c  $< $(LIBS) 

clean:
	rm -f *.out *.o *.optrpt *.exe
