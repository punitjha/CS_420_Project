SHELL=/bin/bash

#To use this makefile type make -f Makefile.win
# Matrix dimensions.

CC=CC

# Compiler optimization level.
OPT_LEVEL=-O3

##########################################
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
##########################################







cc_PAPI_CC=$(CC) -I $(GSL_INCLUDE_PATH) -L $(GSL_LIBRARY_PATH) $(cc_PROG_ARGS)
LIBS= -lpapi 

OBJ = LU_mpi_final.o 

all: MPI.exe

MPI.exe: $(OBJ)
	$(cc_PAPI_CC) -o MPI.exe $^ $(LIBS) -DUSE_PAPI

%.o: %.cpp
	$(cc_PAPI_CC) -c  $< $(LIBS) -DUSE_PAPI 

clean:
	rm -f *.out *.o *.optrpt *.exe

