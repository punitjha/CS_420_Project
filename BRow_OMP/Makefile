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

OBJ = RBlock_OMP.o 

all: ROMP_Block.exe

ROMP_Block.exe: $(OBJ)
	$(cc_PAPI_CC) -o ROMP_Block.exe $^ $(LIBS) -DUSE_PAPI

%.o: %.cpp
	$(cc_PAPI_CC) -c  $< $(LIBS) -DUSE_PAPI 

clean:
	rm -f *.out *.o *.optrpt *.exe

