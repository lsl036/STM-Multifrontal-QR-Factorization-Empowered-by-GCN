#-------------------------------------------------------------------------------
# UMFPACK Makefile for compiling on Unix systems
#-------------------------------------------------------------------------------

INSTALL = $(realpath $(CURDIR)/../..)
AR_TARGET = libsparselu.a

#-------------------------------------------------------------------------------

default: lu

lu:
	make install INSTALL=$(INSTALL)

install: $(AR_TARGET)
	@mkdir -p $(INSTALL)/lib
	rm *.o
	mv -f *.a $(INSTALL)/lib

uninstall: 
	rm $(INSTALL)/lib/$(AR_TARGET)

clean: 
	rm *.o
#-------------------------------------------------------------------------------
CC = gcc
CF = -O2 -fexceptions -fPIC -fopenmp
I = -I../../include/

C = $(CC) $(CF) $(I)

#-------------------------------------------------------------------------------
# include files:
#-------------------------------------------------------------------------------

INC = ../../include/SparseBase_config.h ../../include/amd_internal.h ../../include/amd.h \
		../../include/SparseLU.h ../../include/SparseLU_config.h \
		../../include/SparseLU_internal.h ../../include/SparseLU_function.h
#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

UMF = LU_analyze LU_assemble LU_store_lu LU_base \
	LU_triplet_map_nox LU_triplet_nomap_x LU_triplet_nomap_nox \
	LU_triplet_map_x LU_assemble_fixq LU_store_lu_drop LU_solve \
	LU_ltsolve LU_utsolve LU_lhsolve LU_uhsolve

UMFUSER = SparseLU_base SparseLU_symbolic SparseLU_numeric SparseLU_solve SparseLU_kernel

#-------------------------------------------------------------------------------

OBJ = $(addsuffix .o, $(UMF) $(UMFUSER))
$(OBJ): $(INC)
$(AR_TARGET): $(OBJ)
	ar rv $@ $^
	- ranlib $@

#-------------------------------------------------------------------------------

LU_%.o: LU_%.c 
	$(C) -DDLONG -c $< -o $@

LU_%hsolve.o: LU_%tsolve.c $(INC)
	$(C) -DDLONG -DCONJUGATE_SOLVE -c $< -o $@

LU_triplet_map_x.o: LU_triplet.c
	$(C) -DDLONG -DDO_MAP -DDO_VALUES -c $< -o $@

LU_triplet_map_nox.o: LU_triplet.c 
	$(C) -DDLONG -DDO_MAP -c $< -o $@

LU_triplet_nomap_x.o: LU_triplet.c 
	$(C) -DDLONG -DDO_VALUES -c $< -o $@

LU_triplet_nomap_nox.o: LU_triplet.c 
	$(C) -DDLONG -c $< -o $@

LU_assemble_fixq.o: LU_assemble.c 
	$(C) -DDLONG -DFIXQ -c $< -o $@

LU_store_lu_drop.o: LU_store_lu.c 
	$(C) -DDLONG -DDROP -c $< -o $@

SparseLU_%.o: SparseLU_%.c 
	$(C) -DDLONG -c $< -o $@

SparseLU_wsolve.o: SparseLU_solve.c $(INC)
	$(C) -DDLONG -DWSOLVE -c $< -o $@


