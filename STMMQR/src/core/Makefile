include ../../Makefile.option
INSTALL = $(realpath $(CURDIR)/../..)
AR_TARGET = libsparsecore.a

#-------------------------------------------------------------------------------

default: core

core:
	make install INSTALL=$(INSTALL)

install: $(AR_TARGET)
	@mkdir -p $(INSTALL)/lib
	@mkdir -p $(INSTALL)/include
	rm *.o
	mv -f *.a $(INSTALL)/lib

uninstall: 
	rm $(INSTALL)/lib/$(AR_TARGET)

#-------------------------------------------------------------------------------

#CC = gcc
# CF = -O2 -fexceptions -fPIC -fopenmp
# I = -I../../include

# C = $(CC) $(CF) $(I)

#-------------------------------------------------------------------------------
# Include/ directory contains all include files:
#-------------------------------------------------------------------------------

INC =   ../../include/Sparse.h \
	../../include/SparseCore.h \
	../../include/Sparse_internal.h \
	../../include/Sparse_template.h

DL = SparseCore_matrixops.o SparseCore_change_factor.o SparseCore_common.o \
	SparseCore_matrix_type.o SparseCore_check.o SparseCore_read_write.o SparseCore_norm.o \
	SparseCore_metis.o SparseCore_nesdis.o SparseCore_csymamd.o SparseCore_camd.o SparseCore_ccolamd.o

#-------------------------------------------------------------------------------

OBJ = $(DL)

$(OBJ): $(INC)

$(AR_TARGET): $(OBJ)
	ar rv $@ $^
	- ranlib $@

#-------------------------------------------------------------------------------
# Core Module:
#-------------------------------------------------------------------------------

SparseCore_common.o: SparseCore_common.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseCore_change_factor.o: SparseCore_change_factor.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseCore_matrix_type.o: SparseCore_matrix_type.c SparseCore_t_triplet.c SparseCore_t_dense.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseCore_matrixops.o: SparseCore_matrixops.c SparseCore_t_transpose.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseCore_check.o: SparseCore_check.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseCore_read_write.o: SparseCore_read_write.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseCore_norm.o: SparseCore_norm.c SparseCore_t_sdmult.c
	$(C) $(LONG_FLAGS) -c $< -o $@

#-------------------------
# METIS module :
#-------------------------
SparseCore_metis.o: SparseCore_metis.c
	$(C) $(LONG_FLAGS) $(NMETIS_FLAG) -c $< -o $@

# ---------------------
#	NESDIS MODEULE
# ---------------------
SparseCore_nesdis.o: SparseCore_nesdis.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseCore_csymamd.o: SparseCore_csymamd.c
	$(C) $(LONG_FLAGS) $(INCCCOL) $(INCSuiteBase) -c $< -o $@

SparseCore_camd.o: SparseCore_camd.c
	$(C) $(LONG_FLAGS) $(INCCAMD) $(INCSuiteBase) -c $< -o $@

SparseCore_ccolamd.o: SparseCore_ccolamd.c
	$(C) $(LONG_FLAGS) $(INCCCOL) $(INCSuiteBase) -c $< -o $@

