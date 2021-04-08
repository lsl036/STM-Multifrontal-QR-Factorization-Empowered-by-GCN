#===============================================================================
# hnuchol/Makefile: for compiling the hnuchol library
#===============================================================================
include ../../Makefile.option
INSTALL = $(realpath $(CURDIR)/../..)
AR_TARGET = libsparsechol.a
AR_TARGET_POOL = libsparsechol_pool.a

#-------------------------------------------------------------------------------

default: chol_pool

chol:
	make install INSTALL=$(INSTALL)

chol_pool:
	make install_pool INSTALL=$(INSTALL)

install: $(AR_TARGET)
	@mkdir -p $(INSTALL)/lib
	@mkdir -p $(INSTALL)/include
	rm *.o
	mv -f *.a $(INSTALL)/lib

install_pool: $(AR_TARGET_POOL)
	@mkdir -p $(INSTALL)/lib
	@mkdir -p $(INSTALL)/include
	rm *.o
	mv -f *.a $(INSTALL)/lib

uninstall: 
	rm $(INSTALL)/lib/$(AR_TARGET)

uninstall_pool: 
	rm $(INSTALL)/lib/$(AR_TARGET_POOL)
	
#-------------------------------------------------------------------------------
#CC = gcc -g
# CF = -O2 -fexceptions -fPIC -fopenmp
# I = -I../../include -I../../include/tpsm

# C = $(CC) $(CF) $(I)

INC =  ../../include/SparseChol.h

#-------------------------------------------------------------------------------
# hnuchol library module (long, double)
#-------------------------------------------------------------------------------

LCHOLESKY = SparseChol_analyze.o SparseChol_factorize.o \
	SparseChol_super_numeric.o SparseChol_super_symbolic.o SparseChol_solve.o SparseChol_super_solve.o

LCHOLESKY_POOL = SparseChol_analyze.o SparseChol_factorize.o \
	SparseChol_super_numeric_pool.o SparseChol_super_symbolic.o SparseChol_solve.o SparseChol_super_solve.o

DL = $(LCHOLESKY)

DL_POOL = $(LCHOLESKY_POOL)

#-------------------------------------------------------------------------------

OBJ = $(DL)
OBJ_POOL = $(DL_POOL)

$(OBJ): $(INC)
$(OBJ_POOL): $(INC)
$(OBJ_ATL): $(INC)

$(AR_TARGET): $(OBJ)
	ar rv $@ $^
	- ranlib $@

$(AR_TARGET_POOL): $(OBJ_POOL)
	ar rv $@ $^
	- ranlib $@

#-------------------------------------------------------------------------------
# Cholesky Module:
#-------------------------------------------------------------------------------

SparseChol_solve.o: SparseChol_solve.c SparseChol_t_solve.c SparseChol_t_lsolve.c SparseChol_t_ltsolve.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseChol_super_solve.o: SparseChol_super_solve.c SparseChol_t_super_solve.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseChol_analyze.o: SparseChol_analyze.c SparseChol_t_rowfac.c
	$(C) $(LONG_FLAGS) $(NMETIS_FLAG) $(ANA_METHOD) -c $< -o $@ 

SparseChol_factorize.o: SparseChol_factorize.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseChol_super_symbolic.o: SparseChol_super_symbolic.c
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseChol_super_numeric.o: SparseChol_super_numeric.c 
	$(C) $(LONG_FLAGS) -c $< -o $@

SparseChol_super_numeric_pool.o: SparseChol_super_numeric.c 
	$(C) $(LONG_FLAGS) -DPTHREAD_POOL -c $< -o $@

SparseChol_super_numeric_atl.o: SparseChol_super_numeric.c 
	$(C) $(LONG_FLAGS) -DUSEATL -c $< -o $@