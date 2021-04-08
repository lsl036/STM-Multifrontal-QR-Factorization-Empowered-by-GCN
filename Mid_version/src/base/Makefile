#-------------------------------------------------------------------------------
# HNU_BASE Makefile
#-------------------------------------------------------------------------------

INSTALL = $(realpath $(CURDIR)/../..)
AR_TARGET = libsparsebase.a

#-------------------------------------------------------------------------------

default: base

base:
	make install INSTALL=$(INSTALL)

install: $(AR_TARGET)
	@mkdir -p $(INSTALL)/lib
	@mkdir -p $(INSTALL)/include
	rm *.o
	mv -f *.a $(INSTALL)/lib

uninstall: 
	rm $(INSTALL)/lib/$(AR_TARGET)

#-------------------------------------------------------------------------------

CC = gcc
CF = -O2 -fexceptions -fPIC -fopenmp
I = -I../../include
C = $(CC) $(CF) $(I)

INC = ../../include/SparseBase_config.h ../../include/xerbla.h ../../include/colamd.h \
	 ../../include/amd.h ../../include/amd_internal.h

OBJ = SparseBase_config.o xerbla.o amd.o colamd.o \
		threadbase.o threadpool_manager.o

$(OBJ): $(INC)

$(AR_TARGET): $(OBJ)
	ar rv $@ $^
	- ranlib $@

#-------------------------------------------------------------------------------

SparseBase_config.o: SparseBase_config.c
	$(C) -c $< -o $@

xerbla.o: xerbla.c
	$(C) -c $< -o $@

amd.o: amd.c
	$(C) -DDLONG -c $< -o $@

colamd.o: colamd.c
	$(C) -DDLONG -c $< -o $@

#  ------------------------------
#           线程池
#  ------------------------------

threadbase.o : ./tpsm_base.c ../../include/tpsm_base.h
	gcc -std=c99 -fPIC -O2 $(I) -c $< -o $@ 

threadpool_manager.o: ./tpsm.c ../../include/tpsm.h
	gcc -std=c99 -fPIC -O2 $(I) -c $< -o $@ -D_GNU_SOURCE 