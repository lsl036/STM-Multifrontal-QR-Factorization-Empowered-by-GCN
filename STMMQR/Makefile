 #===============================================================================
# demo/Makefile
#===============================================================================
include ./Makefile.option
#-------------------------------------------------------------------------------
INSTALL = $(realpath $(CURDIR))

INCHdir = ./include

INCLUDE = -I$(INCHdir) -I$(INCHdir)/tpsm

#-------------------------------------------------------------------------------

default: all
#test_qr 
all: library qrtest

library: 
	( cd src/base ; make )
	( cd src/core ; make )
	( cd src/chol ; make )
	( cd src/qr ; make )
	( cd src/qr ; make qrnumalib)

#-------------------------------------------------------------------------------

distclean: clean 
	- cd lib ; rm -r -f libsparsebase.a libsparsecore.a libsparsechol.a libsparsechol_pool.a libsparselu.a libsparseqr.a libsparseqr_numa.a

clean:
	- rm -r -f chol choltest lutest qrtest lqtest qrtest_numa_data

clean_result:
	- cd Results ; rm -r -f *.txt

optset:
	export OMP_NUM_THREADS=32
	export TPSM_ENV_ON=1
	export TPSM_BASIC_SIZE=2000
	ulimit -s unlimited
#-------------------------------------------------------------------------------

#--------------------#
#  使用原来的openblas #
#--------------------#
choltest: ./test/choltest.c
	$(CC) $(CF) $(INCLUDE) -o $@ $< -DPTHREAD_POOL $(LIB)

qrtest: ./test/qrtest.c
	$(CC) $(CF) $(INCLUDE) -o $@ $<  $(LIB)

qrtest_write: ./test/qrtest.c
	$(CC) $(CF) $(INCLUDE) -o $@ $< $(LIB) -Dwrite_graph

qrtest_numa_data: ./test/qrtest.c
	$(CC) $(CF) $(NUMA_FLAGS) $(INCLUDE) -o $@ $<  $(LIB_NUMA)

lqtest: ./test/lqtest.c
	$(CC) $(CF) $(INCLUDE) -o $@ $<  $(LIB)

amdtest: ./test/AMDtest.c
	$(CC) $(CF) $(INCLUDE) -o $@ $<  $(LIB)

saveorder: ./test/save_reorder.c
	$(CC) $(CF) $(INCLUDE) -o $@ $<  $(LIB)

chol_saveorder: ./test/chol_save_order.c
	$(CC) $(CF) $(INCLUDE) -o $@ $<  $(LIB)

update: 
	(cd src/qr; make)
	- (rm qrtest)
	(make qrtest)

DATA_NAME= se
test_lq: 
	@echo $(DATA_NAME)
	./lqtest_default /home/public/DATA_TEST/$(DATA_NAME)/$(DATA_NAME).mtx
	./lqtest_amd /home/public/DATA_TEST/$(DATA_NAME)/$(DATA_NAME).mtx
	./lqtest_metis /home/public/DATA_TEST/$(DATA_NAME)/$(DATA_NAME).mtx

