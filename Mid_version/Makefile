#===============================================================================
# demo/Makefile
#===============================================================================

ARCHDEFS = -DHNU_ARM_LINUX64
PRECISION = -DDOUBLE_PRECISION
#-------------------------------------------------------------------------------
INSTALL = $(realpath $(CURDIR))
LDLIBS = -lm -lrt
LDFLAGS = -L$(INSTALL)/lib

# LIB_ORI =  $(LDFLAGS) -lopenblas -lsparselu -lsparseqr -lsparsechol -lsparsecore -lsparsebase $(LDLIBS)
LIB =  $(LDFLAGS) -lopenblas -lsparselu -lsparseqr -lsparsechol_pool -lsparsecore -lsparsebase $(LDLIBS)

INCHdir = ./include

I = -I$(INCHdir)

CC = gcc -g

CF = -O2 -fexceptions -fPIC -fopenmp
#-------------------------------------------------------------------------------

default: all

all: library choltest qrtest lqtest lutest

library: 
	( cd src/base ; make )
	( cd src/core ; make )
	( cd src/chol ; make )
	( cd src/qr ; make )
	( cd src/lu ; make )

test: chol_test chol_pool_test qr_test qr_pool_test 

#-------------------------------------------------------------------------------

distclean: clean
	- cd lib ; rm -r -f libsparsebase.a libsparsecore.a libsparsechol.a libsparsechol_pool.a libsparselu.a libsparseqr.a

clean:
	- rm -r -f chol choltest lutest qrtest lqtest *.mtx
#-------------------------------------------------------------------------------

#--------------------#
#  使用原来的openblas #
#--------------------#
choltest: ./test/choltest.c
	$(CC) $(CF) $(I) -o $@ $< -DPTHREAD_POOL $(LIB) -lnuma

lutest: ./test/lutest.c
	$(CC) $(CF) $(I) -o $@ $<  $(LIB) -lnuma

qrtest: ./test/qrtest.c
	$(CC) $(CF) $(I) -o $@ $< -DCALTOL -DMULTI $(LIB) -lnuma

lqtest: ./test/lqtest.c
	$(CC) $(CF) $(I) -o $@ $< -DCALTOL -DMULTI $(LIB) -lnuma

update: 
	(cd src/qr; make)
	(rm qrpool)
	(make qrpool)

#--------------------#
#      测试用例       #
#--------------------#

chol_test: chol
	# - ./chol /home/public/DATA_TEST/ pwtk
	# - ./chol /home/public/DATA_TEST/ 1138_bus
	# - ./chol /home/public/DATA_TEST/ 494_bus
	# - ./chol /home/public/DATA_TEST/ 662_bus
	# - ./chol /home/public/DATA_TEST/ 685_bus
	# - ./chol /home/public/DATA_TEST/ bcsstk01
	# - ./chol /home/public/DATA_TEST/ bcsstk02
	# - ./chol /home/public/DATA_TEST/ bcsstk03
	# - ./chol /home/public/DATA_TEST/ bcsstk04
	# - ./chol /home/public/DATA_TEST/ bcsstk05
	# - ./chol /home/public/DATA_TEST/ bcsstk06
	# - ./chol /home/public/DATA_TEST/ bcsstk07
	# - ./chol /home/public/DATA_TEST/ bcsstk08
	# - ./chol /home/public/DATA_TEST/ bcsstk09
	# - ./chol /home/public/DATA_TEST/ bcsstk10
	# - ./chol /home/public/DATA_TEST/ bcsstk11
	# - ./chol /home/public/DATA_TEST/ bcsstk12
	# - ./chol /home/public/DATA_TEST/ bcsstk13
	# - ./chol /home/public/DATA_TEST/ bcsstk14
	# - ./chol /home/public/DATA_TEST/ bcsstk15
	# - ./chol /home/public/DATA_TEST/ bcsstk16
	# - ./chol /home/public/DATA_TEST/ bcsstk17
	# - ./chol /home/public/DATA_TEST/ bcsstk18
	# - ./chol /home/public/DATA_TEST/ bcsstk19


chol_pool_test: chol_pool
	# - ./chol_pool /home/public/DATA_TEST/ pwtk
	# - ./chol_pool /home/public/DATA_TEST/ 1138_bus
	# - ./chol_pool /home/public/DATA_TEST/ 494_bus
	# - ./chol_pool /home/public/DATA_TEST/ 662_bus
	# - ./chol_pool /home/public/DATA_TEST/ 685_bus
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk01
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk02
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk03
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk04
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk05
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk06
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk07
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk08
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk09
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk10
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk11
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk12
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk13
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk14
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk15
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk16
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk17
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk18
	# - ./chol_pool /home/public/DATA_TEST/ bcsstk19


lu_test: lu
	# - ./lu /home/public/DATA_TEST/ tols4000
	# - ./lu /home/public/DATA_TEST/ jan99jac060sc
	# - ./lu /home/public/DATA_TEST/ nnc1374
	# - ./lu /home/public/DATA_TEST/ sx-mathoverflow
	# - ./lu /home/public/DATA_TEST/ psmigr_1
	# - ./lu /home/public/DATA_TEST/ g7jac080sc
	# - ./lu /home/public/DATA_TEST/ Goodwin_071

qrpool_test: qrpool
	# - ./qrpool /home/public/DATA_TEST/ inlet_C
	# - ./qrpool /home/public/DATA_TEST/ LF10_B
	# - ./qrpool /home/public/DATA_TEST/ LF10_C



