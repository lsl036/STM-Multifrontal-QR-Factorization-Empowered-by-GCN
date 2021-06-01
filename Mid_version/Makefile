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
LIB =  $(LDFLAGS) -lsparselu -lsparseqr -lsparsechol_pool -lsparsecore -lsparsebase -llapack -lopenblas $(LDLIBS)

INCHdir = ./include

I = -I$(INCHdir)

CC = gcc -g

CF = -O2 -fexceptions -fPIC -fopenmp
#-------------------------------------------------------------------------------

default: all

all: library qrtest

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
	- rm -r -f qrtest *.mtx *.txt
#-------------------------------------------------------------------------------

#--------------------#
#  openblas #
#--------------------#

qrtest: ./test/qrtest.c
	$(CC) $(CF) $(I) -o $@ $< -DCALTOL -DMULTI $(LIB) -lnuma


update: 
	(cd src/qr; make)
	(rm qrpool)
	(make qrpool)

#--------------------#
#      测试用例       #
#--------------------#

qrpool_test: qrpool
	# - ./qrpool /home/public/DATA_TEST/ inlet_C
	# - ./qrpool /home/public/DATA_TEST/ LF10_B
	# - ./qrpool /home/public/DATA_TEST/ LF10_C



