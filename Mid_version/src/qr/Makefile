CF = -O2 -fPIC -fopenmp

INCHdir = ../../include
I = -I$(INCHdir)
CC = gcc -std=c99 $(CF) $(I)

ARCHDEFS = -DHNU_ARM_LINUX64

PRECISION = -DDOUBLE_PRECISION

default:qrpoollib

# qrlib: SparseQR.o SparseQR_analyze.o SparseQR_factorize.o 
# 	ar crv libhnuqr.a SparseQR.o SparseQR_analyze.o SparseQR_factorize.o 
# 	mv libhnuqr.a ../../lib/ 
# 	rm *.o

# qratllib: SparseQR_atl.o SparseQR_analyze.o SparseQR_factorize_atl.o 
# 	ar crv libSparseQR_atl.a SparseQR.o SparseQR_analyze.o SparseQR_factorize.o 
# 	mv libSparseQR_atl.a ../../lib/ 
# 	rm *.o

# 上层任务开启线程池的并行库
qrpoollib: SparseQR.o SparseQR_analyze.o SparseQR_factorize_pool.o SparseQR_multithreads.o SparseLQ.o
	ar crv libsparseqr.a SparseQR.o SparseQR_analyze.o SparseQR_factorize_pool.o SparseQR_multithreads.o SparseLQ.o
	mv libsparseqr.a ../../lib/ 
	rm *.o

#  ------------------------------
#        三个稀疏分解函数
#  ------------------------------

SparseQR.o: ./SparseQR.c
	$(CC) -c  $<

SparseQR_analyze.o: ./SparseQR_analyze.c
	$(CC) -c  $<

SparseQR_factorize.o: ./SparseQR_factorize.c
	$(CC) -c  $<

SparseQR_factorize_pool.o: ./SparseQR_factorize.c
	$(CC) -o $@ -c -DMULTI $<

SparseQR_multithreads.o: ./SparseQR_multithreads.c
	$(CC) -c  $<

SparseLQ.o: ./SparseLQ.c
	$(CC) -c  $<