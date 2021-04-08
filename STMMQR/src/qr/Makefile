# CF = -O2 -fPIC -fopenmp -DPRINT_TIME

# INCHdir = ../../include
# I = -I$(INCHdir) -I$(INCHdir)/tpsm
#CC = gcc -std=c99 $(CF) $(I) 
include ../../Makefile.option

 default:qrpoollib
#default:qrnumalib 
# 上层任务开启线程池的并行库
qrpoollib: SparseQR.o SparseQR_analyze.o SparseQR_factorize.o SparseQR_multithreads.o SparseLQ.o
	ar crv libsparseqr.a $^
	mv libsparseqr.a ../../lib/ 
	rm *.o

qrnumalib: SparseQR.o SparseQR_analyze.o SparseQR_factorize_numa.o SparseQR_multithreads_numa.o SparseLQ.o
	ar crv libsparseqr_numa.a $^
	mv libsparseqr_numa.a ../../lib/ 
	rm *.o

#  ------------------------------
#        三个稀疏分解函数
#  ------------------------------

SparseQR.o: ./SparseQR.c
	$(C) $(TIME_FLAGS) $(NMETIS_FLAG) $(ANA_METHOD) -c  $< -o $@  

SparseQR_analyze.o: ./SparseQR_analyze.c
	$(C) $(TIME_FLAGS) -c  $< -o $@

SparseQR_factorize.o: ./SparseQR_factorize.c
	$(C) $(TIME_FLAGS) -c  $< -o $@  

SparseQR_factorize_numa.o: ./SparseQR_factorize.c
	$(C) $(TIME_FLAGS) $(NUMA_FLAGS) -c  $< -o $@  

SparseQR_multithreads.o: ./SparseQR_multithreads.c
	$(C) $(TIME_FLAGS) -c  $< -o $@  

SparseQR_multithreads_numa.o: ./SparseQR_multithreads.c
	$(C) $(TIME_FLAGS) $(NUMA_FLAGS) -c  $< -o $@  

SparseLQ.o: ./SparseLQ.c
	$(C) $(TIME_FLAGS) -c  $< -o $@