
getconf GNU_LIBPTHREAD_VERSION
NPTL 2.17

# STM-Multifrontal QR Factorization Code
Our optimization measures are mainly concentrated in the numerical factorization phase.
Typing `make` can generate 5 static lib in *./lib*, and an executable file in root directory named *qrtest*.
Two additional parameters are required to execute the executable file, which are the matrix data set and matrix-ID(Free to fill in). Examples are as follows:
`./qrtest /data/sme3Dc/sme3Dc.mtx 0`
The source code of our streaming task mapping framework is in *./src/base/tpsm_**.
The functions implemented in each subsequent section need to recompile the entire package, type:
```
make distclean
make
```
and then *qrtest* will achieve the result you want.
The library version under the default *Makefile.option* is the STM-Multifrontal QR factorization with optimal performance (qrtest) that we compared with MKL sparse QR function.

## Change the algorithm of fill-in reducing reorder
Enter the *Makefile.option* file and modify the **ANA_METHOD** option to the reordering method you need.
The **DEFAULT** method is used in the code, which is COLAMD algorithm.
When you want to replace the reordering method, after modifying the *Makefile.option*, you need to recompile the entire library in the root directory.
At this time, *QR_time.txt* will be generated in the *./Results* directory to record the matrix-ID, analysis factorization time, numerical factorization time, and error verification.
We recommend to rename and save the QR factorization of a method after compiling it, such as:`mv qrtest qrtest_default`.
Then *qrtest* will run under the reordering method you selected.

## Brute-force method to find the fewest fill-in elements
Enter the *Makefile.option* file and modify the **CF** option, uncomment code : `CF += -Dall_methods_time`, then recompile the entire library in the root directory.
In this method, *Brute-force-fill.txt* will be generated in the *./Results* directory to record the number of filling elements of AMD, COLAMD, METIS and NESDIS method.
*Brute_force_time.txt* will be generated in the *./Results* directory to record the time cost of the brute-force method.
Then *qrtest* will run under the brute-force method to select reordering method.

## Generate features of GCN
Enter the *Makefile.option* file and modify the **CF** option, uncomment code : `CF += -Dwrite_graph`, then recompile the entire library in the root directory.
In this case, the test function will not perform factorization, but will only calculate the attributes we need during the matrix reading process.
*QR_Edge.txt*, *QR_Node.txt* and *QR_exinfo.txt* will be generated in this mode.
*QR_Edge.txt* records graph-ID, row_index, column_index and values.
*QR_Node.txt* records graph-ID, row_index, in-degree, out-degree and external degree.
*QR_exinfo.txt* recodes the properties of entire graph, seen in our paper.
Then *qrtest* will run just for generate GCN input dataset.

## Performance improvement of optimization measures
Change the `CORE` and `cc->SPQR_grain` in *./test/qrtest.c* to adjust the tree-level parallelism.
When `cc->SPQR_grain = 1`, the QR factorization is the single-thread mode in tree-level (still can parallel in BLAS).

We provide an additional code package HNUSparse (Mid-term version of a project), which uses the streaming task scheduling framework and does not implement thread affinity (due to different interfaces and limited time, we did not integrate it into one package).
You can just type `make` at the root directory, then you will get an executable file named `qrtest`.
The specifications of the test function are as follows (ID number is not required):
`./qrtest /data/sme3Dc/sme3Dc.mtx`, which is corresponding to the *Mid* shown in figure.

The data set we used and our test results are summarized in the table *STM-MQR.xlsx* in the root directory.
The *GCNdata_408.txt* save the 408 data sets we selected in the GCN classifier experiment.



## Compared the assembly and packaging time
If you want to view the time of each stack assembling and packaging data in numerical factorization, enter the *Makefile.option* file and modify the **TIME_FLAGS** option.
Uncomment the option `TIME_FLAGS += -DPRINT_PACK_TIME` ( By default, due to too many tasks, we will turn off this option ).
Then recompile the entire library in the root directory.
`make qrtest` can get the executable file without NUMA data affinity.
`make qrtest_numa_data` can get the executable file with NUMA data affinity.
Then running *qrtest* and *qrtest_numa_data* will show the time of assembly and packaging.
Typically,`./qrtest /data/sme3Dc/sme3Dc.mtx 0` and `./qrtest_numa_data /data/sme3Dc/sme3Dc.mtx 0` is the results shown in paper.

NUMA affinity of data can reduce the time of assembly and packaging.
However, due to the opacity of the BLAS thread pool, it is difficult for us to control computationally intensive data migration, and the total time will increase slightly.
