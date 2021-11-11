# STM-Multifrontal-QR-factorization-Empowered-by-GCN
[![DOI](https://zenodo.org/badge/355837565.svg)](https://zenodo.org/badge/latestdoi/355837565)

Multifrontal QR algorithm, which consists of symbolic analysis and numerical factorization, is a high-performance algorithm for orthogonal factorizing sparse matrix.

In this work, a graph convolutional network (GCN) for adaptively selecting the optimal reordering algorithm is proposed in **symbolic analysis**.
Using our GCN adaptive classifier, the average numerical factorization time is reduced by 20.78% compared with the default approach, and the additional memory overhead is only about 4% higher than that of the previous method.

Moreover, for **numerical factorization**, an optimized tasks stream parallel processing strategy is proposed and a more efficient computing task mapping framework for NUMA architecture is adopted in this paper, which we called STM-Multifrontal QR factorization.
On Taishan Server ( 120 processors, 4 NUMA nodes), the performance of our method is equal to or higher than that of MKL library on Intel Xeon 6248 for nearly 80% matrices of the University of Florida Sparse Matrix Collection.

In this repository, we offer an highly optimized Multifrontal QR factorization algorithm in directory *STMMQR* named STM-Multifrontal QR. It is highly recommended to verify our optimization performance on a multi-NUMA multi-processor architecture.
In numerical phase, the task mapping framework with stream processing strategy improves parallelism during factorization and NUMA affinity optimization improves the factorization efficiency of each thread on NUMA architecture.
And in symbolic analysis step, our work combines the graph convolutional network (GCN) to make a self-adapting selection for several current reordering methods, whose code is stored in the directory *GCN_classifier*.

## Directory overview
For building the optimal STM-Multifrontal QR factorization package and verify the conclusions in the paper, please refer README under each directory.

Part of the drawing data shown in the paper is placed in the ***Data*** directory.
The complete code of optimal STM-Multifrontal QR is placed in the ***STMMQR*** directory.
Directory ***GCN_classifier*** contains the code and data for training and validating our GCN classifier, our model is directly modified based on the template of ***PyG package***.

## Notice
Before compiling the **STMMQR** library, we hope you can adjust the system settings by modifying *./STMMQR/include/tpsm/tpsm_sysinfo.h* according to the machine architecture :

```c
//-------CPU------->
#define TPSM_CPUSOCKET 		(2) 		//2 CPU sockets
#define TPSM_NUMANODES 		(4) 		//4 NUMA NODES
#define TPSM_SYSCORES 		(128)  		//128 CORES
```

 Otherwise, the thread-pool function may not execute correctly.
