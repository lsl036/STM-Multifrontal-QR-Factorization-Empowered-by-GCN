# STM-Multifrontal-QR-factorization-combined-with-GCN

Multifrontal QR algorithm, which consists of symbolic analysis and numerical factorization, 
is a high-performance algorithm for orthogonal factorizing sparse matrix.
In this work, a graph convolutional network (GCN) for adaptively selecting the optimal reordering 
algorithm is proposed in symbolic analysis.
Using our GCN adaptive classifier, the average numerical factorization time is reduced by 
20.78% compared with the default approach, and the additional memory overhead is only 
about 4% higher than that of the previous method.
Moreover, for numerical factorization, an optimized tasks stream parallel processing 
strategy is proposed and a more efficient computing task mapping framework for NUMA 
architecture is adopted in this paper, which called STM-Multifrontal QR factorization.
On Kunpeng 920 processors with 4 NUMA nodes, the performance of our method is equal to or higher than 
that of MKL library on Intel Xeon 6248 for nearly 80% matrices of the University of Florida Sparse 
Matrix Collection.

In this repository, we offer an highly optimized Multifrontal QR factorization named STM-Multifrontal QR.
In numerical phase, the task mapping framework with stream processing strategy improves parallelism 
during factorization and NUMA affinity optimization improves the factorization efficiency of 
each thread on NUMA architecture.
And in symbolic analysis step, our work combines the graph convolutional network (GCN) to make 
a self-adapting selection for several current reordering methods.

## Build Procedure
For building the optimal STM-Multifrontal QR factorization package and verify the conclusions in the paper, please refer README under each directory.

The complete code of STM-Multifrontal QR is placed in the *STMMQR* directory, and the *Mid_version* is the intermediate version to verify the optimization effect of task streaming scheduling in the paper.
Directory *GCN_classifier* contains the code and data for training and validating our GCN classifier, our model is directly modified based on the template in the *PyG package*.
