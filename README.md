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

## Reference
Lin S, Yang W, Wang H, et al. STM-multifrontal QR: streaming task mapping multifrontal QR factorization empowered by GCN[C]//Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis. 2021: 1-14.

DOI: https://dl.acm.org/doi/abs/10.1145/3458817.3476199

### BibTeX citation
```
@inproceedings{lin2021stm,
  title={STM-multifrontal QR: streaming task mapping multifrontal QR factorization empowered by GCN},
  author={Lin, Shengle and Yang, Wangdong and Wang, Haotian and Tsai, Qinyun and Li, Kenli},
  booktitle={Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis},
  pages={1--14},
  year={2021}
}
```

## Notice
For details about operations on ***STMMQR*** and ***GCN_classifier***, please see **README** in their respective folders.

Please note that GCN and STM-MQR are **not really combined** due to the lack of accuracy of GCN automatic classifier and the extra cost of attributes generation. And we hope that our work can inspire more researchers to combine graph convolutional network with sparse matrices problems to form more practical and valuable work.

Thank you very much!

