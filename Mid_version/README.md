
## Mid version of STM-Multifrontal QR 
This is the intermediate version code used to verify the optimization effect of the streaming scheduling strategy mentioned in the paper.

`make` can build this static library and and generate an executable file *qrtest*.
The specifications of the test function are as follows: `./qrtest /data/sme3Dc/sme3Dc.mtx`

`make distclean` can remove all compiled products.

**Note**: During the development project, the improved thread pool interface has changed. So hereby put this intermediate version of the code package to verify the optimization effect after adding thread NUMA affinity.
