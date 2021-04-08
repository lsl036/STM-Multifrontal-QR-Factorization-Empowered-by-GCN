
#ifndef SparseQR_DEFINITIONS_H
#define SparseQR_DEFINITIONS_H

/* 排序选项 */
#define QR_ORDERING_FIXED 0
#define QR_ORDERING_NATURAL 1
#define QR_ORDERING_COLAMD 2
#define QR_ORDERING_GIVEN 3       
#define QR_ORDERING_CHOL 4    
#define QR_ORDERING_AMD 5         

#define QR_ORDERING_NESDIS 6

#define QR_ORDERING_DEFAULT 7     
#define QR_ORDERING_BEST 8       
#define QR_ORDERING_BESTAMD 9  
#define QR_ORDERING_METIS 10 
#define QR_ORDERING_ONLYMETIS 11

#define QR_DEFAULT_TOL (-2)       
#define QR_NO_TOL (-1)          

#define QR_QTX 0
#define QR_QX  1
#define QR_XQT 2
#define QR_XQ  3

#define QR_RX_EQUALS_B    0       /* 求解 R*X=B      or X = R\B          */
#define QR_RETX_EQUALS_B  1       /* 求解 R*E'*X=B   or X = E*(R\B)      */
#define QR_RTX_EQUALS_B   2       /* 求解 R'*X=B     or X = R'\B         */
#define QR_RTX_EQUALS_ETB 3       /* 求解 R'*X=E'*B  or X = R'\(E'*B)    */

#endif
