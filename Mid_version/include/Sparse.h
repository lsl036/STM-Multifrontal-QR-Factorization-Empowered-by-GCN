#ifndef SPARSE_H
#define SPARSE_H

/* 假设支持大文件。如果出现问题，则用-DNLARGEFILE 编译 */
#ifndef NLARGEFILE

#undef  _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#undef  _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64

#endif

#include "SparseBase_config.h"
#include "SparseCore.h"
#include "SparseChol.h"

#endif
